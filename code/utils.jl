""" Utilities for logit regression
"""

using MAT

""" One-liner for getting list of files from a dir """
searchdir(path,key) = filter(file->occursin(key,file), readdir(path))


""" Load clean rat behavioral data from .mat files """
function load_rat_behavioral_data(importpath::String)::Tuple{Dict, Dict}
    data = matread(importpath)
    data = data["ratdata"]
    parsed = data["parsed"]
    for k in keys(parsed)
        parsed[k] = parsed[k][1]
    end
    return data, parsed
end


""" bin_bups_trial(pd, rbups, lbups, nbins, timelock, stim_period_lim)
Bins left and right clicks into arrays of length nbins, timelocked to
either the onset or offset of the stimulus.

Parameters:
- pd::Float64: Duration of the trial, in ms
- rbups::Array{Float64}: Right clicktimes, in ms
- lbups::Array{Float64}: Left clicktimes, in ms
- nbins::Int64: Number of bins with which to tile space between 0 and
    stim_period_lim
- timelock::String: What to timelock to. Options:
    - "stimoff": Timelock to stimulus offset. Clicks between
        [stimoffset-stim_period_lim , stimoffset) are binned
    - "stimon": Timelock to the stimulus onset. Clicks between
        [0, stim_period_lim] are binned
- stim_window_len::Int64: Stimulus period window in which to bin clicks

Returns:
- rbinned::Array{Int64}: Binned right clicktimes
- lbinned::Array{Int64}: Binned left clicktimes
- binlims::Array{Int64}: Lower & upper limits for bins
"""
function bin_bups_trial(pd::Float64,
                     rbups::Array{Float64},
                     lbups::Array{Float64},
                     nbins::Int64,
                  timelock::String,
           stim_window_len::Int64)::Tuple{Array{Int64}, Array{Int64}, Array{Int64}}
    rbinned = zeros(Int64, nbins)
    lbinned = zeros(Int64, nbins)
    ## Align rbup/lbup times if necessary, and cut off bups that occur
    ## outside of window specified by timelock and stim_period_lim
    if timelock == "stimoff"
        time_difference = pd - stim_window_len
        rbups = rbups .- time_difference
        lbups = lbups .- time_difference
        rbups = rbups[rbups .> 0]
        lbups = lbups[lbups .> 0]
    elseif timelock == "stimon"
        rbups = rbups[rbups .< stim_window_len]
        lbups = lbups[lbups .< stim_window_len]
    end
    ## Construct time bin limits and bin the data
    binlims = Int.(Array(range(0, stim_window_len, length=nbins+1)))
    for ibin = 1 : nbins
        ## Currently: inclusive start, exclusive end
        rbinned[ibin] = length(findall(b -> b >= binlims[ibin] && b < binlims[ibin+1], rbups))
        lbinned[ibin] = length(findall(b -> b >= binlims[ibin] && b < binlims[ibin+1], lbups))
    end
    if timelock == "stimoff" binlims = binlims .- stim_window_len end
    return rbinned, lbinned, binlims
end


""" construct_logit_expr
Constructs string to be parsed and evaluated as first arg to glm, such
as:
    @formula(y ~ wt_1 + wt_2)
Parameters:
- dep_var::String: Dependent variable for the logit model. e.g., "y"
- regr_prefixes::Array{String}: Prefixes of regressors to the model.
    e.g., ["wtR_", "wtL_"] or ["wt_"]
- nregr::Int64: Number of regressors for each prefix

Returns:
- expr::String: Expression for the logistic regression

Examples:
    > expr = construct_logit_expr("y", ["wtR_", "wtL_"], 2)
    "@formula(y ~ wtR_1 + wtR_2 + wtL_1 + wtL_2)"
    > logit = glm(eval(Meta.parse(expr)), df, Binomial(), LogitLink())
"""
function construct_logit_expr(dep_var::String,
                        regr_prefixes::Array{String},
                                nregr::Int64)::String
    expr = "@formula(" * dep_var * " ~ ";
    for iprefix = 1 : length(regr_prefixes)
        for iregr = 1 : nregr
            if iprefix == length(regr_prefixes) && iregr == nregr
                ## If we're at the end of the expression, end it with ")"
                expr = expr * regr_prefixes[iprefix] * string(iregr) * ")"
            else
                expr = expr * regr_prefixes[iprefix] * string(iregr) * " + "
            end
        end
    end
    return expr
end


""" bd_gr_surface_matrix(rb, lb, gr)
Computes matrix whose axes are bup amounts for left and right,
with magnitude for proportion went right.
Params:
- rb::Array{Int64}: Total right bups per trial
- lb::Array{Int64}: Total left bups per trial
- gr::Array{Int64}: Whether animal went right (1) or left (0)
"""
function bd_gr_surface_matrix(rb::Array{Int64},
                              lb::Array{Int64},
                              gr::Array{Int64})::Array{Float64}
    maxbups = maximum([rb ; lb])
    nclicks = 0 : 1 : maxbups
    matr = zeros(Float64, length(nclicks), length(nclicks))
    for nr = nclicks
        for nl = nclicks
            ## For each (nr, nl) pair, find all trials that match
            matching_trials = (rb .== nr) .& (lb .== nl)
            ## And compute the proportion of these trials where rat
            ## went right
            matr[nr+1,nl+1] = length(findall(gr[matching_trials] .== 1)) / length(gr[matching_trials])
        end
    end
    return matr
end


""" logistic_func(a)
Logistic function, mapping ℜ→[0,1)
"""
function logistic_func(a)
    return 1 / (1 + exp(a))
end


"""
    compute_hh(total_diff, y)
Computes hit history from total trial evidences in total_diff and
responses in y. Returns hit history.
"""
function compute_hh(total_diff, y)
    hh = zeros(length(y))
    for i = 1 : length(total_diff)
        # A right trial
        if total_diff[i] > 0
            if y[i] == 1 hh[i] = 1 # Rat went right
            else         hh[i] = 2 # Rat went left
            end
        # A left trial
        elseif total_diff[i] < 0
            if y[i] == 0 hh[i] = 1 # Rat went left
            else         hh[i] = 2 # Rat went right
            end
        end
    end
    return hh
end


function group_trials_on_evidence(rat_struct, NGROUPS)
    # Compute group ranges based on distribution of
    # trial total evidences
    if abs(maximum(rat_struct["totaldiff"])) > abs(minimum(rat_struct["totaldiff"]))
        rat_struct["grps_ranges"] = Array(range(-maximum(rat_struct["totaldiff"]),
                                                 maximum(rat_struct["totaldiff"]), length=NGROUPS+1))
    else
        rat_struct["grps_ranges"] = Array(range(minimum(rat_struct["totaldiff"]),
                                                abs(minimum(rat_struct["totaldiff"])), length=NGROUPS+1))
    end

    # Assign each trial to a group
    rat_struct["grps_idx"] = []
    for i = 1 : length(rat_struct["grps_ranges"]) - 1
        push!(rat_struct["grps_idx"], findall(x -> x >= rat_struct["grps_ranges"][i] &&
                                                   x <  rat_struct["grps_ranges"][i+1], rat_struct["totaldiff"]))
    end

    # Compute proportion correct of each trial group
    rat_struct["prop_corr"] = zeros(length(rat_struct["grps_idx"]))
    for i = 1 : length(rat_struct["prop_corr"])
        n_corr   = length(findall(rat_struct["hh"][rat_struct["grps_idx"][i]] .== 1))
        n_trials = length(rat_struct["grps_idx"][i])
        rat_struct["prop_corr"][i] = n_corr / n_trials
    end

    # And compute the proportion that went right, for easy psychometric curve
    rat_struct["prop_right"] = zeros(length(rat_struct["grps_idx"]))
    for i = 1 : length(rat_struct["prop_right"])
        n_right  = length(findall(rat_struct["y"][rat_struct["grps_idx"][i]] .== 1))
        n_trials = length(rat_struct["grps_idx"][i])
        rat_struct["prop_right"][i] = n_right / n_trials
    end
    return rat_struct
end
