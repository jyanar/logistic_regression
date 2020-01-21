""" Utilities for logit regression
"""

""" One-liner for getting list of files from a dir
"""
searchdir(path,key) = filter(x->occursin(key,x), readdir(path))


function load_rat_behavioral_data(importpath)
    data = matread(importpath)
    data = data["ratdata"]
    parsed = data["parsed"]
    for k in keys(parsed)
        parsed[k] = parsed[k][1]
    end
    return data, parsed
end


""" construct_logit_expr
Constructs string to be parsed and evaluated as first arg to glm, such
as:
    @formula(y ~ wt_1 + wt_2)
Parameters:
    dep_var: string, dependent variable for the logit model. e.g., "y"
    regr_prefixes: list of strings, prefixes of regressors to the
        model. e.g., ["wtR_", "wtL_"] or ["wt_"]
    nregr: integer, number of regressors for each prefix

Examples:
    > expr = construct_logit_expr("y", ["wtR_", "wtL_"], 2)
    "@formula(y ~ wtR_1 + wtR_2 + wtL_1 + wtL_2)"
    > logit = glm(eval(Meta.parse(expr)), df, Binomial(), LogitLink())
"""
function construct_logit_expr(dep_var, regr_prefixes, nregr)
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
- rb : Vector of length ntrials. Total right bups
- lb : Vector of length ntrials. Total left bups
- gr : Vector of length ntrials. Whether animal went right
"""
function bd_gr_surface_matrix(rb, lb, gr)#, hh)
    maxbups = maximum([rb ; lb])
    nclicks = 0 : 1 : maxbups
    matr = zeros(length(nclicks), length(nclicks)) .+ NaN
    for nr = nclicks
        for nl = nclicks
            ## For each (nr, nl) pair, find all trials that match
            matching_trials = (rb .== nr) .& (lb .== nl)
            ## And compute the proportion of these trials where rat
            ## went right
            matr[Int(nr+1),Int(nl+1)] = length(findall(gr[matching_trials] .== 1)) / length(gr[matching_trials])
        end
    end
    return matr
end


""" read_ratdata(importpath)
Reads in ratdata from a given file containing pbups behavioral data.
"""
function read_ratdata(importpath)
    data = matread(importpath)
    return data["ratdata"]
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
