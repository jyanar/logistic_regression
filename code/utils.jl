""" Utilities for logit regression
"""

"""
    read_ratdata(importpath)

Reads in ratdata from a given file containing pbups behavioral data.
"""
function read_ratdata(importpath)
    data = matread(importpath)
    return data["ratdata"]
end

"""
    logistic_func(a)
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
