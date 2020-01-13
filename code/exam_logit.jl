# Messing around with the final exam stuff from NEU314, as it
# forms a good problem for constructing a simple implementation
# of a logit model.
#
# X - Matrix containing click differences in each 50 ms timebin. Specifically,
#     #right - #left. Thus a trial with mostly right clicks will be positive,
#     whereas a trial with mostly left clicks will be negative.
# y - 1, went right. 0, went left
# hh - 1, correct. 2, incorrect.

## PARAMS #####################################################################

NGROUPS = 9

###############################################################################

using JLD
using GLM
using PyPlot
using Statistics

include("utils.jl")

data = load("data/P2_data.jld")
ratA = data["rat_A"]
ratB = data["rat_B"]

ratA["totaldiff"] = sum(ratA["X"], dims=2)
ratB["totaldiff"] = sum(ratB["X"], dims=2)
ratA["cumdiff"] = cumsum(ratA["X"], dims=2)
ratB["cumdiff"] = cumsum(ratB["X"], dims=2)

## Hit history - whether the trial was correctly answered
ratA["hh"] = compute_hh(ratA["totaldiff"], ratA["y"])
ratB["hh"] = compute_hh(ratB["totaldiff"], ratB["y"])

## Parse trials into groups based on total evidence (#R-#L)
ratA = group_trials_on_evidence(ratA, NGROUPS)
ratB = group_trials_on_evidence(ratB, NGROUPS)

## Compute
