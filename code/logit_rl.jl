""" Fits logistic regression model to data 

TODO:
- Plot with error bars
- Make surface plots like the ones Tyler has
"""

using MAT
using GLM
using PyPlot
using Statistics
using DataFrames
using JLD2, FileIO
include("utils.jl")

################################################################################

cfg = Dict()
# cfg["IMPORTPATH_DATA"] = "data/regrMats_allrats_500-lim_10-bin.jld"
cfg["IMPORTPATH_DATA"] = "data/regrMats_allrats_500msLim_25msBin_timeLockEnd-false_0msOverlap.jld2"
cfg["EXPORTPATH_DATA"] = "data/"

################################################################################

data = load(cfg["IMPORTPATH_DATA"])
data = data["regrMats"]
for irat = 1 : data["nrats"]
    println("Fitting: " * string(irat) * " ...")

    ## Hacky way of constructing expressions needed for the @formula
    ## macro in the glm function from GLM.jl
    expr_d = "@formula(y ~ wt_1"
    for iwt = 2 : data[irat]["nbins"]
        expr_d = expr_d * " + wt_" * string(iwt)
    end
    expr_d = expr_d * ")" 

    expr_rl = "@formula(y ~ wtR_1"
    for iwt = 2 : data[irat]["nbins"]
        expr_rl = expr_rl * " + wtR_" * string(iwt)
    end
    for iwt = 1 : data[irat]["nbins"]
        expr_rl = expr_rl * " + wtL_" * string(iwt)
    end
    expr_rl = expr_rl * ")"

    logit_d = glm(eval(Meta.parse(expr_d)), data[irat]["df_Xd"], Binomial(), LogitLink())
    logit_rl = glm(eval(Meta.parse(expr_rl)), data[irat]["df_Xrl"], Binomial(), LogitLink())
    data[irat]["logit_d"]  = logit_d
    data[irat]["logit_rl"] = logit_rl
end

## Gather weights from all rats into single matrix, as well
## as the 95% CIs
wts_d_allrats  = zeros(data["nbins"], data["nrats"])
wts_rl_allrats = zeros(data["nbins"], 2, data["nrats"])
cis_d_allrats  = zeros(data["nbins"], data["nrats"])
cis_rl_allrats = zeros(data["nbins"], 2, data["nrats"])
for irat = 1 : data["nrats"]
    wts_d_allrats[:,irat] = data[irat]["logit_d"].model.pp.beta0[2:end]
    cis_d_allrats[:,irat] = 1.96 * stderror(data[irat]["logit_d"])[2:end]

    rl_beta0 = data[irat]["logit_rl"].model.pp.beta0[2:end]
    cis_95   = 1.96 * stderror(data[irat]["logit_rl"])[2:end]
    wts_rl_allrats[:,1,irat] = rl_beta0[1:data["nbins"]]
    wts_rl_allrats[:,2,irat] = rl_beta0[data["nbins"]+1:end]
    cis_rl_allrats[:,1,irat] = cis_95[1:data["nbins"]]
    cis_rl_allrats[:,2,irat] = cis_95[data["nbins"]+1:end]
end

## Plot weights from click difference model
cfg = data[1]["cfg"]
if cfg["TIMELOCK_STIM_END"]
    xaxis = range(-cfg["STIM_LENGTH_LIM"], 0, length=data["nbins"])
else
    xaxis = range(0, cfg["STIM_LENGTH_LIM"], 0, length=data["nbins"])
end

figure() ; plot(xaxis, zeros(length(xaxis)), "k--", alpha=0.5)
for irat = 1 : data["nrats"]
    lower = wts_d_allrats[:,irat] - cis_d_allrats[:,irat]
    upper = wts_d_allrats[:,irat] + cis_d_allrats[:,irat]

    plot(xaxis, wts_d_allrats[:,irat])
    fill_between(xaxis, lower, upper, alpha=0.5)
end
xlim([minimum(xaxis), maximum(xaxis)]); ylim([-0.25, 0.25])
title("Logit regression weights, click difference model")

## Plot weights from L/R model
figure() ; plot(xaxis, zeros(length(xaxis)), "k--", alpha=0.5)
for irat = 1 : data["nrats"]
    lowerR = wts_rl_allrats[:,1,irat] - cis_rl_allrats[:,1,irat]
    upperR = wts_rl_allrats[:,1,irat] + cis_rl_allrats[:,1,irat]
    lowerL = wts_rl_allrats[:,2,irat] - cis_rl_allrats[:,2,irat]
    upperL = wts_rl_allrats[:,2,irat] + cis_rl_allrats[:,2,irat]

    plot(xaxis, wts_rl_allrats[:,1,irat], "*-", color="C" * string(irat))
    fill_between(xaxis, lowerR, upperR, alpha=0.5, color="C" * string(irat))

    plot(xaxis, wts_rl_allrats[:,2,irat], "*-", color="C" * string(irat))
    fill_between(xaxis, lowerL, upperL, alpha=0.5, color="C" * string(irat))
end
xlim([minimum(xaxis), maximum(xaxis)]); ylim([-0.25, 0.25])
title("Logit regression weights, L/R model")

