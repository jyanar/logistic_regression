""" Fits logistic regression model to data 

TODO:
- [x] Plot with error bars
- [x] Make surface plots like the ones Tyler has
 - [ ] How does he assign prop right?
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
# cfg["IMPORTPATH_DATA"] = "data/regrMats_allrats_500msLim_25msBin_timeLockEnd-true_0msOverlap.jld2"
cfg["IMPORTPATH_DATA"] = "data/regrMats_allrats_500msLim_25msBin_timeLockEnd-false_0msOverlap.jld2"
cfg["EXPORTPATH_DATA"] = "data/"

################################################################################

data = load(cfg["IMPORTPATH_DATA"])
data = data["regrMats"]
for irat = 1 : data["nrats"]
    println("Fitting: " * string(irat) * " ...")

    ## Construct expressions to pass into the glm function from GLM.jl.
    ## Because it requires explicit variables, we have to use the eval
    ## and Meta.parse functions
    expr_d  = construct_logit_expr("y", ["wt_"], data[irat]["nbins"])
    expr_rl = construct_logit_expr("y", ["wtR_", "wtL_"], data[irat]["nbins"])

    logit_d  = glm(eval(Meta.parse(expr_d)), data[irat]["df_Xd"], Binomial(), LogitLink())
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
cfg = data[1]["cfg"] ; xaxis = data["xaxis"]
figure() ; plot(xaxis, zeros(length(xaxis)), "k--", alpha=0.5)
for irat = 1 : data["nrats"]
    lower = wts_d_allrats[:,irat] - cis_d_allrats[:,irat]
    upper = wts_d_allrats[:,irat] + cis_d_allrats[:,irat]

    plot(xaxis, wts_d_allrats[:,irat])
    fill_between(xaxis, lower, upper, alpha=0.5)
end
xlim([minimum(xaxis), maximum(xaxis)])
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
xlim([minimum(xaxis), maximum(xaxis)])
title("Logit regression weights, L/R model")

function plot_logit_weights(logit, xaxis)
    H = figure();
    plot()

end

