""" Fits logistic regression model to data 
"""

using MAT
using JLD
using GLM
using PyPlot
using Statistics
using DataFrames
include("utils.jl")

################################################################################

cfg = Dict()
# cfg["IMPORTPATH_DATA"] = "data/regrMats_allrats_500-lim_10-bin.jld"
cfg["IMPORTPATH_DATA"] = "data/regrMats_allrats_500msLim_25msBin_0msOverlap.jld"
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

wts_d_allrats  = zeros(data["nbins"], data["nrats"])
wts_rl_allrats = zeros(data["nbins"], 2, data["nrats"])
for irat = 1 : data["nrats"]
    wts_d_allrats[:,irat] = data[irat]["logit_d"].model.pp.beta0[2:end]
    rl_beta0 = data[irat]["logit_rl"].model.pp.beta0
    wts_rl_allrats[:,1,irat] = rl_beta0[2:data["nbins"]+1]
    wts_rl_allrats[:,2,irat] = rl_beta0[data["nbins"]+2:end]
end


