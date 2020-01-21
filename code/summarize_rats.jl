""" Generates figure for each rat summarizing #Lbups/#Rbups decision
landscape, 2D model fits, and time-varying fits for stimoff and stimon
locking.
"""

using MAT
using GLM
using PyPlot
using Statistics
using DataFrames
using JLD2, FileIO
include("utils.jl")

###############################################################################

cfg = Dict()
# cfg["IMPORTPATH_DATA"] = "data/regrMats_allrats_500-lim_10-bin.jld"
cfg["IMPORTPATH_DATA"] = "data/regrMats_allrats_500msLim_25msBin_timeLockEnd-true_0msOverlap.jld2"
# cfg["IMPORTPATH_DATA"] = "data/regrMats_allrats_500msLim_25msBin_timeLockEnd-false_0msOverlap.jld2"
cfg["EXPORTPATH_DATA"] = "data/"

###############################################################################

## Load data and fit logit model to click totals
data = load(cfg["IMPORTPATH_DATA"])
data = data["regrMats"]
for irat = 1 : data["nrats"]
    println("Processing: " * string(irat) * " ...")

    ## Fitting
    ## Construct expressions for glm: "y ~ wtR_1 + wtR_2 + ..."
    expr_d  = construct_logit_expr("y", ["wt_"], data[irat]["nbins"])
    expr_rl = construct_logit_expr("y", ["wtR_", "wtL_"], data[irat]["nbins"])
    logit_d     = glm(eval(Meta.parse(expr_d)), data[irat]["df_Xd"], Binomial(), LogitLink())
    logit_rl    = glm(eval(Meta.parse(expr_rl)), data[irat]["df_Xrl"], Binomial(), LogitLink())
    logit_total = glm(@formula(gr ~ wtR + wtL), data[irat]["df_Xtotal"], Binomial(), LogitLink())
    data[irat]["logit_d"]  = logit_d
    data[irat]["logit_rl"] = logit_rl
    data[irat]["logit_total"] = logit_total

    ## Generate model 2D psychometric surface
    maxbups = maximum([data[irat]["df_Xtotal"].wtL ; data[irat]["df_Xtotal"].wtR])
    nclicks = 0 : 1 : maxbups
    model_matr = zeros(length(nclicks), length(nclicks))
    for nr = nclicks
        for nl = nclicks
            model_matr[Int(nr+1),Int(nl+1)] = predict(logit_total, DataFrame(wtR=Int(nr), wtL=Int(nl)))[1]
        end
    end
    figure(figsize=(16,10))
    subplot(234) ; imshow(model_matr, origin="lower", cmap=get_cmap("RdBu"))
    xlabel("#R Clicks") ; ylabel("#L Clicks") ; colorbar(); clim([0, 1])
    title("2D psychometric, model")

    ## Generate rat 2D psychometric surface
    rat_matr = bd_gr_surface_matrix(data[irat]["df_Xtotal"].wtR, data[irat]["df_Xtotal"].wtL, data[irat]["df_Xtotal"].gr)
    subplot(231) ; imshow(rat_matr, origin="lower", cmap=get_cmap("RdBu"))
    xlabel("#R Clicks") ; ylabel("#L Clicks") ; colorbar(); clim([0, 1])
    title("2D psychometric, rat")

    ## Compute L/R time-varying functions
    ## TODO: should we just compute the stimoff and stimon dfs and keep them in
    ## the same file for easy processing & plotting later?

    ## Generate rate 1D psychometric function (click diff), and model as well
    bd = data[irat]["df_Xtotal"].wtR - data[irat]["df_Xtotal"].wtL
    gr = data[irat]["df_Xtotal"].gr
    df = DataFrame(gr=data[irat]["df_Xtotal"].gr, bd=bd)
    logit_bd = glm(@formula(gr ~ bd), df, Binomial(), LogitLink())
    maxbd = maximum(abs.(bd)) ; xaxis = -maxbd : 1 : maxbd

    propgr = zeros(length(xaxis))
    for i = 1 : length(propgr)
        matching_trials = bd .== xaxis[i]
        propgr[i] = length(findall(gr[matching_trials] .== 1))/length(gr[matching_trials])
    end
    subplot(132)
    plot(xaxis, predict(logit_bd, DataFrame(bd=xaxis)), "k", label="model")
    plot(xaxis, propgr, "k", alpha=0.5, label="rat")
    suptitle(data[irat]["fname"])
end

