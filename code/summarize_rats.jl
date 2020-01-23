""" Generates figure for each rat summarizing #Lbups/#Rbups decision
landscape, 2D model fits, and time-varying fits for stimoff and stimon
locking.
"""

using MAT
using GLM
using PyPlot
using PyCall
using Statistics
using DataFrames
using JLD2, FileIO
include("utils.jl")
# sns = pyimport("seaborn")
# sns.set_style("whitegrid")

###############################################################################

cfg = Dict()
cfg["IMPORTPATH_DATA"] = "data/regrMats_allrats_500msLim_50msBin_0msOverlap.jld2"
cfg["EXPORTPATH_DATA"] = "data/"

###############################################################################

## Load data and fit logit model to click totals
# data = load(cfg["IMPORTPATH_STIMON_DATA"])["regrMats"]
data = load(cfg["IMPORTPATH_DATA"])["regrMats"]
nbins = data["nbins"]
nrats = data["nrats"]
xaxis_stimon  = data["xaxis_stimon"]
xaxis_stimoff = data["xaxis_stimoff"]

for irat = 1 : data["nrats"]
    println("Processing: " * string(irat) * " ...")

    ## Fitting
    ## Construct expressions for glm: "y ~ wtR_1 + wtR_2 + ...", fit to
    ## data locked to the stimulus onset and offset, as well as to the
    ## whole trial bups
    expr_d  = construct_logit_expr("gr", ["wtd_"], data["nbins"])
    expr_rl = construct_logit_expr("gr", ["wtR_", "wtL_"], data["nbins"])
    for timelock in ["stimoff", "stimon"]
        logit_d = glm(eval(Meta.parse(expr_d)), data[irat][timelock]["X"], Binomial(), LogitLink())
        logit_rl = glm(eval(Meta.parse(expr_rl)), data[irat][timelock]["X"], Binomial(), LogitLink())
        logit_total = glm(@formula(gr ~ wtRtot + wtLtot), data[irat][timelock]["X"], Binomial(), LogitLink())
        data[irat][timelock]["logit_d"] = logit_d
        data[irat][timelock]["logit_rl"] = logit_rl
        data[irat][timelock]["logit_total"] = logit_total

        # data[irat][timelock]["lwts_rl"] = 
    end
    data[irat]["wholetrl"]["logit"] = glm(@formula(gr ~ wtRtot + wtLtot),
                                          data[irat]["wholetrl"]["X"], Binomial(), LogitLink())

    ## Construct figure
    figure(figsize=(16,10))

    ## RAT PSYCHOMETRIC SURFACE
    subplot(231)
    rat_decision_surf = bd_gr_surface_matrix(data[irat]["wholetrl"]["X"].wtRtot,
                                             data[irat]["wholetrl"]["X"].wtLtot,
                                             data[irat]["wholetrl"]["X"].gr)
    imshow(rat_decision_surf, origin="lower", cmap=get_cmap("RdBu"))
    xlabel("#R Clicks") ; ylabel("#L Clicks") ; colorbar() ; clim([0, 1])
    title("2d psychometric, rat")

    ## 2D MODEL PSYCHOMETRIC SURFACE
    subplot(234)
    maxbups = maximum([data[irat]["wholetrl"]["X"].wtLtot ;
                       data[irat]["wholetrl"]["X"].wtRtot])
    nclicks = 0 : 1 : maxbups
    model_decision_surf = zeros(length(nclicks), length(nclicks))
    for nr = nclicks
        for nl = nclicks
            model_decision_surf[nl+1,nr+1] = predict(data[irat]["wholetrl"]["logit"], 
                                                     DataFrame(wtRtot=Int(nr),wtLtot=Int(nl)))[1]
        end
    end
    imshow(model_decision_surf, origin="lower", cmap=get_cmap("RdBu"))
    xlabel("#R Clicks") ; ylabel("#L Clicks") ; colorbar(); clim([0, 1])
    title("2D psychometric, model")

    ## TIME-VARYING L/R MODEL, STIM ONSET
    subplot(232)
    # betas = data[irat]["stimon"]["logit_rl"].model.pp.beta0
    # stdm  = stderror(data[irat]["stimon"]["logit_rl"])
    # rwts_stimon = betas[2:1+nbins] ; lwts_stimon = betas[2+nbins:end]
    rwts_stimon = data[irat]["stimon"]["logit_rl"].model.pp.beta0[2:1+nbins]
    lwts_stimon = data[irat]["stimon"]["logit_rl"].model.pp.beta0[2+nbins:end]
    plot(data["xaxis_stimon"], rwts_stimon)
    plot(data["xaxis_stimon"], lwts_stimon)
    plot(data["xaxis_stimon"], zeros(length(data["xaxis_stimon"])), "k--", alpha=0.5)
    grid()
    xlabel("Time relative to stimulus onset [ms]") ; ylabel("logit weights")
    ylims = ylim() ; ylim([-maximum(abs.(ylims)), maximum(abs.(ylims))])
    xlim([data["xaxis_stimon"][1] , data["xaxis_stimon"][end]])

    subplot(235)
    rwts_stimoff = data[irat]["stimoff"]["logit_rl"].model.pp.beta0[2:1+nbins]
    lwts_stimoff = data[irat]["stimoff"]["logit_rl"].model.pp.beta0[2+nbins:end]
    plot(data["xaxis_stimoff"], rwts_stimoff)
    plot(data["xaxis_stimoff"], lwts_stimoff)
    plot(data["xaxis_stimoff"], zeros(length(data["xaxis_stimon"])), "k--", alpha=0.5)
    grid()
    xlabel("Time relative to stimulus offset [ms]") ; ylabel("logit weights")
    ylims = ylim() ; ylim([-maximum(abs.(ylims)), maximum(abs.(ylims))])
    xlim([data["xaxis_stimoff"][1] , data["xaxis_stimoff"][end]])

    ## Generate rate 1D psychometric function (click diff), and model as well
    bd = data[irat]["wholetrl"]["X"].wtRtot - data[irat]["wholetrl"]["X"].wtLtot
    gr = data[irat]["wholetrl"]["X"].gr
    df = DataFrame(gr=data[irat]["wholetrl"]["X"].gr, bd=bd)
    logit_bd = glm(@formula(gr ~ bd), df, Binomial(), LogitLink())
    maxbd = maximum(abs.(bd))
    xaxis = -maxbd : 1 : maxbd

    propgr = zeros(length(xaxis))
    for i = 1 : length(propgr)
        matching_trials = bd .== xaxis[i]
        propgr[i] = length(findall(gr[matching_trials] .== 1))/length(gr[matching_trials])
    end
    subplot(133)
    plot(xaxis, predict(logit_bd, DataFrame(bd=xaxis)), "k", label="model")
    plot(xaxis, propgr, "k", alpha=0.5, label="rat")
    suptitle(data[irat]["fname"])
end
