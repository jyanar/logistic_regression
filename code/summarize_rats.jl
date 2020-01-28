""" Generates figure for each rat summarizing #Lbups/#Rbups decision
landscape, 2D model fits, and time-varying fits for stimoff and stimon
locking.
"""

using MAT
using GLM
using Dates
using Printf
using PyPlot
using PyCall
using Statistics
using DataFrames
using JLD2, FileIO
include("utils.jl")

###############################################################################

cfg = (
## io options
PROGRAM_NAME    = "summarize_rats.jl",
IMPORTPATH_DATA = "data/regrMats_allrats_500msLim_50msBin_0msOverlap.jld2",
EXPORTPATH_DATA = "data/",
SAVE_DATA       = false,

EXPORTPATH_FIGS = "figs/",
SAVE_FIGS       = true,

## analysis and plotting options
PLOT_KERNEL_MAGNITUDE = true  # Whether to plot L/R time-varying kernels'
                              # magnitudes, instead of opposite to one
                              # another.
)

###############################################################################

if cfg.SAVE_FIGS
    exportpath = mkdir_dated(cfg.EXPORTPATH_FIGS, cfg.PROGRAM_NAME)
end

## Load data and fit logit model to click totals
# data = load(cfg["IMPORTPATH_STIMON_DATA"])["regrMats"]
data = load(cfg.IMPORTPATH_DATA)["regrMats"]
nbins = data["nbins"]
nrats = data["nrats"]
xaxis_stimon  = data["xaxis_stimon"]
xaxis_stimoff = data["xaxis_stimoff"]

for irat = 1 : data["nrats"]
    println("Processing: " * string(irat) * " ...")
    ratid = data[irat]["fname"][1:4]

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
    res = get_wts_sterr(data[irat]["stimon"]["logit_rl"], nbins, true)
    if cfg.PLOT_KERNEL_MAGNITUDE lwts = res.lwts .* -1 end
    errorbar(data["xaxis_stimon"], res.rwts, yerr=1.96.*res.rsterr, label="right weights", color="b")
    errorbar(data["xaxis_stimon"], lwts, yerr=1.96.*res.lsterr, label="left weights", color="r")
    plot(data["xaxis_stimon"], zeros(length(data["xaxis_stimon"])), "k--", alpha=0.5)
    grid() ; legend()
    if cfg.PLOT_KERNEL_MAGNITUDE
        ylabel("logit weights magnitude");
    else
        ylabel("logit weights <-left -- right->")
    end
    # xlabel("Time relative to stimulus onset [ms]") ; ylabel("logit weights")
    if !cfg.PLOT_KERNEL_MAGNITUDE
        ylims = ylim() ; ylim([-maximum(abs.(ylims)), maximum(abs.(ylims))])
    end
    xlim([data["xaxis_stimon"][1]-5 , data["xaxis_stimon"][end]+5])
    title("Stim onset-locked logit weights with bias=" * @sprintf "%0.3f+/-%0.3f" res.bias res.berr);

    subplot(235)
    res = get_wts_sterr(data[irat]["stimoff"]["logit_rl"], nbins, true)
    if cfg.PLOT_KERNEL_MAGNITUDE lwts = res.lwts .* -1 end
    errorbar(data["xaxis_stimoff"], res.rwts, yerr=1.96.*res.rsterr, label="right weights", color="b")
    errorbar(data["xaxis_stimoff"], lwts, yerr=1.96.*res.lsterr, label="left weights", color="r")
    plot(data["xaxis_stimoff"], zeros(length(data["xaxis_stimoff"])), "k--", alpha=0.5)
    grid() ; legend()
    if cfg.PLOT_KERNEL_MAGNITUDE
        ylabel("logit weights magnitude");
    else
        ylabel("logit weights <-left -- right->")
    end
    xlabel("Time relative to stimulus offset/onset [ms]")
    if !cfg.PLOT_KERNEL_MAGNITUDE
        ylims = ylim() ; ylim([-maximum(abs.(ylims)), maximum(abs.(ylims))])
    end
    xlim([data["xaxis_stimoff"][1]-5 , data["xaxis_stimoff"][end]+5])
    title("Stim offset-locked logit weights with **bias**=" * @sprintf "%0.3f+/-%0.3f" res.bias res.berr);

    ## Generate rate 1D psychometric function (click diff), and model as well
    subplot(233)
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
    plot(xaxis, predict(logit_bd, DataFrame(bd=xaxis)), "k", label="model")
    plot(xaxis, propgr, "k", alpha=0.5, label="rat")
    legend(); xlabel("Click difference [#R-#L]") ; ylabel("% go right")
    title("1D psychometric")

    ## Performance stacked bar plots
    subplot(236)
    corr_right = sum(Bool.(data[irat]["wholetrl"]["X"].cr) .& Bool.(data[irat]["wholetrl"]["X"].hh))
    corr_left  = sum(Bool.(data[irat]["wholetrl"]["X"].cl) .& Bool.(data[irat]["wholetrl"]["X"].hh))
    incorr_right = sum(Bool.(data[irat]["wholetrl"]["X"].cr) .& .!Bool.(data[irat]["wholetrl"]["X"].hh))
    incorr_left  = sum(Bool.(data[irat]["wholetrl"]["X"].cl) .& .!Bool.(data[irat]["wholetrl"]["X"].hh))
    lapse_right = sum(data[irat]["violtrls"]["cr"])
    lapse_left  = sum(data[irat]["violtrls"]["cl"])

    p1 = bar(1:2, height=[corr_right, corr_left], )
    p2 = bar(1:2, height=[incorr_right, incorr_left], bottom=[corr_right, corr_left])
    p3 = bar(1:2, height=[lapse_right, lapse_left], bottom=[corr_right+incorr_right, corr_left+incorr_left],
        tick_label=["right", "left"])
    legend((p1[1],p2[1],p3[1]), ["corr", "incorr", "viol"])

    tight_layout()
    if cfg.SAVE_FIGS
        savefig(exportpath * "/" * data[irat]["fname"][1:end-8] * ".png", dpi=150);
    end
end
