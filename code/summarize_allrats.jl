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
PROGRAM_NAME    = "summarize_allrats.jl",
IMPORTPATH_DATA = "data/regrMats_logitFits_allrats_500msLim_50msBin_0msOverlap.jld2",
EXPORTPATH_FIGS = "figs/",
SAVE_FIGS       = false,
)

###############################################################################

if cfg.SAVE_FIGS
    exportpath = mkdir_dated(cfg.EXPORTPATH_FIGS, cfg.PROGRAM_NAME)
end

## Load data and fit logit model to click totals
data = load(cfg.IMPORTPATH_DATA)["data"]
nbins = data["nbins"]
nrats = data["nrats"]
xaxis_stimon  = data["xaxis_stimon"]
xaxis_stimoff = data["xaxis_stimoff"]

perf = [sum(data[irat]["wholetrl"]["X"].hh)/data[irat]["ntrls"] for irat = 1 : data["nrats"]]
rlbias = [data[irat]["wholetrl"]["logit"].model.pp.beta0[2]/data[irat]["wholetrl"]["logit"].model.pp.beta0[end] for irat = 1 : data["nrats"]]
bias = [data[irat]["wholetrl"]["logit"].model.pp.beta0[1] for irat = 1 : data["nrats"]]

wts_bd_off = zeros(data["nbins"])
wts_bd_on  = zeros(data["nbins"])
for irat = 1 : length(data["nrats"])
    wts_bd_off = wts_bd_off + data[irat]["stimoff"]["logit_d"].model.pp.beta0[2:end]
    wts_bd_on  = wts_bd_on  + data[irat]["stimon"]["logit_d"].model.pp.beta0[2:end]
end
wts_bd_off = wts_bd_off ./ data["nrats"]
wts_bd_on  = wts_bd_on  ./ data["nrats"]

figure()
plot(xaxis_stimon, wts_bd_on)
plot(xaxis_stimoff, wts_bd_off)
