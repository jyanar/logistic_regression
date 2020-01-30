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
PROGRAM_NAME    = "presentation_fig.jl",
IMPORTPATH_DATA = "data/regrMats_logitFits_allrats_500msLim_50msBin_0msOverlap.jld2",
EXPORTPATH_DATA = "data/",
SAVE_DATA       = false,

EXPORTPATH_FIGS = "figs/",
SAVE_FIGS       = true,

## analysis and plotting options
PLOT_KERNEL_MAGNITUDE = true,  # Whether to plot L/R time-varying kernels'
                               # magnitudes, instead of opposite to one
                               # another.
PLOT_BUPDIFF_KERNEL   = false  # Whether to plot the time-varying click
                               # difference kernel as well
)

###############################################################################

data = load(cfg.IMPORTPATH_DATA)["data"]

## PLOT 1D OF ALL RATS
figure();
for irat = 1 : data["nrats"]
    bd = data[irat]["wholetrl"]["X"].wtRtot - data[irat]["wholetrl"]["X"].wtLtot
    gr = data[irat]["wholetrl"]["X"].gr
    df = DataFrame(gr=data[irat]["wholetrl"]["X"].gr, bd=bd)
    logit_bd = glm(@formula(gr ~ bd), df, Binomial(), LogitLink())
    res = grp_trls_evidence(data[irat]["wholetrl"]["X"], 10)
    xaxis = [(res.grps_ranges[i]+res.grps_ranges[i+1])/2 for i = 1 : length(res.grps_ranges)-1]
    # plot(xaxis, res.prop_gr, "o--", color="k", alpha=0.5)
    plot(xaxis, res.prop_gr, "o--")
    grid() ; xlabel("Click difference [#R - #L]") ; ylabel("proportion go right")
    title("1D psychometric, all rats")
end
savefig("figs/1dpsych_allrats.png", dpi=300)
