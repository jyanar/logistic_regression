""" Preprocessing the pbups data from .mat files into usable JLD files.
Essentially, parses data out into DataFrames that can then easily be
fed into GLM.jl logistic regression models.

TODO
- Timelock to the end or the beginning of the stimulus
- Different time bin sizes
- Can choose whether to exclude trials based on stim length
"""

using MAT
using JLD
using DataFrames
include("utils.jl")

################################################################################

cfg = Dict()
cfg["IMPORTPATH_DATA"] = "data/frozen_noise/"
cfg["EXPORTPATH_DATA"] = "data/"
cfg["STIM_LENGTH_LIM"] = 500
cfg["MSPERSEG"]        = 25
cfg["MSOVERLAP"]       = 0

################################################################################

regrMats = Dict()
filelist = searchdir(cfg["IMPORTPATH_DATA"], "_ALL.mat")
nrats = length(filelist)
nbins = Int(cfg["STIM_LENGTH_LIM"] / cfg["MSPERSEG"])
for irat = 1 : nrats
    println("Running rat: " * filelist[irat] * " ...")

    # Read in data and organize parsed data
    data = matread(cfg["IMPORTPATH_DATA"] * filelist[irat])
    data = data["ratdata"]
    parsed = data["parsed"]
    for k in keys(parsed)
        parsed[k] = parsed[k][1]
    end

    # 1. Exclude the trials where stimulus is less than 500 ms in length
    # Essentially, we're going to look at the 500 ms before the stim ends.
    long_trls = parsed["pd"] .>= cfg["STIM_LENGTH_LIM"] / 1000
    for k in keys(parsed)
        if k == "b"
            parsed[k]["left"]  = parsed["b"]["left"][long_trls]
            parsed[k]["right"] = parsed["b"]["right"][long_trls]
        else
            parsed[k] = parsed[k][long_trls]
        end
    end

    # 2. Organize click times into time bins for L and R.
    ntrls = length(parsed["g"])

    ## Components needed for the logistic regression model
    y  = parsed["gr"]
    Xd = zeros(ntrls, nbins)  # For keeping track of bup difference
    Xr = zeros(ntrls, nbins)  # For keeping track of right bups
    Xl = zeros(ntrls, nbins)  # For keeping track of left bups

    ## Bin clicks from all trials into our regressor matrices
    for itrl = 1 : ntrls
        ## Grab information for this trial
        yi = parsed["gr"][itrl]
        pd = parsed["pd"][itrl] * 1000
        rbups = parsed["b"]["right"][itrl] .* 1000
        lbups = parsed["b"]["left"][itrl] .* 1000

        ## Allocate arrays for holding binned data
        rbinned = zeros(Int(nbins))
        lbinned = zeros(Int(nbins))
        bupdiff = zeros(Int(nbins))

        # Re-express bups as times relative to -500 before the stimulus end
        # Take the difference between stimulus length and 500.
        time_difference = pd - cfg["STIM_LENGTH_LIM"]
        rbups = rbups .- time_difference
        lbups = lbups .- time_difference
        if isa(rbups, Float64) rbups = [rbups] end
        if isa(lbups, Float64) lbups = [lbups] end
        rbups = rbups[rbups .> 0]
        lbups = lbups[lbups .> 0]

        # Construct time bin limits
        binlims = Int.(Array(range(0, 500, length=Int(nbins+1))))
        for ibin = 1 : nbins
            # Can try inclusive start, exclusive end, later try switching this
            rbinned[ibin] = length(findall(b -> b >= binlims[ibin] && b < binlims[ibin+1], rbups))
            lbinned[ibin] = length(findall(b -> b >= binlims[ibin] && b < binlims[ibin+1], lbups))
        end
        Xr[itrl,:] = rbinned
        Xl[itrl,:] = lbinned
        Xd[itrl,:] = rbinned - lbinned
    end

    ## Construct dataframes
    dict_Xrl = Dict("y" => y)
    dict_Xd  = Dict("y" => y)
    for ibin = 1 : nbins
        dict_Xd["wt_" * string(ibin)] = Xd[:,ibin]
        dict_Xrl["wtR_" * string(ibin)] = Xr[:,ibin]
        dict_Xrl["wtL_" * string(ibin)] = Xl[:,ibin]
    end
    df_Xrl = DataFrame(dict_Xrl)
    df_Xd  = DataFrame(dict_Xd)

    regrMats[irat] = Dict("fname" => filelist[irat], 
                            "cfg" => cfg,
                          "nrats" => nrats,
                          "nbins" => nbins,
                         "df_Xrl" => df_Xrl,
                          "df_Xd" => df_Xd);
end
regrMats["nrats"] = nrats
regrMats["nbins"] = nbins

# Save dataframes for next step in logit analysis
save(cfg["EXPORTPATH_DATA"] * "regrMats_allrats_" *
    string(cfg["STIM_LENGTH_LIM"]) * "msLim_" * 
    string(cfg["MSPERSEG"]) * "msBin_" *
    string(cfg["MSOVERLAP"]) * "msOverlap.jld", "regrMats", regrMats)
