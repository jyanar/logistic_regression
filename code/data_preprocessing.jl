""" Preprocessing the pbups data from .mat files into usable JLD files.
Essentially, parses data out into DataFrames that can then easily be
fed into GLM.jl logistic regression models.

TODO
- Timelock to the end or the beginning of the stimulus
- Different time bin sizes
- Can choose whether to exclude trials based on stim length
- Export time bin centers for plotting later
"""

using MAT
using DataFrames
using JLD2, FileIO
include("utils.jl")

################################################################################

cfg = Dict()
cfg["IMPORTPATH_DATA"]   = "data/frozen_noise/"
cfg["EXPORTPATH_DATA"]   = "data/"
cfg["STIM_LENGTH_LIM"]   = 500
cfg["TIMELOCK_STIM_END"] = false
cfg["MSPERSEG"]          = 25
cfg["MSOVERLAP"]         = 0

################################################################################

filelist = searchdir(cfg["IMPORTPATH_DATA"], "_ALL.mat")

nrats = length(filelist)
nbins = Int(cfg["STIM_LENGTH_LIM"] / cfg["MSPERSEG"])

regrMats = Dict()

for irat = 1 : nrats
    regrMats[irat] = Dict()
    println("Running rat: " * filelist[irat] * " ...")

    ## Read in data and organize parsed data
    data = matread(cfg["IMPORTPATH_DATA"] * filelist[irat])
    data = data["ratdata"]
    parsed = data["parsed"]
    for k in keys(parsed)
        parsed[k] = parsed[k][1]
    end

    ## Exclude the trials where stimulus is less than 500 ms in length
    ## Essentially, we're going to look at the 500 ms before the stim ends.
    long_trls = parsed["pd"] .>= cfg["STIM_LENGTH_LIM"] / 1000
    for k in keys(parsed)
        if k == "b"
            parsed[k]["left"]  = parsed["b"]["left"][long_trls]
            parsed[k]["right"] = parsed["b"]["right"][long_trls]
        else
            parsed[k] = parsed[k][long_trls]
        end
    end

    ## Organize click times into time bins for L and R.
    ntrls = length(parsed["g"])

    ## Components needed for the logistic regression model
    y  = parsed["gr"]
    Xd = zeros(ntrls, nbins)  ## For keeping track of bup difference
    Xr = zeros(ntrls, nbins)  ## For keeping track of right bups
    Xl = zeros(ntrls, nbins)  ## For keeping track of left bups

    ## Bin clicks from all trials into our regressor matrices
    for itrl = 1 : ntrls
        ## Grab information for this trial
        yi = parsed["gr"][itrl]
        pd = parsed["pd"][itrl] * 1000
        rbups = parsed["b"]["right"][itrl] .* 1000
        lbups = parsed["b"]["left"][itrl]  .* 1000

        ## Allocate arrays for holding binned data
        rbinned = zeros(Int(nbins))
        lbinned = zeros(Int(nbins))
        bupdiff = zeros(Int(nbins))

        ## Time-lock to the end of the stimulus. Throw out bups that occur
        ## earlier than -cfg["STIM_LENGTH_LIM"]
        if cfg["TIMELOCK_STIM_END"]
            time_difference = pd - cfg["STIM_LENGTH_LIM"]
            rbups = rbups .- time_difference
            lbups = lbups .- time_difference
            if isa(rbups, Float64) rbups = [rbups] end
            if isa(lbups, Float64) lbups = [lbups] end
            rbups = rbups[rbups .> 0]
            lbups = lbups[lbups .> 0]

            ## Construct time bin limits
            binlims = Int.(Array(range(0, cfg["STIM_LENGTH_LIM"], length=Int(nbins+1))))
            for ibin = 1 : nbins
                ## Can try inclusive start, exclusive end, later try switching this
                rbinned[ibin] = length(findall(b -> b >= binlims[ibin] && b < binlims[ibin+1], rbups))
                lbinned[ibin] = length(findall(b -> b >= binlims[ibin] && b < binlims[ibin+1], lbups))
            end
        ## Otherwise, time-lock to the beginning of the stimulus. Throw out
        ## bups that occur later than cfg["STIM_LENGTH_LIM"]
        else
            if isa(rbups, Float64) rbups = [rbups] end
            if isa(lbups, Float64) lbups = [lbups] end
            rbups = rbups[rbups .< cfg["STIM_LENGTH_LIM"]]
            lbups = lbups[lbups .< cfg["STIM_LENGTH_LIM"]]

            ## Construct time bin limits
            binlims = Int.(Array(range(-cfg["STIM_LENGTH_LIM"], 0, length=Int(nbins+1))))
            for ibin = 1 : nbins
                ## Can try inclusive start, exclusive end, later try switching this
                rbinned[ibin] = length(findall(b -> b >= binlims[ibin] && b < binlims[ibin+1], rbups))
                lbinned[ibin] = length(findall(b -> b >= binlims[ibin] && b < binlims[ibin+1], lbups))
            end
            

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

    regrMats[irat]["fname"]   = filelist[irat]
    regrMats[irat]["cfg"]     = cfg
    regrMats[irat]["nrats"]   = nrats
    regrMats[irat]["nbins"]   = nbins
    regrMats[irat]["df_Xrl"]  = df_Xrl
    regrMats[irat]["df_Xd"]   = df_Xd
    regrMats[irat]["binlims"] = binlims
    regrMats[irat]["xaxis"]   = 
end
regrMats["nrats"] = nrats
regrMats["nbins"] = nbins

## Save dataframes for next step in logit analysis
save(cfg["EXPORTPATH_DATA"] * "regrMats_allrats_" *
    string(cfg["STIM_LENGTH_LIM"]) * "msLim_" * 
    string(cfg["MSPERSEG"])        * "msBin_" *
    "timeLockEnd-" * string(cfg["TIMELOCK_STIM_END"]) * "_" *
    string(cfg["MSOVERLAP"]) * "msOverlap.jld2", "regrMats", regrMats)


# function bin_bups_data(cfg, b, pd, nbins)
#     for itrl = 1 : ntrls
#         pd = pd[itrl] * 1000
#         rbups = b["right"][itrl] .* 1000
#         lbups = b["left"][itrl]  .* 1000
#         ## Allocate arrays for holding binned data
#         rbinned = zeros(Int(nbins))
#         lbinned = zeros(Int(nbins))
#         bupdiff = zeros(Int(nbins))
#         ## Re-express bups as times relative to -500 before the stimulus end
#         ## Take the difference between stimulus length and 500
#     end
#     return Xr, Xl, Xd
# end
