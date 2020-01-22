""" Preprocessing the pbups data from .mat files into usable JLD files.
Essentially, parses data out into DataFrames that can then easily be
fed into GLM.jl logistic regression models.

TODO
- Export matrices timelocked to both start and end of stimulus
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
cfg["IMPORTPATH_DATA"] = "data/frozen_noise/"   # Import path to .mat behavior
cfg["EXPORTPATH_DATA"] = "data/"                # Export path for processed data
cfg["STIM_WINDOW_LEN"] = 500                    # Stimulus window length, in ms
cfg["MSPERSEG"]        = 25                     # Bin width, in ms
cfg["MSOVERLAP"]       = 0                      # Bin overlap, in ms

################################################################################

filelist = searchdir(cfg["IMPORTPATH_DATA"], "_ALL.mat")

nrats = length(filelist)
nbins = Int(cfg["STIM_WINDOW_LEN"] / cfg["MSPERSEG"])

regrMats = Dict{Any,Any}("cfg" => cfg,
                       "nrats" => nrats,
                       "nbins" => nbins,
                    "filelist" => filelist)

for irat = 1 : nrats
    println("Running rat: " * filelist[irat] * " ...")

    ## Allocate Dict for rat
    regrMats[irat] = Dict()
    regrMats[irat]["fname"] = filelist[irat]
    regrMats[irat]["stimon"] = Dict()
    regrMats[irat]["stimoff"] = Dict()

    ## Read in data
    data, parsed = load_rat_behavioral_data(cfg["IMPORTPATH_DATA"] * filelist[irat])

    ## Exclude trials where stimulis less than cfg["STIM_WINDOW_LEN"]
    long_trls = parsed["pd"] .>= cfg["STIM_WINDOW_LEN"] / 1000
    for k in keys(parsed)
        if k == "b"
            parsed[k]["left"]  = parsed["b"]["left"][long_trls]
            parsed[k]["right"] = parsed["b"]["right"][long_trls]
        else
            parsed[k] = parsed[k][long_trls]
        end
    end

    ## Generate regressor matrices with bups timelocked to stimoff and stimon
    ntrls = length(parsed["pd"])
    for timelock in ["stimoff", "stimon"]
        Xd = zeros(ntrls, nbins)  ## For keeping track of bup difference
        Xr = zeros(ntrls, nbins)  ## For keeping track of right bups
        Xl = zeros(ntrls, nbins)  ## For keeping track of left bups
        Xtotal = zeros(ntrls, 2)  ## For keeping track of total L/R bups

        ## Iterate through trials, binning and storing clicktimes into
        ## the allocated X matrices
        for itrl = 1 : ntrls
            ## Grab data for this trial
            pd    = parsed["pd"][itrl] * 1000
            rbups = parsed["b"]["right"][itrl] .* 1000
            lbups = parsed["b"]["left"][itrl] .* 1000
            if isa(rbups, Float64) rbups = [rbups] end
            if isa(lbups, Float64) lbups = [lbups] end
            ## Bin left and right clicks aligned to stimulus onset and offset
            rbinned, lbinned, binlims = bin_bups_trial(pd, rbups, lbups, nbins, timelock, cfg["STIM_WINDOW_LEN"])
            Xr[itrl,:] = rbinned
            Xl[itrl,:] = lbinned
            Xd[itrl,:] = rbinned - lbinned
            Xtotal[itrl,1], Xtotal[itrl,2] = sum(rbinned), sum(lbinned)
        end

        ## Construct dataframe
        dict_X = Dict("hh" => parsed["hh"],
                      "gr" => parsed["gr"],
                  "wtRtot" => Xtotal[:,1],
                  "wtLtot" => Xtotal[:,2])
        for ibin = 1 : nbins
            dict_X["wt_"  * string(ibin)] = Xd[:,ibin]
            dict_X["wtR_" * string(ibin)] = Xr[:,ibin]
            dict_X["wtL_" * string(ibin)] = Xl[:,ibin]
        end
        regrMats[irat][timelock] = Dict()
        regrMats[irat][timelock]["X"] = DataFrame(dict_X)
    end

    ## Generate x axis (in the form of centers of time bins) for the
    ## binned data
    regrMats["binlims_stimoff"] = Array(range(-cfg["STIM_WINDOW_LEN"], 0, length=Int(nbins+1)))
    regrMats["xaxis_stimoff"] = [(regrMats["binlims_stimoff"][i]+regrMats["binlims_stimoff"][i+1])/2 for i = 1 : nbins]

    regrMats["binlims_stimon"] = Array(range(0, cfg["STIM_WINDOW_LEN"], length=Int(nbins+1)))
    regrMats["xaxis_stimon"] = [(regrMats["binlims_stimon"][i]+regrMats["binlims_stimon"][i+1])/2 for i = 1 : nbins]
end

## Save dataframes for next step in logit analysis
filename = cfg["EXPORTPATH_DATA"] * "regrMats_allrats_" *
           string(cfg["STIM_WINDOW_LEN"]) * "msLim_" *
           string(cfg["MSPERSEG"])        * "msBin_" *
           string(cfg["MSOVERLAP"]) * "msOverlap.jld2"
println("Saving: " * string(filename))
save(filename, "regrMats", regrMats)
