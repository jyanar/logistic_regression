""" Preprocessing the pbups data from .mat files into usable JLD files
"""

using MAT
using JLD
include("code/utils.jl")

################################################################################

cfg = Dict()
cfg["IMPORTPATH_DATA"] = "data/frozen_noise/K283_classic_40_ALL.mat"

################################################################################

data = matread(cfg["IMPORTPATH_DATA"])
data = data["ratdata"]
