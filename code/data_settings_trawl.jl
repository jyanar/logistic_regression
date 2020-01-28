""" Trawl through settings files in order to find suitable sessions
for adding to the logistic regression analysis.

Should basically, for each rat, spit out the dates that we want to grab
for analysis. Then we can construct bdata statements that let us grab
precisely those.

Want to look at the following:
- PBupsSection_base_freq         : [int,int] if frequency task, [int] if classic task
- RewardsSection_reward_type     : 
- SessionDefinition_active_stage : 
"""

using MAT
include("utils.jl")

###############################################################################

cfg = (
IMPORTPATH = "/jukebox/brody/RATTER/SoloData/Settings/Chuck/K311/",
RATIDS     = ["K311"]
)

###############################################################################

files = searchdir(cfg.IMPORTPATH, "@PBups")

rwrds_active_stage  = zeros(length(files))
rwrds_active_scheme = []
pbups_base_freqs    = zeros(length(files), 2)

for ifile = 1 : length(files)
    data = matread(files[ifile])["saved"]
    if length(data["PBupsSection_base_freq"]) == 2
        pbups_base_freqs[ifile,:] = data["PBupsSection_base_freq"]
    end
    push!(rwrds_active_scheme, data["RewardsSection_reward_type"])
    rwrds_active_stage[ifile] = data["SessionDefinition_active_stage"]
end


















