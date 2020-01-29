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
IMPORTPATH = "/jukebox/brody/RATTER/SoloData/Settings/Chuck/",
RATIDS     = ["K" * string(id) for id = [213, 214, 215, 216, 236, 238, 241, 248,
                                         265, 280, 283, 284, 285, 289, 294, 296,
                                         297, 299, 310, 311, 313, 314, 316, 317,
                                         319, 322, 323, 324, 328, 329, 330, 331,
                                         332, 335, 336, 338, 339, 340, 341]
)

fields = Dict(
        "RewardsSection_reward_type" => "delta clicks",
    "SessionDefinition_active_stage" => 12
)

###############################################################################

"""
1. Run through all files, marking which ones fit the criteria specified
   in fields.
2. Return list of files that fulfill criteria.
"""

for irat = 1 : length(cfg.RATIDS)
    if isdir(cfg.IMPORTPATH * cfg.RATIDS[irat] * "/")
        files = searchdir(cfg.IMPORTPATH * cfg.RATIDS[irat] * "/", "@PBups")
        goodfiles = []
        for ifile = 1 : length(files)
            data = matread(cfg.IMPORTPATH * cfg.RATIDS[irat] * "/" * files[ifile])["saved"]
            if data["RewardsSection_reward_type"] == "delta clicks" && data["SessionDefinition_active_stage"] == 12
                push!(goodfiles, ifile)
            end
        end
        if !isempty(goodfiles)
            ## Find continuous set of sessions
            transitions = diff(goodfiles)
            if length(findall(transitions .> 1)) > 1
                goodfiles = goodfiles[1 : argmax(transitions .> 1)]
            end
            println("First file:  " * files[goodfiles[1]])
            println("Total files: " * length(goodfiles))
            println("Last file:   " * files[goodfiles[end]])
        end
        println("")
    end
end
