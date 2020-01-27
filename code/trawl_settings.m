% Trawls through Chuck's settings in 
% /jukebox/brody/RATTER/SoloData/Settings/Chuck to find suitable
% sessions for the pbups logistic regression analysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg.IMPORTPATH = '/jukebox/brody/RATTER/SoloData/Settings/Chuck/'
cfg.RATIDS = {'K311'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for irat = 1 : length(cfg.RATIDS)
    fprintf('Rat ID: %s\n', cfg.RATIDS{irat});
    files = fullfile(dir(cfg.IMPORTPATH + "*@PBups*"));
    for ifile = 1 : length(files)
        disp(files(ifile));
    end
end


