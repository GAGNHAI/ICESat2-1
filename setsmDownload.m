function setsmDownload(dataset, product, dataver, resolution, local_setsm_dir)
% Download full ArcticDEM or REMA dataset
%
%% Syntax 
% 
%   setsmDownload((dataset, product, dataver, resolution, local_setsm_dir)
%
%% User Input
%   dataset = 'REMA'; % REMA or ArcticDEM
%   product = 'mosaic'; % mosaic or geocell
%   dataver = 1; % REMA[1], ArcticDEM[3]
%   resolution = 8; % mosaic(REMA[8, 100, 200, 1000], ArcticDEM[2, 10, 32, 100, 500, 1000])
%
%% Warning 
%
% REMA: 43 TB for geocells and 1 TB for mosaics 
%
% wget must be added as an alias to your .bash_rc file after adding wget to
% homebrew
%
% first: install wget @ terminal
%       brew install wget
% second: add this to .bash_rc
%       "function _wget() { curl "${1}" -o $(basename "${1}") ; };
%       alias wget='_wget'"
%
%% Author Info
% This function was written by Alex S. Gardner, JPL-Caltech, Oct 2018. 
%
ver = sprintf('v%.1f', dataver);
if resolution < 999
    res = sprintf('%.0fm', resolution);
else
    res = sprintf('%.0fkm', resolution/1E3);
end

index_folder = fullfile(dataset, 'indexes');
unixCall = ['wget -r -N -nH -np -R index.html* --cut-dirs=3 ' fullfile('http://data.pgc.umn.edu/elev/dem/setsm/', index_folder, '/') ' -P '  fullfile(local_setsm_dir, '/')];
unix(unixCall, '-echo')

folder2mirror = fullfile(dataset, product, ver, res);
unixCall = ['wget -r -N -nH -np -R index.html* --cut-dirs=3 ' fullfile('http://data.pgc.umn.edu/elev/dem/setsm/', folder2mirror, '/') ' -P '  fullfile(local_setsm_dir, '/')];
unix(unixCall, '-echo')