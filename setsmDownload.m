function setsmDownload(product, version, localFolder)
% Created by: 
% Alex Gardner, NASA-JPL, California Institute of Technology
% email: alex.s.gardner@jpl.nasa.gov

% User Inputs
% Example 1:
% product = 'ArcticDEM_tiles';
% version = 'v2.0'
% localFolder = '/Volumes/MasterBrain/data';
% 
% Example 2
% product = 'REMA_tiles';
% version = 'v1.0'
% localFolder = '/Volumes/MasterBrain/data';
%

%% ------------------------------------------------------------------------
% NOTES 
%
% !!! for REMA: 43 TB for geocells and 1 TB for mosaics !!!!
%
% wget must be added as an alias to your .bash_rc file after adding wget to
% homebrew
%
% first: install wget @ terminal
%       brew install wget
% second: add this to .bash_rc
%       "function _wget() { curl "${1}" -o $(basename "${1}") ; };
%       alias wget='_wget'"
%% ------------------------------------------------------------------------

switch [product version]
    case ['REMA_tiles' 'v1.0']
        % copy all data from requested folder
        % https://docs.google.com/document/d/1XlSk1wK_KHaYSdKVp3Fq_tJA8gh0v2YqNi8OANj7_IQ/view
        % elevations are w.r.t the WGS84 ellipsoid
        folder2mirror = 'REMA/mosaic/v1.0/8m';
        indexFiles = 'REMA/indexes';
        
    case ['ArcticDEM_tiles' 'v2.0']
        folder2mirror = 'ArcticDEM/mosaic/v2.0';
        indexFiles = 'ArcticDEM/indexes';
    otherwise
        error('unrecognized product and/or version')
end


unixCall = ['wget -r -N -nH -np -R index.html* --cut-dirs=3 ' fullfile('http://data.pgc.umn.edu/elev/dem/setsm/', folder2mirror, '/') ' -P '  fullfile(localFolder)];
unix(unixCall, '-echo')

unixCall = ['wget -r -N -nH -np -R index.html* --cut-dirs=3 ' fullfile('http://data.pgc.umn.edu/elev/dem/setsm/', indexFiles, '/') ' -P '  fullfile(localFolder)];
unix(unixCall, '-echo')