%% Copy featurtracking data to server
%relativePath = 'Code/RIFT/';
%relativePath = 'Code/projects/Greenland/';
%relativePath = 'Code/matlabLibrary/';
%relativePath = 'Code/dataRetrieval';
%relativePath = 'Code/dataRetrieval/Sentinel2/';
%relativePath = 'Code/projects/ICESat2/GlobalAnalysis/';

path2target = fullfile('~', relativePath);
path2master = ['gardnera@bylot.jpl.nasa.gov:' fullfile('/u/bylot0/gardnera/', relativePath)];
unixCall = 'rsync -urltvh --progress -e ssh '; % path2master
unix([unixCall path2master ' ' path2target], '-echo')

path2target = ['gardnera@bylot.jpl.nasa.gov:' fullfile('/u/bylot0/gardnera/', relativePath)];
path2master = fullfile('~', relativePath);
unixCall = 'rsync -urltvh --progress -e ssh '; % path2master
unix([unixCall path2master ' ' path2target], '-echo')



% rsync -urltvh --progress /Volumes/MasterBrain/data/FluxGates/* gardnera@bylot.jpl.nasa.gov:/u/devon-r2/data/FluxGates/