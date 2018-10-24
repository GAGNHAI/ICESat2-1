function [elev, err, day, dem_idx, path2file]  = setsmPointExtract(x,y, dataset, local_setsm_dir, varargin)
% This function extracts point data from SETSM DEMs
%
%% Syntax 
% 
%  [r,c] = setsmPointExtract(x,y, dataset, local_setsm_dir)
%  [r,c] = setsmPointExtract(x,y, dataset, local_setsm_dir, 'Name', value)
%
%% User Input
%
%   x = real world location in x, same coodinate system as setsem dataset
%	y = real world location in y, same coodinate system as setsem dataset
%	dataset = REMA or ArcticDEM
%   location of local SETSM dataset (this is where SETSM will download)
%
% optional name pair inputs
%   no_data_value = [-9999]
%   saveClass = ['single']
%   product = ['mosaic']
%   dataver = 1/3
%   resolution = 8/10
%   download_dataset = [false]
%   uncompress_all = [false]
%% Author Info
% This function was written by Alex S. Gardner, JPL-Caltech, Oct 2018. 

%% set defaults
p = inputParser;
addParameter(p,'no_data_value',-9999);
addParameter(p,'saveClass','single');
addParameter(p,'download_dataset',false);
addParameter(p,'uncompress_all',false);

switch dataset
    case 'REMA'
        addParameter(p,'product','mosaic') % mosaic or geocell
        addParameter(p,'dataver',1) % REMA[1], ArcticDEM[3]
        addParameter(p,'resolution',8) % mosaic(REMA[8, 100, 200, 1000], ArcticDEM[2, 10, 32, 100, 500, 1000])
    case 'ArcticDEM'
        addParameter(p,'product','mosaic') % mosaic or geocell
        addParameter(p,'dataver',3) % REMA[1], ArcticDEM[3]
        addParameter(p,'resolution',10) % mosaic(REMA[8, 100, 200, 1000], ArcticDEM[2, 10, 32, 100, 500, 1000])
    otherwise
        error('unknown product')
end

parse(p,varargin{:});


%% download full dataset [if requested]
if p.Results.download_dataset
    setsmDownload(dataset, p.Results.product, p.Results.dataver, p.Results.resolution, local_setsm_dir)
end

%% construct path 2 index file
switch p.Results.product
    case 'mosaic'
        idx_prod = 'Tile';
    case 'geocell'
        idx_prod = 'Strip';
    otherwise
        error('product not recognized')
end

% find latest release
foo = fullfile(local_setsm_dir, dataset, 'indexes', sprintf('%s_%s_Index_Rel*', dataset, idx_prod));
foo = dir(foo);
rel = nan(size(foo));
for i = 1:length(foo)
    sIdx = strfind(foo(i).name,'_Rel')+4;
    eIdx = strfind(foo(i).name,'.');
    eIdx = eIdx(1)-1;
    rel(i) = str2double(foo(i).name(sIdx:eIdx));
end

% select most recent version of index file
idx_file = fullfile(local_setsm_dir, dataset, 'indexes', sprintf('%s_%s_Index_Rel%.0f', dataset, idx_prod, max(rel)));
    
% look to see if .shp file exists
if exist([idx_file, '.shp'], 'file') ~= 2
    if exist([idx_file, '.zip'], 'file') == 2
        % unzip index file
        unzip([idx_file, '.zip'], fullfile(local_setsm_dir, dataset, 'indexes'))
    else
        error('index file not found')
    end
end

% load index file of polygons
poly = shaperead([idx_file, '.shp']);

% determine file limits and size from info in poly file
foo = [poly.BoundingBox];
XWorldLimits = foo(:,1:2:end)';
YWorldLimits = foo(:,2:2:end)';
PixelSize = vertcat(poly(:).resolution);
RasterSize = [(XWorldLimits(:,2)- XWorldLimits(:,1)), (YWorldLimits(:,2)- YWorldLimits(:,1))] ./ p.Results.resolution;

FileUrl = {poly(:).fileurl};
m = length(poly);
clear('poly')

is_geocell = strcmp(p.Results.product, 'geocell');
if  is_geocell
    % date is contained in metadata (no seperate grid like mosaic)
    DateN = vertcat(poly(:).acquisitio);
    DateN = cast(datenum(DateN,'YYYYMMDD'), saveClass);
end

%%
cut_url_at = sprintf('v%.1f', p.Results.dataver); %url will be trimmed keeping this directory
localdir = fullfile(local_setsm_dir, dataset);
n = length(x);

% convert X and Y to epsg
% if EPSG == 3031
%     % is in polar stereographic projection
%     [x, y] = ...
%         polarstereo_fwd(y, x, ...
%         6378137.0 ,0.08181919, -71, 0);
% else
%     [x, y] = latlon2utm(y, ...
%         x, num2str(zone0), 'WGS84');
% end

% initialize index

idx = false(m,n);

% initialize file paths
path2file = cell([m,1]);

% for each polygon find intersecting point ata
% loop for each polygon
% PARFOR
% low ram option [takes ~2-3s for 170k points, 4 cores]
parfor i = 1:m
    xLims = XWorldLimits(i,:);
    yLims = YWorldLimits(i,:);
    idx(i,:) = x > xLims(1) & x < xLims(2) & y > yLims(1) & y < yLims(2);
    foo = FileUrl{i};
    path2file{i} = foo(strfind(foo, ['/' cut_url_at '/']):end);
    
    if strcmp(dataset, 'ArcticDEM') && p.Results.dataver == 3
        % path is not correct
        path2file{i} = fileparts(path2file{i});
        [~,filePrefx] = fileparts(path2file{i});
        path2file{i} = fullfile(path2file{i}, sprintf('%s_%.0fm_v3.0.tar.gz', filePrefx, p.Results.resolution));
        path2file{i} = strrep(path2file{i},'2m', sprintf('%.0fm',p.Results.resolution));
    end   
end

% ok now we have all points crossed with all files... next up... open files
% and extract data

%% parfor
dt = datenum(2000,1,1);
elev = nan(size(idx), p.Results.saveClass);
err = elev;
day = elev;
dem_idx = zeros(size(err), 'single');

% parfor
for i = 1:length(path2file) % [takes ~23 min for 170k points, 4 cores]
    idx0 = idx(i,:);
    
    if p.Results.uncompress_all || any(idx0) 
        tarfile = fullfile(localdir, p.Results.product, path2file{i});
        untarfold = tarfile(1:(strfind(tarfile, '.tar.gz')-1));
        
        demFile = dir(fullfile(untarfold, '*dem.tif'));
        if ~isempty(demFile)
            demFile = fullfile(untarfold,demFile.name);
        else
            demFile = false;
        end
        
        if exist(tarfile,'file') == 2 && demFile(1) == 0
            fprintf('uncomressing: %s\n', tarfile)
            untar(tarfile, untarfold);
            
            %% do not delete, .tar.gz file is needed for syncing with PGS server
            % delete(tarfile);
            %% -----------------------------------------------------------------
            
            demFile = dir(fullfile(untarfold, '*dem.tif'));
            if ~isempty(demFile)
                demFile = fullfile(untarfold,demFile.name);
            else
                demFile = false;
            end 
            
        elseif exist(tarfile,'file') ~= 2
            fprintf('skipping, %s does not exist \n', tarfile)
            continue
        end
        
        % points that intersect
        elev0 = nan(size(x(idx0)));
        err0 = elev0;
        day0 = elev0;
        dem_idx0 = zeros(size(elev0), 'single');
        [r,c] = worldToIntrinsicArray(x(idx0),y(idx0),XWorldLimits(i,:),YWorldLimits(i,:), RasterSize(i,:));
        
        errFile = dir(fullfile(untarfold, '*err.tif'));
        if ~isempty(errFile)
            errFile = fullfile(untarfold,errFile.name);
        else
            errFile = false;
        end
        
        dayFile = dir(fullfile(untarfold, '*day.tif'));
        if ~isempty(dayFile)
            dayFile = fullfile(untarfold,dayFile.name);
        else
            dayFile = false;
        end
        
        % **** to further spped up code find unique r,c locaitons ****
        for j = 1:length(r)
            elev0(j) = cast(imread(demFile, 'PixelRegion', {[r(j) r(j)], [c(j), c(j)]}), p.Results.saveClass);
            dem_idx0(j) = i;
            if errFile(1) ~= 0
                err0(j) = cast(imread(errFile, 'PixelRegion', {[r(j) r(j)], [c(j), c(j)]}), p.Results.saveClass);
            else
                err0(j) = nan;
            end
            
            if  is_geocell
                day0(j) = DateN(i);
            else
                if dayFile(1) ~= 0
                    foo = cast((imread(dayFile, 'PixelRegion', {[r(j) r(j)], [c(j), c(j)]})), p.Results.saveClass);
                    foo(isnan(foo)) = nan;
                    day0(j) = foo + dt;
                else
                    day0(j) = nan;
                end
            end
        end
        
        % require for parallel loop
        elev1 = elev(i,:);
        err1 = err(i,:);
        day1 = day(i,:);
        dem_idx1 = dem_idx(i,:);
        
        elev1(idx0) = elev0;
        err1(idx0) = err0;
        day1(idx0) = day0;
        dem_idx1(idx0) = dem_idx0;
        
        elev(i,:) = elev1;
        err(i,:) = err1;
        day(i,:) = day1;
        dem_idx(i,:) = dem_idx1;
    end
end

% replace no data value
notValid = elev == p.Results.no_data_value;
elev(notValid) = nan;
err(notValid) = nan;
day(notValid) = nan;
dem_idx(notValid) = nan;

% are there more than one elevation per point?
if any(sum(~isnan(elev)) >1)
    error('more than one elevation point returned')
else
    elev = nanmean(elev);
    err = nanmean(err);
    day = nanmean(day);
    dem_idx = nanmean(dem_idx);
end

if size(x,2) == 1
     elev = elev';
    err = err';
    day = day';
    dem_idx = dem_idx';
end