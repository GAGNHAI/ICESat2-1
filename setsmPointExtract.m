
%trackfile = '/Volumes/MasterBrain/data/ATM/FlightLines/Antarctic/S_OIB_2010.shp';
trackfile = '/Volumes/MasterBrain/data/FluxGates/Antarctic/moholdt_v2/v2_points_shp/GLAD_RESmod1_pts_v2.shp';
%polyfile = '/Volumes/MasterBrain/data/REMA/indexes/REMA_Strip_Index_Rel1/REMA_Strip_Index_Rel1.shp';
polyfile = '/Volumes/MasterBrain/data/REMA/indexes/REMA_Tile_Index_Rel1/REMA_Tile_Index_Rel1.shp';
uncompressAll = false;
cuturlat = 'v1.0'; %url will be trimmed keeping this directory
localdir = '/Volumes/MasterBrain/data/REMA2/mosaic';

localdir = strrep(localdir, '/Volumes/MasterBrain/data/REMA/mosaic','/gardnera/data/REMA/mosaic');

noDataValue = -9999;

updateDEMs= false;

if updateDEMs
    product = 'REMA_tiles';
    version = 'v1.0';
    localFolder = '/Volumes/MasterBrain/data';
    setsmDownload(product, version, localFolder)
end


EPSG = 3031;

% load in point data
pointdata = shaperead(trackfile);
x = [pointdata.X];
y = [pointdata.Y];
h = [pointdata(:).h_bedm2];


n = length(pointdata);
clear pointdata
% load data
poly = shaperead(polyfile);

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
m = length(poly);

idx = false(m,n);

% initialize file paths
path2file = cell([m,1]);
% sep = find(demlocal == '/', 2, 'last');
% if sep(end) == length(demlocal)
%     linkdir = demlocal((sep(1)+1):(sep(2)-1));
%     demlocal = demlocal(1:sep(1)-1);
% else
%     linkdir = demlocal((sep(2)+1):end);
%     demlocal = demlocal(1:sep(2)-1);
% end

% determine file limits and size from info in poly file
foo = [poly.BoundingBox];
XWorldLimits = foo(:,1:2:end)';
YWorldLimits = foo(:,2:2:end)';
PixelSize = vertcat(poly(:).resolution);
RasterSize = [(XWorldLimits(:,2)- XWorldLimits(:,1)), (YWorldLimits(:,2)- YWorldLimits(:,1))] ./ PixelSize;
FileUrl = {poly(:).fileurl};

%DateN = vertcat(poly(:).acquisitio);
%DateN = datenum(DateN,'YYYYMMDD');

clear('poly')
% for each polygon find intersecting point ata
% loop for each polygon
% PARFOR

%% low ram option [takes ~2-3s for 170k points, 4 cores]
parfor i = 1:m
    xLims = XWorldLimits(i,:);
    yLims = YWorldLimits(i,:);
    idx(i,:) = x > xLims(1) & x < xLims(2) & y > yLims(1) & y < yLims(2);
    foo = FileUrl{i};
    path2file{i} = foo(strfind(foo, ['/' cuturlat '/']):end);
end

% ok now we have all points crossed with all files... next up... open files
% and extract data

%% parfor
tic
dt = datenum(2000,1,1);
elev = nan(size(idx));
err = elev;
day = elev;


% parfor
for i = 1:length(path2file) % [takes ~23 min for 170k points, 4 cores]
    idx0 = idx(i,:);
    
    if uncompressAll || any(idx0) 
        tarfile = fullfile(localdir, path2file{i});
        untarfold = tarfile(1:(strfind(tarfile, '.tar.gz')-1));
        [~,b] = fileparts(untarfold);
        
        if exist(tarfile,'file') == 2 && exist(fullfile(untarfold,[b '_dem.tif']), 'file') ~= 2
            untar(tarfile, untarfold);
            
            %% do not delete, .tar.gz file is needed for syncing with PGS server
            % delete(tarfile); 
            %% -----------------------------------------------------------------
           
        elseif exist(fullfile(untarfold,[b '_dem.tif']), 'file') ~= 2
            fprintf('skipping, %s does not exist \n', fullfile(untarfold,[b '_dem.tif']))
            continue
        end
        
        % move file to SSD for testing
        if 1 == 0
      
            return
            move2 = strrep(untarfold, '/Volumes/MasterBrain/data/REMA','/Volumes/MasterBrain/data/REMA2');

            d = dir(fullfile(move2, '*dem.tif'));
            if isempty(d)
            mkdir(move2)
            copyfile(untarfold, move2,'f')
            end
            disp(i)
            continue
        end
        
        % points that intersect
        x0 = x(idx0);
        y0 = y(idx0);
        elev0 = nan(size(x0));
        err0 = elev0;
        day0 = elev0;
        
        [r,c] = worldToIntrinsicArray(x(idx0),y(idx0),XWorldLimits(i,:),YWorldLimits(i,:), RasterSize(i,:));
        
        for j = 1:length(r)
            elev0(j) = imread(fullfile(untarfold,[b '_dem.tif']), 'PixelRegion', {[r(j) r(j)], [c(j), c(j)]});
            err0(j) = imread(fullfile(untarfold,[b '_err.tif']), 'PixelRegion', {[r(j) r(j)], [c(j), c(j)]});
            
            foo = double(imread(fullfile(untarfold,[b '_day.tif']), 'PixelRegion', {[r(j) r(j)], [c(j), c(j)]}));
            foo(isnan(foo)) = nan;
            day0(j) = foo + dt;
        end
        
        elev1 = elev(i,:);
        err1 = err(i,:);
        day1 = day(i,:);
        
        elev1(idx0) = elev0;
        err1(idx0) = err0;
        day1(idx0) = day0;
        
        elev(i,:) = elev1;
        err(i,:) = err1;
        day(i,:) = day1;
    end
end
toc

% replace no data value
notValid = elev == noDataValue;
elev(notValid) = nan;
err(notValid) = nan;
day(notValid) = nan;

% are there more than one elevation per point?
if any(sum(~isnan(elev)) >1)
    error('more than one elevation point returned')
else
    elev = nanmean(elev);
    err = nanmean(err);
    day = nanmean(day);
end

% check data
x = h;
y = elev;
histogram2(x,y,'DisplayStyle','tile','ShowEmptyBins','off', 'XBinLimits', [0 1000], 'YBinLimits', [0 1000], 'BinWidth', 10);

mean(y-x, 'omitnan')
std(y-x, 'omitnan')

plot(x); hold on
plot(y)

