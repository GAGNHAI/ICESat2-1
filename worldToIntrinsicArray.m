function [r,c] = worldToIntrinsicArray(x,y,XWorldLimits,YWorldLimits, RasterSize)

m = length(x);
n = size(XWorldLimits,1);

r = zeros(m,n,'uint16');
c = r;

x0 = reshape(x, [m,1]) * ones(1,n);
y0 = reshape(y, [m,1]) * ones(1,n);

overlaping = x0 >= ones(m,1)* (XWorldLimits(:,1)') & ...
    x0 <= ones(m,1) * (XWorldLimits(:,2)') & ...
    y0 >= ones(m,1)* (YWorldLimits(:,1)') & ...
    y0 <= ones(m,1) * (YWorldLimits(:,2)');

% find just those data that have atleast 1 overalping points
idxRaster = any(overlaping,1);
overlaping = overlaping(:,idxRaster);

XWorldLimits = XWorldLimits(idxRaster',:);
YWorldLimits = YWorldLimits(idxRaster',:);
RasterSize = RasterSize(idxRaster',:);

r0 = zeros(length(x),size(XWorldLimits,1),'uint16');
c0 = r0;

parfor i = 1:size(XWorldLimits,1)
    R = maprefcells(XWorldLimits(i,:), YWorldLimits(i,:), RasterSize(i,:), 'ColumnsStartFrom','north'); 
    [c2, r2] = R.worldToIntrinsic(x(overlaping(:,i)),y(overlaping(:,i)));
    
    r1 = r0(:,i);
    c1 = c0(:,i);
    
    r1(overlaping(:,i)) = r2;
    c1(overlaping(:,i)) = c2;
    
    r0(:,i) = r1;
    c0(:,i) = c1;
end


r(:,idxRaster) = r0;
c(:,idxRaster) = c0;