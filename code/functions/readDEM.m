function [X,Y,Z] = readDEM(file)
% For use with LagTrack

% Read the DEM header
fid_1=fopen(file);
head_dem=textscan(fid_1,'%s',12);
fclose(fid_1);
% ncols_dem=str2double(head_dem{1}{2});
% nrows_dem=str2double(head_dem{1}{4});
xllcorner   = round(str2double(head_dem{1}{6}));
yllcorner   = round(str2double(head_dem{1}{8}));
cellsize    = str2double(head_dem{1}{10});
NODATA      = str2double(head_dem{1}{12});

Z           = dlmread(file,' ',6,0);
Z           = Z(:,1:end-1);
Z(Z==NODATA)= nan;

[X, Y]      = meshgrid(xllcorner:cellsize:xllcorner+(size(Z,2)-1)*cellsize,...
    yllcorner:cellsize:yllcorner+ (size(Z,1)-1)*cellsize);

%figure; surf((X-X(1))*10,(Y-Y(1))*10,flipud(Z)); shading interp;  colormap(gray);axis equal;