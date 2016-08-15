function plot_map(part)   

load(part.path.dem);

tmpF = figure%('Visible', 'off');
tmpA = plot([dem.X(1), dem.X(end)], [dem.Y(1), dem.Y(end)]);
[lonVec, latVec, imag] = plot_google_map('Maptype', 'terrain'); hold on;                           % Plot google background
surface(dem.X, dem.Y, dem.Z./1000, prepare_google_map(dem, lonVec, latVec, imag)); shading flat
delete(tmpA);

xlabel('Longitude');
ylabel('Latitude');
zlabel('Altitude (km asl)');


%     plot_google_map('Maptype', 'terrain'); hold on;
%     surf(dem.X, dem.Y, dem.Z./1000); shading flat; alpha 0.75               % Plot topography
%     colormap(landcolor);
%     c = colorbar;
%     ylabel(c, 'Topography altitude (km asl)');
    freezeColors;
    
    
    [C, htmp] = contour(dem.X, dem.Y, ll2dist(part.vent.lat, part.vent.lon, dem.Y, dem.X));
    
    count = 1;
    contour_mat = [];
    while go == 1
        
        nEl = C(1,count);
        lev = C(2,count);
        tmp = zeros(lev,4);
        tmp(:,1) = C(2,count+1:count+lev)';
        tmp(:,2) = C(1,count+1:count+lev)';
        tmp(:,4) = repmat(nEl,lev,1);
        for i = 2:nEl+1
            dem.Z(dem.X)
        end
    end