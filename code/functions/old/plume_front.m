
%% Plume front
front.s          = logspace(1,6,50); % Spacing, find a way to adapt
front.t(1)       = 0;
front.x(1)       = 0;
front.y(1)       = 1;
front.z(1)       = erup.ht*1e3;
front.lat(1)     = vent.lat;
front.lon(1)     = vent.lon;

[~, front.xI(1)] = min(abs(lon-front.lon(1)));
[~, front.yI(1)] = min(abs(lat-front.lat(1)));
[~, front.tI(1)] = min(abs(time-(erup.time+front.t(1))));
[~, front.zI(1)] = min(abs(alt(front.yI(1), front.xI(1), :, front.tI(1)) - front.z(1)));

% Interpolation
alt_vec          = squeeze(alt(front.yI(1), front.xI(1), :, front.tI(1)));
front.u(1)       = interpn(lat,lon,alt_vec,time,u,front.lat(1),front.lon(1),front.z(1),erup.time+front.t(1)/3600/24);
front.v(1)       = interpn(lat,lon,alt_vec,time,v,front.lat(1),front.lon(1),front.z(1),erup.time+front.t(1)/3600/24);

front.wind(1)    = sqrt(front.u(1)^2 + front.v(1)^2);

[front.ub(1), front.width(1)] = get_gravity_current(front.s(1), erup.ht-vent.ht/1e3, front.wind(1));

front.wxT(1)     = front.x(1);
front.wyT(1)     = front.y(1);
front.wxB(1)     = front.x(1);
front.wyB(1)     = front.y(1);
front.wlatT(1)   = front.lat(1);
front.wlonT(1)   = front.lon(1);
front.wlatB(1)   = front.lat(1);
front.wlonB(1)   = front.lon(1);

front.rho(1)    = 0;

for iF = 2:length(front.s)
    
    dS          = front.s(iF) - front.s(iF-1);
    
    vel_tot     = front.wind(iF-1) + front.ub(iF-1);
    dtF         = dS/vel_tot;
    dir         = atan2d(front.v(iF-1),front.u(iF-1));
    
    front.x(iF) = front.x(iF-1) + dS*cosd(dir);
    front.y(iF) = front.y(iF-1) + dS*sind(dir);
    front.z(iF) = front.z(1);%-front.s(iF-1)*(5e3/6e5);
    front.t(iF) = front.t(iF-1) + dtF;
    
    [front.ub(iF), front.width(iF)]              = get_gravity_current(front.s(iF), erup.ht-vent.ht/1e3, front.wind(iF-1));
    [front.wxT(iF), front.wyT(iF), front.wxB(iF), front.wyB(iF)] = get_plume_edge(dir, front.width(iF), front.x(iF), front.y(iF));
    
    
    % Update geographic coordinates
    [front.lat(iF), front.lon(iF)]         = dist2ll(front.lat(iF-1), front.lon(iF-1), dS*cosd(dir), dS*sind(dir));
    [front.wlatT(iF), front.wlonT(iF)]     = dist2ll(front.lat(iF), front.lon(iF), front.wxT(iF)-front.x(iF), front.wyT(iF)-front.y(iF));
    [front.wlatB(iF), front.wlonB(iF)]     = dist2ll(front.lat(iF), front.lon(iF), front.wxB(iF)-front.x(iF), front.wyB(iF)-front.y(iF));
    % Update indices
    [~, front.xI(iF)] = min(abs(lon-front.lon(iF)));
    [~, front.yI(iF)] = min(abs(lat-front.lat(iF)));
    [~, front.tI(iF)] = min(abs(time-(erup.time + front.t(iF)/3600/24)));
    [~, front.zI(iF)] = min(abs(alt(front.yI(iF), front.xI(iF), :, front.tI(iF)) - front.z(iF)));
 
    % Get new plume velocity
    if strcmp(interpolation, 'none')
        %Simple Scheme without interpolation
        front.u(iF)       = u(front.yI(iF), front.xI(iF), front.zI(iF), front.tI(iF));
        front.v(iF)       = v(front.yI(iF), front.xI(iF), front.zI(iF), front.tI(iF));

    elseif strcmp(interpolation, 'complete')
        % Complete interpolation method (Slow)
        alt_sub            = squeeze(alt(front.yI(iF), front.xI(iF), :, front.tI(iF)));
        front.u(iF)        = interpn(lat,lon,alt_sub,time,u,front.lat(iF),front.lon(iF),front.z(iF),erup.time+front.t(iF)/3600/24);
        front.v(iF)        = interpn(lat,lon,alt_sub,time,v,front.lat(iF),front.lon(iF),front.z(iF),erup.time+front.t(iF)/3600/24);

    elseif strcmp(interpolation, 'subset')
        if front.xI(iF)-subsetI < 1 || front.xI(iF)+subsetI > size(lon,1) || front.yI(iF)-subsetI < 1 || front.yI(iF)+subsetI > size(lat,1) || front.zI(iF)-subsetI < 1 || front.zI(iF)+subsetI > size(alt,3)
            break
        end

        alt_sub   = squeeze(alt(front.yI(iF), front.xI(iF), front.zI(iF)-subsetI:front.zI(iF)+subsetI, front.tI(iF)));
        time_sub  = time (front.tI(iF)-subsetI:front.tI(iF)+subsetI);
        lon_sub   = lon(front.xI(iF)-subsetI:front.xI(iF)+subsetI);
        lat_sub   = lat(front.yI(iF)-subsetI:front.yI(iF)+subsetI);
        u_sub     = u(front.yI(iF)-subsetI:front.yI(iF)+subsetI,front.xI(iF)-subsetI:front.xI(iF)+subsetI,front.zI(iF)-subsetI:front.zI(iF)+subsetI,front.tI(iF)-subsetI:front.tI(iF)+subsetI);
        v_sub     = v(front.yI(iF)-subsetI:front.yI(iF)+subsetI,front.xI(iF)-subsetI:front.xI(iF)+subsetI,front.zI(iF)-subsetI:front.zI(iF)+subsetI,front.tI(iF)-subsetI:front.tI(iF)+subsetI);
        rho_sub   = rhoair(front.yI(iF)-subsetI:front.yI(iF)+subsetI,front.xI(iF)-subsetI:front.xI(iF)+subsetI,front.zI(iF)-subsetI:front.zI(iF)+subsetI,front.tI(iF)-subsetI:front.tI(iF)+subsetI);

        
        
        front.u(iF)        = interpn(lat_sub,lon_sub,alt_sub,time_sub,u_sub,front.lat(iF),front.lon(iF),front.z(iF),erup.time+front.t(iF)/3600/24);
        front.v(iF)        = interpn(lat_sub,lon_sub,alt_sub,time_sub,v_sub,front.lat(iF),front.lon(iF),front.z(iF),erup.time+front.t(iF)/3600/24);
        front.rho(iF)      = interpn(lat_sub,lon_sub,alt_sub,time_sub,rho_sub,front.lat(iF),front.lon(iF),front.z(iF),erup.time+front.t(iF)/3600/24);

    end
        
    front.wind(iF)   = sqrt(front.u(iF)^2 + front.v(iF)^2);
     
end


% tTMP    = [front.t, front.t, front.t];
% latTMP  = [front.lat, front.wlatB, front.wlatT]; 
% lonTMP  = [front.lon, front.wlonB, front.wlonT];
% [Xq,Yq,Tq] = griddata(lonTMP, latTMP, tTMP./3600, dem.X, dem.Y);



% figure; 
% plot(-360+[lon_min, lon_max], [lat_min, lat_max], '.k');                    % Plot frame
% plot_google_map('Maptype', 'terrain'); hold on;                             % Plot google background
% surf(-360+dem.X, dem.Y, dem.Z./1000); shading flat; alpha 0.5               % Plot topography
% colormap(landcolor);
% c = colorbar;
% ylabel(c, 'Topography altitude (km asl)');
% freezeColors;
% [C,h] = contour(-360+dem.X, dem.Y, ll2dist(vent.lat, vent.lon, dem.Y, dem.X));
% set(h, 'LineColor', 'k', 'Linewidth', 1.5);
% clabel(C,h);
% 
% plot3(-360+front.lon, front.lat, front.z./1000, '-r', 'LineWidth', 2); 
% plot3(-360+front.wlonT, front.wlatT, front.z./1000, '-k', 'LineWidth', 2); 
% plot3(-360+front.wlonB, front.wlatB, front.z./1000, '-k', 'LineWidth', 2); 
% 
% % contour3(-360+Xq,Yq,Tq); alpha 0.5

