function part = get_trajectory(P)

% UPDATES
% 2016/02/02 (Seb): Added bearing

% TO DO:
% - RungeKutta

% Release point, relative to the vent, in metres
% x0: Positive when E, negative when W
% y0: Positive when N, negative when S
% z0: Above sea level

%% ________________________________________________________________________
%% Input parameters
%% ________________________________________________________________________

%% Load atmospheric data
atm             = load(P.path.nc); atm = atm.atm;
%% DEM      -> Check whether symetry should be on WIND instead
dem             = load(P.path.dem); dem = dem.dem;                                      % Load DEM
dem.X(dem.X<0)  = 360+(dem.X(dem.X<0));                                     % Symetry latitude

% Set domain extent based on the union of DEM and wind data
lat_min         = min([min(atm.lat), min(dem.Y(:,1))]);
lat_max         = max([max(atm.lat), max(dem.Y(:,1))]);
lon_min         = min([min(atm.lon), min(dem.X(1,:))]);
lon_max         = max([max(atm.lon), max(dem.X(1,:))]);

%% Indices
% Initial particle release position/time
part.x(1)       = P.rel.x;    
part.y(1)       = P.rel.y;
part.t(1)       = P.rel.t;
part.z(1)       = P.vent.alt+P.rel.z;
part.dis(1)     = sqrt(P.rel.x^2 + P.rel.y^2);  % Distance with added ofsets
part.time(1)    = 0;
part.bear(1)    = 0; % Particle bearing

% Calculate lat and lon of particle release
[part.lat(1),part.lon(1)] = dist2ll(P.vent.lat, P.vent.lon, part.x, part.y);

% Find particle indices relative to atmospheric conditions
[~, part.xI(1)] = min(abs(atm.lon - part.lon(1)));
[~, part.yI(1)] = min(abs(atm.lat - part.lat(1)));
[~, part.tI(1)] = min(abs(atm.time - P.date + (P.rel.t/3600/24) ));         % Add time offset to eruption date
[~, part.zI(1)] = min(abs(atm.alt(part.yI(1), part.xI(1), :, part.tI(1)) - part.z(1)));

% Find particle indices relative to DEM
[~, part.xD(1)] = min(abs(dem.X(1,:) - part.lon(1)));
[~, part.yD(1)] = min(abs(dem.Y(:,1) - part.lat(1)));

% Initial atmosperic condidtions
part.uf(1)      = atm.u(part.yI(1), part.xI(1), part.zI(1), part.tI(1));
part.vf(1)      = atm.v(part.yI(1), part.xI(1), part.zI(1), part.tI(1));

%% Particle
part.Fs         = P.part.flat * P.part.elon^1.3;                            % Shape descriptor - Stoke's regime
part.Fn         = P.part.flat^2 * P.part.elon;                              % Shape descriptor - Newton regime
part.Ks         = .5 * (part.Fs^(1/3) + part.Fs^(-1/3));                    % Stoke's drag correction
part.Kn         = 10^(.45 * (-log10(part.Fn))^.99);                         % Newton's drag correction     

% Initial parameters

% If no initial x or y velocity is given, assume the particle is carried by
% the wind
if P.rel.vx == 0
    part.u(1)   = part.uf(1);
else
    part.u(1)   = P.rel.vx;
end
if P.rel.vy == 0
    part.v(1)   = part.vf(1);
else
    part.v(1)   = P.rel.vy;
end

part.w(1)       = P.rel.vz;
part.Re_w(1)    = 0;
part.Cd_w(1)    = 0;
part.tau(1)     = P.part.dens * P.part.diam^2 / (18 * atm.muair(part.yI(1), part.xI(1), part.zI(1), part.tI(1)));    % Particle relaxation time (s)

%% Interpolation
% Initialize counters
int_count1      = 0;                                                        % counter for skip check for density and viscosity
int_count2      = 0;                                                        % counter forskip check for wind velocity


%% ________________________________________________________________________
%% Calculations
%% ________________________________________________________________________

%% Particle trajectory
g               = 9.806;                                                    % Gravity m/s2    
i               = 1;
test_run        = 0;
while test_run == 0
    i = i+1;
    % Interpolation of physical properties of atmosphere     
    if strcmp(P.adv.interp, 'none')
    % Simple indexing method
        denf    = atm.rhoair(part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));   % Fluid density
        visf    = atm.muair( part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));   % fluid viscosity
    
    elseif strcmp(P.adv.interp, 'complete')
    % Complete interpolation method (Slow)
        alt_sub = squeeze(atm.alt(part.yI(i-1), part.xI(i-1), :, part.tI(i-1)));
        denf    = interpn(atm.lat, atm.lon, alt_sub, atm.time, atm.rhoair, part.lat(i-1), part.lon(i-1), part.z(i-1), P.date+part.t(i-1)/3600/24, P.adv.method);
        visf    = interpn(atm.lat, atm.lon, alt_sub, atm.time, atm.muair,  part.lat(i-1), part.lon(i-1), part.z(i-1), P.date+part.t(i-1)/3600/24, P.adv.method);

    elseif strcmp(P.adv.interp, 'subset')
    % Interpolation by subsetting method
        % check if the subsetting range is within in the downloaded data
        if part.xI(i-1)-P.adv.range < 1 || part.xI(i-1)+P.adv.range > size(atm.lon,1) || part.yI(i-1)-P.adv.range < 1 || part.yI(i-1)+P.adv.range > size(atm.lat,1) || part.zI(i-1)-P.adv.range < 1 || part.zI(i-1)+P.adv.range > size(atm.alt,3)
             denf       = atm.rhoair(part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));   % Fluid density
             visf       = atm.muair( part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));   % fluid viscosity
        elseif int_count1==0 || P.adv.skip==0
             int_count1 = int_count1+1;
             alt_vec    = squeeze(atm.alt(part.yI(i-1), part.xI(i-1), part.zI(i-1)-P.adv.range:part.zI(i-1)+P.adv.range, part.tI(i-1)));
             time_sub   = atm.time(part.tI(i-1) - P.adv.range:part.tI(i-1) + P.adv.range);
             lon_sub    = atm.lon(part.xI(i-1)  - P.adv.range:part.xI(i-1) + P.adv.range);
             lat_sub    = atm.lat(part.yI(i-1)  - P.adv.range:part.yI(i-1) + P.adv.range);
             den_sub    = atm.rhoair(part.yI(i-1)-P.adv.range:part.yI(i-1) + P.adv.range, part.xI(i-1)-P.adv.range:part.xI(i-1)+P.adv.range, part.zI(i-1)-P.adv.range:part.zI(i-1)+P.adv.range, part.tI(i-1)-P.adv.range:part.tI(i-1)+P.adv.range);
             vis_sub    = atm.muair(part.yI(i-1) -P.adv.range:part.yI(i-1) + P.adv.range, part.xI(i-1)-P.adv.range:part.xI(i-1)+P.adv.range, part.zI(i-1)-P.adv.range:part.zI(i-1)+P.adv.range, part.tI(i-1)-P.adv.range:part.tI(i-1)+P.adv.range);
             denf       = interpn(lat_sub, lon_sub, alt_vec, time_sub, den_sub, part.lat(i-1), part.lon(i-1), part.z(i-1), P.date+part.t(i-1)/3600/24, P.adv.method);   % Fluid density
             visf       = interpn(lat_sub, lon_sub, alt_vec, time_sub, vis_sub, part.lat(i-1), part.lon(i-1), part.z(i-1), P.date+part.t(i-1)/3600/24, P.adv.method);   % Fluid viscosity
         else
             int_count1  = int_count1+1;
             if int_count1==P.adv.skip
                  int_count1 = 0; 
             end
         end
    end


    % Update particle state
    velr            = sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2+part.w(i-1)^2);	% Particle relative velocity
    
    if sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2)==0
        theta       = 90;
    else
        theta       = atand((part.vf(i-1)-part.v(i-1))/sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2));	% direction of particle in the horizontal plane (x-y)
    end
    if sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2)==0
        Beta        = 90;
    else
        Beta        = atand(part.w(i-1)/sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2)); % direction of particle in the vertical plane (x-z or y-z)
    end
    part.Re(i)      = denf * velr * P.part.diam / visf;                     % Reynolds
    part.Re_S(i)    = part.Re(i) * part.Kn / part.Ks;                       % Ganser Re
    
    if part.Re_S(i) > 3e10^5
        part.Cd(i)  = 0.2;
    else
        part.Cd(i)  = part.Kn * (24 * (1 + .125 * part.Re_S(i)^(2/3)) / part.Re_S(i) + .46 / (1+5330/part.Re_S(i)));  % Drag coef
    end
    
    Fd              = part.Cd(i) * part.Re(i) / (24*part.tau(i-1));         % Total Drag force
    Fd_u            = abs(Fd * cosd(Beta) * cosd(theta));                   % Drag force in x direction
    G_u             = 0.;                                                   % Other forces in x direction 
    Fd_v            = abs(Fd * cosd(Beta) * sind(theta));                   % Drag force in y direction
    G_v             = 0.;                                                   % Other forces in y direction
    Fd_w            = abs(Fd * sind(Beta));                                 % Drag force in z direction
    G_w             = (1-denf/P.part.dens) * -g;                              % Gravity force
    
    if strcmp(P.adv.solution, 'analytical')
    % Analytical solution
        part.u(i)   = part.uf(i-1) + exp(-Fd * P.adv.dt) *(part.u(i-1)-part.uf(i-1)) - G_u * (1/Fd) * (exp(-P.adv.dt*Fd)-1);
        part.x(i)   = part.x(i-1) + (G_u * (1/Fd)+part.uf(i-1) )* P.adv.dt + (1/Fd)* (1-exp(-P.adv.dt*Fd)) * (part.u(i) -part.uf(i-1)- G_u/Fd);

        part.v(i)   = part.vf(i-1) + exp(-Fd * P.adv.dt) *(part.v(i-1)-part.vf(i-1)) - G_v * (1/Fd) * (exp(-P.adv.dt*Fd)-1);
        part.y(i)   = part.y(i-1) + (G_v * (1/Fd)+part.vf(i-1)) * P.adv.dt + (1/Fd)* (1-exp(-P.adv.dt*Fd)) * (part.v(i)-part.vf(i-1) - G_v/Fd);

        part.w(i)   = exp(-Fd * P.adv.dt) *(part.w(i-1)) - G_w * (1/Fd) * (exp(-P.adv.dt*Fd)-1);
        part.z(i)   = part.z(i-1) + G_w * (1/Fd) * P.adv.dt + (1/Fd)* (1-exp(-P.adv.dt*Fd)) * (part.w(i) - G_w/Fd);

    elseif strcmp(P.adv.solution, 'euler')
    % Euler semi-implicit 
        part.u(i)   = ((G_u + Fd_u * part.uf(i-1))* P.adv.dt + part.u(i-1)) / (1+Fd_u * P.adv.dt);
        part.x(i)   = part.x(i-1) + .5 * P.adv.dt * (part.u(i) + part.u(i-1));

        part.v(i)   = ((G_v + Fd_v * part.vf(i-1)) * P.adv.dt + part.v(i-1)) / (1+Fd_v * P.adv.dt);
        part.y(i)   = part.y(i-1) + .5 * P.adv.dt * (part.v(i) + part.v(i-1));

        part.w(i)   = (G_w * P.adv.dt + part.w(i-1)) / (1+Fd_w * P.adv.dt);
        part.z(i)   = part.z(i-1) + .5 * P.adv.dt * (part.w(i) + part.w(i-1));
    
    
    elseif strcmp(P.adv.solution, 'rangekutta')
    % Runge Kutta    
        
    end
    
    
    % Update particle time
    part.t(i) = part.t(i-1) + P.adv.dt;
    
    % Update geographic coordinates
    [part.lat(i), part.lon(i)] = dist2ll(part.lat(i-1), part.lon(i-1), part.x(i)-part.x(i-1), part.y(i)-part.y(i-1));
    
    % Update bearing
    part.bear(i)    = -(atan2d(part.lat(i)-P.vent.lat, part.lon(i)-P.vent.lon)-90);
    if part.bear(i) < 0
        part.bear(i)= 360+part.bear(i);
    end
    
    % Update indices
    [~, part.xI(i)] = min(abs(atm.lon-part.lon(i)));
    [~, part.yI(i)] = min(abs(atm.lat-part.lat(i)));
    [~, part.tI(i)] = min(abs(atm.time-(P.date + part.t(i)/3600/24)));
    [~, part.zI(i)] = min(abs(atm.alt(part.yI(i), part.xI(i), :, part.tI(i)) - part.z(i)));
    part.dis(i)     = sqrt((part.x(i)-part.x(1))^2+(part.y(i)-part.y(1))^2+(part.z(i)-part.z(1))^2);
    part.time(i)    = part.time(i-1)+P.adv.dt;
    
    [~, part.xD(i)] = min(abs(dem.X(1,:)-part.lon(i)));
    [~, part.yD(i)] = min(abs(dem.Y(:,1)-part.lat(i)));
    
    
    % Update tau
    part.tau(i) = P.part.dens * P.part.diam^2 / (18 * atm.muair(part.yI(i), part.xI(i), part.zI(i), part.tI(i)));    % Particle relaxation time (s)
    % Note: I took out the update of tau from the following conditional
    % test so the final tau vector has the same length as the other
    % variables, and so it can be plotted.
    
    
    %% Test conditions
    % Test altitude
    if part.z(i) <= dem.Z(part.yD(i), part.xD(i))
        part.out_msg = sprintf('Particle landed on the domain');
        test_run     = 1;
        % Test domain
    elseif part.lon(i) <= lon_min || part.lon(i) >= lon_max || part.lat(i) <= lat_min || part.lat(i) >= lat_max
        %part.out_msg = sprintf('Particle reached domain border at an altitude of %4.0f m a.s.l.', dem.Z(part.yD(i), part.xD(i)));
        part.out_msg = sprintf('Particle reached the domain border');
        test_run     = 1;
        % Test time
    elseif (P.date + i*P.adv.dt/3600/24) > max(atm.time)+0.25
        %display(sprintf('Particle residence time longer than atmospheric data (currently %2.0f km away from the vent at an altitude of %4.0f m a.s.l.)', sqrt(part.x(i)^2 + part.z(i)^2)/1000, dem.Z(part.yD(i), part.xD(i))));
        part.out_msg = sprintf('Particle residence time longer than atmospheric data');
        test_run     = 1;
    else
                               
        % Get u and v wind coordinates at new location
        if strcmp(P.adv.interp, 'none')
            % simple indexing method
            part.uf(i)  = atm.u(part.yI(i), part.xI(i), part.zI(i), part.tI(i));
            part.vf(i)  = atm.v(part.yI(i), part.xI(i), part.zI(i), part.tI(i));
            
        elseif strcmp(P.adv.interp, 'complete')
            %complete interpolation method (Slow)
            alt_sub      = squeeze(atm.alt(part.yI(i), part.xI(i), :, part.tI(i)));
            part.uf(i)   = interpn(atm.lat, atm.lon, alt_sub, atm.time, atm.u, part.lat(i), part.lon(i), part.z(i), P.date+part.t(i)/3600/24, P.adv.method);
            part.vf(i)   = interpn(atm.lat, atm.lon, alt_sub, atm.time, atm.v, part.lat(i), part.lon(i), part.z(i), P.date+part.t(i)/3600/24, P.adv.method);
            
        elseif strcmp(P.adv.interp, 'subset')
            %interpolation by subsetting method
            % check if the subsetting range is within in the downloaded data
            if part.xI(i)-P.adv.range < 1 || part.xI(i)+P.adv.range > size(atm.lon,1) || part.yI(i)-P.adv.range < 1 || part.yI(i)+P.adv.range > size(atm.lat,1) || part.zI(i)-P.adv.range < 1 || part.zI(i)+P.adv.range > size(atm.alt,3)
                part.uf(i)  = atm.u(part.yI(i), part.xI(i), part.zI(i), part.tI(i));
                part.vf(i)  = atm.v(part.yI(i), part.xI(i), part.zI(i), part.tI(i));
            elseif int_count2==0 || P.adv.skip==0
                int_count2 = int_count2+1;
                alt_vec    = squeeze(atm.alt(part.yI(i), part.xI(i), part.zI(i)-P.adv.range:part.zI(i)+P.adv.range, part.tI(i)));
                time_sub   = atm.time(part.tI(i) - P.adv.range:part.tI(i) + P.adv.range);
                lon_sub    = atm.lon(part.xI(i)  - P.adv.range:part.xI(i) + P.adv.range);
                lat_sub    = atm.lat(part.yI(i) - P.adv.range:part.yI(i) + P.adv.range);
                u_sub      = atm.u(part.yI(i) - P.adv.range:part.yI(i) + P.adv.range, part.xI(i) - P.adv.range:part.xI(i) + P.adv.range, part.zI(i) - P.adv.range:part.zI(i) + P.adv.range, part.tI(i) - P.adv.range:part.tI(i) + P.adv.range);
                v_sub      = atm.v(part.yI(i) - P.adv.range:part.yI(i) + P.adv.range, part.xI(i) - P.adv.range:part.xI(i) + P.adv.range, part.zI(i) - P.adv.range:part.zI(i) + P.adv.range, part.tI(i) - P.adv.range:part.tI(i) + P.adv.range);
                part.uf(i) = interpn(lat_sub, lon_sub, alt_vec, time_sub, u_sub, part.lat(i), part.lon(i), part.z(i), P.date+part.t(i)/3600/24, P.adv.method);
                part.vf(i) = interpn(lat_sub, lon_sub, alt_vec, time_sub, v_sub, part.lat(i), part.lon(i), part.z(i), P.date+part.t(i)/3600/24, P.adv.method);
            else
                int_count2  = int_count2+1;
                if int_count2==P.adv.skip
                    int_count2 = 0;
                end
                part.uf(i)  = part.uf(i-1);
                part.vf(i)  = part.vf(i-1);
            end
        end
    end    
end

% 
% %% Plotting
% figure;
% if strcmp(vent.geo, 'West')
%     plot(-360+[lon_min, lon_max], [lat_min, lat_max], '.k');                    % Plot frame
%     plot_google_map('Maptype', 'terrain'); hold on;                             % Plot google background
%     surf(-360+dem.X, dem.Y, dem.Z./1000); shading flat; alpha 0.5               % Plot topography
%     colormap(landcolor);
%     c = colorbar;
%     ylabel(c, 'Topography altitude (km asl)');
%     freezeColors;
%     [C,h] = contour(-360+dem.X, dem.Y, ll2dist(P.vent.lat, P.vent.lon, dem.Y, dem.X));
%     set(h, 'LineColor', 'k', 'Linewidth', 1.5);
%     clabel(C,h);
%     plot3(-360+part.lon, part.lat, part.z/1000, '.r');                          % Plot trajectory
%     plot3([-360+P.vent.lon,-360+P.vent.lon], [P.vent.lat,P.vent.lat], [P.vent.ht/1000, erup.ht], '--k')
%     plot3(-360+P.vent.lon, P.vent.lat, P.vent.ht/1000, '^k', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
% 
% else
%     
%     plot([min(dem.X(1,:)), max(dem.X(1,:))], [min(dem.Y(:,1)), max(dem.Y(:,1))], '.k');                    % Plot frame
%     plot_google_map('Maptype', 'terrain'); hold on;                             % Plot google background
%     surf(dem.X, dem.Y, dem.Z./1000); shading flat; alpha 0.5               % Plot topography
%     colormap(landcolor);
%     c = colorbar;
%     ylabel(c, 'Topography altitude (km asl)');
%     freezeColors;
%     [C,h] = contour(dem.X, dem.Y, ll2dist(P.vent.lat, P.vent.lon, dem.Y, dem.X));
%     set(h, 'LineColor', 'k', 'Linewidth', 1.5);
%     clabel(C,h);
%     plot3(part.lon, part.lat, part.z/1000, '.r');                          % Plot trajectory
%     plot3([P.vent.lon,P.vent.lon], [P.vent.lat,P.vent.lat], [P.vent.ht/1000, erup.ht], '--k')
%     plot3(P.vent.lon, P.vent.lat, P.vent.ht/1000, '^k', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
% 
% end
% title(sprintf('%4.4f mm, %4.f kg/m^3, fl = %2.1f, el = %2.1f', part.diam*1000, part.dens, part.fl, part.el));






