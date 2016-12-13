function part = LagTrack(x0, y0, z0)

% UPDATES
% 2016/02/02 (Seb): Added bearing

% TO DO:
% Variable to choose interpolation schemen

% Release point, relative to the vent, in metres
% x0: Positive when E, negative when W
% y0: Positive when N, negative when S
% z0: Above sea level

%% ________________________________________________________________________
%% Input parameters
%% ________________________________________________________________________

%% volacano 
Volcano = 'MstHelens_May_1980';
% Volcano = 'Etna_Nov_2013';

%% Input files
if strcmp(Volcano, 'MstHelens_May_1980')
    %% atmospheric file
    ncfile          = 'input/wind/MSH_May_1980/MSH_May_1980.nc';                         % NetCDF file
    
    %% DEM file
    demfile         = 'input/dem/MSH/MSH.mat';                                    % DEM file
    
    %% Vent specifications
    vent.lat        = 46.1912;                                                 % Vent latitude
    vent.lon        = 237.8056;                                                % Vent longitude (deg East)
    vent.ht         = 2549;                                                    % Vent altitude
    vent.geo        = 'West';                                                  % Vent hemisphere

    %% Eruption
    erup.ht         = 16;                                                       % km asl
    erup.time       = datenum([1980,5,18,18,0,0]);                              % Eruption time

elseif strcmp(Volcano, 'Etna_Nov_2013')
    %% atmospheric file
    ncfile          = 'NetCDF_Etna_Nov_2013.nc';                                           % NetCDF file
    
    %% DEM file
    demfile         = 'dem_etna.mat';                                               % DEM file
    
    %% Vent specifications
    vent.lat        = 37.74;                                                 % Vent latitude
    vent.lon        = 15;                                                    % Vent longitude (deg East)
    vent.ht         = 3217;                                                  % Vent altitude
    vent.geo        = 'East';                                                % Vent hemisphere

    %% Eruption
    erup.ht         = 4.217;                                                       % km asl
    erup.time       = datenum([2013,11,23,10,55,0]);                              % Eruption time
end


%% Particle
part.diam       = 125 * 10^-6 ;                                              % Diameter (m)
part.dens       = 1000;                                                     % Density (kg/m3)
part.fl         = 0.7;                                                      % Flatness ( = shortest_length / intermediate_length)
part.el         = 0.7;                                                      % Elongation ( = intermediate_length / longest_length)

%% Interpolation
interpolation   = 'subset';                                                 % Interpolation of atmospheric data, can be either 'none', 'subset' or 'complete'

% Used if interpolation = 'subset'
subsetI         = 1;                                                        % Range of neighbor indexes used for susetting (higher: slower interpolation)
int_skip        = 0;                                                        % Number of time steps to skip interpolation (higher: faster interpolation)
int_method      = 'linear';                                                 % Interpolation method - 'linear', 'neareast', 'pchip', 'cubic' or 'spline'
                                                                            % See help interpn for more details
%% Trajectory solution
solution        = 'Euler';                                                  % Method to solve the trajectory, can be either 'analytical', 'Euler' or 'RungeKutta'
                                                                            
%% Time step
dt              = .1;                                                       % Time step (s)

%% ________________________________________________________________________
%% Pre-processing
%% ________________________________________________________________________

home

%% Atmospheric properties
% Load and pre-process
lat             = double(ncread(ncfile, 'latitude'));                       % Latitude (degrees)
lon             = double(ncread(ncfile, 'longitude'));                      % Longitude (degrees)
level           = ncread(ncfile, 'level');                                  % Pressure level (mb)
time            = datenum(datenum([1900,1,1,0,0,0])+double(ncread(ncfile, 'time'))./24);    % Time (year month day hour min sec)
temp            = permute(ncread(ncfile, 't'), [2,1,3,4]);                  % Temperature (deg K)
alt             = permute(ncread(ncfile, 'z')./9.80665, [2,1,3,4]);         % Altitude (m asl)
humid           = permute(ncread(ncfile, 'r'), [2,1,3,4]);                  % Relative humidity (%)
u               = permute(ncread(ncfile, 'u'), [2,1,3,4]);                  % U wind (ms-1)
v               = permute(ncread(ncfile, 'v'), [2,1,3,4]);                  % V wind (ms-1)
% Air density
pres    = double(repmat(reshape(level,1,1,length(level)), length(lat), length(lon), 1, length(time)));            
tempc           = temp-273.15;                                              % Temperature degC
Psat            = 6.1078 .* 10.^(7.5*tempc./(tempc+237.3)).*100;            % Saturation vapor pressure (Pa)
Pv              = humid./100 .* Psat;                                       % Vapor pressure of water (Pa)
Pd              = pres.*100 - Pv;                                           % Partial pressure dry air
Rd              = 287.058;                                                  % Specific gas constant for dry air (J/kg*K)
Rv              = 461.495;                                                  % Specific gas constant for water vapor (J/kg*K)
rhoair          = Pd./(Rd.*temp) + Pv./(Rv.*temp);                          % Air denstiy (kg m-3);
% Air viscosity
tempR           = (tempc+273.15) .* (9/5);                                  % Temperature (deg Rankine)
mu0             = 0.01827;                                                  % Reference viscosity (centipoise)
t0R             = 524.07;                                                   % Reference temperature (deg Rankine)
C               = 120;                                                      % Sutherland's constant
a               = .555*t0R+C;                                                      
b               = .555.*tempR+C;
muair           = mu0.*(a./b).*(tempR./t0R).^(3/2)./10^3;

%% DEM
load(demfile);                                                              % Load DEM
dem.X(dem.X<0)  = 360+(dem.X(dem.X<0));                                     % Symetry latitude
% Set boundaries based on DEM
lat_min         = min([min(lat), min(dem.Y(:,1))]);
lat_max         = max([max(lat), max(dem.Y(:,1))]);
lon_min         = min([min(lon), min(dem.X(1,:))]);
lon_max         = max([max(lon), max(dem.X(1,:))]);

%% Indices
% Initial particle release position/time
part.x(1)       = x0;    
part.y(1)       = y0;
part.t(1)       = 0;
part.z(1)       = z0;
part.dis(1)     = 0;
part.time(1)    = 0;
part.bear(1)    = 0; % Particle bearing

% Calculate lat and lon of particle release
[part.lat(1),part.lon(1)] = dist2ll(vent.lat, vent.lon, part.x, part.y);

% Time index
[~,erup.timeI]  = min(abs(time-erup.time));                                 % Index of eruption time for atmospheric data

% Find particle indices relative to atmospheric conditions
[~, part.xI(1)] = min(abs(lon-part.lon(1)));
[~, part.yI(1)] = min(abs(lat-part.lat(1)));
[~, part.tI(1)] = min(abs(time-erup.time));
[~, part.zI(1)] = min(abs(alt(part.yI(1), part.xI(1), :, part.tI(1)) - part.z(1)));

% Find particle indices relative to DEM
[~, part.xD(1)] = min(abs(dem.X(1,:)-part.lon(1)));
[~, part.yD(1)] = min(abs(dem.Y(:,1)-part.lat(1)));

% Initial atmosperic condidtions
part.uf(1)       = u(part.yI(1), part.xI(1), part.zI(1), part.tI(1));
part.vf(1)       = v(part.yI(1), part.xI(1), part.zI(1), part.tI(1));

%% Particle
part.Fs         = part.fl * part.el^1.3;                                    % Shape descriptor - Stoke's regime
part.Fn         = part.fl^2 * part.el;                                      % Shape descriptor - Newton regime
part.Ks         = .5*(part.Fs^(1/3) + part.Fs^(-1/3));                      % Stoke's drag correction
part.Kn         = 10^(.45 * (-log10(part.Fn))^.99);                         % Newton's drag correction     
% Initial parameters
part.u(1)       = part.uf(1);
part.v(1)       = part.vf(1);
part.w(1)       = 1e-4;
part.Re_w(1)    = 0;
part.Cd_w(1)    = 0;
part.tau(1)     = part.dens * part.diam^2 / (18 * muair(part.yI(1), part.xI(1), part.zI(1), part.tI(1)));    % Particle relaxation time (s)

%% Interpolation
% Initialize counters
int_count1      = 0;                                                        % counter for int_skip check for density and viscosity
int_count2      = 0;                                                        % counter for int_skip check for wind velocity

%% Printing paramters
Iprint_Skip     = 100;                                                       % number of time steps to skip printing results on the screen
Iprint_count    = 0;                                                        % counter for Iprint_Skip






%% ________________________________________________________________________
%% Calculations
%% ________________________________________________________________________

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



%% Particle trajectory
g               = 9.806;                                                    % Gravity m/s2
%dt              = .5;                                                       % Time step
i               = 1;
test_run        = 0;
count_table     = 1;
while test_run == 0
    i = i+1;
    
%     if i < 100
%         dt = tau/100;
%     elseif i>=100 || 1<1000
%         dt = tau/10;
%     else
%         dt = 10*tau;
%     end
    
% Interpolation of physical properties of atmosphere 
    
    if strcmp(interpolation, 'none')
    % Simple indexing method
        denf = rhoair(part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));  % Fluid density
        visf = muair(part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));   % fluid viscosity
    
    elseif strcmp(interpolation, 'complete')
    % Complete interpolation method (Slow)
     alt_sub = squeeze(alt(part.yI(i-1), part.xI(i-1), :, part.tI(i-1)));
     denf    = interpn(lat,lon,alt_sub,time,rhoair,part.lat(i-1),part.lon(i-1),part.z(i-1),erup.time+part.t(i-1)/3600/24);
     visf    = interpn(lat,lon,alt_sub,time,muair,part.lat(i-1),part.lon(i-1),part.z(i-1),erup.time+part.t(i-1)/3600/24);

    elseif strcmp(interpolation, 'subset')
    % Interpolation by subsetting method
        % check if the subsetting range is within in the downloaded data
        if part.xI(i-1)-subsetI < 1 || part.xI(i-1)+subsetI > size(lon,1) || part.yI(i-1)-subsetI < 1 || part.yI(i-1)+subsetI > size(lat,1) || part.zI(i-1)-subsetI < 1 || part.zI(i-1)+subsetI > size(alt,3)
             denf       = rhoair(part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));  % Fluid density
             visf       = muair(part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));   % fluid viscosity
        elseif int_count1==0 || int_skip==0
             int_count1 = int_count1+1;
             alt_vec    = squeeze(alt(part.yI(i-1), part.xI(i-1), part.zI(i-1)-subsetI:part.zI(i-1)+subsetI, part.tI(i-1)));
             time_sub   = time (part.tI(i-1)-subsetI:part.tI(i-1)+subsetI);
             lon_sub    = lon(part.xI(i-1)-subsetI:part.xI(i-1)+subsetI);
             lat_sub    = lat(part.yI(i-1)-subsetI:part.yI(i-1)+subsetI);
             den_sub    = rhoair(part.yI(i-1)-subsetI:part.yI(i-1)+subsetI,part.xI(i-1)-subsetI:part.xI(i-1)+subsetI,part.zI(i-1)-subsetI:part.zI(i-1)+subsetI,part.tI(i-1)-subsetI:part.tI(i-1)+subsetI);
             vis_sub    = muair(part.yI(i-1)-subsetI:part.yI(i-1)+subsetI,part.xI(i-1)-subsetI:part.xI(i-1)+subsetI,part.zI(i-1)-subsetI:part.zI(i-1)+subsetI,part.tI(i-1)-subsetI:part.tI(i-1)+subsetI);
             denf       = interpn(lat_sub,lon_sub,alt_vec,time_sub,den_sub,part.lat(i-1),part.lon(i-1),part.z(i-1),erup.time+part.t(i-1)/3600/24);   % Fluid density
             visf       = interpn(lat_sub,lon_sub,alt_vec,time_sub,vis_sub,part.lat(i-1),part.lon(i-1),part.z(i-1),erup.time+part.t(i-1)/3600/24);   % Fluid viscosity
         else
             int_count1  = int_count1+1;
             if int_count1==int_skip
                  int_count1 = 0; 
             end
         end
    end
%    visf=18.32e-6;
%     denf = 1.19; % Fluid density
%     visf = 1.828e-5;% fluid viscosity
%     part.tau(i)     = part.dens * part.diam^2 / (18 * visf);    % Particle relaxation time (s)

    % Update particle state
    velr            = sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2+part.w(i-1)^2);     % Particle relative velocity
    w_norm          = part.w(i-1)/velr;
    if sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2)==0
        theta   = 90;
    else
        theta   = atand((part.vf(i-1)-part.v(i-1))/sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2));                       % direction of particle in the horizontal plane (x-y)
    end
    if sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2)==0
        Beta    = 90;
    else
        Beta    = atand(part.w(i-1)/sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2)); % direction of particle in the vertical plane (x-z or y-z)
    end
    part.Re(i)      = denf * velr * part.diam / visf;                       % Reynolds
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
    G_w             = (1-denf/part.dens) * -g;                              % Gravity force
    
    if strcmp(solution, 'analytical')
    % Analytical solution
        part.u(i) = part.uf(i-1) + exp(-Fd * dt) *(part.u(i-1)-part.uf(i-1)) - G_u * (1/Fd) * (exp(-dt*Fd)-1);
        part.x(i) = part.x(i-1) + (G_u * (1/Fd)+part.uf(i-1) )* dt + (1/Fd)* (1-exp(-dt*Fd)) * (part.u(i) -part.uf(i-1)- G_u/Fd);

        part.v(i) = part.vf(i-1) + exp(-Fd * dt) *(part.v(i-1)-part.vf(i-1)) - G_v * (1/Fd) * (exp(-dt*Fd)-1);
        part.y(i) = part.y(i-1) + (G_v * (1/Fd)+part.vf(i-1)) * dt + (1/Fd)* (1-exp(-dt*Fd)) * (part.v(i)-part.vf(i-1) - G_v/Fd);

        part.w(i) = exp(-Fd * dt) *(part.w(i-1)) - G_w * (1/Fd) * (exp(-dt*Fd)-1);
        part.z(i) = part.z(i-1) + G_w * (1/Fd) * dt + (1/Fd)* (1-exp(-dt*Fd)) * (part.w(i) - G_w/Fd);

    elseif strcmp(solution, 'Euler')
        % Euler semi-implicit 
        part.u(i) = ((G_u + Fd_u * part.uf(i-1))* dt + part.u(i-1)) / (1+Fd_u * dt);
        part.x(i) = part.x(i-1) + .5 * dt * (part.u(i) + part.u(i-1));

        part.v(i) = ((G_v + Fd_v * part.vf(i-1)) * dt + part.v(i-1)) / (1+Fd_v * dt);
        part.y(i) = part.y(i-1) + .5 * dt * (part.v(i) + part.v(i-1));

        part.w(i) = (G_w * dt + part.w(i-1)) / (1+Fd_w * dt);
        %test=part.w(i)
        part.z(i) = part.z(i-1) + .5 * dt * (part.w(i) + part.w(i-1));
    
    
    elseif strcmp(solution, 'RangeKutta')
    % Runge Kutta    
        
    end
    
    
    % Update particle time
    part.t(i) = part.t(i-1) + dt;
    
    % Update geographic coordinates
    [part.lat(i), part.lon(i)] = dist2ll(part.lat(i-1), part.lon(i-1), part.x(i)-part.x(i-1), part.y(i)-part.y(i-1));
    
    % Update bearing
    part.bear(i)    = -(atan2d(part.lat(i)-vent.lat, part.lon(i)-vent.lon)-90);
    if part.bear(i) < 0
        part.bear(i) = 360+part.bear(i);
    end
    
    % Update indices
    [~, part.xI(i)] = min(abs(lon-part.lon(i)));
    [~, part.yI(i)] = min(abs(lat-part.lat(i)));
    [~, part.tI(i)] = min(abs(time-(erup.time + part.t(i)/3600/24)));
    [~, part.zI(i)] = min(abs(alt(part.yI(i), part.xI(i), :, part.tI(i)) - part.z(i)));
    part.dis(i)     = sqrt((part.x(i)-part.x(1))^2+(part.y(i)-part.y(1))^2+(part.z(i)-part.z(1))^2);
    part.time(i)    = part.time(i-1)+dt;
    
    [~, part.xD(i)] = min(abs(dem.X(1,:)-part.lon(i)));
    [~, part.yD(i)] = min(abs(dem.Y(:,1)-part.lat(i)));
    
    % Test altitude
    if part.z(i) <= dem.Z(part.yD(i), part.xD(i))
        
        table([{'Landing altitude (m asl)'}; {'Landing distance (km)'}; {'Total flight time (hours)'}; {'Terminal velocity at impact (m/s)'}; {'Maximum terminal velocity (m/s)'}],...
            [dem.Z(part.yD(i), part.xD(i)); sqrt((part.x(end)-x0)^2 + (part.y(end)-y0)^2)/1000; (i*dt)/3600; part.w(end);min(part.w(:))],...
            'VariableNames', {'Parameter' 'Value'});          
        
        test_run        = 1;
        Iprint_count    = 0;      
    % Test domain    
    elseif part.lon(i) <= lon_min || part.lon(i) >= lon_max || part.lat(i) <= lat_min || part.lat(i) >= lat_max        
        display(sprintf('Particle reached domain border at an altitude of %4.0f m a.s.l.', dem.Z(part.yD(i), part.xD(i))));
        test_run        = 1;
        Iprint_count    = 0;
    % Test time
    elseif (erup.time + i*dt/3600/24) > max(time)+0.25        
        display(sprintf('Particle residence time longer than atmospheric data (currently %2.0f km away from the vent at an altitude of %4.0f m a.s.l.)', sqrt(part.x(i)^2 + part.z(i)^2)/1000, dem.Z(part.yD(i), part.xD(i))));
        test_run        = 1;
        Iprint_count    = 0;
    else
        
    % Update tau
    part.tau(i) = part.dens * part.diam^2 / (18 * muair(part.yI(i), part.xI(i), part.zI(i), part.tI(i)));    % Particle relaxation time (s)

    
    % Get u and v wind coordinates at new location
    if strcmp(interpolation, 'none')    
    % simple indexing method
        part.uf(i) = u(part.yI(i), part.xI(i), part.zI(i), part.tI(i));
        part.vf(i) = v(part.yI(i), part.xI(i), part.zI(i), part.tI(i));

    elseif strcmp(interpolation, 'complete')
    %complete interpolation method (Slow)    
        alt_sub     = squeeze(alt(part.yI(i), part.xI(i), :, part.tI(i)));
        part.uf(i)   = interpn(lat,lon,alt_sub,time,u,part.lat(i),part.lon(i),part.z(i),erup.time+part.t(i)/3600/24);
        part.vf(i)   = interpn(lat,lon,alt_sub,time,v,part.lat(i),part.lon(i),part.z(i),erup.time+part.t(i)/3600/24);

    elseif strcmp(interpolation, 'subset')
    %interpolation by subsetting method
        % check if the subsetting range is within in the downloaded data
        if part.xI(i)-subsetI < 1 || part.xI(i)+subsetI > size(lon,1) || part.yI(i)-subsetI < 1 || part.yI(i)+subsetI > size(lat,1) || part.zI(i)-subsetI < 1 || part.zI(i)+subsetI > size(alt,3)
             part.uf(i)  = u(part.yI(i), part.xI(i), part.zI(i), part.tI(i));
             part.vf(i)  = v(part.yI(i), part.xI(i), part.zI(i), part.tI(i));
        elseif int_count2==0 || int_skip==0
             int_count2 = int_count2+1;
             alt_vec    = squeeze(alt(part.yI(i), part.xI(i), part.zI(i)-subsetI:part.zI(i)+subsetI, part.tI(i)));
             time_sub   = time (part.tI(i)-subsetI:part.tI(i)+subsetI);
             lon_sub    = lon(part.xI(i)-subsetI:part.xI(i)+subsetI);
             lat_sub    = lat(part.yI(i)-subsetI:part.yI(i)+subsetI);
             u_sub      = u(part.yI(i)-subsetI:part.yI(i)+subsetI,part.xI(i)-subsetI:part.xI(i)+subsetI,part.zI(i)-subsetI:part.zI(i)+subsetI,part.tI(i)-subsetI:part.tI(i)+subsetI);
             v_sub      = v(part.yI(i)-subsetI:part.yI(i)+subsetI,part.xI(i)-subsetI:part.xI(i)+subsetI,part.zI(i)-subsetI:part.zI(i)+subsetI,part.tI(i)-subsetI:part.tI(i)+subsetI);
             part.uf(i)  = interpn(lat_sub,lon_sub,alt_vec,time_sub,u_sub,part.lat(i),part.lon(i),part.z(i),erup.time+part.t(i)/3600/24);   
             part.vf(i)  = interpn(lat_sub,lon_sub,alt_vec,time_sub,v_sub,part.lat(i),part.lon(i),part.z(i),erup.time+part.t(i)/3600/24);   
        else
             int_count2  = int_count2+1;
             if int_count2==int_skip
                int_count2 = 0; 
             end
             part.uf(i)  = part.uf(i-1);
             part.vf(i)  = part.vf(i-1);
        end
    end
    
        % When u and v conditions are changing, display to the screen
        if Iprint_count==Iprint_Skip
                Iprint_count                    = 0;    
        end         
        if Iprint_count==0 
            
            table_iteration(count_table)    = {sprintf(num2str(i))};
            table_height(count_table)       = {sprintf([num2str(part.z(i), '%4.0f') ' m asl'])};
            table_bear(count_table)         = {sprintf([num2str(part.bear(i), '%3.0f') ' degree'])};
            table_uwind(count_table)        = {sprintf([num2str(part.uf(i-1), '%4.2f') ' m/s'])};
            table_vwind(count_table)        = {sprintf([num2str(part.vf(i-1), '%4.2f'), ' m/s'])};
            table_u(count_table)            = {sprintf([num2str(part.u(i-1), '%4.2f'), ' m/s'])};
            table_v(count_table)            = {sprintf([num2str(part.v(i-1), '%4.2f'), ' m/s'])};
            table_w(count_table)            = {sprintf([num2str(part.w(i-1), '%4.2f'), ' m/s'])};
            table_tau(count_table)          = {sprintf([num2str(part.tau(i), '%4.0e'), ' s'])};
            table_distance(count_table)     = {sprintf([num2str(part.dis(i)/1000, '%6.2f'), ' km'])};
            table_time(count_table)         = {sprintf([num2str(part.time(i), '%4.0e'), ' s'])};
            var_name                        = {'Iteration' 'Distance' 'Height' 'Bearing' 'Uwind' 'Vwind' 'Particle_u' 'Particle_v' 'Terminal_velocity' 'Flight_time' 'Tau'};
            
            home
            display(table(table_iteration', table_distance', table_height', table_bear', table_uwind', table_vwind', table_u', table_v', table_w', table_time', table_tau', 'VariableNames', var_name))
            
            count_table = count_table + 1;
        end
        Iprint_count = Iprint_count + 1;
    end
    
end


%% Plotting
figure;
if strcmp(vent.geo, 'West')
    plot(-360+[lon_min, lon_max], [lat_min, lat_max], '.k');                    % Plot frame
    plot_google_map('Maptype', 'terrain'); hold on;                             % Plot google background
    surf(-360+dem.X, dem.Y, dem.Z./1000); shading flat; alpha 0.5               % Plot topography
    colormap(landcolor);
    c = colorbar;
    ylabel(c, 'Topography altitude (km asl)');
    freezeColors;
    [C,h] = contour(-360+dem.X, dem.Y, ll2dist(vent.lat, vent.lon, dem.Y, dem.X));
    set(h, 'LineColor', 'k', 'Linewidth', 1.5);
    clabel(C,h);
    plot3(-360+part.lon, part.lat, part.z/1000, '.r');                          % Plot trajectory
    plot3([-360+vent.lon,-360+vent.lon], [vent.lat,vent.lat], [vent.ht/1000, erup.ht], '--k')
    plot3(-360+vent.lon, vent.lat, vent.ht/1000, '^k', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

else
    
    plot([min(dem.X(1,:)), max(dem.X(1,:))], [min(dem.Y(:,1)), max(dem.Y(:,1))], '.k');                    % Plot frame
    plot_google_map('Maptype', 'terrain'); hold on;                             % Plot google background
    surf(dem.X, dem.Y, dem.Z./1000); shading flat; alpha 0.5               % Plot topography
    colormap(landcolor);
    c = colorbar;
    ylabel(c, 'Topography altitude (km asl)');
    freezeColors;
    [C,h] = contour(dem.X, dem.Y, ll2dist(vent.lat, vent.lon, dem.Y, dem.X));
    set(h, 'LineColor', 'k', 'Linewidth', 1.5);
    clabel(C,h);
    plot3(part.lon, part.lat, part.z/1000, '.r');                          % Plot trajectory
    plot3([vent.lon,vent.lon], [vent.lat,vent.lat], [vent.ht/1000, erup.ht], '--k')
    plot3(vent.lon, vent.lat, vent.ht/1000, '^k', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

end
title(sprintf('%4.4f mm, %4.f kg/m^3, fl = %2.1f, el = %2.1f', part.diam*1000, part.dens, part.fl, part.el));




function [LAT, LON] = dist2ll(lat,lon,dx,dy)
% lat = old lat
% lon = old lon
% dx = x increm
% dy = y increm
% LAT and LON = new coordinates

R = 6371*10^3;

d = sqrt(dx^2+dy^2);
b = mod(atan2d(dx,dy),360);

LAT = asind(sind(lat) * cosd(d/R) + cosd(lat) * sin(d/R) * cosd(b));
LON = lon + atan2d( sind(b) * sin(d/R) * cosd(lat), cos(d/R) - sind(lat) * sind(LAT));

function D = ll2dist(lat1, lon1,lat2, lon2)

radius=6371;

lat1 = lat1*pi/180;
lat2 = lat2*pi/180;
lon1 = lon1*pi/180;
lon2 = lon2*pi/180;

deltaLat= lat2-lat1;
deltaLon= lon2-lon1;
a       = sin((deltaLat)./2).^2 + cos(lat1).*cos(lat2) .* sin(deltaLon./2).^2;
c       = 2.*atan2(sqrt(a),sqrt(1-a));
D       = radius.*c;    %Haversine distance




