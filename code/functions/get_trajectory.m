function part = get_trajectory(P)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If a single particle is given, convert the structure to a cell
if isstruct(P); P = {P}; end

% Check that all particles point to the same dem/atmospheric file
for iP = 2:length(P)
    if length(unique(cellfun(@getfield, P,repmat({'path'}, size(P)), repmat({'nc'}, size(P)), 'UniformOutput', false))) > 1
        error('All particles must point to the same atmospheric data');
    end
    
    if length(unique(cellfun(@getfield, P,repmat({'path'}, size(P)), repmat({'dem'}, size(P)), 'UniformOutput', false))) > 1
        error('All particles must point to the same DEM');
    end
    
    if length(unique(cellfun(@getfield, P,repmat({'run_name'}, size(P)),  'UniformOutput', false))) > 1
        error('All particles must be part of the same project');
    end
end

% Create output folder if does not exist
if ~exist(['projects', filesep, P{1}.run_name], 'dir')
    mkdir(['projects', filesep, P{1}.run_name]);
end

% Print message
fprintf('_____________________________________________________________________________________\n')
fprintf('LAGTRACK run %s started on %s...\n\n', P{1}.run_name, datestr(now));

% Load input data
atm  = load(P{1}.path.nc); atm = atm.atm;                           % Atmospheric data
dem  = load(P{1}.path.dem); dem = dem.dem;                          % DEM

% Input parameters
if length(P) == 1 % In case one particle is input
    multi= 0;                                                        % If single particle
    part = run_trajectory(P{1}, atm, dem, multi);
else %iscell(P) % In case multiple particles are input    
    multi= 1;                                                        % If multiple particles
    if license('checkout', 'Distrib_Computing_Toolbox') == 1 && length(P)>1 % If parallel toolbox available
        gcp;
        parfor iP = 1:length(P)
            run_trajectory(P{iP}, atm, dem, multi);
        end
    else
        for iP = 1:length(P)
            run_trajectory(P{iP}, atm, dem, multi);
        end    
    end
end

fprintf('LAGTRACK run %s ended on %s...\n\n', P{1}.run_name, datestr(now));
fprintf('_____________________________________________________________________________________\n')

function part = run_trajectory(P, atm, dem, multi)

% In case the grid is not a DEM and standard atmosphere, set the vent
% reference to 0
if strcmp(dem.type, 'GRID') && ~isfield(atm, 'humid')  
    P.vent.lon = 0;
    P.vent.lat = 0;
    P.vent.alt = dem.Z;
end

% Record when the particle was ran
P.timestamp     = now;

% Initialize particle release position/time
part.x(1)       = 0;                                                        % X (m, Positive in E, negative in W)
part.y(1)       = 0;                                                        % Y (m, Positive in N, negative in S)
part.t(1)       = 0;                                                        % Time of eruption
part.z(1)       = P.vent.alt+P.rel.z;                                       % Z (m asl)
part.dis(1)     = 0;                                                        % Distance along the flight path
part.disP(1)    = 0;                                                        % Distance projected on the Earth surface using the Haversine formula
part.bear(1)    = 0;                                                        % Particle bearing (Degree from north)

% Calculate lat and lon of particle release
[part.lat(1),part.lon(1)] = dist2ll(P.vent.lat, P.vent.lon, part.x+P.rel.x, part.y+P.rel.y);
% Distance from release point to the vent
part.dis0       = ll2dist(P.vent.lat, P.vent.lon, part.lat(1), part.lon(1), part.z(1)); 

% Find particle indices relative to atmospheric conditions
[~, part.xI(1)] = min(abs(atm.lon - part.lon(1)));
[~, part.yI(1)] = min(abs(atm.lat - part.lat(1)));
[~, part.tI(1)] = min(abs(atm.time - P.date + (P.rel.t/3600/24) ));         % Add time offset to eruption date
[~, part.zI(1)] = min(abs(atm.alt(part.yI(1), part.xI(1), :, part.tI(1)) - part.z(1)));

% Find particle indices relative to DEM
[~, part.xD(1)] = min(abs(dem.X(1,:) - part.lon(1)));
[~, part.yD(1)] = min(abs(dem.Y(:,1) - part.lat(1)));
% Check if release altitude is above DEM
if part.z(1) < dem.Z(part.yD(1), part.xD(1))
    part.out_msg = sprintf('Release altitude is below the DEM (%4.0f m)', dem.Z(part.yD(1), part.xD(1)));
    error('Release altitude is below the DEM (%4.0f m)', dem.Z(part.yD(1), part.xD(1))),
end
% If run mode is backward, overwrite the initial part.z(i) to DEM value + release elevation
if P.run_mode == 2
    part.z(1)   = dem.Z(part.yD(1), part.xD(1)) + P.rel.z;
    [~, part.zI(1)] = min(abs(atm.alt(part.yI(1), part.xI(1), :, part.tI(1)) - part.z(1)));
end

% Initial atmosperic condidtions
part.uf(1)      = atm.u(part.yI(1), part.xI(1), part.zI(1), part.tI(1));    % Initial u wind (m/s)
part.vf(1)      = atm.v(part.yI(1), part.xI(1), part.zI(1), part.tI(1));    % Initial v wind (m/s)
part.wf(1)      = 0;

% Initial dynamic parameters
part.Fs         = P.part.flat * P.part.elon^1.3;                            % Shape descriptor - Stoke's regime
part.Fn         = P.part.flat^2 * P.part.elon;                              % Shape descriptor - Newton regime
part.Ks         = .5 * (part.Fs^(1/3) + part.Fs^(-1/3));                    % Stoke's drag correction
part.Kn         = 10^(.45 * (-log10(part.Fn))^.99);                         % Newton's drag correction     
part.tau(1)     = P.part.dens * P.part.diam^2 / (18 * atm.muair(part.yI(1), part.xI(1), part.zI(1), part.tI(1)));    % Particle relaxation time (s)
% Mohssen: Are we actually using that?
part.Re_w(1)    = 0;        
part.Cd_w(1)    = 0;


% Initial particle velocities
% If no initial x or y velocity is given, assume the particle is carried by the wind
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

% Control
part.out_msg    = ' ';                                                      % Initial output message
test_run        = 0;                                                        % Control variable

% Constants
earth_radius    = 6371*1e3;                                                 % Earth radius
g               = 9.806;                                                    % Gravity m/s2  
if P.run_mode == 2
    g = -g;
end

% Initialize counters
int_count1      = 0;                                                        % Counter for skip check for density and viscosity
int_count2      = 0;                                                        % Counter forskip check for wind velocity
i               = 1;                                                        % Main iteration counter

% If single particle, print the detail
if multi == 0
    fprintf('Alt.\tDist.\t\tBearing\tw\tlat\tlon\tDEM\tTime\tuwind\tvwind\n')
    fprintf('______\t______\t\t______\t______\t______\t______\t______\t______\t______\t______\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of particle trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while test_run == 0
    i = i+1;
    
    % If single particle, print the detail
    if multi == 0
        if mod(i,500) == 0  % Print the particle status every 500 iterations
            fprintf('%4.0f\t%4.0f\t\t%3.2f\t%3.1f\t%2.2f\t%3.2f\t%4.0f\t%4.0f\t%3.2f\t%3.2f\n', part.z(end), part.disP(end), part.w(end), part.bear(end), part.lat(end), part.lon(end), dem.Z(part.yD(end),part.xD(end)), part.t(end), part.uf(end), part.vf(end))
        end
    end
    
    % Interpolation of physical properties of atmosphere     
    if strcmp(P.adv.interp, 'none')
    % Simple indexing method
        denf    = atm.rhoair(part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));   % Fluid density
        visf    = atm.muair( part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));   % Fluid viscosity
        tempf   = atm.temp( part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));    % Air temperature
    elseif strcmp(P.adv.interp, 'complete')
    % Complete interpolation method (Slow)
        alt_sub = squeeze(atm.alt(part.yI(i-1), part.xI(i-1), :, part.tI(i-1)));
        denf    = interpn(double(atm.lat), double(atm.lon), alt_sub, atm.time, atm.rhoair, part.lat(i-1), part.lon(i-1), part.z(i-1), P.date+part.t(i-1)/3600/24, P.adv.method);
        visf    = interpn(double(atm.lat), double(atm.lon), alt_sub, atm.time, atm.muair,  part.lat(i-1), part.lon(i-1), part.z(i-1), P.date+part.t(i-1)/3600/24, P.adv.method);
        tempf   = interpn(double(atm.lat), double(atm.lon), alt_sub, atm.time, atm.temp,  part.lat(i-1), part.lon(i-1), part.z(i-1), P.date+part.t(i-1)/3600/24, P.adv.method);
        if isnan(denf)
            error('Interpolation returned a Nan. Trying changing the interpolation method.')
        end
        
    elseif strcmp(P.adv.interp, 'subset')
    % Interpolation by subsetting method
        % check if the subsetting range is within in the downloaded data
        if part.xI(i-1)-P.adv.range < 1 || part.xI(i-1)+P.adv.range > size(atm.lon,1) || part.yI(i-1)-P.adv.range < 1 || part.yI(i-1)+P.adv.range > size(atm.lat,1) || part.zI(i-1)-P.adv.range < 1 || part.zI(i-1)+P.adv.range > size(atm.alt,3)
             denf       = atm.rhoair(part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));   % Fluid density
             visf       = atm.muair( part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));   % fluid viscosity
             tempf      = atm.temp( part.yI(i-1), part.xI(i-1), part.zI(i-1), part.tI(i-1));    % Air temperature
        elseif int_count1==0 || P.adv.skip==0
             int_count1 = int_count1+1;
             alt_vec    = squeeze(atm.alt(part.yI(i-1), part.xI(i-1), part.zI(i-1)-P.adv.range:part.zI(i-1)+P.adv.range, part.tI(i-1)));
             time_sub   = atm.time(part.tI(i-1) - P.adv.range:part.tI(i-1) + P.adv.range);
             lon_sub    = atm.lon(part.xI(i-1)  - P.adv.range:part.xI(i-1) + P.adv.range);
             lat_sub    = atm.lat(part.yI(i-1)  - P.adv.range:part.yI(i-1) + P.adv.range);
             den_sub    = atm.rhoair(part.yI(i-1)-P.adv.range:part.yI(i-1) + P.adv.range, part.xI(i-1)-P.adv.range:part.xI(i-1)+P.adv.range, part.zI(i-1)-P.adv.range:part.zI(i-1)+P.adv.range, part.tI(i-1)-P.adv.range:part.tI(i-1)+P.adv.range);
             vis_sub    = atm.muair(part.yI(i-1) -P.adv.range:part.yI(i-1) + P.adv.range, part.xI(i-1)-P.adv.range:part.xI(i-1)+P.adv.range, part.zI(i-1)-P.adv.range:part.zI(i-1)+P.adv.range, part.tI(i-1)-P.adv.range:part.tI(i-1)+P.adv.range);
             temp_sub   = atm.temp(part.yI(i-1) -P.adv.range:part.yI(i-1) + P.adv.range, part.xI(i-1)-P.adv.range:part.xI(i-1)+P.adv.range, part.zI(i-1)-P.adv.range:part.zI(i-1)+P.adv.range, part.tI(i-1)-P.adv.range:part.tI(i-1)+P.adv.range);
             denf       = interpn(lat_sub, lon_sub, alt_vec, time_sub, den_sub, part.lat(i-1), part.lon(i-1), part.z(i-1), P.date+part.t(i-1)/3600/24, P.adv.method);   % Fluid density
             visf       = interpn(lat_sub, lon_sub, alt_vec, time_sub, vis_sub, part.lat(i-1), part.lon(i-1), part.z(i-1), P.date+part.t(i-1)/3600/24, P.adv.method);   % Fluid viscosity
             tempf      = interpn(lat_sub, lon_sub, alt_vec, time_sub, temp_sub, part.lat(i-1), part.lon(i-1), part.z(i-1), P.date+part.t(i-1)/3600/24, P.adv.method);  % Air temperature
         else
             int_count1  = int_count1+1;
             if int_count1==P.adv.skip
                  int_count1 = 0; 
             end
         end
    end

    % Update particle state
    velr            = sqrt((part.u(i-1)-part.uf(i-1))^2 + (part.v(i-1)-part.vf(i-1))^2 + (part.w(i-1)-part.wf(i-1))^2);	% Particle relative velocity

    if sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2)==0
        theta       = 90;
    else
        % theta       = atand((part.vf(i-1) - part.v(i-1)) / sqrt((part.u(i-1) - part.uf(i-1))^2 + (part.v(i-1) - part.vf(i-1))^2));	% direction of particle in the horizontal plane (x-y)
        theta       = atand( (part.v(i-1) - part.vf(i-1)) / (part.u(i-1) - part.uf(i-1)) );	% direction of particle in the horizontal plane (x-y) % Fixed 2020-20-09
    end
    if sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2)==0
        Beta        = 90;
    else
        % Beta        = atand(part.w(i-1)/sqrt((part.u(i-1)-part.uf(i-1))^2+(part.v(i-1)-part.vf(i-1))^2)); % direction of particle in the vertical plane (x-z or y-z)
        Beta        = atand( (part.w(i-1)-part.wf(i-1)) / sqrt( (part.u(i-1)-part.uf(i-1))^2 + (part.v(i-1)-part.vf(i-1))^2) ); % direction of particle in the vertical plane (x-z or y-z) % Fixed 2020-20-09
    end
    
    % Avoid possible divisions by 0
    if ismember(theta, [0,90,180,270]); theta = theta + 1e-6; end
    if ismember(Beta, [0,90,180,270]); Beta = Beta + 1e-6; end
    
    part.Re(i)      = denf * velr * P.part.diam / visf;                     % Reynolds
    part.Re_S(i)    = part.Re(i) * part.Kn / part.Ks;                       % Ganser Re
    part.Mach(i)    = velr / sqrt(1.4*286*273.15+tempf);                    % Mach number
    
    % Drag coefficient
    if part.Re_S(i) > 3e5
        part.Cd(i)  = 0.2;
    else
        part.Cd(i)  = part.Kn * (24 * (1 + .125 * part.Re_S(i)^(2/3)) / part.Re_S(i) + .46 / (1+5330/part.Re_S(i)));  % Drag coef (eq. 34)
    end
    
    % Correct drag coefficient at high mach number following the method of Mastin/Eject!
    if P.adv.mach ~= -1 && part.Mach(i)>P.adv.mach
       part.Cd(i)   = getCdHighMach(part.Mach(i), P.adv.machMeth); 
    end
    
    % If particle within region of reduced drag
    if part.dis(i-1) < P.adv.drag
        part.Cd(i)  = part.Cd(i) * ( part.dis(i-1) / P.adv.drag)^2;
    end
    
    Fd              = part.Cd(i) * part.Re(i) / (24*part.tau(i-1));         % Total Drag force
    
    % Added 2020-03-09 to avoid division by 0 - check with Mohssen
    if Fd == 0
        Fd = 1e-6;
    end
    
    Fd_u            = abs(Fd * cosd(Beta) * cosd(theta));                   % Drag force in x direction
    G_u             = 0.;                                                   % Other forces in x direction 
    Fd_v            = abs(Fd * cosd(Beta) * sind(theta));                   % Drag force in y direction
    G_v             = 0.;                                                   % Other forces in y direction
    Fd_w            = abs(Fd * sind(Beta));                                 % Drag force in z direction
    G_w             = (1-denf/P.part.dens) * -g;                            % Gravity force
    
    if strcmp(P.adv.solution, 'analytical')
    % 2020-03-09: Fixed problem - i.e. the full drag force was considered instead of its components
    % Analytical solution
%         part.u(i)   = part.uf(i-1) + exp(-Fd * P.adv.dt) *(part.u(i-1)-part.uf(i-1)) - G_u * (1/Fd) * (exp(-P.adv.dt*Fd)-1);
%         part.x(i)   = part.x(i-1) + (G_u * (1/Fd)+part.uf(i-1)) * P.adv.dt + (1/Fd)* (1-exp(-P.adv.dt*Fd)) * (part.u(i)-part.uf(i-1)- G_u/Fd);
        part.u(i)   = part.uf(i-1) + exp(-Fd_u * P.adv.dt) *(part.u(i-1)-part.uf(i-1)) - G_u * (1/Fd_u) * (exp(-P.adv.dt*Fd_u)-1);
        part.x(i)   = part.x(i-1) + (G_u * (1/Fd_u)+part.uf(i-1)) * P.adv.dt + (1/Fd_u)* (1-exp(-P.adv.dt*Fd_u)) * (part.u(i)-part.uf(i-1) - G_u/Fd_u);

%         part.v(i)   = part.vf(i-1) + exp(-Fd * P.adv.dt) *(part.v(i-1)-part.vf(i-1)) - G_v * (1/Fd) * (exp(-P.adv.dt*Fd)-1);
%         part.y(i)   = part.y(i-1) + (G_v * (1/Fd)+part.vf(i-1)) * P.adv.dt + (1/Fd)* (1-exp(-P.adv.dt*Fd)) * (part.v(i)-part.vf(i-1)- G_v/Fd);
        part.v(i)   = part.vf(i-1) + exp(-Fd_v * P.adv.dt) *(part.v(i-1)-part.vf(i-1)) - G_v * (1/Fd_v) * (exp(-P.adv.dt*Fd_v)-1);
        part.y(i)   = part.y(i-1) + (G_v * (1/Fd_v)+part.vf(i-1)) * P.adv.dt + (1/Fd_v)* (1-exp(-P.adv.dt*Fd_v)) * (part.v(i)-part.vf(i-1)- G_v/Fd_v);
        
%         part.w(i)   = part.wf(i-1) + exp(-Fd * P.adv.dt) *(part.w(i-1)-part.wf(i-1)) - G_w * (1/Fd) * (exp(-P.adv.dt*Fd)-1);
%         part.z(i)   = part.z(i-1) + (G_w * (1/Fd)+part.wf(i-1)) * P.adv.dt + (1/Fd)* (1-exp(-P.adv.dt*Fd)) * (part.w(i)-part.wf(i-1)- G_w/Fd);        
        part.w(i)   = part.wf(i-1) + exp(-Fd_w * P.adv.dt) *(part.w(i-1)-part.wf(i-1)) - G_w * (1/Fd_w) * (exp(-P.adv.dt*Fd_w)-1);
        part.z(i)   = part.z(i-1) + (G_w * (1/Fd_w)+part.wf(i-1)) * P.adv.dt + (1/Fd_w)* (1-exp(-P.adv.dt*Fd_w)) * (part.w(i)-part.wf(i-1)- G_w/Fd_w);
                      
    elseif strcmp(P.adv.solution, 'euler')
    % Euler semi-implicit 
        part.u(i)   = ((G_u + Fd_u * part.uf(i-1)) * P.adv.dt + part.u(i-1)) / (1+Fd_u * P.adv.dt);
        part.x(i)   = part.x(i-1) + .5 * P.adv.dt * (part.u(i) + part.u(i-1));

        part.v(i)   = ((G_v + Fd_v * part.vf(i-1)) * P.adv.dt + part.v(i-1)) / (1+Fd_v * P.adv.dt);
        part.y(i)   = part.y(i-1) + .5 * P.adv.dt * (part.v(i) + part.v(i-1));

        part.w(i)   = ((G_w + Fd_w * part.wf(i-1)) * P.adv.dt  + part.w(i-1)) / (1+Fd_w * P.adv.dt);
        part.z(i)   = part.z(i-1) + .5 * P.adv.dt * (part.w(i) + part.w(i-1));
        
    end
    
    % Update particle time
    if P.run_mode == 1
        part.t(i) = part.t(i-1) + P.adv.dt;
    else
        part.t(i) = part.t(i-1) - P.adv.dt;
    end
    
    % Update geographic coordinates
    [part.lat(i), part.lon(i)] = dist2ll(part.lat(i-1), part.lon(i-1), part.x(i)-part.x(i-1), part.y(i)-part.y(i-1));
    
    % Update bearing
    part.bear(i) = atan2d( sind(abs(part.lon(i)-part.lon(i-1))) * cosd(part.lat(i)) ,  cosd(part.lat(i-1)) * sind(part.lat(i)) - sind(part.lat(i-1)) * cosd(part.lat(i))  * cosd(abs(part.lon(i)-part.lon(i-1))) );
    if part.bear(i) < 0
        part.bear(i)= 360+part.bear(i);
    end
    
    % Update indices relative to atmospheric parameters
    [~, part.xI(i)] = min(abs(atm.lon-part.lon(i)));
    [~, part.yI(i)] = min(abs(atm.lat-part.lat(i)));
    [~, part.tI(i)] = min(abs(atm.time-(P.date + part.t(i)/3600/24)));
    [~, part.zI(i)] = min(abs(atm.alt(part.yI(i), part.xI(i), :, part.tI(i)) - part.z(i)));
    part.dis(i)     = sqrt((part.x(i)-part.x(1))^2+(part.y(i)-part.y(1))^2+(part.z(i)-part.z(1))^2);
    part.disP(i)    = ll2dist(part.lat(1), part.lon(1), part.lat(i), part.lon(i), earth_radius);
    
    % Update indices relative to DEM
    [~, part.xD(i)] = min(abs(dem.X(1,:)-part.lon(i)));
    [~, part.yD(i)] = min(abs(dem.Y(:,1)-part.lat(i)));
       
    % Update tau
    part.tau(i) = P.part.dens * P.part.diam^2 / (18 * atm.muair(part.yI(i), part.xI(i), part.zI(i), part.tI(i)));    % Particle relaxation time (s)
      
    % Test conditions
    % Test altitude of forward run
    if P.run_mode == 1 && part.z(i) <= dem.Z(part.yD(i), part.xD(i))
        part.out_msg = sprintf('\tParticle %s landed on the domain', P.part.name);
        test_run     = 1;       
    % Test altitude of backward run 
    elseif P.run_mode == 2 && part.z(i) >= P.vent.alt
        part.out_msg = sprintf('\tParticle %s reached the maximum altitude', P.part.name);
        test_run     = 1;       
    % Test domain
    %   1: Both DEM and atmospheric data are defined
%     elseif (isfield(atm, 'humid') && strcmp(dem.type, 'DEM')) && (part.lon(i) <= min([min(atm.lon), min(dem.X(1,:))]) || part.lon(i) >= max([max(atm.lon), max(dem.X(1,:))]) || part.lat(i) <= min([min(atm.lat), min(dem.Y(:,1))]) || part.lat(i) >= max([max(atm.lat), max(dem.Y(:,1))]))    
%         part.out_msg = sprintf('Particle reached the domain border');
%         test_run     = 1;  
    %   2: Only atmospheric data is defined && strcmp(dem.type, 'GRID')) 
    elseif isfield(atm, 'humid') && (part.lon(i) <= min(atm.lon) || part.lon(i) >= max(atm.lon) || part.lat(i) <= min(atm.lat) || part.lat(i) >= max(atm.lat)) 
        part.out_msg = sprintf('\tParticle %s reached the border of atmospheric domain', P.part.name);
        test_run     = 1;    
    %   3: Only DEM is defined (~isfield(atm, 'humid') && 
    elseif strcmp(dem.type, 'DEM') && (part.lon(i) <= min(dem.X(1,:)) || part.lon(i) >= max(dem.X(1,:)) || part.lat(i) <= min(dem.Y(:,1)) || part.lat(i) >= max(dem.Y(:,1))) 
        part.out_msg = sprintf('\tParticle %s reached the border of DEM', P.part.name);
        test_run     = 1;      
    % Test time
    elseif isfield(atm, 'humid') && ((P.date + i*P.adv.dt/3600/24) > max(atm.time)+0.25)
        part.out_msg = sprintf('\tParticle %s residence time is longer than atmospheric data', P.part.name);
        test_run     = 1;
    else                             
        % Get u and v wind coordinates at new location
        if strcmp(P.adv.interp, 'none')
            % simple indexing method
            part.uf(i)   = atm.u(part.yI(i), part.xI(i), part.zI(i), part.tI(i));
            part.vf(i)   = atm.v(part.yI(i), part.xI(i), part.zI(i), part.tI(i));
            part.wf(i)   = 0; % No data for vertical velocity
            
        elseif strcmp(P.adv.interp, 'complete')
            %complete interpolation method (Slow)
            alt_sub      = squeeze(atm.alt(part.yI(i), part.xI(i), :, part.tI(i)));
            part.uf(i)   = interpn(double(atm.lat), double(atm.lon), alt_sub, atm.time, atm.u, part.lat(i), part.lon(i), part.z(i), P.date+part.t(i)/3600/24, P.adv.method);
            part.vf(i)   = interpn(double(atm.lat), double(atm.lon), alt_sub, atm.time, atm.v, part.lat(i), part.lon(i), part.z(i), P.date+part.t(i)/3600/24, P.adv.method);
            part.wf(i)   = 0; % No data for vertical velocity
            
        elseif strcmp(P.adv.interp, 'subset')
            %interpolation by subsetting method
            % check if the subsetting range is within in the downloaded data
            if part.xI(i)-P.adv.range < 1 || part.xI(i)+P.adv.range > size(atm.lon,1) || part.yI(i)-P.adv.range < 1 || part.yI(i)+P.adv.range > size(atm.lat,1) || part.zI(i)-P.adv.range < 1 || part.zI(i)+P.adv.range > size(atm.alt,3)
                part.uf(i) = atm.u(part.yI(i), part.xI(i), part.zI(i), part.tI(i));
                part.vf(i) = atm.v(part.yI(i), part.xI(i), part.zI(i), part.tI(i));
                part.wf(i) = 0; % No data for vertical velocity
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
                part.wf(i) = 0; % No data for vertical velocity
            else
                int_count2 = int_count2+1;
                if int_count2==P.adv.skip
                    int_count2 = 0;
                end
                part.uf(i) = part.uf(i-1);
                part.vf(i) = part.vf(i-1);
                part.wf(i) = 0; % No data for vertical velocity
            end
        end
    end    
end

% Display message
if multi == 0
    fprintf('______\t______\t\t______\t______\t______\t______\t______\t______\t______\t______\n')
end

% Trick to make the uf and vf vectors the same size as everything else
part.uf(i) = part.uf(i-1); 
part.vf(i) = part.vf(i-1);
part.wf(i) = part.wf(i-1);

% Clean output
P.traj = part;
part   = P;

fprintf([P.traj.out_msg, '\n']); 

% Save particle
save(['projects', filesep, part.run_name, filesep, part.part.name, '.mat'], 'part');