function map_part(varargin)
% map_part Plot particles on a 3D map
%   map_part
%       Open the GUI to select particles
%
%   see also plot_part, detail_part.
% 
% This function is part of LagTrack.
% Written by Sebastien Biass & Gholamhossein Bagheri
% GPLv3

% Plotting steps to plot only a subset of the entire trajectory
nSteps = 100;   % Number of steps to plot
nWind  = 10;    % Percentage of nSteps to plot wind vectors

if nargin == 0  % If calling the function independently
    [fl, pth] = uigetfile('projects/*.mat', 'Multiselect', 'on');
    if isnumeric(fl(1)); return; end

    pltData = struct;
    
    for i = 1:length(fl)
        % Check if single or multiple files
        if iscell(fl)
            load(fullfile(pth, fl{i}));
        else
            load(fullfile(pth, fl));
        end
        
        % Check if the traj field exists to validate file
        if ~isfield(part, 'traj')
            errordlg('The selected file is not a valid run')
            return
        end
        
        % Retrieve data
        pltData.(part.part.name) = part;
    end
    
    FG      = figure;
    AX      = axes('Parent', FG);                                   % Set plotting target - i.e. new figure
    fld     = fieldnames(pltData);                                  % Retrieve fieldnames
else % If the function is called through the GUI
    src     = varargin{1};
    pltData = varargin{2};
    fld     = fieldnames(pltData);                                  % Retrieve fieldnames
    
     if nargin == 2
        AX      = findobj(ancestor(src, 'figure'), 'Tag', 'MapAx');    % Set plotting target - i.e. GUI axis
        cla(AX);
        delete(findobj(ancestor(src, 'figure'), 'Tag', 'LegMap'));
    elseif nargin == 3
        FG      = figure;
        AX      = axes('Parent', FG);                                   % Set plotting target - i.e. new figure
    end
end

POS = AX.Position;

% Retrieve the dem for plotting
load(pltData.(fld{1}).path.dem); 

cmap    = lines(length(fld));     % Setup colormap
leg     = cell(length(fld),1);  % Setup legend
legH    = zeros(length(fld),1); % Legend handles
hold(AX, 'on')
for i = 1:length(fld)
    % Indices for plotting only a subset of the trajectory (nSteps)
    if length(pltData.(fld{i}).traj.lon) >= nSteps
        idx = ceil(linspace(1, length(pltData.(fld{i}).traj.lon), nSteps));     
    else
        idx = 1:length(pltData.(fld{i}).traj.lon);
    end
    idx_wind = ceil(linspace(1, length(pltData.(fld{i}).traj.lon), nSteps/nWind));
        
    if strcmp(dem.type, 'DEM')
        vent_lon    = pltData.(fld{i}).vent.lon;
        vent_lat    = pltData.(fld{i}).vent.lat;
        part_x      = pltData.(fld{i}).traj.lon(idx);
        part_y      = pltData.(fld{i}).traj.lat(idx);
        wind_x      = pltData.(fld{i}).traj.lon(idx_wind);
        wind_y      = pltData.(fld{i}).traj.lat(idx_wind);
    else
        vent_lon    = 0;
        vent_lat    = 0;
        part_x      = pltData.(fld{i}).traj.x(idx);
        part_y      = pltData.(fld{i}).traj.y(idx);
        wind_x      = pltData.(fld{i}).traj.x(idx_wind);
        wind_y      = pltData.(fld{i}).traj.y(idx_wind);
    end
    % Plot vent
    if pltData.(fld{i}).run_mode == 1
        plot3(AX, vent_lon, vent_lat, pltData.(fld{i}).vent.alt./1000,...
            '^', 'MarkerSize', 10, 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k',...
            'Tag', ['vent',pltData.(fld{i}).part.name]);
    end
    
    % Plot vertical line
    plot3(AX, repmat(vent_lon,2,1), repmat(vent_lat,2,1), [pltData.(fld{i}).vent.alt./1000; pltData.(fld{i}).traj.z(1)./1000],...
        ':', 'Color', cmap(i,:), 'linewidth',1,'Tag', ['alt',pltData.(fld{i}).part.name]);
      
    % Plot trajectory
    legH(i) = plot3(AX, part_x, part_y, pltData.(fld{i}).traj.z(idx)./1000,...
        '-', 'Color', cmap(i,:), 'linewidth',3,'Tag', ['traj',pltData.(fld{i}).part.name]);
     
    % Plot start and stop
    plot3(AX, [part_x(1), part_x(end)],...
        [part_y(1), part_y(end)],...
        [pltData.(fld{i}).traj.z(1)./1000, pltData.(fld{i}).traj.z(end)./1000],...
        'o', 'MarkerSize', 10, 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k',...
        'Tag', ['end',pltData.(fld{i}).part.name]);
    
    leg{i} = pltData.(fld{i}).part.name;
    
    % If only one particle, plot the wind field
    if length(fld) == 1
    legH(i+1) = quiver3(AX, wind_x, wind_y,pltData.(fld{i}).traj.z(idx_wind)./1000,...
        pltData.(fld{i}).traj.uf(idx_wind),pltData.(fld{i}).traj.vf(idx_wind), zeros(size(pltData.(fld{i}).traj.u(idx_wind))),...
        1e-2, 'Color', 'r');
    leg{i+1} = 'Wind';
    
    legH(i+2) = quiver3(AX, wind_x, wind_y,pltData.(fld{i}).traj.z(idx_wind)./1000,...
        pltData.(fld{i}).traj.u(idx_wind),pltData.(fld{i}).traj.v(idx_wind), pltData.(fld{i}).traj.w(idx_wind),...
        1e-2, 'Color', 'k');
    leg{i+2} = 'Particle';
    
    
    end
    
end
axis tight

% Retrieve the extent of plotted particles
fact = .1; % Factor to extend the map extent
XMin = AX.XLim(1) - fact*abs(AX.XLim(2)-AX.XLim(1));
XMax = AX.XLim(2) + fact*abs(AX.XLim(2)-AX.XLim(1));
YMin = AX.YLim(1) - fact*abs(AX.YLim(2)-AX.YLim(1));
YMax = AX.YLim(2) + fact*abs(AX.YLim(2)-AX.YLim(1));

% Indices of the new extent on the dem
[~, XiMin] = min(abs(dem.X(1,:)-XMin));
[~, XiMax] = min(abs(dem.X(1,:)-XMax));
[~, YiMin] = min(abs(dem.Y(:,1)-YMin));
[~, YiMax] = min(abs(dem.Y(:,1)-YMax));

% If extent covers only one pixel, extend it to at least 3 pixels
if XiMin-XiMax<3 && XiMin>2 && XiMax<size(dem.X,2)-1
    XiMin = XiMin-2; XiMax = XiMax+2;
    XMin  = dem.X(1,XiMin); XMax  = dem.X(1,XiMax); 
end
if YiMin-YiMax<3 && YiMin>2 && YiMax<size(dem.Y,1)-1
    YiMin = YiMin-2; YiMax = YiMax+2;
    YMin  = dem.Y(YiMin,1); YMax  = dem.Y(YiMax,1); 
end

if strcmp(dem.type, 'DEM') % In case the grid is a DEM
    % If the Google Map API exists
    if exist(fullfile('code', 'var', 'google_api.mat'), 'file') == 2
        % Create a temporary figure to retrieve Google background
        f_tmp   = figure('visible', 'off');
        a_tmp   = axes('Parent', f_tmp);
        plot(a_tmp, [XMin, XMax], [YMin, YMax], '.');
        [lonVec, latVec, imag] = plot_google_map('Axis', a_tmp, 'Maptype', 'terrain');
        delete(f_tmp);
        
        if ~isempty(imag)
            % Interpolate for a sharp background
            [Xp, Yp] = meshgrid(linspace(dem.X(1,XiMin), dem.X(1,XiMax), size(imag,2)), linspace(dem.Y(YiMin,1), dem.Y(YiMax,1), size(imag,1)));
            Zp       = interp2(dem.X, dem.Y, dem.Z, Xp, Yp);
            
            % Set topography and corrects ratio
            surface( Xp,Yp,Zp./1000,...
                prepare_google_map(Xp, Yp, lonVec, latVec, imag), 'Parent', AX); % Map the background to the topography
            shading(AX, 'flat'); grid(AX, 'on');
        end
    else
        surf(dem.X(YiMin:YiMax, XiMin:XiMax), dem.Y(YiMin:YiMax,XiMin:XiMax), dem.Z(YiMin:YiMax,XiMin:XiMax)./1e3);
        shading(AX, 'flat'); grid(AX, 'on'); colormap(landcolor);
    end
    % Work on axes
    lat_lon_proportions(AX);

    xlabel('Longitude');
    ylabel('Latitude');
else
    xlabel('X displacement (m)');    ylabel('Y displacement (m)');
end

box(AX, 'on')
zlabel('Altitude (km asl)');
axis(AX, [XMin, XMax, YMin, YMax]);

% Legend
legend(AX, legH, leg, 'Location', 'Best', 'Tag', 'LegMap', 'Interpreter', 'none');
AX.Position = POS;
