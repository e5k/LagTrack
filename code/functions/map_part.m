function map_part(varargin)

nSteps = 100;   % Number of steps to plot

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

% Create a temporary figure to retrieve Google background
f_tmp   = figure('visible', 'off'); 
a_tmp   = axes('Parent', f_tmp);
plot(a_tmp, [dem.X(1), dem.X(end)], [dem.Y(1), dem.Y(end)], '.'); 
[lonVec, latVec, imag] = plot_google_map('Axis', a_tmp, 'Maptype', 'terrain');
delete(f_tmp);

% Set topography and corrects ratio
surface( dem.X, dem.Y, dem.Z./1000, prepare_google_map(dem, lonVec, latVec, imag), 'Parent', AX); % Map the background to the topography
shading(AX, 'flat'); hold(AX, 'on');   grid(AX, 'on'); 
axis(AX, 'tight')
lat_lon_proportions(AX);

for i = 1:length(fld)
    % Plot vent
    plot3(AX, pltData.(fld{i}).vent.lon, pltData.(fld{i}).vent.lat, pltData.(fld{i}).vent.alt./1000,...
        '^', 'MarkerSize', 10, 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k',...
        'Tag', ['vent',pltData.(fld{i}).part.name]);
    
    % Plot vertical line
    plot3(AX, repmat(pltData.(fld{i}).vent.lon,2,1), repmat(pltData.(fld{i}).vent.lat,2,1), [pltData.(fld{i}).vent.alt./1000; pltData.(fld{i}).traj.z(1)./1000],...
        ':', 'Color', cmap(i,:), 'linewidth',1,...
        'Tag', ['alt',pltData.(fld{i}).part.name]);
    
    % Plot trajectory
    idx = ceil(linspace(1, length(pltData.(fld{i}).traj.lon), nSteps));
    legH(i) = plot3(AX, pltData.(fld{i}).traj.lon(idx), pltData.(fld{i}).traj.lat(idx), pltData.(fld{i}).traj.z(idx)./1000,...
        '-', 'Color', cmap(i,:), 'linewidth',3,...
        'Tag', ['traj',pltData.(fld{i}).part.name]);
    
    % Plot start and stop
    plot3(AX, [pltData.(fld{i}).traj.lon(1), pltData.(fld{i}).traj.lon(end)],...
        [pltData.(fld{i}).traj.lat(1), pltData.(fld{i}).traj.lat(end)],...
        [pltData.(fld{i}).traj.z(1)./1000, pltData.(fld{i}).traj.z(end)./1000],...
        'o', 'MarkerSize', 10, 'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k',...
        'Tag', ['end',pltData.(fld{i}).part.name]);
    
    leg{i} = pltData.(fld{i}).part.name;
end
box(AX, 'on')

xlabel('Longitude');
ylabel('Latitude');
zlabel('Altitude (km asl)');

legend(AX, legH, leg, 'Location', 'Best', 'Tag', 'LegMap', 'Interpreter', 'none');
AX.Position = POS;
