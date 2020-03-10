function plot_part(varargin)
% plot_part Plot particle parameters along flight path
%   plot_part
%       Open the GUI to select particles
%
%   see also map_part, detail_part.
% 
% This function is part of LagTrack.
% Written by Sebastien Biass & Gholamhossein Bagheri
% GPLv3


varList     =     {'Time (s)',...    
    'Altitude (m asl)',...
    'Projected distance (m)',...
    'X distance (m)',...
    'Y distance (m)',...
    'Flight distance (m)',...
    'Latitude',...
    'Longitude',...
    'U velocity (m/s)',...
    'V velocity (m/s)',...
    'Z velocity (m/s)',...
    'U wind (m/s)',...
    'V wind (m/s)',...
    'Re',...
    'Ganser Re',...
    'Drag coefficient', ...
    'Relaxation time (s)',...
    'Mach number'};

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
    
    % Choose variables to plot
    [vX,c] = listdlg('PromptString','Select the X variable:',...
                'SelectionMode','single',...
                'ListString',varList);
    if c==0; return; end
    [vY,c] = listdlg('PromptString','Select the Y variable:',...
                'SelectionMode','single',...
                'ListString',varList);
    if c==0; return; end
    
    FG      = figure;
    AX      = axes('Parent', FG);                                   % Set plotting target - i.e. new figure
    fld     = fieldnames(pltData);                                  % Name of particles
    

else % If the function is called through the GUI
    src     = varargin{1};
    pltData = varargin{2};
    fld     = fieldnames(pltData);                                          % Name of particles
    
    if nargin == 2
        AX      = findobj(ancestor(src, 'figure'), 'Tag', 'PlotAx');        % Set plotting target - i.e. GUI axis       
        cla(AX, 'reset');                                                   % Clear axes
        AX.Tag = 'PlotAx';
        delete(findobj(ancestor(src, 'figure'), 'Tag', 'LegPlot'));         % Delete legend
    elseif nargin == 3
        FG      = figure;
        AX      = axes('Parent', FG);                                       % Set plotting target - i.e. new figure
    end
    
    vX      = get(findobj(ancestor(src, 'figure'), 'Tag', 'varX'), 'Value');
    vY      = get(findobj(ancestor(src, 'figure'), 'Tag', 'varY'), 'Value');
        
end



% Retrieve field names
varX    = get_field(varList{vX});
varY    = get_field(varList{vY});

% Plot
hold(AX, 'on');   grid(AX, 'on');
for i = 1:length(fld)    
    plot(AX, pltData.(fld{i}).traj.(varX)(3:end), pltData.(fld{i}).traj.(varY)(3:end), '-', 'linewidth', .75);
    hold on
end
hold(AX, 'off')
grid(AX, 'on')
box(AX, 'on')
axis(AX, 'square')

legend(AX, fld, 'Tag', 'LegPlot', 'interpreter', 'none')
xlabel(AX, varList{vX});
ylabel(AX, varList{vY});


function fldOut = get_field(fldIn)

if strcmp(fldIn, 'Time (s)')
    fldOut = 't';
elseif strcmp(fldIn, 'X distance (m)')
    fldOut = 'x';
elseif strcmp(fldIn, 'Y distance (m)')
    fldOut = 'y';
elseif strcmp(fldIn, 'Altitude (m asl)')
    fldOut = 'z';
elseif strcmp(fldIn, 'Flight distance (m)')
    fldOut = 'dis';
elseif strcmp(fldIn, 'Projected distance (m)')
    fldOut = 'disP';
elseif strcmp(fldIn, 'Latitude')
    fldOut = 'lat';
elseif strcmp(fldIn, 'Longitude')
    fldOut = 'lon';
elseif strcmp(fldIn, 'U velocity (m/s)')
    fldOut = 'u';
elseif strcmp(fldIn, 'V velocity (m/s)')
    fldOut = 'v';
elseif strcmp(fldIn, 'Z velocity (m/s)')
    fldOut = 'w';
elseif strcmp(fldIn, 'U wind (m/s)')
    fldOut = 'uf';
elseif strcmp(fldIn, 'V wind (m/s)')
    fldOut = 'vf';
elseif strcmp(fldIn, 'Re')
    fldOut = 'Re';
elseif strcmp(fldIn, 'Ganser Re')
    fldOut = 'Re_S';
elseif strcmp(fldIn, 'Drag coefficient')
    fldOut = 'Cd';
elseif strcmp(fldIn, 'Relaxation time (s)')
    fldOut = 'tau';
elseif strcmp(fldIn, 'Mach number')
    fldOut = 'Mach';
end

