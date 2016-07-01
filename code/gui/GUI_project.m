function GUI_project
global t

% Define GUI
scr = get(0,'ScreenSize');
w   = 300;
h   = 350;

BGC = get(0,'DefaultUicontrolBackgroundColor');

% Main figure
t.fig = figure(...
    'position', [scr(3)/2-w/2 scr(4)/2-h/2 w h],...
    'Color', BGC,...
    'Resize', 'off',...
    'Toolbar', 'none',...
    'Menubar', 'none',...
    'Name', 'LagTrack',...
    'NumberTitle', 'off');
            
    t.main = uipanel(...
        'parent', t.fig,...
        'units', 'normalized',...
        'position', [.04 .025 .92 .95],...
        'BackgroundColor', BGC,...
        'ForegroundColor', [.0 .0 .4],...
        'HighlightColor', [.0 .0 .4],...
        'BorderType', 'line');
    
                t.name = uicontrol(...
                'parent', t.main,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.025 .84 .95 .135],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Project name',...
                'ToolTipString', 'Name given to main project');
    
                t.atm = uicontrol(...
                'parent', t.main,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.025 .68 .75 .135],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Input atmospheric data',...
                'ToolTipString', 'Full path to the NetCDF file');
            
            
                t.atmP = uicontrol(...
                    'parent', t.main,...
                    'Style', 'pushbutton',...
                    'units', 'normalized',...
                    'position', [.8 .68 .175 .135],...
                    'BackgroundColor', BGC,...
                    'ForegroundColor', [.2 .2 .2],...
                    'String', '...',...
                    'ToolTipString', 'Browse',...
                    'Callback', {@browse_atm});
    
    
                t.dem = uicontrol(...
                'parent', t.main,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.025 .52 .75 .135],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Input DEM',...
                'ToolTipString', 'Full path to SRTM data');
            
            
                t.demP = uicontrol(...
                    'parent', t.main,...
                    'Style', 'pushbutton',...
                    'units', 'normalized',...
                    'position', [.8 .52 .175 .135],...
                    'BackgroundColor', BGC,...
                    'ForegroundColor', [.2 .2 .2],...
                    'String', '...',...
                    'ToolTipString', 'Browse',...
                    'Tag', 'atm',...
                    'Callback', {@browse_dem});
        
                %% Vent
            
                t.lat = uicontrol(...
                'parent', t.main,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.025 .36 .3 .135],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Vent lat.',...
                'ToolTipString', sprintf('Vent latitude (degrees) \nPositive in Northern hemisphere\nNegative in Southern hemisphere'));
    
                t.lon = uicontrol(...
                'parent', t.main,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.35 .36 .3 .135],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Vent lon.',...
                'ToolTipString', sprintf('Vent longitude (degrees)\nPositive in Eastern hemisphere\nNegative in Western hemisphere'));
    
                t.alt = uicontrol(...
                'parent', t.main,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.675 .36 .3 .135],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Vent alt.',...
                'ToolTipString', sprintf('Vent altitude (m asl)'));
            
                %%
        
            
                t.date = uicontrol(...
                'parent', t.main,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.025 .2 .95 .135],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Eruption date',...
                'ToolTipString', 'Eruption date. Search datestr in the Matlab help for more information.');
            
            
                t.loadP = uicontrol(...
                    'parent', t.main,...
                    'Style', 'pushbutton',...
                    'units', 'normalized',...
                    'position', [.025 .025 .45 .15],...
                    'BackgroundColor', BGC,...
                    'ForegroundColor', [.2 .2 .2],...
                    'String', 'Load',...
                    'Callback', {@load_project});
            
            
                t.okP = uicontrol(...
                    'parent', t.main,...
                    'Style', 'pushbutton',...
                    'units', 'normalized',...
                    'position', [.525 .025 .45 .15],...
                    'BackgroundColor', BGC,...
                    'ForegroundColor', [.2 .2 .2],...
                    'String', 'Ok',...
                    'Callback', {@set_project});
                
function browse_atm(~,~)
global t
[fl,pt] = uigetfile('input/wind/*.nc');
set(t.atm, 'String', fullfile(pt,fl));

function browse_dem(~,~)
global t
[fl,pt] = uigetfile('input/dem/*.mat');
set(t.dem, 'String', fullfile(pt,fl));

function set_project(~,~)
global t

% Check if folder already exist
if exist(['projects/', get(t.name, 'String')], 'dir') == 7
    choice = questdlg('A folder with the same name already exists. Do you want to delete it?', ...
	'Project Name', ...
	'No', 'Yes', 'Yes');
    % Handle response
    switch choice
        case 'Yes'
            rmdir(['projects/', get(t.name, 'String')], 's');
        case 'No'
            return
    end
end

% Make output folder
mkdir(['projects/', get(t.name, 'String')])

project.name        = get(t.name, 'String'); 
project.atm         = get(t.atm, 'String'); 
project.dem         = get(t.dem, 'String'); 
project.vent_lat    = get(t.lat, 'String'); 
project.vent_lon    = get(t.lon, 'String'); 
project.vent_alt    = get(t.alt, 'String'); 
project.date   = get(t.date, 'String'); 

save(['projects/', get(t.name, 'String'), '/', get(t.name, 'String'), '.mat'], 'project')

function load_project(~,~)
global t

[fl,pt] = uigetfile('projects/*.mat');
load(fullfile(pt,fl));

set(t.name, 'String',project.name); 
set(t.atm, 'String',project.atm); 
set(t.dem, 'String',project.dem); 
set(t.lat, 'String',project.vent_lat); 
set(t.lon, 'String',project.vent_lon); 
set(t.alt, 'String',project.vent_alt); 
set(t.date, 'String',project.date); 