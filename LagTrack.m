function LagTrack

% Setup figure
BGC     = get(0,'DefaultUicontrolBackgroundColor');
sz      = [700 1000]; % figure size
screenS = get(0,'ScreenSize');
xpos    = ceil((screenS(3)-sz(2))/2);
ypos    = ceil((screenS(4)-sz(1))/2);

% Add folders to search path
addpath(genpath('code/'));

% Check if GUI toolbox is installed
if ~exist('layoutRoot', 'file')==2
    error('The GUI Toolbox App is not installed. Please visit: <a href="https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox">https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox</a>')
end

% Check Matlab version
if verLessThan('matlab','8.4')
    error('You need at least Matlab R2014b to run the GUI, but you can still use separate functions.')
end    
% GUI
f = figure( 'Name', 'LagTrack', 'position',[xpos, ypos, sz(2), sz(1)], 'Toolbar','none', 'Menubar', 'none', 'NumberTitle', 'off', 'Visible', 'off' );
display('Preparing the interface, please wait...')

% Menu
m = uimenu('Label', 'File');
uimenu(m, 'Label', 'Load particle', 'Accelerator', 'O', 'Callback', @load_part);

m2 = uimenu('Label', 'Input parameters');
uimenu(m2, 'Label', 'Download amospheric data', 'Callback', @download_ATM);
uimenu(m2, 'Label', 'Process amospheric data', 'Callback', @process_ATM);
uimenu(m2, 'Label', 'ECMWF - Set API key', 'callback', @writeECMWFAPIKey);
uimenu(m2, 'Label', 'ECMWF - Install library', 'callback', @installECMWFAPI);
uimenu(m2, 'Label', 'Create standard atmosphere', 'Callback', @makeStandardAtm);
uimenu(m2, 'Label', 'Display atmospheric data', 'Callback', @ViewAtm);
uimenu(m2, 'Label', 'Download SRTM DEM', 'Separator', 'on', 'Callback', @download_SRTM);
uimenu(m2, 'Label', 'Process SRTM DEM', 'Separator', 'on', 'Callback', @process_SRTM);

m3 = uimenu('Label', 'Tools');
uimenu(m3, 'Label', 'Get u,v,w velocities', 'Callback', @sphere2cart);
uimenu(m3, 'Label', 'Variable parameters', 'Callback', @vary_param, 'Separator', 'on', 'enable', 'off', 'tag', 'variable');

% Main container
MAIN    = uix.VBoxFlex( 'Parent', f, 'BackgroundColor', BGC, 'Padding', 5 );

    % Top container
    top     = uix.HBox( 'Parent', MAIN , 'BackgroundColor', BGC);
        
    TOP    = uix.HBoxFlex( 'Parent', top, 'BackgroundColor', BGC, 'Padding', 5 );
    TOPL   = uix.BoxPanel( 'Parent', TOP, 'Title', 'Input', 'FontWeight', 'Bold', 'TitleColor', [.2 .2 .2], 'BackgroundColor', BGC,  'Padding', 5 );
    TOPR   = uix.BoxPanel( 'Parent', TOP, 'Title', 'Display', 'FontWeight', 'Bold', 'TitleColor', [.2 .2 .2], 'BackgroundColor', BGC,  'Padding', 5 );
    
        % Top left box
        topL = uix.VBox( 'Parent', TOPL); 
            %Pannels
            topLT = uix.TabPanel( 'Parent', topL, 'Padding', 5 , 'BackgroundColor', BGC, 'Tag', 'TopPanel' );
                topL_PROJ = uix.Panel( 'Parent', topLT );
                topL_PART = uix.Panel( 'Parent', topLT );
                topL_REL  = uix.Panel( 'Parent', topLT );
                %topL_VEL  = uix.Panel( 'Parent', topLT );
                topL_ADV  = uix.Panel( 'Parent', topLT );
                topLT.TabTitles = {'Project', 'Particle', 'Release', 'Advanced'};
                topLT.TabWidth  = 65;
                    
                    % Project
                    topL_proj = uix.Grid( 'Parent', topL_PROJ, 'Padding', 15, 'Spacing', 8, 'BackgroundColor', BGC );
                        topL_proj_nameL     = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'Run name', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_proj_latL      = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'Vent latitude', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_proj_lonL      = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'Vent longitude', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_proj_altL      = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'Vent altitude', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_proj_dateL     = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'Eruption date', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        uix.Empty( 'Parent', topL_proj );
                        topL_proj_atmL     = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'Atmospheric data', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_proj_demL     = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'DEM', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 

                        topL_proj_name      = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'Tooltip', 'Run name', 'Tag', 'name', 'callback', @check_var);
                        topL_proj_lat       = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'Tooltip', 'Vent latitude (negative in southern hemisphere)', 'Tag', 'vent_lat', 'callback', @check_var);
                        topL_proj_lon       = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'Tooltip', 'Vent longitude (negative in western hemisphere)', 'Tag', 'vent_lon', 'callback', @check_var);
                        topL_proj_alt       = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'Tooltip', 'Vent elevation (m asl)', 'Tag', 'vent_alt', 'callback', @check_var);
                        topL_proj_date      = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', datestr(now), 'Tooltip', 'Eruption date (hours in UTC)', 'Tag', 'date', 'callback', @check_var);
                        uix.Empty( 'Parent', topL_proj );
                        topL_proj_atm       = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'Tooltip', 'Path to .nc file of atmospheric data', 'Tag', 'atm', 'callback', @check_var);
                        topL_proj_dem       = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'Tooltip', 'Path to DEM file', 'Tag', 'dem', 'callback', @check_var);

                        uix.Empty( 'Parent', topL_proj );
                        uix.Empty( 'Parent', topL_proj );
                        uix.Empty( 'Parent', topL_proj );
                        uix.Empty( 'Parent', topL_proj );
                        uix.Empty( 'Parent', topL_proj );
                        uix.Empty( 'Parent', topL_proj );
                        topL_proj_atmP      = uicontrol( 'Parent', topL_proj, 'Style', 'Pushbutton', 'String', '...', 'callback', {@set_path, topL_proj_atm, '*.mat', 'input/wind/', 'Load .nc file'});
                        topL_proj_demP      = uicontrol( 'Parent', topL_proj, 'Style', 'Pushbutton', 'String', '...', 'callback', {@set_path, topL_proj_dem, '*.mat', 'input/dem/', 'Load .mat file'});

                        set( topL_proj,  'Heights', [35 35 35 35 35 0 35 35], 'Widths', [120, -1 40] );
                    
                    % Particle  
                    topLT.Selection = 2;
                    topL_part = uix.Grid( 'Parent', topL_PART, 'Padding', 15, 'Spacing', 8, 'BackgroundColor', BGC );
                        topL_part_nameL     = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'String', 'Particle name', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        uix.Empty( 'Parent', topL_part );
                        topL_part_diamL     = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'String', 'Diameter (mm)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_part_densL     = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'String', 'Density (kg/m^3)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_part_flatL     = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'String', 'Flatness', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_part_elonL     = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'String', 'Elongation', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 

                        topL_part_name      = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'Tooltip', 'Particle name', 'Tag', 'part_name', 'callback', @check_var);
                        uix.Empty( 'Parent', topL_part );
                        topL_part_diam      = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'Tooltip', 'Particle diameter (mm)', 'Tag', 'part_diam', 'String', '0.5', 'callback', @check_var);
                        topL_part_dens      = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'Tooltip', 'Density (kg/m^3)', 'Tag', 'part_dens', 'String', '1000', 'callback', @check_var);
                        topL_part_flat      = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'Tooltip', 'Flatness (i.e. ratio of small and intermediate axes)', 'Tag', 'part_flat', 'String', '0.7', 'callback', @check_var);
                        topL_part_elon      = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'Tooltip', 'Elongation (i.e. ration of intermediate and large axes)', 'Tag', 'part_elon', 'String', '0.7', 'callback', @check_var);
                        
                        set( topL_part,  'Heights', [35 0 35 35 35 35], 'Widths', [120, -1] );
                    
                    % Release    
                    topLT.Selection = 3;
                    topL_rel = uix.Grid( 'Parent', topL_REL, 'Padding', 15, 'Spacing', 8, 'BackgroundColor', BGC );                        
                        topL_rel_xL         = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'String', 'X offset (m)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);                         
                        topL_rel_yL         = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'String', 'Y offset (m)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);                         
                        topL_rel_zL         = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'String', 'Altitude (m above vent)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        uix.Empty( 'Parent', topL_rel );
                        topL_rel_tL         = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'String', 'Time offset (s)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);
                        uix.Empty( 'Parent', topL_rel );                        
                        topL_rel_vxL        = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'String', 'Initial X velocity (m/s)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);                         
                        topL_rel_vyL        = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'String', 'Initial Y velocity (m/s)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);                         
                        topL_rel_vzL        = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'String', 'Initial Z velocity (m/s)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                                               
                        topL_rel_x          = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'Tooltip', sprintf('X offset relative to the vent (m)\nPositive towards E, negative towards W'), 'Tag', 'rel_x', 'String', '0', 'callback', @check_var);                     
                        topL_rel_y          = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'Tooltip', sprintf('Y offset relative to the vent (m)\nPositive towards N, negative towards S'), 'Tag', 'rel_y', 'String', '0', 'callback', @check_var);                     
                        topL_rel_z          = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'Tooltip', 'Altitude above vent (m)', 'Tag', 'rel_z', 'String', '0', 'callback', @check_var);  
                        uix.Empty( 'Parent', topL_rel );
                        topL_rel_t          = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'Tooltip', sprintf('Time offset relative to the eruption date (sec)\nPositive in future, negative in past'), 'Tag', 'rel_t', 'String', '0', 'callback', @check_var);
                        uix.Empty( 'Parent', topL_rel );                        
                        topL_rel_vx         = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'Tooltip', sprintf('Initial velocity in X direction (m/s)\nPositive towards E, negative towards W\nEnter -1 to initialize the particle velocity with the u wind component'), 'Tag', 'rel_vx', 'String', '-1', 'callback', @check_var);                     
                        topL_rel_vy         = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'Tooltip', sprintf('Initial velocity in Y direction (m/s)\nPositive towards N, negative towards S\nEnter -1 to initialize the particle velocity with the v wind component'), 'Tag', 'rel_vy', 'String', '-1', 'callback', @check_var);                        
                        topL_rel_vz         = uicontrol( 'Parent', topL_rel, 'Style', 'Edit', 'Tooltip', sprintf('Initial velocity in Z direction (m/s)\nPositive upwards, negative downwards'), 'Tag', 'rel_vz', 'String', '0', 'callback', @check_var);
                       
                        set( topL_rel,  'Heights', [35 35 35 0 35 0 35 35 35], 'Widths', [120, -1] );
                
                    % Advanced 
                    topLT.Selection = 4;
                    topL_adv = uix.Grid( 'Parent', topL_ADV, 'Padding', 15, 'Spacing', 8, 'BackgroundColor', BGC );
                        topL_adv_solL       = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', 'Solution', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        uix.Empty( 'Parent', topL_adv );
                        topL_adv_dTL        = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', 'Time step (s)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_adv_dragL      = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', 'Reduced drag (m)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        uix.Empty( 'Parent', topL_adv );
                        topL_adv_interL     = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', 'Interpolation', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_adv_methodL    = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', 'Method', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_adv_rangeL     = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', 'Range', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_adv_skipL      = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', 'Interpolation', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        
                        
                        %topL_adv_sol        = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', 'Euler', 'Tooltip', sprintf('Solution. Enter Euler, Analytical or RungeKutta'), 'Tag', 'adv_sol', 'callback', @check_var);
                        topL_adv_sol        = uicontrol( 'Parent', topL_adv, 'Style', 'Popupmenu', 'String', {'Euler', 'Analytical'}, 'Tooltip', sprintf('Solution'), 'Tag', 'adv_sol', 'callback', @check_var);
                        uix.Empty( 'Parent', topL_adv );
                        topL_adv_dT         = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', '0.1', 'Tooltip', sprintf('Time step between advancing states (s)'), 'Tag', 'adv_dt', 'callback', @check_var);
                        topL_adv_drag       = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', '0', 'Tooltip', sprintf('Radius of reduced drag (m above vent)'), 'Tag', 'adv_drag', 'callback', @check_var); 
                        uix.Empty( 'Parent', topL_adv );
                        %topL_adv_inter      = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', 'Subset', 'Tooltip', sprintf('Interpolation of atmospheric data, either None, Subset or Complete\n-None: No interpolation\n-Subset: Only specified indices around current points are interpolated\n-Complete: Then entire atmospheric dataset is interpolated at each time step'), 'Tag', 'adv_int','callback', @check_var);
                        topL_adv_inter      = uicontrol( 'Parent', topL_adv, 'Style', 'Popupmenu', 'String', {'None', 'Subset', 'Complete'}, 'Value', 2, 'Tooltip', sprintf('Interpolation of atmospheric data, either None, Subset or Complete\n-None: No interpolation\n-Subset: Only specified indices around current points are interpolated\n-Complete: Then entire atmospheric dataset is interpolated at each time step'), 'Tag', 'adv_int','callback', @check_var);    
                        %topL_adv_method     = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', 'Linear', 'Tooltip', sprintf('Interpolation method, either linear, nearest, pchip, cubic or spline\nSee the Matlab documentation for interpn for more detail'), 'Tag', 'adv_meth', 'callback', @check_var);
                        topL_adv_method     = uicontrol( 'Parent', topL_adv, 'Style', 'Popupmenu', 'String', {'Linear', 'Nearest', 'Pchip', 'Cubic', 'Spline'}, 'Tooltip', sprintf('Interpolation method, either linear, nearest, pchip, cubic or spline\nSee the Matlab documentation for interpn for more detail'), 'Tag', 'adv_meth', 'callback', @check_var);
                        topL_adv_range      = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', '1', 'Tooltip', sprintf('Range of neighbor indexes used for subset interpolation (higher: slower interpolation)'), 'Tag', 'adv_range', 'callback', @check_var);
                        topL_adv_skip       = uicontrol( 'Parent', topL_adv, 'Style', 'Edit', 'String', '0', 'Tooltip', sprintf('Number of time steps to skip between two interpolations (higher: faster interpolation)'), 'Tag', 'adv_skip', 'callback', @check_var);
                        
                        set( topL_adv,  'Heights', [35 0 35 35 0 35 35 35 35], 'Widths', [120, -1] );
                        
                        
                    topLT.Selection = 1;
       
            topLB = uix.HBox('Parent', topL, 'Padding', 5 , 'BackgroundColor', BGC);
            errorE = uicontrol( 'Parent', topLB, 'Style', 'Edit', 'String', '  ', 'HorizontalAlign', 'left', 'Tag','Errmsg');
            uicontrol( 'Parent', topLB, 'Style', 'Pushbutton', 'String', 'Run', 'Enable', 'on', 'Tag', 'run_btn', 'callback', @runIt );
            set( topLB, 'Widths', [-1 120], 'Spacing', 5 );
        
        set(topL, 'Heights', [-1 45]);
        
    set( TOP, 'Widths', [-1 -1], 'Spacing', 5 );

%% Bottom container
    bottom  = uix.BoxPanel( 'Parent', MAIN, 'Title', 'Particles', 'FontWeight', 'Bold', 'TitleColor', [.2 .2 .2], 'BackgroundColor', BGC,  'Padding', 5 );
    
    BOT = uix.HBox( 'Parent', bottom);
    TB  = uitable( 'Parent', BOT, 'Tag', 'DataTable');
    LST = uicontrol('Parent', BOT', 'Style', 'listbox', 'Tag', 'DataList', 'Max', 100, 'Min',0);
    GRD = uix.Grid( 'Parent', BOT, 'Spacing', 5);
        uicontrol( 'Parent', GRD, 'Style', 'Pushbutton', 'String', 'Plot', 'Enable', 'off', 'Tag', 'Bplot', 'Tooltip', 'Update plot in the GUI', 'Callback', @update_plots )
        uicontrol( 'Parent', GRD, 'Style', 'Pushbutton', 'String', 'Clear', 'Enable', 'off', 'Tag', 'Bclear', 'Tooltip', 'Clear selected particles', 'Callback', @clear_part )
        uicontrol( 'Parent', GRD, 'Style', 'Pushbutton', 'String', 'Delete', 'Enable', 'off', 'Tag', 'Bdelete', 'Tooltip', 'Delete selected particles', 'Callback', @delete_part )
        uicontrol( 'Parent', GRD, 'Style', 'Pushbutton', 'String', 'Export', 'Enable', 'off', 'Tag', 'Bexport', 'Tooltip', 'Export figure to new axes', 'Callback', @update_plots )
        uicontrol( 'Parent', GRD, 'Style', 'Pushbutton', 'String', 'Details', 'Enable', 'off', 'Tag', 'Bdetail', 'Callback', @show_details )
        uicontrol( 'Parent', GRD, 'Style', 'Pushbutton', 'String', 'Close', 'Enable', 'on', 'Tag', 'Bclose', 'Callback', 'delete(gcf)' )
    set(GRD, 'Heights', [-1 -1 -1], 'Widths', [-1 -1] );
    set(BOT, 'Widths', [-4 -1 200], 'Spacing', 5 );
    
set(MAIN, 'Heights', [-1 200], 'Spacing', 5 );

% Table data
CName = {'Name',...
    '<html><center>Diameter<br />(mm)</center></html>',...
    '<html><center>Density<br />(kg/m<sup>2</sup>)</center></html>',...
    'Flatness', 'Elongation',...
    '<html><center>X offset<br />(m)</center></html>',...
    '<html><center>Y offset<br />(m)</center></html>',...
    '<html><center>Altitude<br />(m asl)</center></html>',...
    'Time', ...
    '<html><center>Fall time<br />(s)</center></html>',...
    '<html><center>Landing altitude<br />(m asl)</center></html>',...
    '<html><center>Terminal velocity<br />(m/s)</center></html>'};
CForm = {'char','bank', 'bank', 'bank', 'bank', 'bank', 'bank', 'bank', 'char', 'bank', 'bank', 'bank'};
CEdit = [false,false, false, false, false ,false ,false ,false ,false ,false ,false ,false ];

set(TB, 'ColumnName', CName, 'ColumnFormat', CForm, 'ColumnEditable', CEdit, 'RowName', [], 'ColumnWidth', 'auto', 'RowName', {});

%% DISPLAY
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
    'Relaxation time (s)'};

% Setup tabs of map/plt
topRT       = uix.TabPanel( 'Parent', TOPR, 'Padding', 5 , 'BackgroundColor', BGC, 'Tag', 'DispTab');
topR_MAP    = uix.Panel( 'Parent', topRT);  
topR_PLOT   = uix.Panel( 'Parent', topRT );
topRT.TabTitles = {'Map', 'Profiles'};
topRT.TabWidth  = 65;

% Map tab
topR_map    = uix.VBox('Parent', topR_MAP, 'BackgroundColor', BGC, 'Padding', 5, 'Tag', 'MapTab');
    ax_map  = axes('Parent', uicontainer('Parent',topR_map), 'Tag', 'MapAx', 'Box', 'on'); % Map axis

    topR_mapB  = uix.HBox( 'Parent', topR_map );
        uix.Empty('Parent', topR_mapB);
        uicontrol('Parent', topR_mapB, 'Style', 'togglebutton', 'String', '3D rotate', 'Tag', 'Map3D', 'Enable', 'off', 'Callback', {@set_map_mode, ax_map.Position});
        uicontrol('Parent', topR_mapB, 'Style', 'togglebutton', 'String', 'Pan', 'Tag', 'MapPan', 'Enable', 'off', 'Callback', {@set_map_mode, ax_map.Position});
        uicontrol('Parent', topR_mapB, 'Style', 'togglebutton', 'String', 'Zoom', 'Tag', 'MapZoom', 'Enable', 'off', 'Callback', {@set_map_mode, ax_map.Position});
        %uicontrol('Parent', topR_mapB, 'Style', 'togglebutton', 'String', 'Legend', 'Tag', 'MapLegend', 'Enable', 'off', 'Callback', {@set_map_mode, ax_map.Position});
        uix.Empty('Parent', topR_mapB);
        set(topR_mapB, 'Widths', [-1 90 90 90 -1]);
        
    set(topR_map, 'Heights', [-1 40]);
    
% Plot tab    
topR_plot   = uix.VBox('Parent', topR_PLOT, 'BackgroundColor', BGC, 'Padding', 5);
    ax_plot = axes('Parent', uicontainer('Parent',topR_plot), 'Tag','PlotAx', 'Box', 'on'); % Plot axis

    topR_plotB  = uix.HBox( 'Parent', topR_plot );
        topR_varX   = uicontrol('Parent', topR_plotB, 'Style', 'popupmenu', 'String', varList, 'Tag', 'varX', 'Enable', 'off', 'value', 3);
        topR_varY   = uicontrol('Parent', topR_plotB, 'Style', 'popupmenu', 'String', varList, 'Tag', 'varY', 'Enable', 'off', 'value', 2);
        set(topR_plotB, 'Widths', [-1 -1]);
        
    set(topR_plot, 'Heights', [-1 25]);

%% Setup GUI data
part.run_name       = -9999;
part.vent.lat       = -9999;
part.vent.lon       = -9999;
part.vent.alt       = -9999;
part.date           = now;
part.path.nc        = -9999;
part.path.dem       = -9999;
part.part.name      = -9999;
part.part.diam      = 0.5;
part.part.dens      = 1000;
part.part.flat      = 0.7;
part.part.elon      = 0.7;
part.rel.x          = 0;
part.rel.y          = 0;
part.rel.z          = 0;
part.rel.t          = 0;
part.rel.vx         = 0;
part.rel.vy         = 0;
part.rel.vz         = 1e-4;
part.adv.solution   = 'euler';
part.adv.dt         = 0.1;
part.adv.drag       = 0;
part.adv.interp     = 'subset';
part.adv.method     = 'linear';
part.adv.range      = 1;
part.adv.skip       = 0;
part.run_check      = 0;

guidata(f, part)

%% Display interface
set(f, 'Visible', 'on');

disp('Done!');

