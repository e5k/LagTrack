BGC = get(0,'DefaultUicontrolBackgroundColor');
PC  = [.9 .9 .9];
sz = [700 1000]; % figure size
screensize = get(0,'ScreenSize');
xpos = ceil((screensize(3)-sz(2))/2); % center the figure on the
ypos = ceil((screensize(4)-sz(1))/2); % center the figure on the



f = figure( 'Name', 'LagTrack', 'position',[xpos, ypos, sz(2), sz(1)], 'Toolbar','none' );

% Main container
main    = uix.VBoxFlex( 'Parent', f, 'BackgroundColor', BGC, 'Padding', 5 );

    % Top container
    top     = uix.HBox( 'Parent', main , 'BackgroundColor', BGC);
        
        % Top left box
        topL = uix.BoxPanel( 'Parent', top, 'Title', 'Input', 'FontWeight', 'Bold', 'TitleColor', [.2 .2 .2], 'BackgroundColor', BGC,  'Padding', 5 );
        topL2= uix.VBox( 'Parent', topL);
            
            %Pannels
            topLT = uix.TabPanel( 'Parent', topL2, 'Padding', 5 , 'BackgroundColor', BGC );
                topL_PROJ = uix.Panel( 'Parent', topLT );
                topL_PART = uix.Panel( 'Parent', topLT );
                topL_POS  = uix.Panel( 'Parent', topLT );
                topL_VEL  = uix.Panel( 'Parent', topLT );
                topL_OPT  = uix.Panel( 'Parent', topLT );
                topLT.TabTitles = {'Project', 'Particle', 'Position', 'Velocity','Option'};
                
                
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
                        topL_proj_atmP      = uicontrol( 'Parent', topL_proj, 'Style', 'Pushbutton', 'String', '...', 'callback', {@set_path, topL_proj_atm, '*.nc', 'input/wind/', 'Load .nc file'});
                        topL_proj_demP      = uicontrol( 'Parent', topL_proj, 'Style', 'Pushbutton', 'String', '...', 'callback', {@set_path, topL_proj_dem, '*.mat', 'input/dem/', 'Load .mat file'});

                        set( topL_proj,  'Heights', [35 35 35 35 35 0 35 35], 'Widths', [120, -1 40] );
                    
                    topLT.Selection = 2;
                    topL_part = uix.Grid( 'Parent', topL_PART, 'Padding', 15, 'Spacing', 8, 'BackgroundColor', BGC );
                        topL_part_nameL     = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'String', 'Particle name', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        uix.Empty( 'Parent', topL_part );
                        topL_part_diamL     = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'String', 'Diameter (mm)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_part_densL     = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'String', 'Density (kg/m^3)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_part_flatL     = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'String', 'Flatness', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
                        topL_part_elonL     = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'String', 'Elongation', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 

                        topL_part_name      = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'Tooltip', 'Particle name', 'Tag', 'name', 'callback', @check_var);
                        uix.Empty( 'Parent', topL_part );
                        topL_part_diam      = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'Tooltip', 'Particle diameter (mm)', 'Tag', 'name', 'callback', @check_var);
                        topL_part_dens      = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'Tooltip', 'Density (kg/m^3)', 'Tag', 'name', 'callback', @check_var);
                        topL_part_flat      = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'Tooltip', 'Flatness (i.e. ratio of small and intermediate axes)', 'Tag', 'name', 'callback', @check_var);
                        topL_part_elon      = uicontrol( 'Parent', topL_part, 'Style', 'Edit', 'Tooltip', 'Elongation (i.e. ration of intermediate and large axes)', 'Tag', 'name', 'callback', @check_var);
                        
                        set( topL_part,  'Heights', [35 0 35 35 35 35], 'Widths', [120, -1
                            
                        ] );
                        
                        
                
                
                    topL_pos = uix.Grid( 'Parent', topL_POS, 'Padding', 15, 'Spacing', 8, 'BackgroundColor', BGC );
                
                
                    topL_vel = uix.Grid( 'Parent', topL_VEL, 'Padding', 15, 'Spacing', 8, 'BackgroundColor', BGC );
                
                
                    topL_opt = uix.Grid( 'Parent', topL_OPT, 'Padding', 15, 'Spacing', 8, 'BackgroundColor', BGC );
                
                
                

                
            
                
                
            
            topLB = uix.HBox('Parent', topL2, 'Padding', 5 , 'BackgroundColor', BGC);
            %uix.HButtonBox( 'Parent', topL2, 'Padding', 5 , 'BackgroundColor', BGC   );
            errorE = uicontrol( 'Parent', topLB, 'Style', 'Edit', 'String', '  ', 'HorizontalAlign', 'left', 'Tag','Errmsg');
            uicontrol( 'Parent', topLB, 'Style', 'Pushbutton', 'String', 'Run' );
            set( topLB, 'Widths', [-1 120], 'Spacing', 5 );
            %set( topLB, 'ButtonSize', [100 45], 'Spacing', 5 );
        
        set(topL2, 'Heights', [-1 45]);
        
        % Top right box
        topR = uix.BoxPanel( 'Parent', top, 'Title', 'Display', 'FontWeight', 'Bold', 'TitleColor', [.2 .2 .2], 'BackgroundColor', BGC,  'Padding', 5 );

    set( top, 'Widths', [-1 -1], 'Spacing', 5 );

    % Bottom container
    bottom  = uix.BoxPanel( 'Parent', main, 'Title', 'Input', 'FontWeight', 'Bold', 'TitleColor', [.2 .2 .2], 'BackgroundColor', BGC,  'Padding', 5 );

set(main, 'Heights', [-1 200], 'Spacing', 5 );










% jEditbox = findjobj(a);
% jEditbox.setBorder([]); % or: set(jEditbox,'Border',[])
% 




%topL = uix.TabPanel( 'Parent', top, 'Padding', 5 , 'BackgroundColor', BGC  );
% topL_proj = uicontrol( 'Parent', topL );
% topL_erup = uicontrol( 'Parent', topL );
% topL_part = uicontrol( 'Parent', topL );
% topL.TabTitles = {'Project', 'Blue', 'Green'};

%uicontrol( 'Parent', top, 'Background', 'r' )
%uicontrol( 'Parent', top, 'Background', 'b' )