BGC = get(0,'DefaultUicontrolBackgroundColor');
PC  = [.9 .9 .9];
sz = [700 1000]; % figure size
screensize = get(0,'ScreenSize');
xpos = ceil((screensize(3)-sz(2))/2); % center the figure on the
ypos = ceil((screensize(4)-sz(1))/2); % center the figure on the



f = figure( 'Name', 'LagTrack', 'position',[xpos, ypos, sz(2), sz(1)] );

% Main container
main    = uix.VBoxFlex( 'Parent', f, 'BackgroundColor', BGC );

    % Top container
    top     = uix.HBox( 'Parent', main , 'BackgroundColor', BGC);
        
        % Top left box
        topL = uix.BoxPanel( 'Parent', top, 'Title', 'Input', 'FontWeight', 'Bold', 'TitleColor', [.2 .2 .2], 'BackgroundColor', BGC,  'Padding', 5 );
        topL2= uix.VBox( 'Parent', topL);

            topLT = uix.TabPanel( 'Parent', topL2, 'Padding', 5 , 'BackgroundColor', BGC );
            
                topL_proj = uix.Grid( 'Parent', topLT, 'Padding', 15, 'Spacing', 15, 'BackgroundColor', PC );
                    topL_proj_nameL     = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'Run name', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left'); remove_frame(topL_proj_nameL  );
                    topL_proj_latL      = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'Latitude', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left'); remove_frame(topL_proj_latL  );
                    topL_proj_lonL      = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'Longitude', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left'); remove_frame(topL_proj_lonL );
                    topL_proj_dateL     = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'Eruption date', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left'); remove_frame(topL_proj_dateL );
                    uix.Empty( 'Parent', topL_proj )
                    topL_proj_atmL     = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'Atmospheric data', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left'); remove_frame(topL_proj_atmL );
                    topL_proj_demL     = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', 'DEM', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left'); remove_frame(topL_proj_demL );

                    topL_proj_name      = uicontrol( 'Parent', topL_proj, 'Style', 'Edit');
                    topL_proj_lat       = uicontrol( 'Parent', topL_proj, 'Style', 'Edit');
                    topL_proj_lon       = uicontrol( 'Parent', topL_proj, 'Style', 'Edit');
                    topL_proj_date      = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', datestr(now));
                    uix.Empty( 'Parent', topL_proj )
                    topL_proj_atm       = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', datestr(now));
                    topL_proj_dem       = uicontrol( 'Parent', topL_proj, 'Style', 'Edit', 'String', datestr(now));

                    uix.Empty( 'Parent', topL_proj );
                    uix.Empty( 'Parent', topL_proj );
                    uix.Empty( 'Parent', topL_proj );
                    uix.Empty( 'Parent', topL_proj );
                    uix.Empty( 'Parent', topL_proj );
                    topL_proj_atmP      = uicontrol( 'Parent', topL_proj, 'Style', 'Pushbutton', 'String', '...');
                    topL_proj_demP      = uicontrol( 'Parent', topL_proj, 'Style', 'Pushbutton', 'String', '...');
                    
                    set( topL_proj,  'Heights', [35 35 35 35 0 35 35], 'Widths', [-2, -4 -1] );
                
            topL_erup = uicontrol( 'Parent', topLT );
            topL_part = uicontrol( 'Parent', topLT );
            topL_pos  = uicontrol( 'Parent', topLT );
            topL_vel  = uicontrol( 'Parent', topLT );
            topL_opr  = uicontrol( 'Parent', topLT );
            topLT.TabTitles = {'Project', 'Eruption', 'Particle', 'Position', 'Velocity','Option'};
            
            topLB = uix.HButtonBox( 'Parent', topL2, 'Padding', 5 , 'BackgroundColor', BGC   );
            uicontrol( 'Parent', topLB, 'String', 'Run' );
            set( topLB, 'ButtonSize', [100 45], 'Spacing', 5 );
        
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