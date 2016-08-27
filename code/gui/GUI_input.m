function GUI_input
global t

% Define GUI
scr = get(0,'ScreenSize');
w   = 400;
h   = 500;

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

t.menu0 = uimenu(t.fig, 'Label', 'ECMWF');
        t.m01 = uimenu(t.menu0, 'Label', 'Set ECMWF API key', 'callback', 'writeECMWFAPIKey');
        t.m02 = uimenu(t.menu0, 'Label', 'Install ECMWF libraries', 'callback', 'installECMWFAPI');

% t.menu1 = uimenu(t.fig, 'Label', 'ECMWF');
%         t.m11 = uimenu(t.menu0, 'Label', 'Set ECMWF API key', 'callback', 'writeECMWFAPIKey');
%         t.m12 = uimenu(t.menu0, 'Label', 'Install ECMWF libraries', 'callback', 'installECMWFAPI');
            
        t.name = uipanel(...
            'parent', t.fig,...
            'units', 'normalized',...
            'position', [.04 .8 .92 .18],...
            'title', 'Name',...
            'BackgroundColor', BGC,...
            'ForegroundColor', [.0 .0 .4],...
            'HighlightColor', [.0 .0 .4],...
            'BorderType', 'line');   
            
            t.nameE = uicontrol(...
                'parent', t.name,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.2 .3 .6 .4],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Input name',...
                'ToolTipString', 'Name given to the input datasets');


        t.extent = uipanel(...
            'parent', t.fig,...
            'units', 'normalized',...
            'position', [.04 .5 .92 .28],...
            'title', 'Extent',...
            'BackgroundColor', BGC,...
            'ForegroundColor', [.0 .0 .4],...
            'HighlightColor', [.0 .0 .4],...
            'BorderType', 'line');  
        
            
            t.extentN = uicontrol(...
                'parent', t.extent,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.3 .725 .4 .24],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Latitude North',...
                'ToolTipString', 'Positive in Northen hemisphere, negative in Southern');
        
            
            t.extentS = uicontrol(...
                'parent', t.extent,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.3 .075 .4 .24],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Latitude South',...
                'ToolTipString', 'Positive in Northen hemisphere, negative in Southern');
        
            
            t.extentW = uicontrol(...
                'parent', t.extent,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.05 .41 .4 .24],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Longitude West',...
                'ToolTipString', 'Positive in Eastern hemisphere, negative in Western');
        
            
            t.extentE = uicontrol(...
                'parent', t.extent,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.55 .41 .4 .24],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Longitude East',...
                'ToolTipString', 'Positive in Eastern hemisphere, negative in Western');
        
        
        
        t.wind = uipanel(...
            'parent', t.fig,...
            'units', 'normalized',...
            'position', [.04 .215 .92 .265],...
            'title', 'Era-Interim',...
            'BackgroundColor', BGC,...
            'ForegroundColor', [.0 .0 .4],...
            'HighlightColor', [.0 .0 .4],...
            'BorderType', 'line');  
        
            t.windMS = uicontrol(...
                'parent', t.wind,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.05 .55 .225 .275],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Start month',...
                'ToolTipString', 'Example: 04 for April');
        
            t.windME = uicontrol(...
                'parent', t.wind,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.325 .55 .225 .275],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'End month',...
                'ToolTipString', 'To download only one month, enter Start month = End month');
        
            t.windYS = uicontrol(...
                'parent', t.wind,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.05 .175 .225 .275],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Start year',...
                'ToolTipString', 'Example: 2003');
        
            t.windYE = uicontrol(...
                'parent', t.wind,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.325 .175 .225 .275],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'End year',...
                'ToolTipString', 'Example: 2003');
            
            t.windD = uicontrol(...
                'parent', t.wind,...
                'Style', 'pushbutton',...
                'units', 'normalized',...
                'position', [.65 .1 .25 .35],...
                'BackgroundColor', BGC,...
                'ForegroundColor', [.2 .2 .2],...
                'String', 'ECMWF',...
                'Tag', 'Interim',...
                'Tooltip', 'ECMWF ERA-Interim',...
                'Callback', {@download});
            
            t.windD2 = uicontrol(...
                'parent', t.wind,...
                'Style', 'pushbutton',...
                'units', 'normalized',...
                'position', [.65 .55 .25 .35],...
                'BackgroundColor', BGC,...
                'ForegroundColor', [.2 .2 .2],...
                'String', 'NOAA',...
                'Tag', 'Reanalysis',...
                'Tooltip', 'NOAA NCEP-DOE Reanalysis 2',...
                'Callback', {@download});
        
        
        t.dem = uipanel(...
            'parent', t.fig,...
            'units', 'normalized',...
            'position', [.04 .025 .92 .17],...
            'title', 'SRTM',...
            'BackgroundColor', BGC,...
            'ForegroundColor', [.0 .0 .4],...
            'HighlightColor', [.0 .0 .4],...
            'BorderType', 'line');  
        
        
            
            t.demR = uicontrol(...
                'parent', t.dem,...
                'style', 'edit',...
                'unit', 'normalized',...
                'position', [.05 .275 .5 .45],...
                'HorizontalAlignment', 'center',...
                'ForegroundColor', [.2 .2 .2],...
                'BackgroundColor', BGC,...
                'String', 'Grid resolution (m)',...
                'ToolTipString', sprintf('Grid resolution after linear interpolation.\nThe original resolution of the SRTM dataset is 90 m'));
            
            t.demD = uicontrol(...
                'parent', t.dem,...
                'Style', 'pushbutton',...
                'units', 'normalized',...
                'position', [.65 .175 .25 .65],...
                'BackgroundColor', BGC,...
                'ForegroundColor', [.2 .2 .2],...
                'String', 'Download',...
                'Tag', 'SRTM',...
                'Callback', {@download});
            
    function download(source,~)
    global t    
    
    if strcmp(get(source, 'Tag'), 'SRTM')
        downloadSRTM(str2double(get(t.extentS, 'String')), str2double(get(t.extentN, 'String')), str2double(get(t.extentW, 'String')), str2double(get(t.extentE, 'String')), str2double(get(t.demR, 'String')), get(t.nameE, 'String'));
    elseif strcmp(get(source, 'Tag'), 'Interim')
        dwind_ECMWF(str2double(get(t.extentS, 'String')), str2double(get(t.extentN, 'String')), str2double(get(t.extentW, 'String')), str2double(get(t.extentE, 'String')), str2double(get(t.windYS, 'String')), str2double(get(t.windYE, 'String')), str2double(get(t.windMS, 'String')), str2double(get(t.windME, 'String')), get(t.nameE, 'String'), 'Interim')       
    elseif strcmp(get(source, 'Tag'), 'NOAA')
        dwind_ECMWF(str2double(get(t.extentS, 'String')), str2double(get(t.extentN, 'String')), str2double(get(t.extentW, 'String')), str2double(get(t.extentE, 'String')), str2double(get(t.windYS, 'String')), str2double(get(t.windYE, 'String')), str2double(get(t.windMS, 'String')), str2double(get(t.windME, 'String')), get(t.nameE, 'String'), 'Interim')
    end