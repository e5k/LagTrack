function vary_param(varargin)

% Retrieve data
part = guidata(ancestor(varargin{1}, 'Figure'));

% Setup figure
BGC     = get(0,'DefaultUicontrolBackgroundColor');
sz      = [600 500]; % figure size
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
f    = figure( 'Name', 'LagTrack', 'position',[xpos, ypos, sz(2), sz(1)], 'Toolbar','none', 'Menubar', 'none', 'NumberTitle', 'off' );
MAIN = uix.VBoxFlex( 'Parent', f, 'BackgroundColor', BGC, 'Padding', 10 );

main = uix.VBox( 'Parent', MAIN );
    TOP_txt = 'This tool is design to run particles with initial conditions defined as ranges instead of single values. For each parameter, define a minimum and a maximum value and the number of intervals in the min-max range.';
    TOP = uix.Panel('Parent', main, 'BackgroundColor', BGC, 'Padding', 5);
        uicontrol( 'Parent', TOP, 'Style', 'text', 'String', TOP_txt, 'BackgroundColor', BGC, 'HorizontalAlign', 'Left'); 
    ROW = uix.Panel('Parent', main, 'BackgroundColor', BGC, 'Padding', 5);
        row = uix.VBox( 'Parent', ROW );
    BOT = uix.Panel('Parent', main, 'BackgroundColor', BGC, 'Padding', 5);
        bot = uix.HButtonBox( 'Parent', BOT );
set( main, 'Heights', [75 -1 75], 'Spacing', 5 );


% Columns
c1 = uix.HBox( 'Parent', row);
c2 = uix.HBox( 'Parent', row);
c3 = uix.HBox( 'Parent', row);
c4 = uix.HBox( 'Parent', row);
c5 = uix.HBox( 'Parent', row);
c6 = uix.HBox( 'Parent', row);
c7 = uix.HBox( 'Parent', row);
c8 = uix.HBox( 'Parent', row);
c9 = uix.HBox( 'Parent', row);
c10 = uix.HBox( 'Parent', row);
c11 = uix.HBox( 'Parent', row);

set( row, 'Heights', [-1 -1 -1 -1 -1 -0.5 -1 -1 -1 -1 -0.5], 'Spacing', 5 );

    % Row 1
    uix.Empty( 'Parent', c1);
    uicontrol( 'Parent', c1, 'Style', 'Edit', 'String', 'Minimum', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);
    uicontrol( 'Parent', c1, 'Style', 'Edit', 'String', 'Maximum', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
    uicontrol( 'Parent', c1, 'Style', 'Edit', 'String', 'Steps', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame); 
    set(c1, 'Widths', [-2 -1 -1 -1], 'spacing', 5);
    
    % Row 2
    uicontrol( 'Parent', c2, 'Style', 'Edit', 'String', 'Diameter (mm)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);
    uicontrol( 'Parent', c2, 'Style', 'Edit', 'String', num2str(part.part.diam*1e3),  'HorizontalAlign', 'Left', 'tag', 'diamMin');
    uicontrol( 'Parent', c2, 'Style', 'Edit', 'String', num2str(part.part.diam*1e3),  'HorizontalAlign', 'Left', 'tag', 'diamMax');
    uicontrol( 'Parent', c2, 'Style', 'Edit', 'String', num2str(1),  'HorizontalAlign', 'Left', 'tag', 'diamInt');
    set(c2, 'Widths', [-2 -1 -1 -1], 'spacing', 5);
    
    % Row 3
    uicontrol( 'Parent', c3, 'Style', 'Edit', 'String', 'Density (kg/m3)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);
    uicontrol( 'Parent', c3, 'Style', 'Edit', 'String', num2str(part.part.dens), 'HorizontalAlign', 'Left', 'tag', 'densMin');
    uicontrol( 'Parent', c3, 'Style', 'Edit', 'String', num2str(part.part.dens), 'HorizontalAlign', 'Left', 'tag', 'densMax');
    uicontrol( 'Parent', c3, 'Style', 'Edit', 'String', num2str(1), 'HorizontalAlign', 'Left', 'tag', 'densInt');
    set(c3, 'Widths', [-2 -1 -1 -1], 'spacing', 5);
    
    % Row 4
    uicontrol( 'Parent', c4, 'Style', 'Edit', 'String', 'Flatness', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);
    uicontrol( 'Parent', c4, 'Style', 'Edit', 'String', num2str(part.part.flat),  'HorizontalAlign', 'Left', 'tag', 'flatMin');
    uicontrol( 'Parent', c4, 'Style', 'Edit', 'String', num2str(part.part.flat),  'HorizontalAlign', 'Left', 'tag', 'flatMax');
    uicontrol( 'Parent', c4, 'Style', 'Edit', 'String', num2str(1),  'HorizontalAlign', 'Left', 'tag', 'flatInt');
    set(c4, 'Widths', [-2 -1 -1 -1], 'spacing', 5);
    
    % Row 5
    uicontrol( 'Parent', c5, 'Style', 'Edit', 'String', 'Elongation', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);
    uicontrol( 'Parent', c5, 'Style', 'Edit', 'String', num2str(part.part.elon),  'HorizontalAlign', 'Left', 'tag', 'elonMin');
    uicontrol( 'Parent', c5, 'Style', 'Edit', 'String', num2str(part.part.elon),  'HorizontalAlign', 'Left', 'tag', 'elonMax');
    uicontrol( 'Parent', c5, 'Style', 'Edit', 'String', num2str(1),  'HorizontalAlign', 'Left', 'tag', 'elonInt');
    set(c5, 'Widths', [-2 -1 -1 -1], 'spacing', 5);
    
    % Row 6
    uix.Empty( 'Parent', c6);
    uix.Empty( 'Parent', c6);
    uix.Empty( 'Parent', c6);
    uix.Empty( 'Parent', c6);
    set(c6, 'Widths', [-2 -1 -1 -1], 'spacing', 5);
    
    % Row 7
    uicontrol( 'Parent', c7, 'Style', 'Edit', 'String', 'X offset (m)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);
    uicontrol( 'Parent', c7, 'Style', 'Edit', 'String', num2str(part.rel.x), 'HorizontalAlign', 'Left', 'tag', 'xoffMin');
    uicontrol( 'Parent', c7, 'Style', 'Edit', 'String', num2str(part.rel.x), 'HorizontalAlign', 'Left', 'tag', 'xoffMax');
    uicontrol( 'Parent', c7, 'Style', 'Edit', 'String', num2str(1), 'HorizontalAlign', 'Left', 'tag', 'xoffInt');
    set(c7, 'Widths', [-2 -1 -1 -1], 'spacing', 5);
    
    % Row 8
    uicontrol( 'Parent', c8, 'Style', 'Edit', 'String', 'Y offset (m)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);
    uicontrol( 'Parent', c8, 'Style', 'Edit', 'String', num2str(part.rel.y), 'HorizontalAlign', 'Left', 'tag', 'yoffMin');
    uicontrol( 'Parent', c8, 'Style', 'Edit', 'String', num2str(part.rel.y), 'HorizontalAlign', 'Left', 'tag', 'yoffMax');
    uicontrol( 'Parent', c8, 'Style', 'Edit', 'String', num2str(1), 'HorizontalAlign', 'Left', 'tag', 'yoffInt');
    set(c8, 'Widths', [-2 -1 -1 -1], 'spacing', 5);
    
    % Row 9
    uicontrol( 'Parent', c9, 'Style', 'Edit', 'String', 'Altitude (m above vent)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);
    uicontrol( 'Parent', c9, 'Style', 'Edit', 'String', num2str(part.rel.z), 'HorizontalAlign', 'Left', 'tag', 'zoffMin');
    uicontrol( 'Parent', c9, 'Style', 'Edit', 'String', num2str(part.rel.z), 'HorizontalAlign', 'Left', 'tag', 'zoffMax');
    uicontrol( 'Parent', c9, 'Style', 'Edit', 'String', num2str(1), 'HorizontalAlign', 'Left', 'tag', 'zoffInt');
    set(c9, 'Widths', [-2 -1 -1 -1], 'spacing', 5);
    
    % Row 10
    uicontrol( 'Parent', c10, 'Style', 'Edit', 'String', 'Time offset (s)', 'BackgroundColor', BGC, 'Enable', 'off', 'HorizontalAlign', 'Left', 'CreateFcn', @remove_frame);
    uicontrol( 'Parent', c10, 'Style', 'Edit', 'String', num2str(part.rel.t),  'HorizontalAlign', 'Left', 'tag', 'toffMin');
    uicontrol( 'Parent', c10, 'Style', 'Edit', 'String', num2str(part.rel.t),  'HorizontalAlign', 'Left', 'tag', 'toffMax');
    uicontrol( 'Parent', c10, 'Style', 'Edit', 'String', num2str(1),  'HorizontalAlign', 'Left', 'tag', 'toffInt');
    set(c10, 'Widths', [-2 -1 -1 -1], 'spacing', 5);
    
    % Row 11
    uix.Empty( 'Parent', c11);
    uix.Empty( 'Parent', c11);
    uix.Empty( 'Parent', c11);
    uix.Empty( 'Parent', c11);
    set(c11, 'Widths', [-2 -1 -1 -1], 'spacing', 5);

% Bottom
uix.Empty( 'Parent', bot);
uix.Empty( 'Parent', bot);
uicontrol( 'Parent', bot, 'String', 'Run', 'callback', {@set_part, part, varargin{1}} );
set( bot, 'ButtonSize', [130 35], 'Spacing', 5 );


function set_part(src, ~, part, SRC)

diamVec = linspace( str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'diamMin'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'diamMax'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'diamInt'), 'String')))./1e3;
densVec = linspace( str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'densMin'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'densMax'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'densInt'), 'String')));
flatVec = linspace( str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'flatMin'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'flatMax'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'flatInt'), 'String')));
elonVec = linspace( str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'elonMin'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'elonMax'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'elonInt'), 'String')));

xoffVec = linspace( str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'xoffMin'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'xoffMax'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'xoffInt'), 'String')));
yoffVec = linspace( str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'yoffMin'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'yoffMax'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'yoffInt'), 'String')));
zoffVec = linspace( str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'zoffMin'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'zoffMax'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'zoffInt'), 'String')));
toffVec = linspace( str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'toffMin'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'toffMax'), 'String')), str2double(get(findobj(ancestor(src, 'Figure'), 'tag', 'toffInt'), 'String')));

nPart   = length(diamVec)*length(densVec)*length(flatVec)*length(elonVec)*length(xoffVec)*length(yoffVec)*length(zoffVec)*length(toffVec); 
PART    = cell(nPart,1);

countP  = 1;
for idiam = 1:length(diamVec)
    for idens = 1:length(densVec)
        for iflat = 1:length(flatVec)
            for ielon = 1:length(elonVec)
                for ixoff = 1:length(xoffVec)
                    for iyoff = 1:length(yoffVec)
                        for izoff = 1:length(zoffVec)
                            for itoff = 1:length(toffVec)
                                PART{countP} = part;
                                PART{countP}.part.name = [PART{countP}.part.name, '_', num2str(countP)];
                                PART{countP}.part.diam = diamVec(idiam);
                                PART{countP}.part.dens = densVec(idens);
                                PART{countP}.part.flat = flatVec(iflat);
                                PART{countP}.part.elon = elonVec(ielon);
                                PART{countP}.rel.x     = xoffVec(ixoff);
                                PART{countP}.rel.y     = yoffVec(iyoff);
                                PART{countP}.rel.z     = zoffVec(izoff);
                                PART{countP}.rel.t     = toffVec(itoff);
                                countP = countP + 1;
                                zoffVec(izoff);
                            end
                        end
                    end
                end
            end
        end
    end
end
delete(ancestor(src, 'Figure'));
runIt(SRC,' ',PART)