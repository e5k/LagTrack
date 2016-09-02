function ViewAtm(varargin)

% Retrieve path to atmospheric file
if nargin == 0
    [fl, pth] = uigetfile('*.mat', 'Select atmospheric file');
    PTH = fullfile(pth,fl);
    if fl == 0; return; end
else
    src = varargin{1};
    if isempty(get(findobj(ancestor(src, 'figure'), 'Tag', 'atm'), 'String'))        
        [fl, pth] = uigetfile('*.mat', 'Select atmospheric file');
        PTH = fullfile(pth,fl);
        if fl == 0; return; end
    else
        PTH = get(findobj(ancestor(src, 'figure'), 'Tag', 'atm'), 'String');
    end
end

% Load atmospheric file
atm     = load(PTH); atm = atm.atm;
varList = {'Wind velocity', 'U wind', 'V wind', 'Temperature', 'Relative humidity', 'Air density', 'Air viscosity'};
levList = cellstr(num2str(atm.level));
timList = cellstr(datestr(atm.time));

% Set initial time
if nargin == 0
    tIstart = 1;
else
    if isempty(datenum(get(findobj(ancestor(src, 'figure'), 'Tag', 'time'))))
        tIstart = 1;
    else
        [~,tIstart] = min(abs(atm.time - datenum(get(findobj(ancestor(src, 'figure'), 'Tag', 'time'), 'String'))));
    end
end

% Setup figure
BGC = get(0,'DefaultUicontrolBackgroundColor');
sz = [600 600]; % figure size
screensize = get(0,'ScreenSize');
xpos = ceil((screensize(3)-sz(2))/2); % center the figure on the
ypos = ceil((screensize(4)-sz(1))/2); % center the figure on the

f = figure( 'Name', 'Display atmospheric data', 'position',[xpos, ypos, sz(2), sz(1)], 'Toolbar','none', 'Menubar', 'none', 'NumberTitle', 'off', 'Visible', 'on' );

% Main V box
main = uix.VBox( 'Parent', f , 'Padding', 10, 'BackgroundColor', BGC);

        axes('Parent', uicontainer('Parent', main), 'Tag', 'AtmAX');
        %xlabel('Longitude'); ylabel('Latitude');
        %quiver(AX, atm.lon, atm.lat, squeeze(atm.u(:,:,1,tIstart)), squeeze(atm.v(:,:,1,tIstart)))
        
        % Secondary bottom container
        bot = uix.VBox( 'Parent', main, 'Padding', 5, 'BackgroundColor', BGC );
            % Top/bottom
            botT = uix.HBox( 'Parent', bot, 'Padding', 5, 'BackgroundColor', BGC );
                uicontrol('Parent', botT, 'Style', 'popupmenu', 'String', varList, 'Tag', 'varList', 'Tooltip', 'Parameter to plot', 'Callback', {@PlotAtm, atm});
                uicontrol('Parent', botT, 'Style', 'popupmenu', 'String', levList, 'Tag', 'levList', 'Tooltip', 'Geopotential height (mb)', 'Callback', {@PlotAtm, atm}); 
                uicontrol('Parent', botT, 'Style', 'popupmenu', 'String', timList, 'Tag', 'timList', 'Value', tIstart, 'Tooltip', 'Time', 'Callback', {@PlotAtm, atm});
            set(botT, 'Widths', [-1 -1 -1],'Spacing', 5 );
            
            % Midle/bottom
            botM = uix.HBox( 'Parent', bot, 'Padding', 5, 'BackgroundColor', BGC );
                uicontrol( 'Parent', botM, 'Style', 'Pushbutton', 'String', '<<', 'Callback', {@changeT, atm});
                uicontrol( 'Parent', botM, 'Style', 'Pushbutton', 'String', '<', 'Callback', {@changeT, atm});
                uicontrol( 'Parent', botM, 'Style', 'Edit', 'String', 'Run name', 'BackgroundColor', BGC,  'HorizontalAlign', 'Center', 'Tag', 'DateInfo', 'Enable', 'off'); 
                uicontrol( 'Parent', botM, 'Style', 'Pushbutton', 'String', '>', 'Callback', {@changeT, atm});
                uicontrol( 'Parent', botM, 'Style', 'Pushbutton', 'String', '>>', 'Callback', {@changeT, atm});
            set(botM, 'Widths', [50 50 -1 50 50],'Spacing', 5 );
        
            % Bottom/bottom
            botB = uix.HBox( 'Parent', bot, 'Padding', 5, 'BackgroundColor', BGC );
                uicontrol( 'Parent', botB, 'Style', 'Edit', 'String', 'Run name', 'BackgroundColor', BGC,  'HorizontalAlign', 'Left', 'Tag', 'info'); 
                uicontrol( 'Parent', botB, 'Style', 'Pushbutton', 'String', 'Close', 'Callback', 'delete(gcf)');%, 'Enable', 'off', 'Tag', 'Bsave' )

            set(botB, 'Widths', [-1 105],'Spacing', 5 );            
        set(bot, 'Heights', [-1 -1.5 -1.5] );        
set( main, 'Heights', [-1 130 ], 'Spacing', 5 );

PlotAtm(f,'',atm);


function changeT(src, ~, atm)

% Time
timL = get(findobj(ancestor(src, 'figure'), 'Tag', 'timList'), 'String');
timS = get(findobj(ancestor(src, 'figure'), 'Tag', 'timList'), 'Value');

if strcmp(src.String, '<')
    step = -1;
elseif strcmp(src.String, '<<')
    step = -10;
elseif strcmp(src.String, '>>')
    step = 10;
elseif strcmp(src.String, '>')
    step = 1;
end

if (timS+step) < 1
    idx = 1;
elseif (step+timS) > length(timL)
    idx = length(timL);
else
    idx = timS+step;
end

set(findobj(ancestor(src, 'figure'), 'Tag', 'timList'), 'Value', idx);
PlotAtm(src,'',atm);