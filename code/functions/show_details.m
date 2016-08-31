function show_details(varargin)

if nargin == 0  % If calling the function independently
    [fl, pth] = uigetfile('*.mat');
    if isempty(fl)
        return
    end
    
    load(fullfile(pth, fl));
    if ~isfield(part, 'traj')
        errordlg('The selected file is not a valid run')
        return
    end
    
    
else % If the function is called through the GUI
    src     = varargin{1};
    tmpS    = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'String');
    tmpV    = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'Value');
    APDTA   = getappdata(ancestor(src, 'figure'));  
    
    if length(tmpV) > 1
        errordlg('Please select only one particle')
        return
    elseif isempty(tmpV)
        errordlg('Please select one particle');
        return
    else
        part = APDTA.pltData.(char(tmpS(tmpV)));
    end
end

header = {'<html><center> Time <br/> (s) </center></html>',...
    '<html><center> X distance <br/> (m) </center></html>',...
    '<html><center> Y distance <br/> (m) </center></html>',...
    '<html><center> Altitude <br/> (m asl) </center></html>',...
    '<html><center> Distance <br/> (m) </center></html>',...
    '<html><center> Latitude <br/> (Deg) </center></html>',...
    '<html><center> Longitude <br/> (Deg) </center></html>',...
    '<html><center> U velocity <br/> (m/s) </center></html>',...
    '<html><center> V velocity <br/> (m/s) </center></html>',...
    '<html><center> Z velocity <br/> (m/s) </center></html>',...
    '<html><center> U wind <br/> (m/s) </center></html>',...
    '<html><center> V wind <br/> (m/s) </center></html>',...
    '<html><center> Re <br/>  </center></html>',...
    '<html><center> Ganser Re <br/>  </center></html>',...
    '<html><center> Drag coef. <br/>  </center></html>',...
    '<html><center> Relaxation time <br/> (s) </center></html>',...    
    };


data = [part.traj.t(1:end-1)',...
    part.traj.x(1:end-1)',...
    part.traj.y(1:end-1)',...
    part.traj.z(1:end-1)',...
    part.traj.dis(1:end-1)',...
    part.traj.lat(1:end-1)',...
    part.traj.lon(1:end-1)',...
    part.traj.u(1:end-1)',...
    part.traj.v(1:end-1)',...
    part.traj.w(1:end-1)',...
    part.traj.uf(1:end)',...
    part.traj.vf(1:end)',...
    part.traj.Re(1:end-1)',...
    part.traj.Re_S(1:end-1)',...
    part.traj.Cd(1:end-1)',...
    part.traj.tau(1:end-1)',...
    ];


CForm = repmat({'bank'}, 1, 16);
CEdit = false(1, 16);

%% GUI

%BGC         = get(0,'DefaultUicontrolBackgroundColor');
%PC          = [.9 .9 .9];
sz          = [700 1000]; % figure size
screensize  = get(0,'ScreenSize');
xpos        = ceil((screensize(3)-sz(2))/2); % center the figure on the
ypos        = ceil((screensize(4)-sz(1))/2); % center the figure on the

f           = figure( 'Name', [part.run_name, ' - ', part.part.name], 'position',[xpos, ypos, sz(2), sz(1)], 'Toolbar','none', 'Menubar', 'none', 'NumberTitle', 'off');

main        = uix.VBox('Parent', f, 'Padding', 15, 'Spacing', 5);

    TB  = uitable( 'Parent', main, 'Tag', 'DataTable');
    BOT = uix.HBox('Parent', main, 'Spacing', 5); 
        
        uix.Empty( 'Parent', BOT );
        uicontrol( 'Parent', BOT, 'Style', 'Pushbutton', 'String', 'Export', 'callback', {@export_part, part, data});
        uicontrol( 'Parent', BOT, 'Style', 'Pushbutton', 'String', 'Close'); %, 'callback', {@set_path, topL_proj_atm, '*.mat', 'input/wind/', 'Load .nc file'});

    set(BOT, 'Widths', [-1 100 100]);
set(main, 'Heights', [-1 50]);

set(TB, 'Data', data, 'ColumnName', header, 'ColumnFormat', CForm, 'ColumnEditable', CEdit, 'ColumnWidth', 'auto');


function export_part(~, ~, part, data)

%[fl, dr] = uiputfile({'*.txt', 'Tab-delimited text file (*.txt)'}, 'Save as')

dr = uigetdir('.', 'Choose target folder');

if isempty(dr)
    return
end

% Details
fl1 = [dr, filesep, part.run_name, '_', part.part.name, '_input.txt'];

fid1 = fopen(fl1, 'w');
fprintf(fid1, 'Run name:\t%s\n', part.run_name);
fprintf(fid1, 'Particle name:\t%s\n\n', part.part.name);

fprintf(fid1, 'Atmpspheric data:\t%s\n', part.path.nc);
fprintf(fid1, 'DEM:\t%s\n\n', part.path.dem);

fprintf(fid1, 'Vent latitude:\t%2.2f\n', part.vent.lat);
fprintf(fid1, 'Vent longitude:\t%2.2f\n', part.vent.lon);
fprintf(fid1, 'Vent altitude:\t%2.2f\n', part.vent.alt);
fprintf(fid1, 'Eruption date:\t%s\n\n', datestr(part.date));

fprintf(fid1, 'Diameter (mm):\t%2.4f\n', part.part.diam*1e2);
fprintf(fid1, 'Density (kg/m3):\t%.2f\n', part.part.dens);
fprintf(fid1, 'Flatness:\t%.2f\n', part.part.flat);
fprintf(fid1, 'Elongation:\t%.2f\n\n', part.part.elon);

fprintf(fid1, 'X offset (m):\t%4.0f\n', part.rel.x);
fprintf(fid1, 'Y offset (m):\t%4.0f\n', part.rel.y);
fprintf(fid1, 'Z offset (m):\t%4.0f\n', part.rel.z);
fprintf(fid1, 'Time offset (s):\t%4.0f\n', part.rel.t);
fprintf(fid1, 'Initial U velocity (m/s):\t%4.2f\n', part.rel.vx);
fprintf(fid1, 'Initial V velocity (m/s):\t%4.2f\n', part.rel.vy);
fprintf(fid1, 'Initial Z velocity (m/s):\t%4.2f\n\n', part.rel.vz);

fprintf(fid1, 'Solution:\t%s\n', part.adv.solution);

fprintf(fid1, 'Time step:\t%.2f\n\n', part.adv.dt);
fprintf(fid1, 'Radius of reduced drag:\t%.2f\n\n', part.adv.drag);

fprintf(fid1, 'Interpolation scheme:\t%s\n', part.adv.interp);
fprintf(fid1, 'Interpolation method:\t%s\n', part.adv.method);
fprintf(fid1, 'Interpolation range:\t%i\n', part.adv.range);
fprintf(fid1, 'Interpolation skip:\t%i\n', part.adv.skip);



fclose(fid1);

% Trajectory
fl2 = [dr, filesep, part.run_name, '_', part.part.name, '_trajectory.txt'];
fid2 = fopen(fl2, 'w');

fprintf(fid2, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n',...
    'Iteration',...
    'Time (s)',...
    'X distance (m)',...
    'Y distance (m)',...
    'Altitude (m asl)',...
    'Distance (m)',...
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
    'Relaxation time (s)');

for i = 1:size(data,1)
    fprintf(fid2, '%i\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t\n',...
        i, data(i,1), data(i,2), data(i,3), data(i,4), data(i,5), data(i,6), data(i,7), data(i,8), data(i,9), data(i,10), data(i,11), data(i,12), data(i,13), data(i,14), data(i,15), data(i,16));
end
fclose(fid2);