function check_var(src, ~)

jEdit = findjobj(src);  % Retieve Java object
part  = guidata(src);   % Retrieve GUI data

%% Project pannels
if strcmp(src.Tag, 'name')
    tmp = src.String;
    err = 'Enter a valid run name';
    if isempty(tmp); change_frame(jEdit,src,0,err); part.run_name = -9999; else; change_frame(jEdit,src,1,' '); part.run_name = tmp; end
    
elseif strcmp(src.Tag, 'vent_lat')
    tmp = str2double(src.String);
    err = 'Latitude should be between +90 and -90';
    if isnan(tmp) || tmp>90 || tmp<-90; change_frame(jEdit,src,0,err); part.vent.lat = -9999; else; change_frame(jEdit,src,1,' '); part.vent.lat = tmp; end

elseif strcmp(src.Tag, 'vent_lon')  
    tmp = str2double(src.String);
    err = 'Check the vent longitude';
    if isnan(tmp); change_frame(jEdit,src,0,err); part.vent.lon = -9999; else; change_frame(jEdit,src,1,' '); part.vent.lon = tmp; end

elseif strcmp(src.Tag, 'vent_alt')
    tmp = str2double(src.String);
    err = 'Check the vent altitude';
    if isnan(tmp) || tmp<0; change_frame(jEdit,src,0,err); part.vent.alt = -9999;  else; change_frame(jEdit,src,1,' '); part.vent.alt = tmp;  end
    
elseif strcmp(src.Tag, 'date')
    tmp = src.String;
    err = 'Enter a valid date (see datestr in Matlab help)';
    if ~isempty(regexp(tmp,'\d\d-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-\d\d\d\d \d\d:\d\d:\d\d', 'once')) ; change_frame(jEdit,src,1,' '); part.date = datenum(tmp); else; change_frame(jEdit,src,0,err); part.date = -9999; end
    
elseif strcmp(src.Tag, 'atm')   
    tmp = src.String;
    err = 'The specified NecCDF file does not exist';
    if exist(tmp, 'file') ; change_frame(jEdit,src,1,' '); part.path.nc = tmp; else; change_frame(jEdit,src,0,err); part.path.nc = -9999; end
    
elseif strcmp(src.Tag, 'dem')
    tmp = src.String;
    err = 'The specified DEM file does not exist';
    if exist(tmp, 'file') ; change_frame(jEdit,src,1,' '); part.path.dem = tmp; else; change_frame(jEdit,src,0,err); part.path.dem = -9999; end

    
%% Part pannels    
elseif strcmp(src.Tag, 'part_name')
    tmp = src.String;
    err = 'Enter a particle name';
    err1 = 'Invalid particle name';
    if isempty(tmp) ; change_frame(jEdit,src,0,err); part.part.name = -9999; elseif isnumeric(str2double(tmp(1))) && ~isnan(str2double(tmp(1))); change_frame(jEdit,src,0,err1); part.part.name = -9999;  else; change_frame(jEdit,src,1,' '); part.part.name = tmp; end

elseif strcmp(src.Tag, 'part_diam')
    tmp = str2double(src.String);
    err = 'Check the particle diameter';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); part.part.diam = -9999; else; change_frame(jEdit,src,1,' '); part.part.diam = tmp/1e3; end % Convert diameter to m

elseif strcmp(src.Tag, 'part_dens')
    tmp = str2double(src.String);
    err = 'Check the particle density';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); part.part.dens = -9999; else; change_frame(jEdit,src,1,' '); part.part.dens = tmp; end

elseif strcmp(src.Tag, 'part_flat')
    tmp = str2double(src.String);
    err = 'Check the particle flatness';
    if isnan(tmp) || tmp<0 || tmp>1; change_frame(jEdit,src,0,err); part.part.flat = -9999; else; change_frame(jEdit,src,1,' '); part.part.flat = tmp; end

elseif strcmp(src.Tag, 'part_elon')
    tmp = str2double(src.String);
    err = 'Check the particle elongation';
    if isnan(tmp) || tmp<0 || tmp>1; change_frame(jEdit,src,0,err); part.part.elon = -9999; else; change_frame(jEdit,src,1,' '); part.part.elon = tmp; end
    

%% Release pannels 
elseif strcmp(src.Tag, 'rel_x')
    tmp = str2double(src.String);
    err = 'Check the x offset';
    if isnan(tmp); change_frame(jEdit,src,0,err); part.rel.x = -9999; else; change_frame(jEdit,src,1,' '); part.rel.x = tmp; end    
    
elseif strcmp(src.Tag, 'rel_y')
    tmp = str2double(src.String);
    err = 'Check the y offset';
    if isnan(tmp); change_frame(jEdit,src,0,err); part.rel.y = -9999; else; change_frame(jEdit,src,1,' '); part.rel.y = tmp; end 
    
elseif strcmp(src.Tag, 'rel_z')
    tmp = str2double(src.String);
    err = 'Check the release altitude';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); part.rel.z = -9999; else; change_frame(jEdit,src,1,' '); part.rel.z = tmp; end    
    
elseif strcmp(src.Tag, 'rel_t')
    tmp = str2double(src.String);
    err = 'Check the time offset';
    if isnan(tmp); change_frame(jEdit,src,0,err); part.rel.t = -9999; else; change_frame(jEdit,src,1,' '); part.rel.t = tmp/(3600*24); end   
    
elseif strcmp(src.Tag, 'rel_vx')
    tmp = str2double(src.String);
    err = 'Check the x velocity';
    if isnan(tmp); change_frame(jEdit,src,0,err); part.rel.vx = -9999; else; change_frame(jEdit,src,1,' '); part.rel.vx = tmp; end    
    
elseif strcmp(src.Tag, 'rel_vy')
    tmp = str2double(src.String);
    err = 'Check the y velocity';
    if isnan(tmp); change_frame(jEdit,src,0,err); part.rel.vy = -9999; else; change_frame(jEdit,src,1,' '); part.rel.vy = tmp; end 
    
elseif strcmp(src.Tag, 'rel_vz')
    tmp = str2double(src.String);
    err = 'Check the z velocity';
    if isnan(tmp); change_frame(jEdit,src,0,err); part.rel.vz = -9999; else; change_frame(jEdit,src,1,' '); if tmp == 0; tmp = 1e-4; end; part.rel.vz = tmp; end % Check if Vz = 0
    
    
%% Advanced pannel
elseif strcmp(src.Tag, 'adv_sol')
    tmp = lower(src.String{src.Value});
    err = 'Solution is Euler, Analytical or RungeKutta';
    if isempty(regexp(tmp, '(euler|analytical|rungekutta)', 'once')); change_frame(jEdit,src,0,err); part.adv.solution = -9999; else; change_frame(jEdit,src,1,' '); part.adv.solution = tmp; end 
    
elseif strcmp(src.Tag, 'adv_dt')
    tmp = str2double(src.String);
    err = 'Check the time step';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); part.adv.dt = -9999; else; change_frame(jEdit,src,1,' '); part.adv.dt = tmp; end    
    
elseif strcmp(src.Tag, 'adv_drag')
    tmp = str2double(src.String);
    err = 'Check the radius of reduced drag';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); part.adv.drag = -9999; else; change_frame(jEdit,src,1,' '); part.adv.drag = tmp; end    
    
elseif strcmp(src.Tag, 'adv_int')
    tmp = lower(src.String{src.Value});
    err = 'Interpolation is None, Subset or Complete';
    if isempty(regexp(tmp, '(none|subset|complete)', 'once')); change_frame(jEdit,src,0,err); part.adv.interp = -9999; else; change_frame(jEdit,src,1,' '); part.adv.interp = tmp; end 

elseif strcmp(src.Tag, 'adv_meth')
    tmp = lower(src.String{src.Value});
    err = 'Check interpolation method (interpn)';
    if isempty(regexp(tmp, '(linear|nearest|pchip|cubic|spline)', 'once')); change_frame(jEdit,src,0,err); part.adm.method = -9999; else; change_frame(jEdit,src,1,' '); part.adv.method = tmp; end 
    
elseif strcmp(src.Tag, 'adv_range')
    tmp = str2double(src.String);
    err = 'Check the interpolation range';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); part.adv.range = -9999; else; change_frame(jEdit,src,1,' '); part.adv.range = tmp; end        
    
elseif strcmp(src.Tag, 'adv_skip')
    tmp = str2double(src.String);
    err = 'Check the time skip interval';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); part.adv.skip = -9999; else; change_frame(jEdit,src,1,' '); part.adv.skip = tmp; end    
    
end

%% Update GUI
% Check if using standard atmosphere AND standard grid
if ~isempty(regexp(get(findobj(ancestor(src, 'figure'), 'Tag', 'dem'), 'String'),'_STD.mat', 'once')) && ~isempty(regexp(get(findobj(ancestor(src, 'figure'), 'Tag', 'atm'), 'String'),'_STD.mat', 'once'))
    set(findobj(ancestor(src, 'figure'), 'Tag', 'vent_lat'), 'Enable', 'off'); 
    set(findobj(ancestor(src, 'figure'), 'Tag', 'vent_lon'), 'Enable', 'off');
    set(findobj(ancestor(src, 'figure'), 'Tag', 'vent_alt'), 'Enable', 'off');
    gui_check = 1;
else
    set(findobj(ancestor(src, 'figure'), 'Tag', 'vent_lat'), 'Enable', 'on'); 
    set(findobj(ancestor(src, 'figure'), 'Tag', 'vent_lon'), 'Enable', 'on');
    set(findobj(ancestor(src, 'figure'), 'Tag', 'vent_alt'), 'Enable', 'on');
    gui_check = 0;
end

% Check to enable Run button
if ischar(part.run_name) && ...
        (gui_check == 1 || (part.vent.lat ~= -9999 && part.vent.lon ~= -9999 && part.vent.alt~= -9999)) && ...
        part.date ~= -9999 && ...
        ischar(part.path.nc) && ischar(part.path.dem) && ...
        ischar(part.part.name) && part.part.diam ~= -9999 && part.part.dens ~= -9999 && part.part.flat ~= -9999 && part.part.elon ~= -9999 && ...
        part.rel.x ~= -9999 && part.rel.y ~= -9999 && part.rel.z ~= -9999 && part.rel.t ~= -9999 && part.rel.vx ~= -9999 && part.rel.vy ~= -9999 && part.rel.vz ~= -9999 &&   ...
        ischar(part.adv.solution) && ischar(part.adv.interp) && ischar(part.adv.method) && ...
        part.adv.dt ~= -9999 && part.adv.drag ~= -9999 && part.adv.range ~= -9999 && part.adv.skip ~= -9999
    
    %set(findobj(ancestor(src, 'figure'), 'Tag', 'run_btn'), 'Enable', 'on');
    part.run_check = 1;
else
    %set(findobj(ancestor(src, 'figure'), 'Tag', 'run_btn'), 'Enable', 'off');
    part.run_check = 0;
end

% Update GUI data
guidata(src, part);

function change_frame(jEdit,src, typ,errmsg)
if typ == 0
    jEdit.Border = javax.swing.border.LineBorder(java.awt.Color(1,0,0),1,false);
    else
    jEdit.Border = javax.swing.border.LineBorder(java.awt.Color(.65,.65,.65),1,false);
end
set(findobj(ancestor(src, 'figure'), 'Tag', 'Errmsg'), 'String', errmsg);
