function check_var(src, ~)
jEdit = findjobj(src);

%% Retrieve GUI data
part = guidata(src);

%% Project pannels
if strcmp(src.Tag, 'name')
    tmp = src.String;
    part.run_name = tmp;
elseif strcmp(src.Tag, 'vent_lat')
    tmp = str2double(src.String);
    err = 'Check the vent latitude';
    if isnan(tmp) || tmp>90 || tmp<-90 ; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.vent.lat = tmp; end
    
elseif strcmp(src.Tag, 'vent_lon')  
    tmp = str2double(src.String);
    err = 'Check the vent longitude';
    if isnan(tmp) || tmp>180 || tmp<-180 ; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.vent.lon = tmp; end
    
elseif strcmp(src.Tag, 'vent_alt')
    tmp = str2double(src.String);
    err = 'Check the vent altitude';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.vent.alt = tmp; end
    
elseif strcmp(src.Tag, 'date')
    tmp = src.String;
    err = 'Enter a valid date (see datestr in Matlab help)';
    if ~isempty(regexp(tmp,'\d\d-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)-\d\d\d\d \d\d:\d\d:\d\d', 'once')) ; change_frame(jEdit,src,1,' '); part.date = tmp; else change_frame(jEdit,src,0,err); end
    
elseif strcmp(src.Tag, 'atm')   
    tmp = src.String;
    err = 'The specified NecCDF file does not exist';
    if exist(tmp, 'file') ; change_frame(jEdit,src,1,' '); part.path.nc = tmp; else change_frame(jEdit,src,0,err); end
    
elseif strcmp(src.Tag, 'dem')
    tmp = src.String;
    err = 'The specified DEM file does not exist';
    if exist(tmp, 'file') ; change_frame(jEdit,src,1,' '); part.path.dem = tmp; else change_frame(jEdit,src,0,err); end

    
%% Part pannels    
elseif strcmp(src.Tag, 'part_name')
    tmp = src.String;
    err = 'Enter a particle name';
    if isempty(tmp) ; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.part.diam = tmp; end

elseif strcmp(src.Tag, 'part_diam')
    tmp = str2double(src.String);
    err = 'Check the particle diameter';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.part.diam = tmp; end

elseif strcmp(src.Tag, 'part_dens')
    tmp = str2double(src.String);
    err = 'Check the particle density';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.part.dens = tmp; end

elseif strcmp(src.Tag, 'part_flat')
    tmp = str2double(src.String);
    err = 'Check the particle flatness';
    if isnan(tmp) || tmp<0 || tmp>1; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.part.flat = tmp; end

elseif strcmp(src.Tag, 'part_elon')
    tmp = str2double(src.String);
    err = 'Check the particle elongation';
    if isnan(tmp) || tmp<0 || tmp>1; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.part.elon = tmp; end
    

%% Release pannels 
elseif strcmp(src.Tag, 'rel_x')
    tmp = str2double(src.String);
    err = 'Check the x offset';
    if isnan(tmp); change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.rel.x = tmp; end    
    
elseif strcmp(src.Tag, 'rel_y')
    tmp = str2double(src.String);
    err = 'Check the y offset';
    if isnan(tmp); change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.rel.y = tmp; end 
    
elseif strcmp(src.Tag, 'rel_z')
    tmp = str2double(src.String);
    err = 'Check the release altitude';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.rel.z = tmp; end    
    
elseif strcmp(src.Tag, 'rel_t')
    tmp = str2double(src.String);
    err = 'Check the time offset';
    if isnan(tmp); change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.rel.t = tmp; end    
    
elseif strcmp(src.Tag, 'rel_vx')
    tmp = str2double(src.String);
    err = 'Check the x velocity';
    if isnan(tmp); change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.rel.vx = tmp; end    
    
elseif strcmp(src.Tag, 'rel_vy')
    tmp = str2double(src.String);
    err = 'Check the y velocity';
    if isnan(tmp); change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.rel.vy = tmp; end 
    
elseif strcmp(src.Tag, 'rel_vz')
    tmp = str2double(src.String);
    err = 'Check the z velocity';
    if isnan(tmp); change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); if tmp == 0; tmp = 1e-4; end; part.rel.vz = tmp; end % Check if Vz = 0
    
    
%% Advanced pannel
elseif strcmp(src.Tag, 'adv_sol')
    tmp = lower(src.String);
    err = 'Solution is Euler, Analytical or RungeKutta';
    if isempty(regexp(tmp, '(euler|analytical|rungekutta)', 'once')); change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.adv.solution = tmp; end 
    
elseif strcmp(src.Tag, 'adv_dt')
    tmp = str2double(src.String);
    err = 'Check the time step';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.adv.dt = tmp; end    
    
elseif strcmp(src.Tag, 'adv_drag')
    tmp = str2double(src.String);
    err = 'Check the radius of reduced drag';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.adv_drag = tmp; end    
    
elseif strcmp(src.Tag, 'adv_int')
    tmp = lower(src.String);
    err = 'Interpolation is None, Subset or Complete';
    if isempty(regexp(tmp, '(none|subset|complete)', 'once')); change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.adv.interp = tmp; end 

elseif strcmp(src.Tag, 'adv_meth')
    tmp = lower(src.String);
    err = 'Check interpolation method (interpn)';
    if isempty(regexp(tmp, '(linear|nearest|pchip|cubic|spline)', 'once')); change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.adv.method = tmp; end 
    
elseif strcmp(src.Tag, 'adv_range')
    tmp = str2double(src.String);
    err = 'Check the interpolation range';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.adv.range = tmp; end        
    
elseif strcmp(src.Tag, 'adv_skip')
    tmp = str2double(src.String);
    err = 'Check the time skip interval';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); else change_frame(jEdit,src,1,' '); part.adv_skip = tmp; end    
    
end





%% Update GUI data
guidata(src, part);



%%



function change_frame(jEdit,src, typ,errmsg)
if typ == 0
    jEdit.Border = javax.swing.border.LineBorder(java.awt.Color(1,0,0),1,false);
    else
    jEdit.Border = javax.swing.border.LineBorder(java.awt.Color(.65,.65,.65),1,false);
end
set(findobj(ancestor(src, 'figure'), 'Tag', 'Errmsg'), 'String', errmsg);



% 
% if strcmp(typ, 'fl')
%     
% elseif strcmp(typ, 'num')
%     
% elseif strcmp(typ, 'str')
%     go = 0;
% end
% 
% 
% if go == 0
% 
%     jEditbox.Border = javax.swing.border.LineBorder(java.awt.Color(1,0,0),1,false);
% else
%     jEditbox.Border = javax.swing.border.LineBorder(java.awt.Color(.171,.173,.179),1,false);
% end