function check_var(src, evt)
jEdit = findjobj(src);

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
    if isnan(tmp) || tmp>180 || tmp<-180 ; change_frame(jEdit,src,0,err); else change_frame(jEdit,1); part.vent.lon = tmp; end
    
elseif strcmp(src.Tag, 'vent_alt')
    tmp = str2double(src.String);
    err = 'Check the vent altitude';
    if isnan(tmp) || tmp<0 ; change_frame(jEdit,src,0,err); else change_frame(jEdit,1); part.vent.alt = tmp; end
    
elseif strcmp(src.Tag, 'date')
    tmp = src.String;
    err = 'Enter a valid date (see datestr in Matlab help)';
    if isdatetime(tmp) ; change_frame(jEdit,src,1,' '); part.date = tmp; else change_frame(jEdit,src,0,err); end
    
elseif strcmp(src.Tag, 'atm')   
    tmp = src.String;
    err = 'The specified NecCDF file does not exist';
    if exist(tmp, 'file') ; change_frame(jEdit,src,1,' '); part.path.nc = tmp; else change_frame(jEdit,src,0,err); end
    
elseif strcmp(src.Tag, 'dem')
    tmp = src.String;
    err = 'The specified DEM file does not exist';
    if exist(tmp, 'file') ; change_frame(jEdit,src,1,' '); part.path.dem = tmp; else change_frame(jEdit,src,0,err); end
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