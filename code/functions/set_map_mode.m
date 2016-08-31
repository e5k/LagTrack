% Changes mode of the 3D map
function set_map_mode(src, ~, position)

if src.Value == 1 && strcmp(src.Tag, 'Map3D')
    set(findobj(ancestor(src, 'figure'), 'Tag', 'MapZoom'), 'Value',0);
    set(findobj(ancestor(src, 'figure'), 'Tag', 'MapPan'), 'Value',0);
    rotate3d on
elseif src.Value == 0 && strcmp(src.Tag, 'Map3D')
    rotate3d off
elseif src.Value == 1 && strcmp(src.Tag, 'MapPan')    
    set(findobj(ancestor(src, 'figure'), 'Tag', 'Map3D'), 'Value',0);
    set(findobj(ancestor(src, 'figure'), 'Tag', 'MapZoom'), 'Value',0);
    pan on
elseif src.Value == 0 && strcmp(src.Tag, 'MapPan')   
    pan off
elseif src.Value == 1 && strcmp(src.Tag, 'MapZoom')
    set(findobj(ancestor(src, 'figure'), 'Tag', 'MapPan'), 'Value',0);
    set(findobj(ancestor(src, 'figure'), 'Tag', 'Map3D'), 'Value',0);
    zoom on
elseif src.Value == 0 && strcmp(src.Tag, 'MapZoom')   
    zoom off        
elseif strcmp(src.Tag, 'MapLegend')
    AX      = findobj(ancestor(src, 'figure'), 'Tag', 'MapAx'); 
    legend(AX, 'toggle')
end

function btn(src,~)
a