% Disable or enable UI

function enableUI(src, state)

set(findobj(ancestor(src, 'figure'), 'Tag', 'Bplot'), 'Enable', state);
set(findobj(ancestor(src, 'figure'), 'Tag', 'Bclear'), 'Enable', state);
set(findobj(ancestor(src, 'figure'), 'Tag', 'Bexport'), 'Enable', state);
set(findobj(ancestor(src, 'figure'), 'Tag', 'Bdetail'), 'Enable', state);
set(findobj(ancestor(src, 'figure'), 'Tag', 'Bdelete'), 'Enable', state);

set(findobj(ancestor(src, 'figure'), 'Tag', 'varX'), 'Enable', state);
set(findobj(ancestor(src, 'figure'), 'Tag', 'varY'), 'Enable', state);

set(findobj(ancestor(src, 'figure'), 'Tag', 'Map3D'), 'Enable', state);
set(findobj(ancestor(src, 'figure'), 'Tag', 'MapPan'), 'Enable', state);
set(findobj(ancestor(src, 'figure'), 'Tag', 'MapZoom'), 'Enable', state);
set(findobj(ancestor(src, 'figure'), 'Tag', 'MapLegend'), 'Enable', state);