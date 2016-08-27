% Update plots

function update_plots(src,~)

APDTA   = getappdata(ancestor(src, 'figure'));
pltData = APDTA.pltData;

if strcmp(src.String, 'Plot')       % Plot in GUI
    if get(findobj(ancestor(src, 'figure'), 'Tag', 'DispTab'), 'Selection') == 1
        set(findobj(ancestor(src, 'figure'), 'Tag', 'Errmsg'), 'String', 'Retrieving basemap, please wait...');
        map_part(src, pltData);
        set(findobj(ancestor(src, 'figure'), 'Tag', 'Errmsg'), 'String', '');
    else
        plot_part(src, pltData);
    end
elseif strcmp(src.String, 'Export') % Plot in new figure
    if get(findobj(ancestor(src, 'figure'), 'Tag', 'DispTab'), 'Selection') == 1
        set(findobj(ancestor(src, 'figure'), 'Tag', 'Errmsg'), 'String', 'Retrieving basemap, please wait...');
        map_part(src, pltData, 1);
        set(findobj(ancestor(src, 'figure'), 'Tag', 'Errmsg'), 'String', '');
    else
        plot_part(src, pltData, 1);
    end
end
    
