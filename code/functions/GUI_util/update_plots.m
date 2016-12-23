% Update plots

function update_plots(src,~)

APDTA   = getappdata(ancestor(src, 'figure'));
pltData = APDTA.pltData;

% Get table data from GUI
List    = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'String');
ListV   = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'Value');

if isempty(ListV)
    errordlg('Select at least one particle to plot');
    return
end

for i = 1:length(ListV)
    toPlot.(List{ListV(i)}) = pltData.(List{ListV(i)});
end

if strcmp(src.String, 'Plot')       % Plot in GUI
    if get(findobj(ancestor(src, 'figure'), 'Tag', 'DispTab'), 'Selection') == 1
        set(findobj(ancestor(src, 'figure'), 'Tag', 'Errmsg'), 'String', 'Retrieving basemap, please wait...');
        map_part(src, toPlot);
        set(findobj(ancestor(src, 'figure'), 'Tag', 'Errmsg'), 'String', '');
    else
        plot_part(src, toPlot);
    end
elseif strcmp(src.String, 'Export') % Plot in new figure
    if get(findobj(ancestor(src, 'figure'), 'Tag', 'DispTab'), 'Selection') == 1
        set(findobj(ancestor(src, 'figure'), 'Tag', 'Errmsg'), 'String', 'Retrieving basemap, please wait...');
        map_part(src, toPlot, 1);
        set(findobj(ancestor(src, 'figure'), 'Tag', 'Errmsg'), 'String', '');
    else
        plot_part(src, toPlot, 1);
    end
end
    
