% Update plots

function update_plots(src,~)

APDTA   = getappdata(ancestor(src, 'figure'));
pltData = APDTA.pltData;

if strcmp(src.String, 'Plot')       % Plot in GUI
    map_part(src, pltData);
    plot_part(src, pltData);
elseif strcmp(src.String, 'Export') % Plot in new figure
    map_part(src, pltData, 1);
    plot_part(src, pltData, 1);
end
    