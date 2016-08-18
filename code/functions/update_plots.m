% Update plots

function update_plots(src,~)

APDTA   = getappdata(ancestor(src, 'figure'));
pltData = APDTA.pltData;

map_part(src, pltData);
plot_part(src, pltData);