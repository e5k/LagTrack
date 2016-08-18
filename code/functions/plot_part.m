function plot_part(src, pltData)

AX      = findobj(ancestor(src, 'figure'), 'Tag', 'PlotAx');     % Set plotting target - i.e. GUI axis
fld     = fieldnames(pltData);                                  % Retrieve fieldnames

delete(AX.Children);
