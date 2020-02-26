% Adapts some of the GUI's parameters as a function of the run mode (i.e.
% forward vs backward

function check_run_mode(src, ~)

% If forward
if src.Value == 1
    x_str = 'Vent latitude';
    y_str = 'Vent longitude';
    z_str = 'Vent altitude';
    
    x_ttp = sprintf('Vent latitude (negative in southern hemisphere)');
    y_ttp = sprintf('Vent longitude (negative in western hemisphere)');
    z_ttp = sprintf('Vent elevation (m asl)');
else
    x_str = 'Landing latitude';
    y_str = 'Landing longitude';
    z_str = 'Maximum altitude';
    
    x_ttp = sprintf('Latitude of impact point (negative in southern hemisphere)');
    y_ttp = sprintf('Longitude of impact point (negative in western hemisphere)');
    z_ttp = sprintf('Maximum altitude reached by the particle (m asl)');
end

% Update labels
x_tmpL = findobj(ancestor(src, 'figure'), 'Tag', 'vent_latL');
x_tmpL.String = x_str;
remove_frame(x_tmpL)
y_tmpL = findobj(ancestor(src, 'figure'), 'Tag', 'vent_lonL');
y_tmpL.String = y_str;
remove_frame(y_tmpL)
z_tmpL = findobj(ancestor(src, 'figure'), 'Tag', 'vent_altL');
z_tmpL.String = z_str;
remove_frame(z_tmpL)

% Update tooltips
x_tmpL = findobj(ancestor(src, 'figure'), 'Tag', 'vent_lat');
x_tmpL.TooltipString = x_ttp;
y_tmpL = findobj(ancestor(src, 'figure'), 'Tag', 'vent_lon');
y_tmpL.TooltipString = y_ttp;
z_tmpL = findobj(ancestor(src, 'figure'), 'Tag', 'vent_alt');
z_tmpL.TooltipString = z_ttp;