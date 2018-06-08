function plot_ATM(src, ~, atm)

AX      = findobj(ancestor(src, 'figure'), 'Type', 'axes'); % Retrieve axes
part    = guidata(src);                                     % Retrieve stored particle
cla(AX);    % Clear axes

% Variable to plot
varL = get(findobj(ancestor(src, 'figure'), 'Tag', 'varList'), 'String');
varS = get(findobj(ancestor(src, 'figure'), 'Tag', 'varList'), 'Value');
varT = varL{varS};

% Time
timL = get(findobj(ancestor(src, 'figure'), 'Tag', 'timList'), 'String');
timS = get(findobj(ancestor(src, 'figure'), 'Tag', 'timList'), 'Value');

% Altitude
levS = get(findobj(ancestor(src, 'figure'), 'Tag', 'levList'), 'Value');
 
% Plot
if strcmp(varT, 'Wind velocity')
    quiver(AX, atm.lon, atm.lat, squeeze(atm.u(:,:,levS,timS)), squeeze(atm.v(:,:,levS,timS)));
    axis equal
    if length(atm.lon)>1 && length(atm.lat)>1
        axis([min(atm.lon), max(atm.lon), min(atm.lat), max(atm.lat)]);
    end
else
    if strcmp(varT, 'U wind')
        var = 'u';
        lab = 'U wind (m/s)';
    elseif strcmp(varT, 'V wind')
        var = 'v';
        lab = 'V wind (m/s)';
    elseif strcmp(varT, 'Temperature')
        var = 'temp';
        lab = 'Temperature (deg K)';
    elseif strcmp(varT, 'Relative humidity')
        var = 'humid';
        lab = 'Relative humidity (%)';
    elseif strcmp(varT, 'Air density')
        var = 'rhoair';
        lab = 'Density (kg/m3)';
    elseif strcmp(varT, 'Air viscosity')
        var = 'muair';
        lab = 'Viscosity';
    end
    
    tmp = atm.(var);
    
    pcolor(AX, atm.lon, atm.lat, squeeze(tmp(:,:,levS,timS)));
    c = colorbar;
    ylabel(c, lab)   
    axis(AX, 'equal'); axis(AX, 'tight');
    caxis([min(min(min(tmp(:,:,levS,:)))) max(max(max(tmp(:,:,levS,:))))]);
end

% Plot options
xlabel(AX, 'Longitude');
ylabel(AX, 'Latitude');

plot_google_map('maptype','terrain');
AX.Layer = 'top';

% If the vent is defined in the GUI, plot it
if ~isempty(part) 
    if part.vent.lon > -9999
        hold(AX, 'on');
        plot(AX, part.vent.lon, part.vent.lat, '^', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
    end
end

% Update info
set(findobj(ancestor(src, 'figure'), 'Tag', 'info'), 'String', ['Altitude: ', num2str(round(mean(mean(squeeze(atm.alt(:,:,levS,timS)))))), ' m asl']);
set(findobj(ancestor(src, 'figure'), 'Tag', 'DateInfo'), 'String', timL{timS});
