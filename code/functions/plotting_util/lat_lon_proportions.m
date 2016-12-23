function lat_lon_proportions(h)
%LAT_LON_PROPORTIONS Proportions a lat/lon bounded map
%
%   LAT_LON_PROPORTIONS Scales a plot in latitude & longitude axis to 
%       meters. It will compress the x axis by cos(latitude) in order to
%       reflect the relationship between a degree latitude and a degree
%       longitude at the center of the map. The major assumption here is
%       sperical Earth.
%
% By: A Weaver, April 2004. Slightly tweaked by J Sullivan, August 2011

%Grab the axis limits
if nargin > 0
    ax = axis(h);
else
    ax = axis;
end

% Calculate the distances along the axis
z_dist = (h.ZLim(2) - h.ZLim(1));
y_dist = ll2dist(ax(3), ax(1), ax(4), ax(1), 6371*1e3)/1e3; %diff(ax(1:2));
x_dist = ll2dist(ax(3), ax(1), ax(3), ax(2), 6371*1e3)/1e3; %diff(ax(3:4));

% Calculate the ratio of axes
xRatio = x_dist/min([x_dist y_dist]);
yRatio = y_dist/min([x_dist y_dist]);
zRatio = min([x_dist y_dist])/z_dist *  min([x_dist y_dist z_dist])/max([x_dist y_dist z_dist])*100;
dar    = [xRatio yRatio zRatio];

set(gca, 'DataAspectRatio',dar);