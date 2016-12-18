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

%Calculate the distances along the axis
y_dist = ll2dist(ax(3), ax(1), ax(4), ax(1)); %diff(ax(1:2));
x_dist = ll2dist(ax(3), ax(1), ax(3), ax(2)); %diff(ax(3:4));
z_dist = h.ZLim(2)*1000;

%Adjust the aspect ratio
c_adj = cosd(mean(ax(3:4)));
dar = [1 c_adj z_dist/mean([x_dist y_dist])/10];
pbar = [x_dist*c_adj/y_dist 1 1];
set(gca, 'DataAspectRatio',dar,'PlotBoxAspectRatio',pbar);
