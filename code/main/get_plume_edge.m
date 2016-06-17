function [xt, yt, xb, yb] = get_plume_edge(theta, w, x0, y0)
w  = w/2;

xb = x0 + sign(cosd(theta)) * (w*sind(theta));
yb = y0 - sign(cosd(theta)) * (w*cosd(theta));
xt = x0 - sign(cosd(theta)) * (w*sind(theta));
yt = y0 + sign(cosd(theta)) * (w*cosd(theta));
