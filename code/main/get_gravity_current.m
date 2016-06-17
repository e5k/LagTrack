function [ub, width] = get_gravity_current(s_avg,H, wind)
% H = km above vent
% x = distance in meters

%x = 0:100:1000000;

%Q = (H./.287).^(1/.19);

Q = 1.95e10;

lam = .8;
eps = 3.9;
N   = .02;

% Velocity due to gravity current
ub  = sqrt(lam*N*Q./(eps.*s_avg));

% Plume width
a       = wind / sqrt(lam*N*Q./(eps));
width   = 2*s_avg / (1+a * sqrt(s_avg));


% figure; 
% 
% semilogy(x./1000,ub1, '-k.');
% hold on
% semilogy(x./1000,ub2, '-r.');