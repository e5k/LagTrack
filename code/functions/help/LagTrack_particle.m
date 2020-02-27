function LagTrack_particle
%
% MODEL THE PARTICLE DISPERSAL WITHOUT THE GUI
% Add main functions to the path:
%   >> addpath(genpath('code/')
% Load default particle
%   >> load default_part
% A variable called part is now created
%
% PARTICLE VARIABLES
% General:
%   part.run_name           Run name, i.e. folder name in projects/ containing the output particles
%   part.run_mode           Defines if the code is run forward (1) or backward (2)
%   part.date               Eruption date, should be within the time extent of atmospheric data
%
% Vent:
%   part.vent.lat           Vent latitude (decimal degree, negative in S hemisphere)
%   part.vent.lon           Vent longitude (decimal degree, negative in W hemisphere)
%   part.vent.alt           Vent elevation (m asl)
%
% Particle
%   part.part.name          Particle name
%   part.part.diam          Diameter (m)
%   part.part.dens          Density (kg/m3)
%   part.part.fl            Flatness
%   part.part.el            Elongation
%
% Path to input files
%   part.path.dem           Path to the .mat file containing atmospheric data
%   part.path.nc            Path to the .mat file containing the calculation grid
%
% Release properties
%   part.rel.x              x displacement (m, negative towards W)
%   part.rel.y              y displacement (m, negative towards S)
%   part.rel.z              Release altitude (m above vent)   
%   part.rel.t              t offset relative to part.date (s)
%   part.rel.vx             Initial vx velocity (m/s). If 0, vx is initialized with the wind velocity
%   part.rel.vy             Initial vy velocity (m/s). If 0, vy is initialized with the wind velocity
%   part.rel.vz             Initial vz velocity (m/s).
%
% Advanced options
%   part.adv.solution       Solution to particle motion - either 'euler' or 'analytical'
%   part.adv.dt             Integration time step (s)
%   part.adv.drag           Radius of region of reduced drag (m, sensus Mastin 2001)
%   part.adv.interp         Interpolation mode for atmospheric data
%                               - 'none'        No interpolation
%                               - 'subset'      Interpolate only a subset
%                               - 'complete'    Interpolate the entire dataset (slower)
%   part.adv.method         Interpolation technique (see Matlab help for interpn)
%                           'linear', 'nearest', 'pchip', 'cubic', 'spline'
%   part.adv.range          Range of interpolation when part.adv.interp = 'subset'
%   part.adv.skip           Number of time steps to skips between interpolations. 
%                           E.g. if part.adv.dt = 0.1 and part.adv.skip = 600, a
%                           new interpolation will be performed every minute
%
% See also LagTrack_functions


