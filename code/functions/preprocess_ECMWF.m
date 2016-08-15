% Extract ECMWF data into a Matlab matrix

function preprocess_ECMWF(filename)

ncfile = ['input/wind/', filename, filesep, filename, '.nc'];

atm.lat             = double(ncread(ncfile, 'latitude'));                       % Latitude (degrees)
atm.lon             = double(ncread(ncfile, 'longitude'));                      % Longitude (degrees)
atm.level           = ncread(ncfile, 'level');                                  % Pressure level (mb)
atm.time            = datenum(datenum([1900,1,1,0,0,0])+double(ncread(ncfile, 'time'))./24);    % Time (year month day hour min sec)
atm.temp            = permute(ncread(ncfile, 't'), [2,1,3,4]);                  % Temperature (deg K)
atm.alt             = permute(ncread(ncfile, 'z')./9.80665, [2,1,3,4]);         % Altitude (m asl)
atm.humid           = permute(ncread(ncfile, 'r'), [2,1,3,4]);                  % Relative humidity (%)
atm.u               = permute(ncread(ncfile, 'u'), [2,1,3,4]);                  % U wind (ms-1)
atm.v               = permute(ncread(ncfile, 'v'), [2,1,3,4]);                  % V wind (ms-1)
% Air density
pres    = double(repmat(reshape(atm.level,1,1,length(atm.level)), length(atm.lat), length(atm.lon), 1, length(atm.time)));            
tempc           = atm.temp-273.15;                                          % Temperature degC
Psat            = 6.1078 .* 10.^(7.5*tempc./(tempc+237.3)).*100;            % Saturation vapor pressure (Pa)
Pv              = atm.humid./100 .* Psat;                                   % Vapor pressure of water (Pa)
Pd              = pres.*100 - Pv;                                           % Partial pressure dry air
Rd              = 287.058;                                                  % Specific gas constant for dry air (J/kg*K)
Rv              = 461.495;                                                  % Specific gas constant for water vapor (J/kg*K)
atm.rhoair      = Pd./(Rd.*atm.temp) + Pv./(Rv.*atm.temp);                  % Air denstiy (kg m-3);
% Air viscosity
tempR           = (tempc+273.15) .* (9/5);                                  % Temperature (deg Rankine)
mu0             = 0.01827;                                                  % Reference viscosity (centipoise)
t0R             = 524.07;                                                   % Reference temperature (deg Rankine)
C               = 120;                                                      % Sutherland's constant
a               = .555*t0R+C;                                                      
b               = .555.*tempR+C;
atm.muair       = mu0.*(a./b).*(tempR./t0R).^(3/2)./10^3;


save(['input/wind/', filename, filesep, filename, '.mat'], 'atm');