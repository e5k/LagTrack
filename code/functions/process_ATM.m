function process_ATM(filename, dataset, varargin)
% PREPROCESS_ATM Converts NetCDF data into a Matlab matrix.
%   PREPROCESS_ATM(filename, dataset) Loads NetCDF files located in
%       input/filename/ and save a .mat file in the same folder.
%       The dataset is either 'Interim' for ECMWF Era-Interim or 
%       'Reanalysis2' for NOAA Reanalysis 2 database.
%
%   See also writeECMWFAPIKey, download_ATM.

% This function is part of LagTrack.
% Written by Sebastien Biass & Gholamhossein Bagheri
% GPLv3

% This is setup in case we ever figure out how to interpolate a 4D matrix
% if isempty(varargin)
% elseif ~isempty(varargin) && nargin==5  
%     resX = varargin{1};
%     resP = varargin{2};
%     resT = varargin{3};
%     intM = 'linear';
% elseif ~isempty(varargin) >0 && nargin==6  
%     resX = varargin{1};
%     resP = varargin{2};
%     resT = varargin{3};
%     intM = varargin{4};
%     if isempty(regexp(tmp, '(linear|nearest|pchip|cubic|spline)', 'once'))
%         error('Wrong interpolation method. Choose linear, nearrest, pchip, cubic or spline')
%     end
% else
%     error('Wrong number of input arguments')
% end

if ~isempty(regexp(dataset, 'Reanalysis', 'once'))
    if length(varargin) ~= 8
        error('Wrong number of arguments. NOAA data requires to input:\nlat_min, lat_max, lon_min, lon_max, year_min (yyyy), year_max (yyyy), month_min (mm), month_max (mm)'); 
    else
        lat_min     = varargin{1};
        lat_max     = varargin{2};
        lon_min     = varargin{3};
        lon_max     = varargin{4};
        year_min    = varargin{5};
        year_max    = varargin{6};
        month_min   = varargin{7};
        month_max   = varargin{8};
    end
else

end


display('Processing atmospheric data...')
if strcmp(dataset, 'Interim')
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

else
        
    if strcmp(dataset, 'Reanalysis1')
        source_dir = 'input/wind/_Reanalysis1_Rawdata/';
    elseif strcmp(dataset, 'Reanalysis2')
        source_dir = 'input/wind/_Reanalysis2_Rawdata/';
    else
        error('Unknown dataset requested')
    end
    
    varList = {'hgt', 'uwnd', 'vwnd', 'rhum', 'air'};
    strList = {'alt', 'u', 'v', 'humid', 'temp'};
 
    for iY = year_min:year_max
        for iV = 1:length(varList)
            fprintf('\tReading file %s\n', [varList{iV}, '.', num2str(iY), '.nc'])
            ncfile = [source_dir, varList{iV}, '.', num2str(iY), '.nc'];
            % If reading the first file, define storage matrix
            if iV == 1 && iY == year_min
                LT          = sort(ncread(ncfile, 'lat')); % Make sure it is in inreasing order
                LN          = sort(ncread(ncfile, 'lon')); % Make sure it is in inreasing order
                [~,latImin] = min(abs(LT - lat_min));
                [~,latImax] = min(abs(LT - lat_max));                
                [~,lonImin] = min(abs(LN - lon_min));
                [~,lonImax] = min(abs(LN - lon_max));
                atm.lat     = double(LT(latImin:latImax));
                atm.lon     = double(LN(lonImin:lonImax));
                atm.level   = double(ncread(ncfile, 'level'));
                atm.time    = datenum([year_min,month_min,1,0,0,0]):0.25:datenum([year_max,month_max,eomday(year_max,month_max),18,0,0]);
                atm.temp    = zeros(length(atm.lat), length(atm.lon), length(atm.level), length(atm.time));
                atm.alt     = zeros(size(atm.temp));
                atm.humid   = zeros(size(atm.temp));
                atm.u       = zeros(size(atm.temp));
                atm.v       = zeros(size(atm.temp));
            end
            % Retrieve the time indices
            time_tmp    = datenum(datenum([1800,1,1,0,0,0]) + double(ncread(ncfile, 'time'))./24)';    % Time (year month day hour min sec)
            timeInc     = ismember(time_tmp, atm.time);
            timeIatm    = ismember(atm.time, time_tmp);
            
            nc_tmp     = permute(double(ncread(ncfile, varList{iV})), [2,1,3,4]);   % Read nc data
            
            % Problem here: Humidity data for Reanalysis1 only goes to a
            % pressure level of 300 mb (~10 km). I fill what is higher with
            % the last value of humidity
            if strcmp(dataset, 'Reanalysis1') && iV == 4
                atm.(strList{iV})(1:length(atm.lat), 1:length(atm.lon), 1:8, timeIatm) = nc_tmp(latImin:latImax, lonImin:lonImax, :, timeInc);
                atm.(strList{iV})(1:length(atm.lat), 1:length(atm.lon), 9:length(atm.level), timeIatm) = ...
                    repmat(nc_tmp(latImin:latImax, lonImin:lonImax, 8, timeInc), 1,1,9,1);
            else
                atm.(strList{iV})(1:length(atm.lat), 1:length(atm.lon), :, timeIatm) = nc_tmp(latImin:latImax, lonImin:lonImax, :, timeInc);
            end
        end
    end
end

% Express degrees E (i.e. >180) into negative values in W hemisphere
atm.lon(atm.lon>180) = -360+atm.lon(atm.lon>180) ;

% Air density
pres            = double(repmat(reshape(atm.level,1,1,length(atm.level)), length(atm.lat), length(atm.lon), 1, length(atm.time)));
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
atm.muair       = mu0.*(a./b).*(tempR./t0R).^(3/2)./10^3;                   % Dynamic viscosity (Pa S)


% %% Interpolation
% if nargin > 2
%     latV = min(atm.lat):resX:max(atm.lat);
%     lonV = min(atm.lon):resX:max(atm.lon);
%     altV = min(atm.level):resP:max(atm.level);
%     timV = min(atm.time):resT/24:max(atm.time);
%     
%     a
%   denf    = interpn(atm.lat, atm.lon, alt_sub, atm.time, atm.rhoair, part.lat(i-1), part.lon(i-1), part.z(i-1), P.date+part.t(i-1)/3600/24, P.adv.method);
%     
%   latT = repmat(atm.lat, 1, size(atm.rhoair,2), size(atm.rhoair,3), size(atm.rhoair,4));
%   lonT = repmat(atm.lon', size(atm.rhoair,1), 1, size(atm.rhoair,3), size(atm.rhoair,4));
%   timT = repmat(reshape(atm.time,1,1,1,length(atm.time)), size(atm.rhoair,1), size(atm.rhoair,2), size(atm.rhoair,3), 1);
%   
%   latM = repmat(flipud(latV'), 1, length(lonV), length(altV), length(timV)); %latM = double(latM);
%   lonM = repmat(lonV, length(latV), 1, length(altV), length(timV)); %lonM = double(lonM);
%   altM = repmat(reshape(altV,1,1,length(altV)), length(latV), length(lonV), 1, length(timV)); %altM = double(altM);
%   timM = repmat(reshape(timV,1,1,1,length(timV)), length(latV), length(lonV), length(altV),1); %timM = double(timM);
%   
%     denf    = interpn(latT, lonT, single(atm.alt), timT, single(atm.rhoair), latM, lonM, altM, timM, intM);
% 
%     denf    = interpn(atm.lat, atm.lon, single(atm.alt), atm.time, single(atm.rhoair), latM, lonM, altM, timM, intM);
% end

display('Saving...')
save(['input/wind/', filename, filesep, filename, '.mat'], 'atm');
display('Done!')