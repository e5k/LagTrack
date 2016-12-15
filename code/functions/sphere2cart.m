function varargout = sphere2cart(varargin)
% Bearing:      Degree from north
% Inclination:  Degree from vertical
% Velocity:     Initial velocity (m/s)

% check number of input parameters
if nargin == 0 || nargin == 2
    answer      = inputdlg({'Bearing (degrees N)', 'Inclination (degrees from vertical)', 'Ejection velocity (m/s)'}, ' ', 1);
    if isempty(answer)
        return
    end
    bearing     = str2double(answer{1});
    inclination = str2double(answer{2});
    velocity    = str2double(answer{3});    
elseif nargin == 3   
    bearing     = varargin{1};
    inclination = varargin{2};
    velocity    = varargin{3};
else
    error('Wrong number of input arguments, should be either 0 or 3');
end

% check number of output parameters
if nargout > 3
    error('Too many output parameters')
end

% Check that input are valid
if bearing<0 || bearing>360
    error('Bearing should be comprised between 0 and 360');
elseif inclination<0 || inclination>90
    error('Inclination should be comprised between 0 and 90');
elseif velocity<0
    error('Velocity should be positive')
end
    
% Convert geographical bearing into trigonometric system
if bearing >90 && bearing < 360
    azimuth = 360-bearing+90;
elseif bearing>=0 && bearing <= 90
    azimuth = bearing*(-1)+90;
end
azimuth   = deg2rad(azimuth);
elevation = deg2rad(90-inclination);

[a(1),a(2),a(3)] = sph2cart(azimuth, elevation, velocity);

if nargin == 0
    fprintf('U velocity: %.2f m/s\n', a(1));
    fprintf('V velocity: %.2f m/s\n', a(2));
    fprintf('W velocity: %.2f m/s\n', a(3));
    
elseif nargin == 2
    src = varargin{1};
    set(findobj(ancestor(src, 'figure'), 'Tag', 'TopPanel'), 'Selection', 3);
    set(findobj(ancestor(src, 'figure'), 'Tag', 'rel_vx'), 'String', num2str(a(1)));
    set(findobj(ancestor(src, 'figure'), 'Tag', 'rel_vy'), 'String', num2str(a(2)));
    set(findobj(ancestor(src, 'figure'), 'Tag', 'rel_vz'), 'String', num2str(a(3)));
end

if nargout > 0
    for i = 1:nargout
        varargout{i} = a(i);
    end
end