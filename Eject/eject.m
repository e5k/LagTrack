function [ traj ] = eject( thetadeg,vi,diam )
 
%% This function calculates the trajectory of a ballistic projetile during a volcanic eruption.

%The eject function takes ejection angle, initial velocity and diameter of ballistic as input parameters.
%It has been defined like this in order to be called within a matlab routine
%that performs probability analyses using random values of the above listed
%three parameters.

%% Other input parameters (set within function)
% adapt these to specific situation as needed

    dragred = 100; %zone of reduced drag (empirical value)
    elev = 0; %elevation of ejection point
    lapse = 6.5; %thermal lapse rate
    rad = diam / 2; %radius of ballistic projectiles
    rhor = 2400; %density of ballistic projectiles
    w = 0; %wind velocity
    xi = -1; %vertical distance between landing point and ejection point 
    TzeroC = 14; %Temperature at sealevel

    g=9.80665; %gravity
    erad=6370800; %earth's radius (meters)
    rair=286.98; %air density
    
 %% Note: Only negative values (<0) are allowed for xi.
 % In a separate script, the trajectory of the projectiles is compared to actual topography to determine the landing point.

%% Check whether xi is negative
    if xi >= 0 
        disp('Please enter a negative value for xi.')
        return
    end 

%% Initialise variables

theta=degtorad(thetadeg); %convert theta to radians
TzeroK = TzeroC + 273.15; %convert temperature to Kelvin
ttotcd0 = (-vi * sin(theta) - sqrt(vi ^ 2 * (sin(theta) ^ 2) - 2 * g * xi))/ (-g); %estimate time of travel, assuming no drag
dt=ttotcd0/500; %length of timesteps
x=zeros(500,1); %vector for storage of horizontal positions from vent
z=zeros(500,1); %vector for storage of heights above vent
vx=zeros(500,1);%vector for storage of x-components of velocity
vz=zeros(500,1);%vector for storage of z-components of velocity
CD=zeros(500,1);
vx(1)=vi*cos(theta); %initial vx
vz(1)=vi*sin(theta); %initial vz
v=vi; %initial velocity
time=0; %initial time (seconds)
zmax=0; %maximum point in block's trajectory
po=101300 * ((TzeroK - elev * lapse / 1000)/ TzeroK) ^ (-g / (rair ^ 2 * lapse / 1000)); %air pressure at vent
mass = 8 * rad ^ 3 * rhor; %mass of clast calculated based on the option "high cube" (as defined in the Visual Basic code by Mastin 2001)
area = 4 * rad * rad; %area of clast calculated based on the option "high cube" (as defined in the Visual Basic code by Mastin 2001)

%% timesteps

xnow=x(1);
znow=z(1);
vxnow=vx(1);
vznow=vz(1);
cdnow = nan;
i=1;

while znow >= xi %run until flightpath crosses the landing plane (xi metres below ejection) 
    
    %store current positions and velocities in vectors
    x(i)=xnow; 
    z(i)=znow;
    vx(i)=vxnow;
    vz(i)=vznow;
    CD(i) = cdnow;

    temp=TzeroK - (znow + elev) * lapse / 1000;  %temperature at znow
    visc = (0.000172 * (390 / (temp + 117)) *(temp / 273) ^ 1.5) / 10; %air viscosity at znow
    pressure = po *((TzeroK - (znow + elev) * lapse / 1000) / TzeroK) ^(g / (rair * lapse / 1000)); %air pressure at znow
    rhoa = pressure / (rair * temp); %air density at znow
    c = 20.116 * sqrt(temp); %sound speed at znow
    vwind = sqrt(vznow ^ 2 + (vxnow - w) ^ 2); %velocity relative to wind
    mach = vwind / c;                            %Mach number
    cd=draghicube(mach);                         % drag coefficient, calculated with the draghicube function
    
    %calculate drag reduction in reduced-drag zone:
    if dragred>0 && sqrt(xnow^2 + znow^2) < dragred
        cdnow = cd * (sqrt(xnow^2 + znow^2)/dragred)^2;
    else
        cdnow=cd;
    end
    
    %% Runge Kutta to determine new position and velocities
   
    RK=rk(xnow, znow, vxnow, vznow,w,dt,rhoa,cdnow,area,mass,g,rhor);
    
    %update position and velocity
    xnow=RK(1);
    znow=RK(2);
    vxnow=RK(3);
    vznow=RK(4);
    
    %update i so that next value is stored in the next position of storage
    %vector
    i=i+1;
    
    %update time
    time = time + dt;

end
% 
% So now the elevation of the clast is less than xi.
% % Interpolate to find the exact value of x when the block hits the earth
% 
%         xfinal = x(i-1) + (xnow - x(i-1)) * (xi - z(i-1)) / (znow - z(i-1));     %final x
%         tfinal = time + dt * (xi - z(i-1)) / (znow - z(i-1));    %final time (s)
%         zfinal = xi;                                         %final z       
%         vxfinal = vx(i-1) + (vxnow - vx(i-1)) * (xi - z(i-1)) / (znow - z(i-1)); %final vx
%         vzfinal = vz(i-1) + (vznow - vz(i-1)) * (xi - z(i-1)) / (znow - z(i-1)); %final vz
%         vfinal = sqrt(vxfinal ^ 2 + vzfinal ^ 2);                         %final velocity
% 
%    
% store the final values        
% x(i)=xfinal;
% z(i)=zfinal;
% 
% delete unneeded entries in storage vector:
% 
% x=x(1:i);
% z=z(1:i);
% 
% Output: the trajectory (in 2D) of the ballistic ejectile. Add values to the output vector as
% desired.
% traj=[x z]; 

idx = x>0;
traj = [x(idx), z(idx), CD(idx)];

end

