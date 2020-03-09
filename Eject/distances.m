function [ distances ] = distances( x,z_ballistic,direction )

%% This function calculates the flight distance of a ballistic clast by comparing its flight trajectory with the topography of the volcano.

%note that x and z_ballistic are matrices where every row represents the
%trajectory of one clast - i.e. x(i,:) is clast i. 
%direction is a vector, with one entry for every clast

%% Storage vectors 

distances_north=zeros(length(direction),1);
distances_south=zeros(length(direction),1);
distances_east=zeros(length(direction),1);
distances_west=zeros(length(direction),1);

%% Definition of topography functions

%% Import topography data (here: example for Ruapehu, elevation profiles for 4 directions: east, west, north and south of the Crater

[~, ~, raw] = xlsread('.../Topography.xlsx','Ruapehu'); %adjust the file path to your storage system!
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
Topography = table;

%% Allocate imported array to column variable names
x_east = data(:,1);
topography_east = data(:,2);
x_north = data(:,3);
topography_north = data(:,4);
x_west = data(:,5);
topography_west = data(:,6);
x_south = data(:,7);
topography_south = data(:,8);

%% Clear temporary variables
clearvars data raw R;

%% For Ruapehu's eastern flank: Topography_east

% convert topography such that the vent elevation is at 0:
topography_east=topography_east-2380;

%fitted function
Topography_east=fit(x_east,topography_east,'linearinterp');
 

%% For Ruapehu's western flank: Topography_west

% convert topography such that the lake is at 0:
topography_west=topography_west-2380;

%fitted function
Topography_west=fit(x_west,topography_west,'linearinterp');


%% For Ruapehu's northern flank: Topography_north

% convert topography such that the lake is at 0:
topography_north=topography_north-2380;

%fitted function
Topography_north=fit(x_north,topography_north,'linearinterp');


%% For Ruapehu's southern flank: Topography_south

% convert topography such that the lake is at 0:
topography_south=topography_south-2380;

%fitted function
Topography_south=fit(x_south,topography_south,'linearinterp');
 

%% Calculate the flight distance by comparing the flight trajectory against topography 

for k=1:length(direction)
    
    if direction(k)==1 %first check which direction the clast is flying (1 means "north")
        % Calculate Distance for northern direction:

        for i=1:length(z_ballistic(k,:))
            if z_ballistic(k,i)<Topography_north(x(k,i))
                x_impact_north=x(k,i-1)+((x(k,i)-x(k,i-1))/2); 
                break
            end
        end

        distances_north(k)=x_impact_north;
    end 


    if direction(k)==2 %first check which direction the clast is flying (2 means "south")
        % Calculate Distance for southern direction:

        for i=1:length(z_ballistic(k,:))
            if z_ballistic(k,i)<Topography_south(x(k,i))
                x_impact_south=x(k,i-1)+((x(k,i)-x(k,i-1))/2);
                break
            end
        end

        distances_south(k)=x_impact_south;
    end 

    if direction(k)==3 %first check which direction the clast is flying (3 means "east")
        
        % Calculate Distance for eastern direction:

        for i=1:length(z_ballistic(k,:))
            if z_ballistic(k,i)<Topography_east(x(k,i))
                x_impact_east=x(k,i-1)+((x(k,i)-x(k,i-1))/2);
                break
            end
        end

        distances_east(k)=x_impact_east;
        
    end

    if direction(k)==4 %first check which direction the clast is flying (4 means "west")
        % Calculate Distance for western direction:

        for i=1:length(z_ballistic(k,:))
            if z_ballistic(k,i)<Topography_west(x(k,i))
                x_impact_west=x(k,i-1)+((x(k,i)-x(k,i-1))/2);
                break
            end
        end

        distances_west(k)=x_impact_west;
    end


end

%% Output: flight distances that take the topography into account
distances=[distances_north distances_south distances_east distances_west];

end
