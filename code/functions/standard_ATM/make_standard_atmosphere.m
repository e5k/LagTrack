a = xlsread('1976_standard_ATM.xlsx');
rhoair = zeros(1,1, size(a,1), 2);
rhoair(:,:,:,1) = reshape(a(:,2), 1,1,size(a,1),1);
rhoair(:,:,:,2) = reshape(a(:,2), 1,1,size(a,1),1);

muair = zeros(1,1, size(a,1), 2);
muair(:,:,:,1) = reshape(a(:,3), 1,1,size(a,1),1);
muair(:,:,:,2) = reshape(a(:,3), 1,1,size(a,1),1);

alt = zeros(1,1, size(a,1), 2);
alt(:,:,:,1) = reshape(a(:,1), 1,1,size(a,1),1);
alt(:,:,:,2) = reshape(a(:,1), 1,1,size(a,1),1);

u = zeros(1,1, size(a,1), 2);
u(:,:,:,1) = zeros(1,1,size(a,1),1);
u(:,:,:,2) = zeros(1,1,size(a,1),1);

v = zeros(1,1, size(a,1), 2);
v(:,:,:,1) = zeros(1,1,size(a,1),1);
v(:,:,:,2) = zeros(1,1,size(a,1),1);

time = [datenum([1948,1,1]); now];

clear atm
atm.lat = lat;
atm.lat = lat;
atm.lon = lon;
atm.time = time;
atm.alt = alt;
atm.u = u;
atm.v = v;
atm.rhoair = rhoair;
atm.muair = muair;