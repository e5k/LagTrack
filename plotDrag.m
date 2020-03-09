noC = load('/Users/sebastien/Documents/Codes/LagTrack/projects/standard/balEU_noCritical.mat');
noC = noC.part;
C = load('/Users/sebastien/Documents/Codes/LagTrack/projects/standard/balEU_Critical.mat');
C = C.part;

figure 
subplot(2,1,1)
plot(C.traj.y, C.traj.z)
hold on
plot(noC.traj.y, noC.traj.z)
plot(ej(:,1), ej(:,2)); 
legend({'Critical Cd = 0.2', 'No Funny stuff', 'Eject'})
xlabel('x (m)')
ylabel('elevation (m)')

subplot(2,1,2)
plot(C.traj.y, C.traj.Cd)
hold on
plot(noC.traj.y, noC.traj.Cd)
plot(ej(:,1), ej(:,3)); 
legend({'Critical Cd = 0.2', 'No Funny stuff', 'Eject'})
xlabel('x (m)')
ylabel('C_d')
