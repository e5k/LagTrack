load '/Users/Seb/Documents/WORK/Codes/LagTrack/projects/MSH/validation.mat'

part.run_name = 'MSH_valid';
diam = [707,500,354,250,177,125,88,63];
P = cell(length(diam),1);

for i = 1:length(P)
   part_tmp = part;
   part_tmp.part.name = [part.part.name, '_', num2str(round(-log2(diam(i)/1e3),1), '%1.1f'), 'phi'];
   part_tmp.part.diam = diam(i)/1e6;
   P{i} = part_tmp;
end

%get_trajectory(P);

pthPart     = '/Users/Seb/Documents/WORK/Codes/LagTrack/projects/MSH_valid/';
pthGIS      = '/Users/Seb/Dropbox/Project_current/LagTrack/Hysplit validation/';

parts = {'707um', '500um', '354um', '250um', '177um', '125um', '88um', '63um'};
cmap = lines(length(parts));

figure; hold on

for i = 1:length(diam)
    load([pthPart, 'validation_', num2str(round(-log2(diam(i)/1e3),1), '%1.1f'), 'phi.mat']);
    
    fl1 = dir([pthGIS, parts{i}, '/gis_*']);
    fl2 = dir([pthGIS, parts{i},'/', fl1.name, '/GIS_DEP*.shp']);
    ctr = shaperead([pthGIS, parts{i}, '/', fl1.name, '/',fl2(end).name]);
    
    idx = floor(linspace(1, length(part.traj.x), 100));
    
    Pl(i) = plot3(part.traj.lon(idx), part.traj.lat(idx), part.traj.z(idx),'-', 'color', cmap(i,:));
    Z(i) = patch(ctr(1).X, ctr(1).Y, cmap(i,:))%zeros(size(ctr(1).Y)), 'FaceVertexCData',cmap(i,:),'FaceColor','flat');
end
plot_google_map
legend(Pl,parts)
box on
grid on