function update_table(src, part)

Tdata = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataTable'), 'Data');
Ldata = get(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'String');

data = {part.part.name, ...
    part.part.diam*1e3,...
    part.part.dens,...
    part.part.flat,...
    part.part.elon,...
    part.rel.x,...
    part.rel.y,...
    part.vent.alt+part.rel.z,...
    datestr(part.date+part.rel.t/3600/24),...
    part.traj.t(end),...
    part.traj.z(end),...
    part.traj.w(end)
    };


if isempty(Tdata)
    Tdata = data;
    Ldata = {part.part.name};
else
    Tdata(size(Tdata,1)+1,:) = data;
    Ldata{size(Ldata,1)+1,1} = part.part.name;
end

set(findobj(ancestor(src, 'figure'), 'Tag', 'DataTable'), 'Data', Tdata);
set(findobj(ancestor(src, 'figure'), 'Tag', 'DataList'), 'String', Ldata);