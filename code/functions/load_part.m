% Pre-process and run
function load_part(src, ~)

[fl,pth] = uigetfile('projects/*.mat');

if isempty(fl)
    return
    
else
    load([pth, filesep, fl])
    APDTA = getappdata(ancestor(src, 'figure'));    

    if isfield(APDTA, 'pltData') && ~isempty(find(strcmp(fieldnames(APDTA.pltData), part.part.name)==1,1))
        errordlg('This particle is already loaded');
        return
    end
end

%load([pth, filesep, fl])

% Fill up GUI
set(findobj(ancestor(src, 'figure'), 'Tag', 'name'), 'String', part.run_name);
set(findobj(ancestor(src, 'figure'), 'Tag', 'vent_lat'), 'String', num2str(part.vent.lat));
set(findobj(ancestor(src, 'figure'), 'Tag', 'vent_lon'), 'String', num2str(part.vent.lon));
set(findobj(ancestor(src, 'figure'), 'Tag', 'vent_alt'), 'String', num2str(part.vent.alt));
set(findobj(ancestor(src, 'figure'), 'Tag', 'date'), 'String', datestr(part.date));
set(findobj(ancestor(src, 'figure'), 'Tag', 'atm'), 'String', part.path.nc);
set(findobj(ancestor(src, 'figure'), 'Tag', 'dem'), 'String', part.path.dem);

set(findobj(ancestor(src, 'figure'), 'Tag', 'part_name'), 'String', part.part.name);
set(findobj(ancestor(src, 'figure'), 'Tag', 'part_diam'), 'String', num2str(part.part.diam*1e2));
set(findobj(ancestor(src, 'figure'), 'Tag', 'part_dens'), 'String', num2str(part.part.dens));
set(findobj(ancestor(src, 'figure'), 'Tag', 'part_flat'), 'String', num2str(part.part.flat));
set(findobj(ancestor(src, 'figure'), 'Tag', 'part_elon'), 'String', num2str(part.part.elon));

set(findobj(ancestor(src, 'figure'), 'Tag', 'rel_x'), 'String', num2str(part.rel.x));
set(findobj(ancestor(src, 'figure'), 'Tag', 'rel_y'), 'String', num2str(part.rel.y));
set(findobj(ancestor(src, 'figure'), 'Tag', 'rel_z'), 'String', num2str(part.rel.z));
set(findobj(ancestor(src, 'figure'), 'Tag', 'rel_t'), 'String', num2str(part.rel.t));
set(findobj(ancestor(src, 'figure'), 'Tag', 'rel_vx'), 'String', num2str(part.rel.vx));
set(findobj(ancestor(src, 'figure'), 'Tag', 'rel_vy'), 'String', num2str(part.rel.vy));
set(findobj(ancestor(src, 'figure'), 'Tag', 'rel_vz'), 'String', num2str(part.rel.vz));

set(findobj(ancestor(src, 'figure'), 'Tag', 'adv_sol'), 'String', part.adv.solution);
set(findobj(ancestor(src, 'figure'), 'Tag', 'adv_dt'), 'String', num2str(part.adv.dt));
set(findobj(ancestor(src, 'figure'), 'Tag', 'adv_drag'), 'String', num2str(part.adv.drag));
set(findobj(ancestor(src, 'figure'), 'Tag', 'adv_int'), 'String', part.adv.interp);
set(findobj(ancestor(src, 'figure'), 'Tag', 'adv_meth'), 'String', part.adv.method);
set(findobj(ancestor(src, 'figure'), 'Tag', 'adv_range'), 'String', num2str(part.adv.range));
set(findobj(ancestor(src, 'figure'), 'Tag', 'adv_skip'), 'String', num2str(part.adv.skip));

set(findobj(ancestor(src, 'figure'), 'Tag', 'run_btn'), 'Enable', 'on');

% Update GUI data
update_table(src, part);
guidata(src, part);

% Set data to plot
APDTA.pltData.(part.part.name) = part;
setappdata(ancestor(src, 'figure'), 'pltData',APDTA.pltData);

% Enable buttons
enableUI(src, 'on');