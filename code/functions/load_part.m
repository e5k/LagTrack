% Load particle to the GUI

function load_part(src, ~)

[fl,pth] = uigetfile('projects/*.mat', 'Multiselect', 'on');

if iscell(fl)                                       % If multiple particles were selected
    for i = 1:length(fl)
        load([pth, filesep, fl{i}]);
        APDTA = getappdata(ancestor(src, 'figure'));
        preprocess_part(part, APDTA, src);
    end
elseif ~iscell(fl) && fl>0                          % If one single particle was selected
    load([pth, filesep, fl]);
    APDTA = getappdata(ancestor(src, 'figure'));
    preprocess_part(part, APDTA, src);
else                                                % If the user cancelled
    return
end


function preprocess_part(part, APDTA, src)
% Fill up GUI
set(findobj(ancestor(src, 'figure'), 'Tag', 'name'), 'String', part.run_name);
set(findobj(ancestor(src, 'figure'), 'Tag', 'vent_lat'), 'String', num2str(part.vent.lat));
set(findobj(ancestor(src, 'figure'), 'Tag', 'vent_lon'), 'String', num2str(part.vent.lon));
set(findobj(ancestor(src, 'figure'), 'Tag', 'vent_alt'), 'String', num2str(part.vent.alt));
set(findobj(ancestor(src, 'figure'), 'Tag', 'date'), 'String', datestr(part.date));
set(findobj(ancestor(src, 'figure'), 'Tag', 'atm'), 'String', part.path.nc);
set(findobj(ancestor(src, 'figure'), 'Tag', 'dem'), 'String', part.path.dem);

set(findobj(ancestor(src, 'figure'), 'Tag', 'part_name'), 'String', part.part.name);
set(findobj(ancestor(src, 'figure'), 'Tag', 'part_diam'), 'String', num2str(part.part.diam*1e3));
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

if strcmp(part.adv.solution, 'analytical'); V = 2; elseif strcmp(part.adv.solution, 'euler'); V = 1; end
set(findobj(ancestor(src, 'figure'), 'Tag', 'adv_sol'), 'Value', V);
set(findobj(ancestor(src, 'figure'), 'Tag', 'adv_dt'), 'String', num2str(part.adv.dt));
set(findobj(ancestor(src, 'figure'), 'Tag', 'adv_drag'), 'String', num2str(part.adv.drag));
if strcmp(part.adv.solution, 'none'); V = 1; elseif strcmp(part.adv.solution, 'subset'); V = 2; elseif strcmp(part.adv.solution, 'complete'); V = 3; end
set(findobj(ancestor(src, 'figure'), 'Tag', 'adv_int'), 'Value', V);
if strcmp(part.adv.solution, 'linear'); V = 1; elseif strcmp(part.adv.solution, 'nearest'); V = 2; elseif strcmp(part.adv.solution, 'pchip'); V = 3; elseif strcmp(part.adv.solution, 'cubic'); V = 4; elseif strcmp(part.adv.solution, 'spline'); V = 5; end
set(findobj(ancestor(src, 'figure'), 'Tag', 'adv_meth'), 'Value', V);
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