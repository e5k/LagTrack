% Used to set path in GUI
function set_path(~, ~, O, filter, pthI, ttl)
[fl, pth] = uigetfile([pthI, filter], ttl);
set(O, 'String', fullfile(pth, fl));
check_var(O);