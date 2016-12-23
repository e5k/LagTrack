% Used to set path in GUI
function set_path(~, ~, src, filter, pthI, ttl)
[fl, pth] = uigetfile([pthI, filter], ttl);
set(src, 'String', fullfile(pth, fl));
check_var(src);