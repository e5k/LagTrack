function set_path(src, evt, O, filter, pthI, ttl)
[fl, pth] = uigetfile([pthI, filter], ttl);
set(O, 'String', fullfile(pth, fl));