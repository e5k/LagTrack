% Remove frame of editbox using JAVA
function remove_frame(src,evt)
jEditbox = findjobj(src);
jEditbox.setBorder([]); % or: set(jEditbox,'Border',[])