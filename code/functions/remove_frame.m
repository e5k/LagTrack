% Remove frame of editbox using JAVA
function remove_frame(objH)
jEditbox = findjobj(objH);
jEditbox.setBorder([]); % or: set(jEditbox,'Border',[])