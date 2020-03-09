% Remove frame of editbox using JAVA
function remove_frame(src,~)

jEditbox = findjobj(src);
try 
    jEditbox.setBorder([]); % or: set(jEditbox,'Border',[])
catch
    
end