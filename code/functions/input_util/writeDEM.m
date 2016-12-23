function writeDEM(out_name, xllcorner, yllcorner, Z, cellsize)
Z(isnan(Z)) = -9999;

fid_3=fopen(out_name,'w');
fprintf(fid_3,'%s\n',['ncols         ' '6001'],['nrows         ' '6001'],...
    ['xllcorner     ' num2str(xllcorner)],['yllcorner     ' num2str(yllcorner)],...
    ['cellsize      ' num2str(cellsize)],['NODATA_value  ' num2str(-9999)]);
fclose(fid_3);

dlmwrite(out_name,Z,'-append','delimiter',' ')