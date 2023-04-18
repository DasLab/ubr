function output_labels_csv(datafile,good_idx,r_norm,ids,sequences,BLANK_OUT5,BLANK_OUT3);

fid = fopen(datafile,'w');
fprintf(fid,'eterna_id,sequence');
Nres = length(sequences{1});
for k = 1:Nres;
    fprintf(fid,',SHAPE_%04d',k);
end
fprintf(fid,'\n');
for i = good_idx'
    fprintf(fid,'%d,%s', ids(i),sequences{i});
    for k = 1:Nres;
        val = r_norm(i,k);
        if (k <= BLANK_OUT5 | k > (Nres - BLANK_OUT3 ) ) val = NaN; end;
        fprintf(fid,',%7.3f',val);
    end
    fprintf(fid,'\n');
end
fclose(fid);
fprintf('Outputted heading and %d rows to %s\n',length(good_idx),datafile);
