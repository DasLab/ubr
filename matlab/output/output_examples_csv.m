function output_examples_csv(datafile,good_idx,ids,sequences,mfe_tags,mfe_structures,mfe_tags_to_output);
% output_examples_csv(datafile,good_idx,ids,sequences,mfe_tags,mfe_structures,mfe_tags_to_output);

if isempty(good_idx); good_idx = [1:length(sequences)]; end;
if size(good_idx,2)>size(good_idx,1); good_idx = good_idx'; end;
fid = fopen( datafile,'w');
fprintf(fid,'eterna_id,sequence');
for k = 1:length(mfe_tags_to_output)
    fprintf(fid,',%s',mfe_tags_to_output{k});
end
fprintf(fid,'\n');

for i = good_idx'
    fprintf(fid,'%d,%s', ids(i),sequences{i});
    for k = 1:length(mfe_tags_to_output)
        m = find(strcmp(mfe_tags,mfe_tags_to_output{k}));
        fprintf(fid,',%s',mfe_structures{m}{i});
    end
    fprintf(fid,'\n');
end
fclose(fid);
fprintf('Outputted heading and %d rows to %s\n',length(good_idx),datafile);
