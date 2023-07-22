function output_labels_csv(datafile,good_idx,r_norm,ids,sequences,BLANK_OUT5,BLANK_OUT3,experiment_type,dataset_name, r_norm_err);
% output_labels_csv(datafile,good_idx,r_norm,ids,sequences,BLANK_OUT5,BLANK_OUT3,experiment_type,dataset_name,r_norm_err);
%
%
if isempty(good_idx); good_idx = [1:size(r_norm,1)]; end;
if size(good_idx,2)>size(good_idx,1); good_idx = good_idx'; end;

fid = fopen(datafile,'w');
fprintf(fid,'eterna_id,sequence');

if exist( 'experiment_type','var'); fprintf(fid,',experiment_type'); end
if exist( 'dataset_name','var'); fprintf(fid,',dataset_name'); end
Nres = length(sequences{1});
for k = 1:Nres;
    fprintf(fid,',reactivity_%04d',k);
end
if exist('r_norm_err') & ~isempty(r_norm_err)
    for k = 1:Nres;
        fprintf(fid,',reactivity_error_%04d',k);
    end
end
fprintf(fid,'\n');

for i = good_idx'
    fprintf(fid,'%d,%s', ids(i),sequences{i});
    if exist( 'experiment_type','var'); fprintf(fid,',%s',experiment_type); end
    if exist( 'dataset_name','var'); fprintf(fid,',%s',dataset_name); end
    for k = 1:Nres;
        val = r_norm(i,k);
        if (k <= BLANK_OUT5 | k > (Nres - BLANK_OUT3 ) ) val = NaN; end;
        fprintf(fid,',%7.3f',val);
    end
    if exist('r_norm_err') & ~isempty(r_norm_err)
        for k = 1:Nres;
            val = r_norm_err(i,k);
            if (k <= BLANK_OUT5 | k > (Nres - BLANK_OUT3 ) ) val = NaN; end;
            fprintf(fid,',%7.3f',val);
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);
fprintf('Outputted heading and %d rows to %s\n',length(good_idx),datafile);
