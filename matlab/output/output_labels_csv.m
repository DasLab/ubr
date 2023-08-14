function output_labels_csv(datafile,good_idx,r_norm,ids,sequences,BLANK_OUT5,BLANK_OUT3,experiment_type,dataset_name, r_norm_err, reads, SN_filter);
% output_labels_csv(datafile,good_idx,r_norm,ids,sequences,BLANK_OUT5,BLANK_OUT3,experiment_type,dataset_name [,r_norm_err, reads, SN_filters]);
%
% Create .csv file of UBR for simple i/o to machine learning efforst.
%
% Inputs
%  datafile = name of data file
%  good_idx = list of integers, indices in r_norm for which to output data
%  r_norm   = [Ndesigns x Nres] normalized reactivity data
%  ids      = array of Ndesigns integers (or cell of Ndesigns) strings with ids
%  sequence = cell of Ndesigns strings with RNA sequences.
%  BLANK_OUT5 = [integer] number of residues to blank (NaN) at 5' end.
%  BLANK_OUT3 = [integer] number of residues to blank (NaN) at 3' end.
%  experiment_type = chemical probing type, e.g. '1M7_MaP',
%         'DMS_MaP','2A3_MaP'
%  dataset_name = arbitrary human readable name that will allow trackback
%          to original experimental files and analysis, e.g., RHEt1_Novoseq_PK50_Twist_multimemberclasses_train
%  r_norm_err = [Optional] [Ndesigns x Nres] errors on normalized
%               reactivity data
%  reads = [Optional] number of reads contributing to profile. 
%  SN_filter = [Optional] array of Ndesigns logicals (booleans) -- whether
%                  design passes S/N filter.
%
% (C) R. Das, Stanford University & HHMI, 2023

if isempty(good_idx); good_idx = [1:size(r_norm,1)]; end;
if size(good_idx,2)>size(good_idx,1); good_idx = good_idx'; end;

% in original use case, ids were integers for Eterna ID, but 
% generalize to strings (e.g., MD5 hashes of sequences) 
id_label = 'id';
if isnumeric( ids )
    id_strings = {};
    for i = 1:length(ids); id_strings{i} = sprintf('%d',ids(i)); end;
    ids = id_strings;
    id_label = 'eterna_id';
end

[filepath,filename] = fileparts(datafile);
if ~exist(filepath,'dir'); mkdir(filepath); end;

fid = fopen(datafile,'w');
fprintf(fid,'%s,sequence',id_label);

if exist( 'experiment_type','var'); fprintf(fid,',experiment_type'); end
if exist( 'dataset_name','var'); fprintf(fid,',dataset_name'); end

if exist( 'reads','var') & length(reads)>0; fprintf(fid,',reads'); end
if exist( 'SN_filter','var') & length(SN_filter)>0; fprintf(fid,',SN_filter'); end


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
    fprintf(fid,'%s,%s', ids{i},sequences{i});
    if exist( 'experiment_type','var'); fprintf(fid,',%s',experiment_type); end
    if exist( 'dataset_name','var'); fprintf(fid,',%s',dataset_name); end
    if exist( 'reads','var') & length(reads)>0; fprintf(fid,',%d',reads(i)); end
    if exist( 'SN_filter','var') & length(SN_filter)>0; fprintf(fid,',%d',SN_filter(i)); end
    for k = 1:Nres;
        val = r_norm(i,k);
        if (k <= BLANK_OUT5 | k > (Nres - BLANK_OUT3 ) ) val = NaN; end;
        fprintf(fid,',%.3f',val);
    end
    if exist('r_norm_err') & ~isempty(r_norm_err)
        for k = 1:Nres;
            val = r_norm_err(i,k);
            if (k <= BLANK_OUT5 | k > (Nres - BLANK_OUT3 ) ) val = NaN; end;
            fprintf(fid,',%.3f',val);
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);
fprintf('Outputted heading and %d rows to %s\n',length(good_idx),datafile);
