function output_labels_csv(datafile,good_idx,r_norm,ids,sequences,BLANK_OUT5,BLANK_OUT3,experiment_type,dataset_name, r_norm_err, reads, SN_filter, signal_to_noise);
% output_labels_csv(datafile,good_idx,r_norm,ids,sequences,BLANK_OUT5,BLANK_OUT3,experiment_type,dataset_name [,r_norm_err, reads, SN_filters, signal_to_noise]);
%
% Create .csv file of UBR for simple i/o to machine learning efforts.
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
%  signal_to_noise = [Optional] numerical signal to noise estimate
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

t = table();

t.(id_label) = ids(good_idx)';
t.sequence = sequences(good_idx)';
if exist( 'experiment_type','var'); t.experiment_type = repmat(experiment_type,length(good_idx),1); end
if exist( 'dataset_name','var'); t.dataset_name = repmat(dataset_name,length(good_idx),1); end
if exist( 'reads','var') & length(reads)>0; t.reads = reads(good_idx);  end
if exist( 'signal_to_noise','var') & length(signal_to_noise)>0; t.signal_to_noise = strtrim(cellstr(num2str(signal_to_noise(good_idx),'%.3f')));  end
if exist( 'SN_filter','var') & length(SN_filter)>0; t.SN_filter = SN_filter(good_idx);  end

Nres = length(sequences{1});
r_norm = single(r_norm);
r_norm(:,1:BLANK_OUT5) = NaN;
r_norm(:,(Nres-BLANK_OUT3+1):Nres) = NaN;
for k = 1:Nres;
    label = sprintf('reactivity_%04d',k);
    t.(label) = strtrim(cellstr(num2str(r_norm(good_idx,k),'%.3f')));
end
if exist('r_norm_err') & ~isempty(r_norm_err)
    r_norm_err = single(r_norm_err);
    r_norm_err(:,1:BLANK_OUT5) = NaN;
    r_norm_err(:,(Nres-BLANK_OUT3+1):Nres) = NaN;
    for k = 1:Nres;
        label = sprintf('reactivity_error_%04d',k);
        t.(label) = strtrim(cellstr(num2str(r_norm_err(good_idx,k),'%.3f')));
    end
end
writetable(t,datafile);
strip_nan(datafile, 0);
fprintf('Outputted heading and %d rows to %s\n',length(good_idx),datafile);
