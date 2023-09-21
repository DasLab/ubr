function rdat = output_rdat_from_ubr_analysis( filename, name, good_idx, r_norm, r_norm_err, sequences, structures, BLANK_OUT5, BLANK_OUT3, comments, annotations, extra_data_annotations, reads );
%
% rdat = output_rdat_from_ubr_analysis( filename, name, good_idx, r_norm, r_norm_err, sequences, structures, BLANK_OUT5, BLANK_OUT3, comments, annotations, extra_data_annotations, reads );
%
%  filename   = name of file to output to disk, e.g., 'RHEt1_PK50_Twist.rdat'
%  name       = will show up in NAME field of RDAT -- short description of
%                RNA library, e.g., 'PK50'
%  good_idx        = indices of which designs to output. Give [] for all. 
%  r_norm     = [Ndesigns x Nres] Matrix of data
%  r_norm_err = [Ndesigns x Nres] Matrix of data errors
%  sequences  = {cell of strings of length Nres} sequences of all RNA's
%  structures  = {cell of strings of length Nres} structures of all RNA's.
%                    If empty, no structures will be output in
%                    DATA_ANNOTATION
%  comments   = comments to put in RDAT
%  annotations = cell of tags -- any annotations to include in ANNOTATION. Make sure to
%                 include tags like {'modifier:1M7','chemical:HEPES:50mM(pH8.0)','chemical:MgCl2:10mM','temperature:24C'}
%  extra_data_annotations = {cell of Ndesign cells} extra annotations to put
%                   in DATA_ANNOTATION. Give [] or {} if no extra
%                   annotations.
%  reads = [Ndesigns] integers for number of reads contributing to profile.
%
% Adapted from an analogous output wrapper from MAPseeker.
%
% (C) R. Das, 2023, HHMI & Stanford University..
%

if nargin == 0; help( mfilename ); return; end;
if ~exist('extra_data_annotations','var') ; extra_data_annotations = {}; end;
if ~exist( 'output_workspace_to_rdat_file', 'file' )
  fprintf( '\nCannot find function output_workspace_to_rdat_file() ...\nWill not output RDAT: %s.\n', filename );
  fprintf( 'Install RDATkit and include   rdatkit/matlab_scripts/   in MATLAB path.\n\n' );
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reformat data to cell.
reactivity     = [];
reactivity_err = [];
data_annotations = {};
rdat = [];

Ndesigns = size( r_norm, 1 );
if length( sequences ) ~= Ndesigns;   fprintf( 'Number of sequences does not match Ndesigns in r_norm\n' ); return; end

Nres = size(r_norm,2);
%if length( sequences{1} ) ~= Nres;   fprintf( 'Sequence length does not match Nres in r_norm\n' ); return; end

if isempty(good_idx); good_idx = [1 : size( r_norm, 1)]'; end;

which_res = [(1+BLANK_OUT5):(Nres-BLANK_OUT3)];

count = 0;
for j = squeeze(good_idx)'

    count = count + 1;
    reactivity(:,count) = r_norm(j,which_res)';
    reactivity_err(:,count) = r_norm_err(j,which_res)';

    SN_ratio(count)   = ubr_estimate_signal_to_noise_ratio( reactivity(:,count), reactivity_err(:,count) );
    SN_classification = classify_signal_to_noise_ratio( SN_ratio(count) );

    modifier = '';
    ID = '';
    project_name = '';
    design_name = '';
    tag_cols = {};
  
    data_annotation = {};
    data_annotation  = [data_annotation,  ['sequence:',sequences{j}] ];
    if exist( 'structures','var') & length(structures)>=j; data_annotation  = [data_annotation,  ['structure:',structures{j}] ]; end;
    if exist( 'reads','var') & length(reads)>=j; data_annotation = [data_annotation, sprintf('reads:%d',reads(j))]; end

    data_annotation = [data_annotation, ['signal_to_noise:',SN_classification,':',num2str(SN_ratio(count),'%8.3f') ] ];
    if length(extra_data_annotations) > 0;  data_annotation  = [data_annotation,  extra_data_annotations{j} ]; end;
    data_annotations{count} = data_annotation;
    
end

sequence  = repmat('X',1,Nres);
structure = repmat('.',1,Nres);
offset = 0; 
seqpos = which_res;

dirname = fileparts( filename );
if length(dirname)>0 & ~exist(dirname,'dir'); mkdir( dirname); end;

rdat = output_workspace_to_rdat_file( filename, name, sequence, offset, seqpos, reactivity, ...
				      structure, ...
				      annotations, data_annotations, ...
				      reactivity_err, ...
				      [],[],[], comments );

fprintf('Outputted %d profiles into %s.\n\n',size(reactivity,2),filename);

