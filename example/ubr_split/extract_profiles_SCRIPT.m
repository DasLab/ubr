%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filedir = './';
shape_nomod_idx = { ...
    {'RTB008_Twist_PK50_1M7','RTB012_Twist_PK50_nomod'},...
    {'RTB010_CustomArray_PK50_1M7','RTB014_CustomArray_PK50_nomod'},...
    };
sequence_file = 'pseudoknot50_puzzle_11318423.tsv.RNA_sequences.fa';
structure_csv_file = '';
d = quick_look_ubr(filedir,sequence_file,shape_nomod_idx,structure_csv_file);
save workspace_PK50_example.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output RDAT block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_prefix = 'OpenKnotPilot_PK50_RH_Et1_MiSeq';
name = 'OpenKnotPilot_PK50';
comments = {...
    'from data: Miseq_2023-03-25_RH_Et1_PK50_PK90',...
    'Mutational profiling',...
    };
condition_comments{1} = {'Library ordered from Twist, error-prone PCR'};
condition_comments{2} = {'Library ordered from CustomArray'};

annotations = {'chemical:MgCl2:10mM','temperature:24C','chemical:bicine:150mM(pH8.5)','modifier:1M7','processing:RNAFramework-v2.8.4','processing:spread_deletions','processing:backgroundSubtraction','reverse_transcriptase:SSII_Mn','experimentType:StandardState'};
condition_annotations = {};
%condition_annotations{1} = {'modifier:DMS','reverse_transcriptase:Marathon'};
%condition_annotations{2} = {'modifier:2A3','reverse_transcriptase:SSII_Mn'};

output_rdats_from_ubr_analysis( d, out_prefix, name, [], comments, annotations, condition_comments, condition_annotations );
save workspace_PK50_example.mat;




