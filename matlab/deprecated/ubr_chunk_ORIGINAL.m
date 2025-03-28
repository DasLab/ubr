function all_chunk_idx = ubr_chunk( filedir, sequence_file, chunk_size );
% prepare_chunks( filedir, sequence_file, chunk_size );
%
% Prepare directories CHUNKS/0001, CHUNKS/0002, etc.
%  which subdivide large sequence file & ubr output .txt files 
%  into chunks.
%
% Inputs:
%  filedir       = Directory name with merged UBR output (files like
%       RTB002.muts.txt, RTB002.coverage.txt, etc.)
%  sequence_file = FASTA file describing all sequences in library. Ideally with tab
%       delimited headers containing id, title, author.
%  chunk_size    = [integer] number of sequences per chunk (Default 10000)
%
% Output:
%  all_chunk_idx = cell of num_chunks arrays, each of which contains
%                   integer indices (chunk_idx) of which sequences ended up
%                   in which chunk.
%
% (C) R. Das, Stanford University and HHMI, 2023

if ~exist( 'chunk_size','var'); chunk_size = 10000; end;

%% Read in sequences
fasta = fastaread( sequence_file );

%% create subdirs
Nseq = length(fasta);
num_chunks = floor(Nseq/chunk_size);
chunk_dir = 'CHUNKS';
if ~exist(chunk_dir,'dir'); mkdir(chunk_dir); end;
for i = 1:num_chunks
    subdir = sprintf('%s/%04d',chunk_dir,i);
    if ~exist(subdir,'dir'); mkdir(subdir); end;
    subdir = sprintf('%s/%04d/raw_counts/',chunk_dir,i);
    if ~exist(subdir,'dir'); mkdir(subdir); end;
end

%% Chunk the sequences
all_chunk_idx = prepare_sequence_chunks(fasta,sequence_file,chunk_dir,chunk_size,num_chunks);

%% figure out tags for data files
x = dir([filedir,'/*.muts.txt*']);
for i = 1:length(x); tags{i} = strrep(strrep(x(i).name,'.muts.txt',''),'.gz',''); end
tags = unique(tags,'stable');

%% Do the chunking of data files
prepare_data_chunks(filedir,chunk_dir,chunk_size,num_chunks,tags,'muts','');
prepare_data_chunks(filedir,chunk_dir,chunk_size,num_chunks,tags,'coverage','');

mut_types = {'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del'};
for k = 1:length(mut_types)
    mut_type = mut_types{k};
    prepare_data_chunks(filedir,chunk_dir,chunk_size,num_chunks,tags,mut_type,'raw_counts/');
end


fprintf('\nCreate %d chunk directories with up to %d lines each, in CHUNKS.\n',num_chunks, chunk_size);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prepare_data_chunks(filedir,chunk_dir,chunk_size,num_chunks,tags,filetype,subdir);
for i = 1:length(tags)
    tag = tags{i};
    split_completed = 1;
    for n = 1:num_chunks
        filename = sprintf('%s/%04d/%s%s.%s.txt',chunk_dir,n,subdir,tag,filetype);
        if ~exist(filename,'file'); split_completed = 0; break; end;
    end
    if ~split_completed; % do the split
        filename = sprintf('%s.%s.txt',tag, filetype);
        original_file = sprintf('%s/%s%s',filedir,subdir,filename);
        if ~exist(original_file,'file')
            original_file_gz = [original_file,'.gz'];
            assert( exist(original_file_gz,'file'));
            command = sprintf('gunzip %s');
            fprintf('%s\n',command);
            system(command);
        end
        command = sprintf('split -d -l %d -a 4 %s %s/%s.',chunk_size,original_file,chunk_dir,filename);
        fprintf('%s\n',command);
        system(command);
        for n = 1:num_chunks
            src_file = sprintf('%s/%s.%s.txt.%04d',chunk_dir,tag,filetype,n-1);
            dest_file = sprintf('%s/%04d/%s%s.%s.txt',chunk_dir,n,subdir,tag,filetype);
            command = sprintf('mv %s %s',src_file,dest_file);
            fprintf('%s\n',command);
            system(command);            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function all_chunk_idx = prepare_sequence_chunks(fasta,sequence_file,chunk_dir,chunk_size,num_chunks);
Nseq = length(fasta);
all_chunk_idx = {};
for n = 1:num_chunks
    [~,seq_file_basename,ext] = fileparts(sequence_file);
    filename = sprintf('%s/%04d/%s%s',chunk_dir,n,seq_file_basename,ext);
    if exist(filename,'file'); delete(filename); end;
    chunk_idx_min = (n-1)*chunk_size + 1;
    chunk_idx_max = min(n*chunk_size,Nseq);
    chunk_idx = [chunk_idx_min:chunk_idx_max];
    all_chunk_idx{n} = chunk_idx;
    fastawrite(filename,fasta(chunk_idx));
    fprintf('Wrote %d sequences to %s\n',length(chunk_idx),filename);
end

