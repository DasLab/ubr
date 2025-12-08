function [m,c,rc,tags] = read_ubr_output( filedir, tags, read_raw_counts, IGNORE_HDF5, seq_range );
% [m,c] = read_ubr_output( filedir, tags, read_raw_counts, IGNORE_HDF5, seq_range );
% Read in data from ubr (ultraplex-bowtie2-RNAframework pipeline)
%
% Inputs
%  filedir = Directory name with RTB002.hdf5, etc. (or, if legacy, 
%              RTB002.muts.txt, RTB002.coverage.txt, etc.)
%  tags = cell of strings with tags {'RTB002', 'RTB004',...}; 
%        If [], will infer from names of .hdf5 (or .muts.txt) files
%  read_raw_counts = read in raw count files (default 0 for dirs 
%                        with .muts.txt).
%  IGNORE_HDF5 = ignore .hdf5 (new format) in favor of .muts.txt (legacy
%                format) Default: 0.
%  seq_range    =  [start_idx, end_idx] two integers with range of sequences
%                   to read in. If oneinteger is specified, assumed to be maximum number of 
%                   sequences to read in. Default [], read all sequences.
% Outputs
%  m = [Ndesigns x Nres x Ntags] mutation counts 
%  c = [Ndesigns x Nres x Ntags] coverage counts 
%  rc = [Ndesigns x Nres x Nmuttypes x Ntags] coverage counts, divided by
%            mutation type (comes out of RNA-framework with -orc flag)
%  tags = tags, could have been inferred from files in directory
%
% (C) R. Das, HHMI/Stanford University 2023.

m = []; c = []; rc = [];

if ~exist( 'tags', 'var') tags = []; end;
if ~exist( 'read_raw_counts','var') read_raw_counts = 0; end;
if ~exist( 'IGNORE_HDF5', 'var') IGNORE_HDF5 = 0; end;
if ~exist( 'seq_range', 'var') seq_range = []; end;
if length(seq_range)==1; seq_range = [1 seq_range(1)]; end;

if ~exist( filedir, 'dir');
    fprintf('Could not find directory %s\n',filedir);
    return;
end

x = dir([filedir,'/*.hdf5']);
hdf5_format = (length(x)>0);
if isempty(tags)
    if ~hdf5_format | IGNORE_HDF5
        % legacy -- used integers to save space.
        m = uint32.empty();
        c = uint32.empty();
        rc = uint32.empty();
        x = dir([filedir,'/*.muts.txt*']);
        for i = 1:length(x); tags{i} = strrep(strrep(x(i).name,'.muts.txt',''),'.gz',''); end
        hdf5_format = 0;
    else
        for i = 1:length(x); tags{i} = strrep(x(i).name,'.hdf5',''); end
    end
    tags = unique(tags,'stable');
end


nts = 'ACGT';
mut_types = {'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del'};

if hdf5_format    
    for n = 1:length(tags)
        hdf5_file = [filedir,'/',tags{n},'.hdf5'];
        fprintf('Reading from ... %s\n',hdf5_file)
        [muts,ins] = load_hdf5_file(hdf5_file, tags{n}, seq_range );
        c(:,:,n) = squeeze(sum(sum(muts,3),4));
        for q = 1:12
            rc(:,:,q,n) = muts(:,:,find(nts==mut_types{q}(1)),find(nts==mut_types{q}(2)));
        end
        rc(:,:,strcmp(mut_types,'del'),n) = squeeze(sum(muts(:,:,:,5),3));
        m(:,:,n) = sum(rc(:,:,~strcmp(mut_types,'ins'),n),3);
        if read_raw_counts
            rc(:,:,strcmp(mut_types,'ins'),n) = squeeze(sum(sum(ins,3),4));
        end
    end
    return;
end

% if not early return from HDF5 block, assume legacy format.
for i = 1:length(tags)
    tag = tags{i};
    fprintf('Reading from ... %s%s\n',filedir,tag)
    m(:,:,i) = read_table_file( [filedir,'/',tag,'.muts.txt']);
    c(:,:,i) = read_table_file( [filedir,'/',tag,'.coverage.txt']);
end

if read_raw_counts
    rc = zeros(size(m,1),size(m,2),length(mut_types),size(m,3),'uint32');
    for i = 1:length(tags)
        tag = tags{i};
        fprintf('Reading raw counts from ... %s%s\n',filedir,tag)
        for k = 1:length(mut_types)
            mut_type = mut_types{k};
            rc(:,:,k,i) = read_table_file( [filedir,'/raw_counts/',tag,'.',mut_type,'.txt']);
        end
    end
end

if ~isempty(seq_range)
    chunk = [seq_range(1):seq_range(2)];
    m = m(chunk,:,:);
    c = c(chunk,:,:);
    rc = rc(chunk,:,:,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [muts,ins] = load_hdf5_file(hdf5_file, tag, seq_range);
hdf5_name = h5info(hdf5_file).Groups(1).Name;
if endsWith(hdf5_name,tag) % legacy cmuts format
    muts = get_hdf5_data(hdf5_file, [hdf5_name,'/mutations'], seq_range);
    ins  = get_hdf5_data(hdf5_file, [hdf5_name,'/insertions'], seq_range);
    return;
end

%if newer cmuts format, recursively figure out name of group
group = h5info(hdf5_file).Groups(1);
while ~isempty(group.Groups)
    group = group.Groups(1);
end
dataset_name = [group.Name,'/',group.Datasets.Name];
counts = get_hdf5_data( hdf5_file, dataset_name, seq_range);
muts = counts(:,:,:,1:5);
ins  = counts(:,:,:,6);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = get_hdf5_data( hdf5_file, dataset_name, seq_range )
hdf5_size = h5info(hdf5_file,dataset_name).Dataspace.Size;
Ndim = length(hdf5_size);
Nseqs = hdf5_size(end);
if ~isempty(seq_range) 
    seq_range(2) = min( seq_range(2), Nseqs );
    start = [ones(1,Ndim-1), seq_range(1)];
    count = [hdf5_size(1:end-1), seq_range(2)-seq_range(1)+1];
    d = h5read(hdf5_file,dataset_name,start,count); 
else
    d = h5read(hdf5_file,dataset_name);
end
d = permute( d, [Ndim:-1:1]);


