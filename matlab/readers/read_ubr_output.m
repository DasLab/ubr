function [m,c,rc,tags] = read_ubr_output( filedir, tags, read_raw_counts );
% [m,c] = read_ubr_output( filedir, tags );
% Read in data from ubr (ultraplex-bowtie2-RNAframework pipeline)
%
% Inputs
%  filedir = Directory name with RTB002.muts.txt, RTB002.coverage.txt, etc.
%  tags = cell of strings with tags {'RTB002', 'RTB004',...}; 
%        If [], will infer from names of .muts.txt files
%  read_raw_counts = read in raw count files (default 0).
%
% Outputs
%  m = [Ndesigns x Nres x Ntags] mutation counts 
%  c = [Ndesigns x Nres x Ntags] coverage counts 
%  rc = [Ndesigns x Nres x Nmuttypes x Ntags] coverage counts, divided by
%            mutation type (comes out of RNA-framework with -orc flag)
%  tags = tags, could have been inferred from files in directory
%
% (C) R. Das, HHMI/Stanford University 2023.
m = uint32.empty();
c = uint32.empty();
rc = uint32.empty();

if ~exist( 'tags', 'var') tags = []; end;
if ~exist( 'read_raw_counts','var') read_raw_counts = 0; end;
if ~exist( filedir, 'dir');
    fprintf('Could not find directory %s\n',filedir);
    return;
end
if isempty(tags)
    x = dir([filedir,'/*.muts.txt*']);
    for i = 1:length(x); tags{i} = strrep(strrep(x(i).name,'.muts.txt',''),'.gz',''); end
    tags = unique(tags,'stable');
end

for i = 1:length(tags)
    tag = tags{i};
    fprintf('Reading from ... %s%s\n',filedir,tag)
    m(:,:,i) = load_file( [filedir,tag,'.muts.txt']);
    c(:,:,i) = load_file( [filedir,tag,'.coverage.txt']);
end

if read_raw_counts
    mut_types = {'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del'};
    rc = zeros(size(m,1),size(m,2),length(mut_types),size(m,3),'uint32');
    for i = 1:length(tags)
        tag = tags{i};
        fprintf('Reading raw counts from ... %s%s\n',filedir,tag)
        for k = 1:length(mut_types)
            mut_type = mut_types{k};
            rc(:,:,k,i) = load_file( [filedir,'raw_counts/',tag,'.',mut_type,'.txt']);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = load_file( filename );
gzip_file = [filename,'.gz'];
if ~exist(filename,'file')
    assert(exist(gzip_file,'file'));
    %fprintf( 'Unzipping %s... \n',gzip_file);
    gunzip( gzip_file );
end

t = readtable(filename);
m = table2array(t,'unt32');

if exist(gzip_file,'file')
    delete(filename); % to save space.
end

