function m = read_table_file( filename );
% read_table_file( filename );
%
% Input
%   filename = txt or txt.gz with counts as integers
% Output
%   m = table in matlab format.
%
% (C) R. Das, Stanford University, HHMI 2023-25
gzip_file = [filename,'.gz'];
if ~exist(filename,'file')
    %fprintf( 'Unzipping %s... \n',gzip_file);
    if ~exist(gzip_file,'file') fprintf( '\n\nDid not find %s or %s! \n\n',filename,gzip_file); end;
    assert(exist(gzip_file,'file'));
    gunzip( gzip_file );
end

t = readtable(filename);
m = table2array(t,'unt32');

if exist(gzip_file,'file')
    delete(filename); % to save space.
end
