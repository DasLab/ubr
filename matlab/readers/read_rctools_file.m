function  [m,c, seq_names, sequences] = read_rctools_file( filename );
%[m,c, seq_names, sequences] = read_rctools_file( filename );
% Read files output by rf_rctools from RNAframework
%
% (C) R. Das, HHMI/Stanford University 2023.

fid = fopen( filename );
line = fgetl( fid );
i = 0;
while ~feof(fid)
    i = i + 1;
    seq_names{i} = line;

    line = fgetl( fid );
    sequences{i} = line;

    line = fgetl( fid );
    m(i,:) = str2num( line );

    line = fgetl( fid );
    c(i,:) = str2num( line );

    line = fgetl( fid );
    line = fgetl( fid );
end

fclose( fid );


