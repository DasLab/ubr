function prepare_data_chunks(filedir,chunk_dir,chunk_size,all_chunk_idx,tags,filetype,subdir);
% prepare_data_chunks(filedir,chunk_dir,chunk_size,all_chunk_idx,tags,filetype,subdir);
%
% Inputs
%  filedir   = input file dir with *muts.txt files, etc.
%  chunk_dir = directory in which chunks will be placed 
%  chunk_size = size of chunks (to use split); give 0 to use all_chunk_idx.
%  all_chunk_idx = cell of arrays of indices of which sequences go into
%                       which chunk
%  tags = cell of strings like 'RTB000_Marathon_Bicine_3pct_DMS' that are
%          the prefix of the data files.
%  filetype = which file type to run ('muts', 'coverage', 'GA', etc.)
%  subdir   = subdirectory of filedir that holds the file (either "" or
%                 "raw_counts")
%
% (C) R. Das, Stanford University & HHMI, 2023

num_chunks = length(all_chunk_idx);
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
            command = sprintf('gunzip %s',original_file_gz);
            fprintf('%s\n',command);
            system(command);
        end
        if chunk_size > 0
            % use Unix split -- fast but limited to even splits!
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
        else
            t = readtable(original_file);
            %x = load(original_file);
            for n = 1:num_chunks
                filename = sprintf('%s/%04d/%s%s.%s.txt',chunk_dir,n,subdir,tag,filetype);
                fprintf('Creating %s\n',filename);
                writetable(t(all_chunk_idx{n},:),filename,'WriteVariableNames',0);
                %y = x(all_chunk_idx{n},:);
                %save(filename,'y','-ascii');
            end
        end
    end
end


