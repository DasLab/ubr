function d = read_ubr_hdf5( filename );
% d = read_ubr_hdf5( filename );
%
% Reads HDF5 file into MATLAB workspace.
%
% Inputs:
%  filename = name of  HDF5-formatted file.
%
% Output:
%  d        = MATLAB workspace 
%
% (C) R. Das, Stanford University & HHMI 2024

d = struct();
ds = h5info(filename).Datasets;
for i = 1:length(ds)
    fieldName = ds(i).Name;
    % if strcmp(fieldName,'Attributes'); 
    %     for j = 1:length( ds(i).Attributes )
    %         attrName = ds(i).Attributes(j).Name;
    %         d = setfield(d,attrName,h5readatt(filename,'/Attributes',attrName));
    %     end
    %     continue; 
    % end;
    data = h5read(filename,['/',fieldName]);
    d = setfield(d,fieldName, permute(data,length(size(data)):-1:1) );
end

% may be empty:
% ds = h5info(filename).Attributes;
% for i = 1:length(ds)
%     fieldName = ds(i).Name;
%     d = setfield(d,fieldName, h5readatt(filename,'/',fieldName));
% end




