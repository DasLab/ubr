function data = h5load(filename)
% Specify your HDF5 file name 
%filename = 'yourfile.hdf5'; 
% Get info about the HDF5 file 
info = h5info(filename); 
% Preallocate a cell array to hold your data 
data = struct(); 


for k = 1:length(info.Datasets)
    dataset_name = info.Datasets(k).Name;
    dataset_path =  ['/',dataset_name];
    % Read the dataset and store it in the corresponding struct
    data.(matlab.lang.makeValidName(dataset_name)) = h5read(filename, dataset_path);
end


% Loop through each group 
for i = 1:length(info.Groups) 
    group_name = info.Groups(i).Name; 
    data.(matlab.lang.makeValidName(strrep(group_name,'/',''))) = struct(); % Create a struct for the group 

    for k = 1:length(info.Groups(i).Datasets)
        dataset_name = info.Groups(i).Datasets(k).Name;
        dataset_path = fullfile(group_name, dataset_name);
        % Read the dataset and store it in the corresponding struct
        data.(matlab.lang.makeValidName(strrep(group_name,'/',''))).(matlab.lang.makeValidName(dataset_name)) = h5read(filename, dataset_path);
    end


    % Loop through each subgroup within the group 
    for j = 1:length(info.Groups(i).Groups) 
        subgroup_name = info.Groups(i).Groups(j).Name; 
        data.(matlab.lang.makeValidName(strrep(group_name,'/',''))).(matlab.lang.makeValidName(subgroup_name)) = struct(); % Create a struct for the subgroup          
        % Loop through each dataset within the subgroup 
        for k = 1:length(info.Groups(i).Groups(j).Datasets) 
            dataset_name = info.Groups(i).Groups(j).Datasets(k).Name; 
            dataset_path = fullfile(subgroup_name, dataset_name); 
            % Read the dataset and store it in the corresponding struct 
            data.(matlab.lang.makeValidName(strrep(group_name,'/',''))).(matlab.lang.makeValidName(subgroup_name)).(matlab.lang.makeValidName(dataset_name)) = h5read(filename, dataset_path); 
        end 
    end 
end 