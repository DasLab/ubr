function d = read_ubr_hdf5_data( filename, name, fields, seq_range)
d = struct();
no_seq_range_fields = {'tags','conditions','norm_idx','norm_val',...
    'shape_nomod_idx','coverage_matrix','mut_rate_matrix','rfcount_mut_rate_profiles'};
no_seq_range_fields = [no_seq_range_fields, {'authors','titles','id_strings','ids'}];
ds = h5info(filename,['/',name]).Datasets;
for i = 1:length(ds)
    fieldName = ds(i).Name;
    field_tag = ['/',name,'/',fieldName];
    if ~isempty( fields ) & ~any(strcmp(fields,fieldName)) continue; end;
    if ~isempty(seq_range)
        hdf5_size = h5info(filename,field_tag).Dataspace.Size;
        Ndim = length(hdf5_size);
        if all(hdf5_size == 1) | any(strcmp(no_seq_range_fields,fieldName));
            start = ones(1,Ndim);
            count = hdf5_size;
        elseif Ndim == 2 & hdf5_size(2) == 1 % array
            Nseqs = hdf5_size(1);
            seq_range_max = min( seq_range(2), Nseqs );
            start = [seq_range(1), ones(1,Ndim-1)];
            count = [ seq_range_max-seq_range(1)+1, hdf5_size(2:end)];
        else % sequences ,headers, etc.
            Nseqs = hdf5_size(end);
            seq_range_max = min( seq_range(2), Nseqs );
            start = [ones(1,Ndim-1), seq_range(1)];
            count = [hdf5_size(1:end-1), seq_range_max-seq_range(1)+1];
        end
        data = h5read(filename,field_tag,start,count);
    else
        data = h5read(filename,field_tag);
    end
    d = setfield(d,fieldName, permute(data,length(size(data)):-1:1) );
end

