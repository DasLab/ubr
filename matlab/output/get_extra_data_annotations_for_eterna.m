function extra_data_annotations = get_extra_data_annotations_for_eterna( ids, titles, authors, extra_info_structs )
% extra_data_annotations = get_extra_data_annotations_for_eterna( ids, titles, authors, extra_info_structs )
%
% Inputs
%  ids = (list of Ndesigns numbers) Eterna ids 
%  titles = (cell of Ndesigns strings) titles 
%  authors = (cell of Ndesigns strings)  Eterna authors 
%  extra_info_structs = Optional (cell of structs) all fields will be reformatted as
%                   annotation strings
%
% Output
%  extra_data_annotations = (cell of Ndesigns cells) cells of tag strings
%
%
% (C) R. Das, 2023, HHMI & Stanford University..
%
if ~exist( 'extra_info_structs','var'); extra_info_structs = []; end;
if nargin == 0; help( mfilename ); return; end;

Ndesigns = length(ids);
assert( length(ids)  == Ndesigns );
assert( length(titles)  == Ndesigns );
assert( length(authors) == Ndesigns );

extra_data_annotations = {};
for i = 1:Ndesigns
    extra_data_annotation = {};
    extra_data_annotation = [extra_data_annotation, {sprintf('Eterna:id:%d',ids(i))} ];
    extra_data_annotation = [extra_data_annotation, {sprintf('Eterna:design_name:%s',titles{i})} ];
    extra_data_annotation = [extra_data_annotation, {sprintf('Eterna:author:%s',authors{i})} ];

    if i <= length(extra_info_structs);    
        extra_data_annotation = [extra_data_annotation, convert_struct_to_annotations( extra_info_structs{i}, 'Eterna' ) ];
    end

    extra_data_annotations{i} = extra_data_annotation;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

