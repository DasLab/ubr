function  view_sublibraries(d, view_tags, print_stuff )
% view_sublibraries(d, view_tags, print_stuff )
%
% Inputs
%  d = output data struct from quick_look_ubr
%  view_tags = cell of strings which will be searched for in header
%               to define each heatmap
%  print_stuff = print out .png's into Heatmaps/ folder. [default 0]
%
% (C) R. Das, Stanford University & HHMI, 2023

if ~exist('print_stuff','var'); print_stuff = 0; end;
set(figure(9),'name','Sublibrary heatmap');
for n = 1:length(view_tags)
    view_tag = view_tags{n};
    extra_tag = '';

    idx = find(contains(d.headers, view_tag));
    if length(idx)>5000
        idx = idx(1:5000);
        extra_tag='_first5k';
    end

    fprintf('Found tag %s in %d sequences with min_idx %d and max_idx %d\n',view_tag,length(idx),min(idx),max(idx));
    make_library_heat_map( d.r_norm, idx, d.structure_map, d.headers, d.BLANK_OUT5, d.BLANK_OUT3, d.conditions);
    a = get(gcf,'Children');
    for i = 1:length(a)
        a(i).Title.String = [{tag};a(i).Title.String];
    end
    if ~print_stuff;  
        pause;
    else
        figdir = 'Heatmaps';
        if ~exist(figdir,'dir'); mkdir(figdir); end;
        output_figure_name = sprintf('%s/Sublibrary_%s.png',figdir,cleanup([view_tag,extra_tag]));
        fprintf('Creating: %s\n',output_figure_name);
        print(output_figure_name,'-dpng');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function name = cleanup(name);
name = strip(name,'_');
name = strrep(name,',','');
name = strrep(name,'/','_');
name = strrep(name,' ','_');
name = strrep(name,'(','');
name = strrep(name,')','');
name = strrep(name,sprintf('\t'),' ');
name = strip(name);


