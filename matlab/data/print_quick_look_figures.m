function print_quick_look_figures( filedir )
% print_quick_look_figures()
% print_quick_look_figures( filedir )
%
% Go through figures prepared by quick_look_ubr() and print .png's.
% 
% Inputs [Optional]
%  filedir = directory into which to place the Figures/ subdirectory
%              (default is current directory)
%
% (C) R. Das, HHMI/Stanford University 2023.
if ~exist( 'filedir','var'); filedir = './'; end;
figdir = sprintf('%s/Figures',filedir);
if ~exist(figdir,'dir'); mkdir(figdir); end;
for i = 1:8
    toggle_to_figure(i);
    name = get(gcf,'Name');
    number = get(gcf,'Number');
    if isempty(get(gcf,'Children')); close(number); continue; end;
    output_figure_name = sprintf('%s/Figure%d_%s.png',figdir,number,cleanup(name));
    fprintf('Creating: %s\n',output_figure_name);
    print(output_figure_name,'-dpng');
end

function name = cleanup(name);
name = strrep(name,',','');
name = strrep(name,'/','_');
name = strrep(name,' ','_');
name = strrep(name,'(','');
name = strrep(name,')','');
name = strrep(name,sprintf('\t'),' ');
name = strip(name);
