function make_library_heat_map( r_norm, good_idx, structure_map, headers, BLANK_OUT5, BLANK_OUT3, tags_conditions);
% make_library_heat_map( r_norm, good_idx, structure_map, headers, BLANK_OUT5, BLANK_OUT3, tags_conditions);
% make_library_heat_map( d, good_idx );
%
% Inputs
%  r_norm = [Ndesign x Nres x Nconditions] Reactivity matrix, normalized.
%             If Nconditions = 1, will make heatmap side by side with
%                   structure map
%             If Nconditions > 1, will tile the heatmaps
%  good_idx = [list of integers] index of designs for which to show heatmap
%  structure_map = [Ndesign x Nres] 0/1 map
%            of paired/unpaired for each predicted structure. Give [] to
%            not show it.
%  headers = cell of Ndesign strings describing each design (titles for
%  plot)
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues 
%  tags_conditions = cell of Nconditions strings with display tags for each
%                          condition
%
% (C) R. Das, HHMI/Stanford, 2023

if isstruct(r_norm); d = r_norm; r_norm = d.r_norm; end;
if ~exist( 'good_idx','var') | isempty(good_idx), good_idx = [1:size(r_norm,1)]; end;
if ~exist( 'structure_map','var' ) & exist('d','var'); structure_map = d.structure_map; end;
if ~exist( 'structure_map','var' ); structure_map = []; end
if ~exist( 'headers','var' ) & exist('d','var'); headers = d.headers; end
if ~exist( 'BLANK_OUT5','var' ) & exist('d','var'); BLANK_OUT5 = d.BLANK_OUT5; end
if ~exist( 'BLANK_OUT3','var' ) & exist('d','var'); BLANK_OUT3 = d.BLANK_OUT3; end
if ~exist( 'tags_conditions','var' ) & exist('d','var'); tags_conditions = d.conditions; end

Nconditions = size(r_norm,3);
N = size(r_norm,2);
Nseq = length(good_idx);
Nplots = Nconditions;
if ~isempty(structure_map) Nplots = Nconditions+1; end;

set(gcf,'Position',[100 100 max(50+Nconditions*150,800) 200+0.6*Nseq]); clf;
set(gcf,'color','white')

plot_width = 0.65/Nplots;
for i = 1:Nconditions
    %subplot(1,Nplots,i)
    %set(gca,'Position',[0.5 0.1 0.2 0.8]);
    axes('Position',[0.3+(i-1)*plot_width, 0.1, plot_width, 0.8]);

    imagesc(r_norm(good_idx,:,i),[-2 2])
    colormap([0.7 0.7 0.7; redwhiteblue(-1,1)])
    set(gca, 'TickLabelInterpreter','none','fontweight','bold','ticklength',[0.01,0.25],'tickdir','out' );
    if i == 1
        for n = 1:length(good_idx); header = headers{good_idx(n)}; labels{n} = header(1:min(50,length(header)));end;
        stride = 1; 
        if length(good_idx)>100; stride = floor(length(good_idx)/100); end;
        which_idx =[1:stride:length(good_idx)];
        set(gca,'ytick',which_idx,'yticklabel', labels(which_idx));
    else
        set(gca,'ytick',[]);
    end
    if length(good_idx)< 100; make_lines_horizontal; end
    xlim([BLANK_OUT5 N-BLANK_OUT3])
    title( 'SHAPE data')
    if exist( 'tags_conditions', 'var');
        if Nconditions>4
            h = title( strsplit(tags_conditions{i},'_'),'interp','none' );
        else
            h = title( tags_conditions{i},'interp','none' );
        end
    end;
    xlabel('Position');
    box off
end

if ~isempty(structure_map)
    axes('Position',[0.3+Nconditions*plot_width 0.1 plot_width 0.8])
    structure_map(:,1:BLANK_OUT5) = 0;
    structure_map(:,end-BLANK_OUT3+1:end) = 0;
    imagesc(0.02*structure_map(good_idx,:),[-0.05 0.05])
    colormap([0.7 0.7 0.7; redwhiteblue(-0.1,0.1)])
    set(gca,'ytick',[1:length(good_idx)],'yticklabel', [],'fontweight','bold','ticklength',[0.01,0.25],'tickdir','out')
    set(gca,'yticklabel',[])
    xlim([BLANK_OUT5 N-BLANK_OUT3])
    xlabel('Position');
    %colorbar()
    if length(good_idx)< 100; make_lines_horizontal; end
    title( 'Predicted (unpaired)')
    box off
end