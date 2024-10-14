function show_mut_type_analysis( mut_rate_matrix,rfcount_mut_rate_profiles,tags,labels,BLANK_OUT5, BLANK_OUT3)
% run_mut_type_analysis( m,c,rc,tags,labels,BLANK_OUT5, BLANK_OUT3)
%
% Inputs
%   m  = [Ndesigns x Nres x Ntags] mutation counts 
%   c  = [Ndesigns x Nres x Ntags] coverage counts 
%   rc = [Ndesigns x Nres x Nmuttypes x Ntags] coverage counts, divided by
%            mutation type (comes out of RNA-framework with -orc flag)
%   tags = {cell of Ntags strings} tags for each condition. Ex:
%                               {'RTB000_1M7','RTB002_nomod'};
%   labels =   {cell of Ntags strings} human-readable labels for each
%            condition. Default: [], which triggers labels = tags.
%   BLANK_OUT5 = Ignore this number of 5' residues in final summary
%   BLANK_OUT3 = Ignore this number of 3' residues in final summary  
%
if ~exist('labels','var') | isempty(labels); labels = tags; end;
if ~exist('BLANK_OUT5','var'); BLANK_OUT5 = 0; end;
if ~exist('BLANK_OUT3','var'); BLANK_OUT3 = 0; end;

toggle_to_figure(8);
set(gcf,'color','white','pos',[ 517   182   806 485+100*length(tags)],'name','Mutation type analysis, position-wise')
clf
mut_types = {'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del'};
for i = 1:length(tags)
    subplot( length(tags),1,i); 
    imagesc(mut_rate_matrix(:,:,i)',[0 0.01])
    set(gca,'ytick',[1:length(mut_types)],'YTickLabels',mut_types,'tickdir','out')
    title( labels{i},'interp','none');
    colorbar();
end
xlabel( 'Position')

%% Mutation rates (by MutType)
toggle_to_figure(7);
set(gcf,'color','white','pos',[1000  1109 613 228],'name','Mutation type analysis, table summary');
clf

mutpos = [1+BLANK_OUT5:(size(mut_rate_matrix,1)-BLANK_OUT3)];
mean_mut_rate_by_type = reshape(mean(mut_rate_matrix(mutpos,:,:),1),size(mut_rate_matrix,2),size(mut_rate_matrix,3));
imagesc( mean_mut_rate_by_type',[0 0.005]);
set(gca,'xtick',[1:length(mut_types)],'xticklabel',mut_types,'ytick',[1:length(labels)],'YTickLabel',labels,'tickdir','out',...
    'TickLabelInterpreter','none');
colorbar();
make_lines_horizontal;
xlabel( 'Mutation Type' )

%% Table summary
% strictmuts = sum(mean_mut_rate_by_type(1:12,:),1);
% inserts = mean_mut_rate_by_type(13,:);
% dels = mean_mut_rate_by_type(14,:);
% strictmuts_dels = sum(mean_mut_rate_by_type([1:12,14],:),1);
% 
% rfcount_mut_rate = mean(rfcount_mut_rate_profiles(mutpos,:),1);
% 
% fprintf('\n\n')
% for i = 1:length(tags)
%     if strcmp(tags{i},labels{i}) labels{i} = ''; end; 
%     fprintf( 'strictmuts:%8.5f ins:%8.5f del:%8.5f strictmuts+dels:%8.5f rf-count:%8.5f  %s %s\n',...
%         strictmuts(i),inserts(i),dels(i),strictmuts_dels(i),rfcount_mut_rate(i),tags{i},labels{i})
% end