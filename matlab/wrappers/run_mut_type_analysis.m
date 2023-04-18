function run_mut_type_analysis( m,c,rc,tags,conditions,BLANK_OUT5, BLANK_OUT3)

set(figure(7),'color','white','pos',[ 517   182   806   685])
clf
mut_types = {'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del'};
mut_rate_matrix = [];
for i = 1:length(tags)
    if length(tags)<= 10; subplot( length(tags),1,i); end;
    mut_count_matrix = squeeze(sum(rc(:,:,:,i)));
    coverage_matrix = repmat(squeeze(sum(c(:,:,i),1)),length(mut_types),1)';
    mut_rate_matrix(:,:,i) = mut_count_matrix./coverage_matrix;
    imagesc(mut_rate_matrix(:,:,i)',[0 0.01])
    set(gca,'ytick',[1:length(mut_types)],'YTickLabels',mut_types,'tickdir','out')
    title( conditions{i},'interp','none');
    if length(tags)> 10; pause; end;
    colorbar();
end
xlabel( 'Position')

%% Mutation rates (by MutType)
set(figure(6),'color','white','pos',[1000        1109         613         228]);
clf

mutpos = [BLANK_OUT5:(size(m,2)-BLANK_OUT3)];
mean_mut_rate_by_type = squeeze(mean(mut_rate_matrix(mutpos,:,:),1));
imagesc( mean_mut_rate_by_type',[0 0.005]);
set(gca,'xtick',[1:length(mut_types)],'xticklabel',mut_types,'ytick',[1:length(conditions)],'YTickLabel',conditions,'tickdir','out',...
    'TickLabelInterpreter','none');
colorbar();
make_lines_horizontal;
xlabel( 'Mutation Type' )

%% Table summary
strictmuts = sum(mean_mut_rate_by_type(1:12,:),1);
inserts = mean_mut_rate_by_type(13,:);
dels = mean_mut_rate_by_type(14,:);
strictmuts_dels = sum(mean_mut_rate_by_type([1:12,14],:),1);

rfcount_mut_rate_profiles = squeeze(sum(m,1)./sum(c,1));
rfcount_mut_rate = mean(rfcount_mut_rate_profiles(mutpos,:),1);
fprintf('\n\n')
for i = 1:length(tags)
    fprintf( 'strictmuts:%f ins:%f del:%f strictmuts+dels:%f rf-count:%f %s %s\n',...
        strictmuts(i),inserts(i),dels(i),strictmuts_dels(i),rfcount_mut_rate(i),tags{i},conditions{i})
end