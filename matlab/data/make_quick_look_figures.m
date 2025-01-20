function make_quick_look_figures( d, f )
% make_quick_look_figures( d, f )

if isfield(d,'f'); f = d.f; end;
if ~exist( 'f','var'); f = []; end;

tic
fprintf('\nCreating figures (provide no_figures in options to skip)...\n')
if ~isempty(f)
    toggle_to_figure(1);
    Ntags = size(f,3);
    set(gcf,'color','white','Position',[100   100  max(50+Ntags*100,800)   526],'name','fraction reacted, first 500 designs')
    Nplot = min(size(f,1),500);
    for i = 1:Ntags
        subplot(1,Ntags,i)
        imagesc(sum(f(1:Nplot,:,i),3),[0 0.1]);
        colormap(1-gray(100));
        title(strsplit(d.tags{i},'_'),'interpreter','none')
        ylim([0.5 Nplot+0.5])
        colorbar('Location','SouthOutside');
    end
end

%% Make plot of signal-to-noise vs. coverage
toggle_to_figure(2);
set(gcf,'color','white','position',[450   991   355   322],'name','S/N vs. reads'); clf
clf
for k = 1:size(d.reads,2)
    semilogx(d.reads(:,k),d.signal_to_noise(:,k),'.'); hold on
    xlabel('Reads (signal+background channels)')
    ylabel('Signal/noise');
end
hold off
legend(d.conditions,'Interpreter','none');

%% Histogram of signal to noise
toggle_to_figure(3);
set(gcf,'color','white','position',[694   960   400   600],'name','Mean S/N')
clf
for i = 1:size(d.r_norm,3)
    subplot(size(d.r_norm,3),1,i);
    s2n = d.signal_to_noise(:,i); s2n( find(isnan(s2n)) ) = 0.0;
    histogram( s2n )
    xlabel( 'Signal/noise' );
    cols = strsplit(what(d.filedir).path,'/');
    dirname = strjoin(cols(end-1:end),'/');
    title( [d.conditions{i},' Mean signal/noise = ',num2str(mean(s2n))],'interp','none');
end

%% Make heat map, up to 500 with high signal to noise 
toggle_to_figure(4);
set(gcf,'color','white','name','first designs with good S/N (up to 500)')
clf
good_idx = find( d.signal_to_noise(:,end)>=1.0 & d.reads(:,end) > 100);
Nplot = min(length(good_idx),500);
good_idx = good_idx(1:Nplot);
make_library_heat_map( d.r_norm, good_idx, d.structure_map, d.headers, d.BLANK_OUT5, d.BLANK_OUT3, d.conditions);

%% Make heat map, up to 10000 with high signal to noise
toggle_to_figure(5);
good_idx = find( d.signal_to_noise(:,end)>=1.0 & d.reads(:,end) > 100);
if length(good_idx)>500
    set(gcf,'color','white','name','first designs with good S/N (up to 10000)')
    clf
    Nplot = min(length(good_idx),10000);
    good_idx = good_idx(1:Nplot);
    make_library_heat_map( d.r_norm, good_idx, d.structure_map, d.headers, d.BLANK_OUT5, d.BLANK_OUT3, d.conditions);
else
   close(5);
end


%% Take a close look at one of the constructs with high apparent signal to noise
toggle_to_figure(6);
if ~isempty(f)
    set(gcf,'color','white','position',[599   477   560   420],'name','Top S/N design')
    clf
    [~,idx] = max(sum(d.signal_to_noise,2));
    for i = 1:length(d.shape_nomod_idx)
        subplot(length(d.shape_nomod_idx),1,i);
        cidx = d.shape_nomod_idx{i}; %(end:-1:1);
        plot( squeeze(f(idx,:,cidx)),'LineWidth',2 )
        h = legend( d.tags{cidx});
        set(h,'interp','none');
    end
    h=title(d.headers(idx));
    set(h,'interpreter','none')
end

%% Look through each of the conditions - mutational profiles (mean over designs)
if ~isempty(d.mut_rate_matrix)
    show_mut_type_analysis(d.mut_rate_matrix,d.rfcount_mut_rate_profiles,d.tags,d.tags,d.BLANK_OUT5, d.BLANK_OUT3);
end
toc
