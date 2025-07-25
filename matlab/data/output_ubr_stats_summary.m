function t = output_ubr_stats_summary( d )
% output_ubr_stats_summary( d )
%
% Input:
%
%  d = workspace output by quick_look_ubr()
%
% (C) R. Das, HHMI & Stanford, 2024

if strcmp(d.filedir,'./') d.filedir = pwd(); end;
dir_parts = strsplit(d.filedir,'/');
dirname = strjoin(dir_parts( max(end-2,1): end ),'/');
fprintf('\n%s\n',dirname)
tot_reads = sum(d.reads(:));
fprintf('Statistics over %d sequences with length %d and total reads %d (%s reads):\n',size(d.r_norm,1),size(d.r_norm,2),tot_reads,human_readable(tot_reads));
padlen = max(10,max( cellfun(@length,d.conditions)));

if isfield(d,'mut_rate_matrix')

    mutpos = [1+d.BLANK_OUT5:(size(d.mut_rate_matrix,1)-d.BLANK_OUT3)];
    mean_mut_rate_by_type = reshape(nanmean(d.mut_rate_matrix(mutpos,:,:),1),size(d.mut_rate_matrix,2),size(d.mut_rate_matrix,3));
    strictmuts = nansum(mean_mut_rate_by_type(1:12,:),1);
    inserts = mean_mut_rate_by_type(13,:);
    dels = mean_mut_rate_by_type(14,:);
    strictmuts_dels = nansum(mean_mut_rate_by_type([1:12,14],:),1);
    rfcount_mut_rate = nanmean(d.rfcount_mut_rate_profiles(mutpos,:),1);
 
    padlen = max(10,max( cellfun(@length,d.tags)));
    fprintf( '%s %11s %9s %8s %6s | %7s %7s %7s %7s %7s\n',pad('Expt tag',padlen,'left'),'reads','mean','median','mn2md',...
        'mut','ins','del','mut+del','rfcount');
    for i = 1:length(d.tags) %7s
        reads(i) = sum(d.coverage(:,i));
        meanreads(i) = nanmean(double(d.coverage(:,i)));
        medianreads(i) = floor(nanmedian(double(d.coverage(:,i))));
        mn2med(i) = nanmean(double(d.coverage(:,i)))/nanmedian(double(d.coverage(:,i)));
        fprintf( '%s %11d %9.1f %8d %6.3f | %7.4f %7.4f %7.4f %7.4f %7.4f\n',...
            pad(d.tags{i},padlen,'left'),...
            uint64(reads(i)), ...
            meanreads(i),...
            medianreads(i),...
            mn2med(i),...
            strictmuts(i),...
            inserts(i),...
            dels(i),...
            strictmuts_dels(i),...
            rfcount_mut_rate(i));
    end
end

fprintf('\n\n');
fprintf( '%s %11s %8s %9s %6s | %7s %7s | %7s %6s %7s\n',pad('Condition',padlen,'left'),'reads','mean','median','mn2md','normval','mnreact','mednorm','s2n','fracSN1' );
reads  = []; meanreads = []; medianreads = []; mn2med = [];
norm_val = d.norm_val;
for i = 1:length(d.conditions)
    conditions{i} = d.conditions{i};
    reads(i) = sum(d.reads(:,i));
    meanreads(i) = nanmean(double(d.reads(:,i)));
    medianreads(i) = floor(nanmedian(double(d.reads(:,i))));
    mn2med(i) = nanmean(double(d.reads(:,i)))/nanmedian(double(d.reads(:,i)));

    vals = d.r_norm(:,(1+d.BLANK_OUT5):(size(d.r_norm,2)-d.BLANK_OUT3),i);
    meanreact(i) = nanmean(vals(:))*d.norm_val(:,i);
 
    idx = find(d.signal_to_noise(:,i)>1.0);
    medreact_norm(i) = nanmedian( reshape(vals(idx,:),1,[]) );
    meanSN(i) = nanmean(d.signal_to_noise(:,i));
    fracSN1(i) = sum(d.signal_to_noise(:,i)>1.0)/size(d.signal_to_noise,1);
    fprintf( '%s %11d %9.1f %8d %6.3f | %7.5f %7.5f | %7.5f %6.3f %7.3f\n',...
        pad(d.conditions{i},padlen,'left'),...
        uint64(reads(i)) , ...
        meanreads(i),...
        medianreads(i),...
        mn2med(i),...
        norm_val(i),...
        meanreact(i),...
        medreact_norm(i),...
        meanSN(i),...
        fracSN1(i) );
end

conditions = conditions';
reads = reads';
meanreads = meanreads';
medianreads = medianreads';
mn2med = mn2med';
norm_val = norm_val';
meanreact = meanreact';
medreact_norm = medreact_norm';
meanSN = meanSN';
fracSN1 = fracSN1';
t = table(conditions,reads,meanreads,medianreads,mn2med,norm_val,meanreact,medreact_norm,meanSN,fracSN1);



function r = human_readable(reads)
s = {'', 'k', 'M', 'B', 'T'}; 
if reads == 0; r = '0'; return; end
e = floor(log(reads)/log(1000));
r = sprintf(['%.1f',s{e+1}], (reads/(1000^floor(e))));

