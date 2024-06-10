function t = output_ubr_stats_summary( d )
% output_ubr_stats_summary( d )
%
% (C) R. Das, HHMI & Stanford, 2024

[~,dirname] = fileparts(pwd);
if ~strcmp(d.filedir,'./') dirname = fileparts(d.filedir); end;
fprintf('\n%s\n',dirname)
fprintf('Statistics over %d sequences with length %d:\n',size(d.r_norm,1),size(d.r_norm,2));
padlen = max(50,max( cellfun(@length,d.conditions)));
fprintf( '%s %9s %8s %6s %8s %7s %7s %7s\n',pad('Condition',padlen,'left'),'reads','mean','median','mn2md','s2n','mnreact','fracSN1' );
for i = 1:length(d.conditions) %7s
    vals = d.r_norm(:,(1+d.BLANK_OUT5):(size(d.r_norm,2)-d.BLANK_OUT3),i);
    conditions{i} = d.conditions{i};
    reads(i) = sum(d.reads(:,i));
    meanreads(i) = nanmean(double(d.reads(:,i)));
    medianreads(i) = floor(nanmedian(double(d.reads(:,i))));
    mn2med(i) = nanmean(double(d.reads(:,i)))/nanmedian(double(d.reads(:,i)));
    meanSN(i) = nanmean(d.signal_to_noise(:,i));
    meanreact(i) = nanmean(vals(:))*d.norm_val(:,i);
    fracSN1(i) = sum(d.signal_to_noise(:,i)>1.0)/size(d.signal_to_noise,1);
    fprintf( '%s %9d %8.1f %6d %8.3f %7.3f %7.5f %7.3f\n',...
        pad(d.conditions{i},padlen,'left'),...
        reads(i) , ...
        meanreads(i),...
        medianreads(i),...
        mn2med(i),...
        meanSN(i),...
        meanreact(i),...
        fracSN1(i) );
end

conditions = conditions';
reads = reads';
meanreads = meanreads';
medianreads = medianreads';
mn2med = mn2med';
meanSN = meanSN';
meanreact = meanreact';
fracSN1 = fracSN1';
t = table(conditions,reads,meanreads,medianreads,mn2med,meanSN,meanreact);

