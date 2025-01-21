function output_reads(d,tag)
%  output_reads(d,tag)
%
% Output coverage and reads vectors from UBR workspace.
%
% Inputs
%  d   = UBR matlab workspace with fields 'coverage' and 'reads'
%  tag = tag to use in output files 
%
% (C) R. Das, HHMI/Stanford University 2025.

if nargin<2; help(mfilename); return; end;
   
if ~exist('data_sets','dir'); mkdir( 'data_sets'); end;

t = table();
t.header = d.headers';
t.sequence = d.sequences';
for i = 1:length(d.tags);
    t.(d.tags{i}) = d.coverage(:,i);
end
outfile = sprintf('data_sets/%s_coverage.csv',tag);
writetable(t,outfile);
fprintf('Outputted %d sequences for %d columns to %s.\n',height(t),size(d.coverage,2),outfile);


t = table();
t.header = d.headers';
t.sequence = d.sequences';
for i = 1:length(d.conditions);
    t.(d.conditions{i}) = d.reads(:,i);
end
outfile = sprintf('data_sets/%s_reads.csv',tag);
writetable(t,outfile);
fprintf('Outputted %d sequences for %d columns to %s.\n',height(t),size(d.reads,2),outfile);





