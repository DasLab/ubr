function good_idx = figure_out_idx_for_normalization( total_coverage, MAX_NORM_SEQS );
% good_idx = figure_out_idx_for_normalization( total_coverage );
%
% Input:
%  total_coverage = [Nseq] integers, list of coverages for sequences
%  MAX_NORM_SEQS  = integer, maximum number of indices to return, based on
%                     sorting by coverage (Default 1000)
%
% Output:
%  good_idx = list of indices of sequences to use for normalization,
%              based on heuristic and COVERAGE_CUTOFF of 10000.
%
% (C) R. Das, Stanford & HHMI

if ~exist('MAX_NORM_SEQS','var'); MAX_NORM_SEQS = 1000; end

COVERAGE_CUTOFF = 10000; % stringent cutoff.
good_idx = find( total_coverage > COVERAGE_CUTOFF );
if length(good_idx) < 10
    COVERAGE_CUTOFF = 1000; % less stringent cutoff.
    good_idx = find( total_coverage > COVERAGE_CUTOFF );
end
if length(good_idx) < 10
    COVERAGE_CUTOFF = 100; % not stringent
    good_idx = find( total_coverage > COVERAGE_CUTOFF );
end
if length(good_idx) < 10
    COVERAGE_CUTOFF = 0; % not stringent
    good_idx = [1:length(total_coverage)];
end

% Too many sequences is leading to out of memory issues downstream!
if length(good_idx) > MAX_NORM_SEQS
     %[~,sortidx] = sort(total_coverage,'descend');
     %good_idx = sortidx(1:MAX_NORM_SEQS);
    stride = floor(length(good_idx)/MAX_NORM_SEQS);
    good_idx = good_idx(stride:stride:end); % just pick out a subset of sequences.
end

fprintf( 'For normalization, using %d sequences with min coverage %d, passing a total coverage cutoff of %d \n',...
    length(good_idx),min(total_coverage(good_idx)),COVERAGE_CUTOFF)
