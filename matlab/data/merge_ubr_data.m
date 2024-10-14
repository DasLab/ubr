function d_out = merge_ubr_data( d, cidx )
% d_out = merge_ubr_data( d, cidx )
%
% Merges average reactivities (appropriately weighted by inverse error
%  squared) from different 'conditions' which, e.g., might have 
%  sublibraries rebalanced in different proportions.
%
% Inputs
%  d = MATLAB data structure with fields like r_norm, reads, etc. (as
%        output from quick_look_ubr)
%  cidx = which condition indices to merge  (e.g., [1,2], or can be
%          cell of index sets, e.g., {[1,3],[2,4]})
%
% Output
%  d_out = MATLAB data structure with merged
%
% (C) R. Das, Stanford University, HHMI, 2024

if ~iscell(cidx) cidx_all = {cidx}; 
else cidx_all = cidx; end

d_out = d;
d_out.r_norm = [];
d_out.r_norm_err = [];
d_out.conditions = {};
d_out.reads = [];
d_out.signal_to_noise = [];
d_out.coverage = []; 
d_out.tags = []; 
d_out.shape_nomod_idx = {};
d_out.norm_val = [];

for n = 1:length(cidx_all)
    cidx = cidx_all{n};
    d_out.r_norm(:,:,n) = nansum( d.r_norm(:,:,cidx)./d.r_norm_err(:,:,cidx).^2, 3 )./ nansum(1./d.r_norm_err(:,:,cidx).^2, 3);
    d_out.r_norm_err(:,:,n) = sqrt(1./nansum(1./ d.r_norm_err(:,:,cidx).^2,3));
    d_out.conditions{n} = strjoin(d.conditions(cidx),'_');
    d_out.reads(:,n) = sum(d.reads(:,cidx),2);
    d_out.signal_to_noise(:,n) = get_signal_to_noise(d_out.r_norm(:,:,n),d_out.r_norm_err(:,:,n),d.BLANK_OUT5,d.BLANK_OUT3, d.sequences);
    d_out.norm_val(n) = mean(d.norm_val(cidx));

    sh_idx = []; nm_idx = [];
    for i = cidx
        sh_idx = [sh_idx, d.shape_nomod_idx{i}(1)];
        nm_idx = [nm_idx, d.shape_nomod_idx{i}(2)];
    end
    d_out.coverage = [d_out.coverage, sum(d.coverage(:,sh_idx),2), sum(d.coverage(:,nm_idx),2) ];
    d_out.tags = [d_out.tags, { strjoin(d.tags( sh_idx ),'_'), strjoin(d.tags( nm_idx ),'_')} ];
    d_out.shape_nomod_idx = [d_out.shape_nomod_idx, length(d_out.tags) + [-1,0]];
    if isfield(d,'mut_rate_matrix');
        d_out.mut_rate_matrix(:,:,2*n-1) = mean( d.mut_rate_matrix(:,:,sh_idx), 3 );
        d_out.mut_rate_matrix(:,:,2*n  ) = mean( d.mut_rate_matrix(:,:,nm_idx), 3 );
    end
    if isfield(d,'rfcount_mut_rate_profiles');    
        d_out.rfcount_mut_rate_profiles(:,2*n-1) = mean( d.rfcount_mut_rate_profiles(:,sh_idx), 2 );
        d_out.rfcount_mut_rate_profiles(:,2*n  ) = mean( d.rfcount_mut_rate_profiles(:,nm_idx), 2 );
    end
end
d_out.total_coverage = sum(d_out.coverage,2);

