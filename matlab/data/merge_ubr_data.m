function d_out = merge_ubr_data( d, cidx, RENORM )
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
%  RENORM = average reactivities after "un"-normalizing; and then
%            renormalize at end (Default 1). Set to 0 for legacy
%            calculation which just averages norm_val.
%
% Output
%  d_out = MATLAB data structure with merged
%
% (C) R. Das, Stanford University, HHMI, 2024

if ~iscell(cidx) cidx_all = {cidx}; 
else cidx_all = cidx; end
if ~exist('RENORM','var'); RENORM = 1; end;

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

if isfield(d,'coverage_matrix');
    d_out.coverage_matrix = [];
    d_out.mut_rate_matrix = [];
    d_out.rfcount_mut_rate_profiles = [];
    for i = 1:size(d.coverage_matrix,2)
        coverage_matrix_tmp = repmat(reshape( d.coverage_matrix(:,i),[],1), 1, size(d.mut_rate_matrix,2));
        mut_count(:,:,i) = d.mut_rate_matrix(:,:,i) .* coverage_matrix_tmp;
    end
    rfcount = d.rfcount_mut_rate_profiles .* d.coverage_matrix;
else
    warning('Did not find field coverage_matrix so mut_rate_matrix will not be exact.');
end
if ~RENORM; warning('Not renormalizing, so norm_val will not be exact.'); end; 

for n = 1:length(cidx_all)
    cidx = cidx_all{n};
    if RENORM
        % may need to do these calculations 'in place' to avoid memory issues
        for q = 1:length(cidx)
            r_cidx(:,:,q)     = d.r_norm(:,:,cidx(q))     .* d.norm_val(cidx(q));
             r_cidx_err(:,:,q) = d.r_norm_err(:,:,cidx(q)) .* d.norm_val(cidx(q));
        end
        %r_cidx     = tensorprod(d.r_norm(:,:,cidx),diag(d.norm_val(cidx)),3,1);
        %r_cidx_err = tensorprod(d.r_norm_err(:,:,cidx),diag(d.norm_val(cidx)),3,1);

        % normalization occurs below:
        r(:,:,n) = nansum( r_cidx./r_cidx_err.^2, 3 )./ nansum(1./r_cidx_err.^2, 3);
        r_err(:,:,n) = sqrt(1./nansum(1./r_cidx_err.^2,3));
    else
        % legacy, without "un"-normalizing -- slightly different numbers.
        d_out.r_norm(:,:,n) = nansum( d.r_norm(:,:,cidx)./d.r_norm_err(:,:,cidx).^2, 3 )./ nansum(1./d.r_norm_err(:,:,cidx).^2, 3);
        d_out.r_norm_err(:,:,n) = sqrt(1./nansum(1./ d.r_norm_err(:,:,cidx).^2,3));
    end

    d_out.conditions{n} = strjoin(d.conditions(cidx),'_');
    d_out.reads(:,n) = sum(d.reads(:,cidx),2);
    d_out.norm_val(n) = mean(d.norm_val(cidx));

    sh_idx = []; nm_idx = [];
    for i = cidx
        sh_idx = [sh_idx, d.shape_nomod_idx{i}(1)];
        if length(d.shape_nomod_idx{i})>1;  nm_idx = [nm_idx, d.shape_nomod_idx{i}(2)]; else ; nm_idx = NaN; end;
    end

    if ~any( isnan( nm_idx) ) 
        d_out.coverage = [d_out.coverage, sum(d.coverage(:,sh_idx),2), sum(d.coverage(:,nm_idx),2) ];
        d_out.tags = [d_out.tags, { strjoin(d.tags( sh_idx ),'_'), strjoin(d.tags( nm_idx ),'_')} ];
        d_out.shape_nomod_idx = [d_out.shape_nomod_idx, length(d_out.tags) + [-1,0]];

        if isfield(d,'coverage_matrix')
            d_out.coverage_matrix(:,2*n-1) = sum( d.coverage_matrix(:,sh_idx),2 );
            d_out.coverage_matrix(:,2*n  ) = sum( d.coverage_matrix(:,nm_idx),2 );

            d_out.mut_rate_matrix(:,:,2*n-1) = sum( mut_count(:,:,sh_idx), 3 )./ d_out.coverage_matrix(:,2*n-1) ;
            d_out.mut_rate_matrix(:,:,2*n  ) = sum( mut_count(:,:,nm_idx), 3 )./ d_out.coverage_matrix(:,2*n  ) ;

            d_out.rfcount_mut_rate_profiles(:,2*n-1) = sum( rfcount(:,sh_idx), 2 )./ d_out.coverage_matrix(:,2*n-1) ;
            d_out.rfcount_mut_rate_profiles(:,2*n  ) = sum( rfcount(:,nm_idx), 2 )./ d_out.coverage_matrix(:,2*n  ) ;
        else % legacy -- in workspaces where we did not save coverage_matrix, just take mean.
            if isfield(d,'mut_rate_matrix');
                d_out.mut_rate_matrix(:,:,2*n-1) = mean( d.mut_rate_matrix(:,:,sh_idx), 3 );
                d_out.mut_rate_matrix(:,:,2*n  ) = mean( d.mut_rate_matrix(:,:,nm_idx), 3 );
            end
            if isfield(d,'rfcount_mut_rate_profiles');
                d_out.rfcount_mut_rate_profiles(:,2*n-1) = mean( d.rfcount_mut_rate_profiles(:,sh_idx), 2 );
                d_out.rfcount_mut_rate_profiles(:,2*n  ) = mean( d.rfcount_mut_rate_profiles(:,nm_idx), 2 );
            end
        end
    else
        assert( all(isnan(nm_idx))); % require nomods to be all defined or all undefined
        d_out.coverage = [d_out.coverage, sum(d.coverage(:,sh_idx),2)];
        d_out.tags = [d_out.tags, { strjoin(d.tags( sh_idx ),'_') } ];
        d_out.shape_nomod_idx = [d_out.shape_nomod_idx, length(d_out.tags) ];

        if isfield(d,'coverage_matrix')
            d_out.coverage_matrix(:,n) = sum( d.coverage_matrix(:,sh_idx),2 );
            d_out.mut_rate_matrix(:,:,n) = sum( mut_count(:,:,sh_idx), 3 )./ d_out.coverage_matrix(:,n) ;
            d_out.rfcount_mut_rate_profiles(:,n) = sum( rfcount(:,sh_idx), 2 )./ d_out.coverage_matrix(:,n) ;
        else % legacy -- in workspaces where we did not save coverage_matrix, just take mean.
            if isfield(d,'mut_rate_matrix');
                d_out.mut_rate_matrix(:,:,n) = mean( d.mut_rate_matrix(:,:,sh_idx), 3 );
            end
            if isfield(d,'rfcount_mut_rate_profiles');
                d_out.rfcount_mut_rate_profiles(:,n) = mean( d.rfcount_mut_rate_profiles(:,sh_idx), 2 );
            end
        end
    end


end

d_out.total_coverage = sum(d_out.coverage,2);

if RENORM
    norm_idx = figure_out_idx_for_normalization( d_out.total_coverage );
    if length(norm_idx)>0;
        [d_out.r_norm, d_out.r_norm_err,~,d_out.norm_val] = normalize_reactivity(r,r_err,norm_idx,d_out.BLANK_OUT5, d_out.BLANK_OUT3, d_out.conditions, [], d_out.sequences );
    end
    d_out.norm_idx = norm_idx;
end

for n = 1:length(cidx_all)
    for n = 1:length(cidx_all)
        d_out.signal_to_noise(:,n) = get_signal_to_noise(d_out.r_norm(:,:,n),d_out.r_norm_err(:,:,n),d.BLANK_OUT5,d.BLANK_OUT3, d.sequences);
    end
end




