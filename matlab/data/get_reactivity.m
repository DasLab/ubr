function [r,r_err,f,f_err,coverage,signal_to_noise,r_nomod] = get_reactivity(rc,c,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3,sequences,options )
% [r,r_err,f,f_err,coverage,signal_to_noise,r_nomod] = get_reactivity(rc,c,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3,sequences,options)
%
% Main data processing routine for UBR. Includes step where deletions,
% which are assumed to be collected as counts placed at 3' end of same-nt
% stretches, are spread out through same-nt stretch based on counts in 'strict' 
% mutations (i.e. positions that are mutations but not indels).
%
% Note that errors are now defined so that the minimum corresponds to a
%  single count in each channel, analogous to MAPseeker.
%
% Inputs
%  rc = [Ndesigns x Nres x Ntags x Nmuttypes] raw counts (Ntags = Number of different
%                  experimental conditions). Last dimension should have size
%                  14 and correspond to: AC,AG,AT,CA,CG,CT,GA,GC,GT,TA,TC,TG,ins,del.
%                  Deletions are assumed to be collected at right side
%                  (most 3' nucleotide) of same-nt stretches, and will be
%                  spread out.
%  c  = [Ndesigns x Nres x Ntags] coverage counts
%  shape_nomod_idx = cell of Nmodconditions pairs of indices into experimental conditions
%            in order [shape_i,nomod_i]. Ex: {[3,1],[4,2]}.
%  BLANK_OUT5 = Ignore this number of 5' residues in signal_to_noise
%                  calculation.
%  BLANK_OUT3 = Ignore this number of 3' residues in signal_to_noise
%                  calculation.
%  sequences = (cell of Ndesigns strings) RNA sequences. If not provided or
%                if empty, will not spread out deletions at homopolymer
%                stretches.
%  options    = cell of strings specifying possible options. (Default: {})
%                Current option strings are:
%                 'no_spread_deletions': for same-nucleotide stretches,
%                   deletions are ambigous and here, by default, are spread across the
%                   stretch according to the mutation signal. Supply
%                   'no_spread_deletions' to suppress this procedure.
%                 'no_GA': Exclude G to A mutations
%                 'only_GA': Only count G to A mutations (e.g. DMS at N7)
%                 'only_AG': Only count A to G mutations (e.g. deaminase)%                
%
% Outputs
%  r     = [Ndesign x Nres x Nmodconditions] Reactivity matrix, as fraction at each position that leads to mutation, background subtracted.
%  r_err = [Ndesign x Nres x Nmodconditions] Error of reactivity matrix, as fraction at each position that leads to mutation.
%  f     = [Ndesign x Nres x Ntags] Fraction at each position that leads to mutation.
%  f_err = [Ndesign x Nres x Ntags] Error ofraction at each position that leads to mutation.
%  coverage = [Ndesign x Ntags] Total reads for each design under each
%               experimental condition
%  signal_to_noise = [Ndesign x Nmodconditions] Estimated signal-to-noise for each
%             design under each experimental modification condition. [Mean of signal
%             divided by mean of error across design]
%  r_nomod = [Ndesign x Nres x Nmodconditions] No mod values subtracted out to get r
%
% (C) Rhiju Das, Stanford University & HHMI, 2023.
%

spread_deletions = ~any(strcmp(options,'no_spread_deletions'));
no_GA = any(strcmp(options,'no_GA'));
only_GA = any(strcmp(options,'only_GA'));
only_AG = any(strcmp(options,'only_AG'));

mut_types = {'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del'};
Nmuttypes = length(mut_types);
if length(size(rc))~=4 & size(rc,3)~=Nmuttypes;
    m = rc;
    [r,r_err,f,f_err,coverage,signal_to_noise,r_nomod] = get_reactivity_OLD(double(m),double(c),shape_nomod_idx,BLANK_OUT5,BLANK_OUT3);
    fprintf(['\n\nWARNING! Get_reactivity requires raw counts (''rc'') as first input:\n[Ndesigns x Nres x Ntags x Nmuttypes] raw counts (Ntags = Number of different',...
        'experimental conditions).\nLast dimension should have size',...
        '14 and correspond to: AC,AG,AT,CA,CG,CT,GA,GC,GT,TA,TC,TG,ins,del.\n\n',...
        'For now, running get_reactivity_OLD, but this will be deprecated in the future.\n\n']);
    return;
else
    assert( size(rc,3) == Nmuttypes );
end

tic
fprintf('Compiling reactivity...\n')
c_reshape = repmat(reshape(c,size(c,1),size(c,2),1,size(c,3)),1,1,Nmuttypes,1);
frc = single(rc)./single(c_reshape);
frc_err = sqrt(single(rc))./single(c_reshape); 
f_pseudocount_err = 1./single(c); % pseudocount -- will be a floor on the error.

mut_del_idx = [1:12,14];
strictmut_idx = [1:12];
del_idx = [14];
if no_GA 
    mut_del_idx   = [1:6,8:12,14];
    strictmut_idx = [1:6,8:12];
elseif only_GA
    fprintf('Only counting mutations of G to A\n')
    mut_del_idx   = [7];
    strictmut_idx = [7];
    del_idx = [];
    frc = put_nan_except_at_nt(frc,sequences,'G');
elseif only_AG
    fprintf('Only counting mutations of A to G\n')
    mut_del_idx   = [2];
    strictmut_idx = [2];
    del_idx = [];
    frc = put_nan_except_at_nt(frc,sequences,'A');
end

f = reshape(sum(frc(:,:,mut_del_idx,:),3) ,size(frc,1),size(frc,2),size(frc,4));
f_err = reshape(sqrt(sum(frc_err(:,:,mut_del_idx,:).^2,3)) ,size(frc,1),size(frc,2),size(frc,4));

coverage = reshape(max(c,[],2),size(c,1),size(c,3));

%% Background subtract, though do keep track of nomod.
for i = 1:length(shape_nomod_idx)
    ci = shape_nomod_idx{i};
    if length(ci)>1
        rsub(:,:,:,i)     = frc(:,:,:,ci(1)) - frc(:,:,:,ci(2));
        rsub_err(:,:,:,i) = sqrt(frc_err(:,:,:,ci(1)).^2 + frc_err(:,:,:,ci(2)).^2);
        r_nomod(:,:,:,i)  = frc(:,:,:,ci(2));
        rsub_pseudocount_err(:,:,i) = sqrt(f_pseudocount_err(:,:,ci(1)).^2 + f_pseudocount_err(:,:,ci(2)).^2);
    else
        assert(length(ci)==1);
        rsub(:,:,:,i)     = frc(:,:,:,ci(1));
        rsub_err(:,:,:,i) = frc_err(:,:,:,ci(1));
        r_nomod(:,:,:,i)  = 0*rsub(:,:,:,i);
        rsub_pseudocount_err(:,:,i) = f_pseudocount_err(:,:,ci(1));
    end
end
rsub_strictmut     = reshape(sum(rsub(:,:,strictmut_idx,:),3),size(rsub,1),size(rsub,2),size(rsub,4));
rsub_strictmut_err = reshape(sqrt(sum(rsub_err(:,:,strictmut_idx,:).^2,3)),size(rsub,1),size(rsub,2),size(rsub,4));

rsub_del     = reshape(rsub(:,:,del_idx,:), size(rsub,1),size(rsub,2),[]);
rsub_del_err = reshape(rsub_err(:,:,del_idx,:), size(rsub,1),size(rsub,2),[]);

r_nomod = reshape(sum(r_nomod(:,:,mut_del_idx,:),3),size(r_nomod,1),size(r_nomod,2),size(r_nomod,4));
toc

if ~exist('sequences') | isempty(sequences)
    fprintf(['\n\nWARNING! Get_reactivity requires sequences to spread out deletions at same-nt stretches.\nPlease provide sequences as last input.\nFor now, proceeding without spread out of deletion.\n\n']);
    spread_deletions = 0;
end
if ~spread_deletions
    [r, r_err, signal_to_noise] = get_r_from_strictmut_del( rsub_strictmut, rsub_del, rsub_strictmut_err, rsub_del_err, BLANK_OUT5, BLANK_OUT3, rsub_pseudocount_err, sequences);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare new deletion profiles, 'spread' out at same-nt stretches on mutation spectrum
tic
fprintf('Spreading out deletions...\n')
assert(length(sequences)==size(rsub,1));
rsub_del_spread = rsub_del*0;
rsub_del_spread_err = rsub_del_err*0;
for i = 1:size(rsub_del,1);
    if mod(i,1000)==0; fprintf('Doing %d out of %d...\n',i,size(rsub_del,1)); end;
    for m = 1:size(rsub_del,3);
        sequence = sequences{i};
        del_profile = squeeze(rsub_del(i,:,m));
        del_profile_err = squeeze(rsub_del_err(i,:,m));
        del_profile2 = del_profile*0;
        del_profile2_err = del_profile_err*0;
        strictmut_profile = squeeze(rsub_strictmut(i,:,m));
        strictmut_profile_err = squeeze(rsub_strictmut_err(i,:,m));
        pos = 1;
        for n = 1:length(sequence)+1
            if n == length(sequence)+1 | sequence(n)~=sequence(pos) % end of a same-nt stretch
                % Fill in beginning of same-nt stretch with any values that
                % have accumulated there. In principle this should be 0 in
                % rf-count, but some counts end up there, perhaps due to
                % use of collapse_count option in RNA framework.  
                del_profile2(pos:n-2) = del_profile(pos:n-2);
                del_profile2_err(pos:n-2) = del_profile_err(pos:n-2);
                % Now handle the counts on last nucleotide in same-nt
                % stretch.
                % clip negative strict mutation counts to zero:
                weights = max(strictmut_profile(pos:n-1),0);
                weights_err = strictmut_profile_err(pos:n-1);
                reweight = weights/sum(weights);
                % The (1-weights) below comes out of the derivative of
                %  x/(x+y) with respect to x.
                reweight_err = weights_err.*(1-reweight)/sum(weights);
                % uniform weighting if strict mutations have inf/nan/0 values:
                if any(isnan(reweight) | isinf(reweight)) reweight = ones(1,(n-pos))/(n-pos); reweight_err = 0*reweight; end;
                del_profile2(pos:n-1) = del_profile2(pos:n-1) + del_profile(n-1) .* reweight;
                del_profile2_err(pos:n-1) = sqrt( del_profile2_err(pos:n-1).^2 + (del_profile_err(n-1) .* reweight).^2 + ...
                    (del_profile(n-1)*reweight_err).^2);
                pos = n;
            end
        end
        % Important: counts should be conserved:
        if all(~isnan(del_profile)); 
            if( abs(sum(del_profile2)-sum(del_profile )) > 1e-3 )
            %if( abs(sum(del_profile2(1+BLANK_OUT5:end-BLANK_OUT3))-sum(del_profile(1+BLANK_OUT5:end-BLANK_OUT3))) > 1e-3 )
                warning( 'spread out deletions may not have worked for profile %d\n',i);
            end
        end;
        rsub_del_spread(i,:,m) = del_profile2;
        rsub_del_spread_err(i,:,m) = del_profile2_err;
    end
end

[r, r_err, signal_to_noise] = get_r_from_strictmut_del( rsub_strictmut, rsub_del_spread, rsub_strictmut_err, rsub_del_spread_err, BLANK_OUT5, BLANK_OUT3, rsub_pseudocount_err, sequences);

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   [r, r_err, signal_to_noise] = get_r_from_strictmut_del( rsub_strictmut, rsub_del, rsub_strictmut_err, rsub_del_err, BLANK_OUT5, BLANK_OUT3, rsub_pseudocount_err, sequences);
if size(rsub_del,3) == 0; rsub_del = 0*rsub_strictmut; rsub_del_err = 0*rsub_strictmut; end;
r     = rsub_strictmut + rsub_del;
r_err = sqrt(rsub_strictmut_err.^2 + rsub_del_err.^2 + rsub_pseudocount_err.^2);
%r_err = r_err .* sqrt(1-r); % closer to correct formula w/ high reactivity.
signal_to_noise = get_signal_to_noise(r,r_err,BLANK_OUT5,BLANK_OUT3, sequences);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function frc = put_nan_except_at_nt(frc,sequences,nt);
% Provide nt as character like 'A' or 'G'
non_nt = find(char(sequences) ~= nt);
for i = 1:size(frc,3);
    for j = 1:size(frc,4);
        tmp = frc(:,:,i,j);
        tmp(non_nt) = NaN;
        frc(:,:,i,j) = tmp;
    end
end


