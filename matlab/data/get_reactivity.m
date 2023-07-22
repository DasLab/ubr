function [r,r_err,f,f_err,coverage,signal_to_noise,r_nomod] = get_reactivity(rc,c,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3,sequences)
% [r,r_err,f,f_err,coverage,signal_to_noise,r_nomod] = get_reactivity(rc,c,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3,sequences)
%
% Main data processing routine for UBR. Includes step where deletions,
% which are assumed to be collected as counts placed at 3' end of same-nt
% stretches, are spread out through same-nt stretch based on counts in 'strict' 
% mutations (i.e. positions that are mutations but not indels).
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
mut_types = {'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del'};
Nmuttypes = length(mut_types);
if length(size(rc))~=4;
    m = rc;
    [r,r_err,f,f_err,coverage,signal_to_noise,r_nomod] = get_reactivity_OLD(m,c,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3);
    fprintf(['\n\nWARNING! Get_reactivity requires raw counts (''rc'') as first input:\n[Ndesigns x Nres x Ntags x Nmuttypes] raw counts (Ntags = Number of different',...
        'experimental conditions).\nLast dimension should have size',...
        '14 and correspond to: AC,AG,AT,CA,CG,CT,GA,GC,GT,TA,TC,TG,ins,del.\n\n',...
        'For now, running get_reactivity_OLD, but this will be deprecated in the future.\n\n']);
    return;
else
    assert( size(rc,3) == Nmuttypes );
end

c_reshape = repmat(reshape(c,size(c,1),size(c,2),1,size(c,3)),1,1,Nmuttypes,1);
f = rc./c_reshape;
f_err = sqrt(2*rc)./c_reshape; % RNAframework typically spits out 2, not 1, for each paired read!

coverage = squeeze(max(c,[],2));

%% Background subtract, though do keep track of nomod.
for i = 1:length(shape_nomod_idx)
    ci = shape_nomod_idx{i};
    if length(ci)>1
        rsub(:,:,:,i)     = f(:,:,:,ci(1)) - f(:,:,:,ci(2));
        rsub_err(:,:,:,i) = sqrt(f_err(:,:,:,ci(1)).^2 + f_err(:,:,:,ci(2)).^2);
        r_nomod(:,:,:,i)  = f(:,:,:,ci(2));
    else
        assert(length(ci)==1);
        rsub(:,:,:,i)     = f(:,:,:,ci(1));
        rsub_err(:,:,:,i) = f_err(:,:,:,ci(1));
        r_nomod(:,:,:,i)  = 0*rsub(:,:,:,i);
    end
end
rsub_strictmut     = squeeze(sum(rsub(:,:,1:12,:),3));
rsub_strictmut_err = squeeze(sqrt(sum(rsub_err(:,:,1:12,:).^2,3)));

rsub_del     = squeeze(rsub(:,:,14,:));
rsub_del_err = squeeze(rsub_err(:,:,14,:));

r_nomod = squeeze(sum(r_nomod(:,:,[1:12 14],:),3));

if ~exist('sequences') | isempty(sequences)
    fprintf(['\n\nWARNING! Get_reactivity requires sequences to spread out deletions at same-nt stretches.\nPlease provide sequences as last input.\nFor now, proceeding without spread out of deletion.\n\n']);
    [r, r_err, signal_to_noise] = get_r_from_strictmut_del( rsub_strictmut, rsub_del, rsub_strictmut_err, rsub_del_err, BLANK_OUT5, BLANK_OUT3);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare new deletion profiles, 'spread' out at same-nt stretches on mutation spectrum
assert(length(sequences)==size(rsub,1));
rsub_del_spread = rsub_del*0;
rsub_del_spread_err = rsub_del_err*0;
for i = 1:size(rsub_del,1);
    for m = 1:size(rsub_del,3);
        sequence = sequences{i};
        del_profile = squeeze(rsub_del(i,:,m));
        del_profile_err = squeeze(rsub_del_err(i,:,m));
        del_profile2 = del_profile*0;
        del_profile2_err = del_profile_err*0;
        strictmut_profile = squeeze(rsub_strictmut(i,:,m));
        strictmut_profile_err = squeeze(rsub_strictmut_err(i,:,m));
        pos = 1;
        for n = 1:length(sequence)
            if sequence(n)~=sequence(pos) % end of a same-nt stretch
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
        if all(~isnan(del_profile)) assert( abs(sum(del_profile2)-sum(del_profile )) < 1e-3 ); end;
        rsub_del_spread(i,:,m) = del_profile2;
        rsub_del_spread_err(i,:,m) = del_profile2_err;
    end
end

[r, r_err, signal_to_noise] = get_r_from_strictmut_del( rsub_strictmut, rsub_del_spread, rsub_strictmut_err, rsub_del_spread_err, BLANK_OUT5, BLANK_OUT3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   [r, r_err, signal_to_noise] = get_r_from_strictmut_del( rsub_strictmut, rsub_del, rsub_strictmut_err, rsub_del_err, BLANK_OUT5, BLANK_OUT3);

r     = rsub_strictmut + rsub_del;
r_err = sqrt(rsub_strictmut_err.^2 + rsub_del_err.^2);
signal_to_noise = get_signal_to_noise(r,r_err,BLANK_OUT5,BLANK_OUT3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  signal_to_noise = get_signal_to_noise(r,r_err,BLANK_OUT5,BLANK_OUT3);
N = size(r,2);
good_seqpos = [(BLANK_OUT5+1):(N-BLANK_OUT3)];
for i = 1:size(r,3)
    for m = 1:size(r,1)
        signal_to_noise(m,i) = ubr_estimate_signal_to_noise_ratio(r(m,good_seqpos,i)',r_err(m,good_seqpos,i)');
    end
end
