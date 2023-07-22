function [r,r_err,f,f_err,coverage,signal_to_noise,r_nomod] = get_reactivity_OLD(m,c,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3,nomod_error_scalefactor);
% [r,r_err,f,f_err,signal_to_noise] = get_reactivity(m,c,shape_nomod_idx,BLANK_OUT5,BLANK_OUT3);
%
% Inputs
%  m = [Ndesigns x Nres x Ntags] mutation counts  (Ntags = Number of different
%                  experimental conditions)
%  c = [Ndesigns x Nres x Ntags] coverage counts 
%  shape_nomod_idx = cell of Nmodconditions pairs of indices into experimental conditions 
%            in order [shape_idx,nomod_idx]. Ex: {[3,1],[4,2]}.
%  BLANK_OUT5 = Ignore this number of 5' residues in signal_to_noise
%                  calculation.
%  BLANK_OUT3 = Ignore this number of 3' residues in signal_to_noise
%                  calculation. 
%
% Outputs
%  r     = [Ndesign x Nres x Nmodconditions] Reactivity matrix, as fraction at each position that leads to mutation, background subtracted.
%  r_err = [Ndesign x Nres x Nmodconditions] Error of reactivity matrix, as fraction at each position that leads to mutation.
%  f     = [Ndesign x Nres x Ntags] Fraction at each position that leads to mutation.
%  f_err = [Ndesign x Nres x Ntags] Error ofraction at each position that leads to mutation.
%  coverage = [Ndesign x Ntags] Total reads for each design under each
%          experimental condition
%  signal_to_noise = [Ndesign x Nmodconditions] Estimated signal-to-noise for each
%  design under each experimental modification condition. [Mean of signal
%  divided by mean of error across design]
%  r_nomod = [Ndesign x Nres x Nmodconditions] No mod values subtracted out to get r
% 
% (C) R. Das, HHMI/Stanford University 2023.

if ~exist('nomod_error_scalefactor','var') nomod_error_scalefactor = 0; end;

f = m./c;
f_err = sqrt(2*m)./c; % RNAframework typically spits out 2, not 1, for each paired read!

if length(shape_nomod_idx)==0; 
    shape_nomod_idx = {};
    for n = 1:size(m,3); shape_nomod_idx{n} = [n]; end;
end;
r = [];
N = size(f,2);
good_seqpos = [(BLANK_OUT5+1):(N-BLANK_OUT3)];
for i = 1:length(shape_nomod_idx)
    if length(shape_nomod_idx{i}) == 2 & shape_nomod_idx{i}(2)>0
        shape_idx = shape_nomod_idx{i}(1);
        nomod_idx = shape_nomod_idx{i}(2);
        r(:,:,i) = f(:,:,shape_idx)-f(:,:,nomod_idx);
        r_err(:,:,i) = sqrt( f_err(:,:,shape_idx).^2 + f_err(:,:,nomod_idx).^2 + (nomod_error_scalefactor*f(:,:,nomod_idx)).^2);
        r_nomod(:,:,i) = f(:,:,nomod_idx);
    else
        shape_idx = shape_nomod_idx{i};
        r(:,:,i) = f(:,:,shape_idx);
        r_err(:,:,i) = f_err(:,:,shape_idx);
        r_nomod(:,:,i) = 0*f(:,:,shape_idx);
    end
    for m = 1:size(r,1)
        signal_to_noise(m,i) = ubr_estimate_signal_to_noise_ratio(r(m,good_seqpos,i)',r_err(m,good_seqpos,i)'); %mean(abs(r(:,good_seqpos,i)),2)./mean(r_err(:,good_seqpos,i),2);
    end
end
coverage = squeeze(max(c,[],2));
