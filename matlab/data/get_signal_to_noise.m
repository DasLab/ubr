function  signal_to_noise = get_signal_to_noise(r,r_err,BLANK_OUT5,BLANK_OUT3, sequences);
% signal_to_noise = get_signal_to_noise(r,r_err,BLANK_OUT5,BLANK_OUT3, sequences);
%
%
% Inputs
%  r     = [Ndesign x Nres x Nmodconditions] Reactivity matrix, as fraction at each position that leads to mutation, background subtracted.
%  r_err = [Ndesign x Nres x Nmodconditions] Error of reactivity matrix, as fraction at each position that leads to mutation.
%  BLANK_OUT5 = Ignore this number of 5' residues in signal_to_noise
%                  calculation.
%  BLANK_OUT3 = Ignore this number of 3' residues in signal_to_noise
%                  calculation.
%  sequences = (cell of Ndesigns strings) RNA sequences. If not provided or
%                if empty, will not spread out deletions at homopolymer
% Output
%  signal_to_noise = [Ndesign x Nmodconditions] Estimated signal-to-noise for each
%             design under each experimental modification condition. [Mean of signal
%             divided by mean of error across design]
%
% (C) R. Das, Stanford & HHMI, 2024

for i = 1:size(r,3)
    for m = 1:size(r,1)
        if exist('sequences','var') & ~isempty(sequences), N = length(sequences{m}); else N = size(r,2); end;
        good_seqpos = [(BLANK_OUT5+1):(N-BLANK_OUT3)];
        signal_to_noise(m,i) = ubr_estimate_signal_to_noise_ratio(r(m,good_seqpos,i)',r_err(m,good_seqpos,i)');
    end
end
