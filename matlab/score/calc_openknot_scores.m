function [eterna_scores,crossed_pair_scores, best_cc, structures] = calc_openknot_scores( r_norm, mfe_structures, good_idx, mfe_tags, BLANK_OUT3, BLANK_OUT5, headers);
% [eterna_scores,crossed_pair_scores, best_cc, structures] = calc_openknot_scores( r_norm, mfe_structures, good_idx, mfe_tags, BLANK_OUT3, BLANK_OUT5, headers);
%
%
%  r_norm = [Ndesign x Nres] Reactivity matrix, normalized.
%               normalized to go from 0 to 1 (~90th percentile)
%  mfe_structures = [Npackages x Ndesign] cell of cell of strings of predicted structures
%  good_idx = [Nidx] index of designs for which to show heatmap
%  mfe_tags = cell of string, name of each package
%  BLANK_OUT5 = gray out this number of 5' residues
%  BLANK_OUT3 = gray out this number of 3' residues
% headers = cell of Ndesign strings describing each design (titles for
%  plot)
%
% (C) Rhiju Das, Stanford, HHMI, 2023

eterna_scores = [];
crossed_pair_scores = [];
best_cc = [];
structures = {};

for idx = good_idx'
    [eterna_score, crossed_pair_score, cc, structure ] = calc_openknot_score( r_norm, mfe_structures, idx, mfe_tags, BLANK_OUT3, BLANK_OUT5, headers );

    eterna_scores = [eterna_scores, eterna_score];
    crossed_pair_scores = [crossed_pair_scores, crossed_pair_score];
    best_cc = [best_cc, cc];
    structures = [structures, {structure}];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eterna_classic_score, crossed_pair_score, best_cc, structure ] = calc_openknot_score( r_norm, mfe_structures, idx, mfe_tags, BLANK_OUT3, BLANK_OUT5, headers )

data = r_norm(idx,:,1);

% remove singlet base pairs
for n = 1:length(mfe_tags)
    structure = mfe_structures{n}{idx};
    mfe_structures{n}{idx} = sanitize_structure( structure, 1);
end
mfe_structure_map = get_mfe_structure_map( mfe_structures, idx );


[cc,cc_sort_idx] = get_corr_coeff( r_norm, mfe_structure_map, idx, mfe_tags, BLANK_OUT3, BLANK_OUT5);

% eterna score classic:
for n = 1:length(mfe_tags)
    pred = squeeze(mfe_structure_map(idx,:,n));
    eterna_classic_scores(n) = calc_eterna_score_classic( data, pred, BLANK_OUT5, BLANK_OUT3);
    crossed_pair_scores(n) = calc_crossed_pair_score(data, mfe_structures{n}{idx}, BLANK_OUT5, BLANK_OUT3 );
end
plot( cc, eterna_classic_scores,'.')
xlabel('Correlation coefficient');
ylabel( 'Eterna score (classic)');
title(headers{idx},'interp','none');

% best_cc = cc(cc_sort_idx(end));
% best_model_idx = find( best_cc == cc);
% 
% eterna_classic_score = max(eterna_classic_scores(best_model_idx));
% best_model_idx = find( best_cc == cc & eterna_classic_score == eterna_classic_scores);
% 
% crossed_pair_score = max(crossed_pair_scores(best_model_idx));
% best_model_idx = find( best_cc == cc & eterna_classic_score == eterna_classic_scores & crossed_pair_scores == crossed_pair_score);

% use eterna score instead of correlation coefficient...
eterna_classic_score = max(eterna_classic_scores);
best_model_idx = find( eterna_classic_score == eterna_classic_scores);

crossed_pair_score = max(crossed_pair_scores(best_model_idx));
best_model_idx = find( eterna_classic_score == eterna_classic_scores & crossed_pair_scores == crossed_pair_score);

structure = mfe_structures{best_model_idx(1)}{idx};
best_cc = cc(best_model_idx(1));

% try ensembling any models that are close, based on correlation
% coefficient
% cc_band = 0.05;
% best_model_idx = find( cc >= (best_cc-cc_band) );
% [~,sortidx] = sort(cc(best_model_idx));
% best_model_idx = best_model_idx( sortidx(end:-1:1) );
% eterna_classic_score = mean( eterna_classic_scores( best_model_idx) );
% crossed_pair_score   = mean( crossed_pair_scores( best_model_idx) );


% try ensembling any models that are close, based on Eterna Score
eterna_score_band = 5;
best_model_idx = find( eterna_classic_scores >= (eterna_classic_score - eterna_score_band) );
[~,sortidx] = sort(eterna_classic_scores(best_model_idx));
best_model_idx = best_model_idx( sortidx(end:-1:1) );
%eterna_classic_score = mean( eterna_classic_scores( best_model_idx) );
crossed_pair_score   = mean( crossed_pair_scores( best_model_idx) );


for n = best_model_idx;
    fprintf( '%30s %s\n',mfe_tags{n},mfe_structures{n}{idx});
end

fprintf( '\nBest model for %s.  Correlation coefficient: %4.2f  Eterna score classic: %5.2f  Crossed pair score: %5.2f\n',...
    headers{idx},best_cc,...
    eterna_classic_score,...
    crossed_pair_score);

%for n = best_model_idx;
%    fprintf( '%30s %s\n',mfe_tags{n},mfe_structures{n}{idx});
%end





