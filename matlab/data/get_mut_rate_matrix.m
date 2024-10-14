function [mut_rate_matrix, rfcount_mut_rate_profiles] = get_mut_rate_matrix( m,c,rc )
% [mut_rate_matrix, rfcount_mut_rate_profiles] = get_mut_rate_matrix( m,c,rc )
%
% Inputs
%   m  = [Ndesigns x Nres x Ntags] mutation counts 
%   c  = [Ndesigns x Nres x Ntags] coverage counts 
%   rc = [Ndesigns x Nres x Nmuttypes x Ntags] coverage counts, divided by
%            mutation type (comes out of RNA-framework with -orc flag)
%
% Output
%   mut_rate_matrix = [Nres x Nmuttypes x Ntags] average mutation
%      frequency per residue with each of 14 mutation types, across each tag
%      (separately barcoded condition, including separate nomod and mod).
%   rfcount_mut_rate_profiles = [Nres x Ntags ] average mutation frequency 
%      in rf-count muts.txt file.
%
% (C) R. Das, HHMI, Stanford 2024.

% Keep here for reference:
% mut_types = {'AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG','ins','del'};

mut_rate_matrix = [];
for i = 1:size(c,3)
    mut_count_matrix = reshape(sum(rc(:,:,:,i),1),size(rc,2),size(rc,3));
    coverage_matrix = repmat(reshape(sum(c(:,:,i),1),size(c,2),1),1,size(rc,3));
    mut_rate_matrix(:,:,i) = mut_count_matrix./coverage_matrix;
end

rfcount_mut_rate_profiles = reshape(sum(m,1)./sum(c,1),size(m,2),size(m,3),size(m,4));
