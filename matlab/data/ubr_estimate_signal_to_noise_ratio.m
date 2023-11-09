function  SN_ratio = ubr_estimate_signal_to_noise_ratio( signal, noise );
%   SN_ratio = ubr_estimate_signal_to_noise_ratio( signal, noise );
%
%  Calculate mean(signal) / mean( noise) over all points with noise>0
%   Ignore  first and last such points.
%
% Inputs
%  signal = [Nres x Nprofiles] data
%  noise  = [Nres x Nprofiles] errors on data
%
% Adapted from a function in MAPseeker (CloudLab data analysis)
%
% (C) R. Das, Stanford 2010-

if nargin==0; help( mfilename ); return; end;

all_ratio = 0;
for i = 1:size( signal, 2 )
  good_points = find( noise(:,i) > 0 & ~isnan(signal) & ~isnan(noise) );

  if length( good_points ) > 3

    good_points = good_points(2:end-1);
    
    S = mean( signal(good_points,i) );
    N = mean( noise( good_points,i) );  
    if ( N > 0 ) all_ratio(i) = S/N; end;
  
  end

end

SN_ratio = mean( all_ratio );

