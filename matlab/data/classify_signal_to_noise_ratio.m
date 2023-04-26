function  classification = classify_signal_to_noise_ratio( SN_ratio );
%  classification = classify_signal_to_noise_ratio( SN_ratio );

if nargin==0; help( mfilename ); return; end;

classification = 'N/A';

if ( SN_ratio == 0 )
  classification = 'N/A';
elseif ( SN_ratio  <  1 )
  classification = 'weak';
elseif ( SN_ratio < 5 )
  classification = 'medium';
elseif ( SN_ratio < 10 )
  classification = 'good';
else
  classification = 'strong';
end