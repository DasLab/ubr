function strip_nan(outfile, verbose)
% strip_nan(outfile, verbose)
%
% Replace NaN with blank ('') for readability, in place.
%
% (C) R. Das, Stanford & HHMI, 2023
if ~exist('verbose'); verbose = 1; end;
if verbose; fprintf('\nReplacing NaNs with blanks in %s...\n',outfile); end;
command = sprintf('sed -I "" "s/NaN//g"  %s',outfile);
if verbose; fprintf([command,'\n']); end;
system(command);
