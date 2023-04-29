function structure_map = get_mfe_structure_map( mfe_structures );
% structure_map = get_mfe_structure_map( mfe_structures );
% (C) R. Das, HHMI, Stanford, 2023
if ~exist( 'get_mfe_structure_map', 'file' )
  fprintf( 'ERROR! You need to install OpenKnotScore from https://github.com/eternagame/OpenKnotScore, and add OpenKnotScore/matlab and subdirectories to your path.' );
  return;
end

structure_map = get_structure_map( mfe_structures );

fprintf( 'WARNING! get_mfe_structure_map has been deprecated in favor of structure_map in the OpenKnotScore repository at https://github.com/eternagame/OpenKnotScore\n' );
