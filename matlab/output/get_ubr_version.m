function version_tag = get_ubr_version();
% version_tag = get_ubr_version();
%
% Read version tag from ubr_version.txt
%
% (C) R. Das, HHMI/Stanford University 2023.

dirname = fileparts(which(mfilename));
dirname_ubr_version = [dirname,'/../../'];
ubr_version = [dirname_ubr_version,'ubr_version.txt'];
assert(exist(ubr_version),'file');

fid = fopen(ubr_version,'r');
version_tag = fgetl(fid);
fclose(fid);




