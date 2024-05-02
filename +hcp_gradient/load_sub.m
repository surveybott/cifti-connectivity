function sub = load_sub(varargin)
% inputs
fpath = fileparts(fileparts(mfilename('fullpath')));
p = inputParser;
p.addParameter('filepath',fullfile(fpath,'data'));
p.addParameter('filename','sub_n120.csv');
p.parse(varargin{:});
inputs = p.Results;
% load
file = fullfile(inputs.filepath,inputs.filename);
fid = fopen(file,'r');
sub = textscan(fid,'%d\n');
fclose(fid);
sub = double(sub{:});
end