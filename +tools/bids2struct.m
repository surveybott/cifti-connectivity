% BIDS2STRUCT extracts metadata from one or more BIDS filenames
%
%  SurveyBott 2019

function s = bids2struct(in,varargin)
% inputs
p = inputParser;
p.addRequired('in',@(x) ischar(x) || iscellstr(x) || isstruct(x));
p.addParameter('file',false,@(x) islogical(x) || isnumeric(x)); % validate file?
p.parse(in,varargin{:});
inputs = p.Results;
if ischar(in)
   in = {in}; 
end
if isstruct(in) && all(isfield(in,{'folder','name'}))
    in = arrayfun(@(x) fullfile(x.folder,x.name),in,'UniformOutput',0);
elseif isstruct(in)
    error('struct input must be from ''dir'' function (has name and folder fields)');
end
% validate files
if inputs.file 
    valid = cellfun(@(x) exist(x,'file'),in);
    if ~all(valid)
        error('%d/%d file(s) invalid',sum(~valid),numel(valid));
    end
end

% filename
for i=1:numel(in)
    s(i).file = in{i};
    [s(i).path,s(i).name,s(i).ext] = fileparts(s(i).file); 
    if strcmp(s(i).ext,'.gz') && numel(s(i).name) >= 4 && strcmp(s(i).name(end-3:end),'.nii')
        s(i).ext = '.nii.gz';
        [~,s(i).name] = fileparts(s(i).name);
    end
    
    % labels
    dash = regexp(s(i).name,'-');
    uscore = regexp(s(i).name,'_');
    for j=1:numel(dash)
        idx = uscore(find(uscore < dash(j),1,'last'));
        if isempty(idx)
            idx = 1;
        else
            idx = idx + 1;
        end
        field = s(i).name(idx:dash(j)-1);
        value = s(i).name(dash(j)+1:uscore(find(uscore > dash(j),1,'first'))-1);
        s(i).(field) = value;
        value = str2double(value);
        if ~isnan(value)
            s(i).(field) = value;
        end
    end
    
    % suffix
    if uscore(end) > dash(end)
        s(i).suffix = s(i).name(uscore(end)+1:end);
    end
end
end
%#ok<*AGROW>