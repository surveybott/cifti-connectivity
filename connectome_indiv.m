function conn = connectome_indiv(indivFolder,varargin)
p = inputParser;
p.addRequired('indivFolder',@isfolder);
p.addParameter('matField','conn');
p.addParameter('prefix','');
p.parse(indivFolder,varargin{:});
inputs = p.Results;
% get mat files, sub and conditions
mat = dir(fullfile(indivFolder,'*.mat'));
sub = cell(size(mat));
cond = cell(size(mat));
for i=1:numel(mat)
    tmp = strsplit(regexprep(mat(i).name,'.mat',''),'_');
    sub{i} = tmp{1};
    cond{i} = tmp{2};
    mat(i).sub = sub{i};
    mat(i).cond = cond{i};
end
sub = unique(sub);
cond = unique(cond);
% create condition-level table of each subject's connectome
conn = [];
% load first connectome to get size and lower triangle indices
indiv = load(fullfile(mat(1).folder,mat(1).name));
idx = tril(true(size(indiv.out.(inputs.matField))),-1);
for i=1:numel(cond)
    t = nan(sum(idx(:)),numel(sub));
    parfor j=1:numel(sub)
        tmp = mat(strcmp({mat.sub},sub{j}) & strcmp({mat.cond},cond{i}));
        tmp_indiv = load(fullfile(tmp.folder,tmp.name));
        t(:,j) = tmp_indiv.out.(inputs.matField)(idx);
    end
    conn.(cond{i}) = array2table(t,'VariableNames',sub);
    if ~isempty(inputs.prefix)
        writetable(conn.(cond{i}),sprintf('%s_%s_n%d.csv',inputs.prefix,lower(cond{i}),numel(sub)),'WriteVariableNames',1);
    end
end
end
