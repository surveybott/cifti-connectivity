% GROUP_PIPELINE - finds functional cifti and confound data (assumed fmriprep derivatives) and runs `individual_gradient.m` on each subject
%       Denoising options:
%           - Voxel-wise regression (`confounds` parameter list corresponds to header of *_confounds.tsv fmriprep file
%               - setting `icaAroma`=true will automatically add noise components to the above regression
%           - `bandpass` filtering 
%           - `scrubThresh` and `scrubBefore`/`scrubAfter` to excise high motion volumes based on framewise displacement (fd)
%
% Requires the following packages to be downloaded and added to the MATLAB path:
%   - BrainSpace:     https://brainspace.readthedocs.io/en/latest/
%   - cifti-matlab:   https://github.com/Washington-University/cifti-matlab
%
% SurveyBott 2020

function sub = gradient_pipeline(varargin)
% setup inputs
p = inputParser;
p.addParameter('dataDir','/scratch/st-tv01-1/hbn/bids/derivatives',@isfolder);
p.addParameter('dataSuffix','*_den-91k_bold.dtseries.nii'); % suffix of cifti files to use
p.addParameter('outDir',''); % directory outputs (.mat, etc.) will be saved to
p.addParameter('outSuffix','');
p.addParameter('gradient',false); % run gradient decomposition or NOT
p.addParameter('includeSub',{});    % list of subjects to include (default (empty) includes all)
p.addParameter('includeTasks', {}); % same as above, but for tasks
% denoising options
p.addParameter('parc','schaefer');
p.addParameter('res',1000);
p.addParameter('confounds',{'a_comp_cor_00', 'a_comp_cor_01','a_comp_cor_02','a_comp_cor_03','a_comp_cor_04','a_comp_cor_05', ...
    'rot_x','rot_y','rot_z','trans_x','trans_y','trans_z', 'framewise_displacement'},@iscellstr);
p.addParameter('icaAroma',true);
p.addParameter('bandpass',[0.008 0.1],@(x) isnumeric(x) && numel(x)==2);
p.addParameter('scrubThresh',0.3,@isnumeric);
p.addParameter('scrubBefore',0,@isnumeric);
p.addParameter('scrubAfter',1,@isnumeric);
% misc 
p.addParameter('saveIndiv',true); 
p.addParameter('saveConn',true);
p.addParameter('runMaxVols',[]); % maximum number of volumes per run
p.addParameter('increment',[]); % # of volumes to successively add (multiple gradients)
p.parse(varargin{:});
inputs = p.Results;

%% create 'sub' struct - find motion and cifti files
% get cifti files, extract info, include based on inputs
cifti = tools.bids2struct(dir(fullfile(inputs.dataDir,'sub-*','**',inputs.dataSuffix)));
if ~isempty(inputs.includeSub)
    cifti = cifti(ismember({cifti.sub},inputs.includeSub));
end
if ~isempty(inputs.tasks)
    tasks = inputs.tasks;
else
    tasks = unique({cifti.task});
end
cifti = cifti(ismember({cifti.task},tasks));

% get confound files, remove those without
endings = {'desc-confounds_timeseries.tsv','desc-MELODIC_mixing.tsv','AROMAnoiseICs.csv'};
names = {'confounds','melodic','aroma'};
rm = false(size(cifti));
for i=1:numel(cifti)
    prefix = strsplit(cifti(i).name,'_space');
    cifti(i).prefix = prefix{1};
    for j=1:numel(endings)
        file = fullfile(cifti(i).path,sprintf('%s_%s',cifti(i).prefix,endings{j}));
        if exist(file,'file')
            cifti(i).(names{j}) = file;
        else
            cifti(i).(names{j}) = [];
            if inputs.icaAroma || strcmp(names{j},'confounds')
                rm(i) = true;
                fprintf('No %s\n',file)
            end
        end
    end
end
cifti(rm) = [];

% check subjects have all tasks
sub = unique({cifti.sub});
rm = false(size(sub));
for i=1:numel(sub)
    tmp = cifti(ismember({cifti.sub},sub{i}));
    if numel(unique({tmp.task})) ~= numel(tasks)
        rm(i) = true;
    end
end
fprintf('%d/%d subjects have all %d tasks\n',sum(~rm),numel(sub),numel(tasks));

%% run peach subject (in parallel)
if ~isempty(inputs.outDir) && ~isfolder(inputs.outDir)
    mkdir(inputs.outDir);
end
parfor i=1:numel(sub)
    fprintf('\t%d/%d\tsub-%s\n',i,numel(sub),sub{i});
    % run (and save)
    try
        out = individual_gradient(cifti(ismember({cifti.sub},sub{i})),'icaAroma',inputs.icaAroma,'confounds',inputs.confounds,'bandpass',inputs.bandpass,...
            'scrubThresh',inputs.scrubThresh,'scrubBefore',inputs.scrubBefore,'scrubAfter',inputs.scrubAfter,...
            'parc',inputs.parc,'res',inputs.res,'gradient',inputs.gradient,'increment',inputs.increment);
        if ~isempty(inputs.outDir) && inputs.saveIndiv
            filename = fullfile(inputs.outDir,sprintf('%s%s.mat',sub{i},inputs.outSuffix));
            save_gradient(out,filename,inputs.saveConn);
        end
    catch
        fprintf('ERROR\t%s\n',sub{i});
    end
end

end

function save_gradient(out,filename,saveConn)
if ~saveConn
   out.conn = [];
end
save(filename,'out');
end
