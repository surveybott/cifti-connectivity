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

function sub = group_pipeline(varargin)
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
p.addParameter('scrubThresh',0.3,@isnumeric); % set to 0 or Inf to disable
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
    if isfield(cifti,'sub') && isnumeric(cifti(i).sub)
        cifti(i).sub = num2str(cifti(i).sub);
    end
    if isfield(cifti,'ses') && isnumeric(cifti(i).ses)
        cifti(i).ses = num2str(cifti(i).ses);
    end
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
    fprintf('%d/%d\tsub-%s',i,numel(sub),sub{i});
    sub_cifti = cifti(ismember({cifti.sub},sub{i}));
    if isfield(sub_cifti, 'ses')
        % loop over sessions
        sub_ses = unique({sub_cifti.ses});
        for j=1:numel(sub_ses)
            fprintf('\n\t\ttses-%s',sub_ses{j});
            sub_ses_cifti = sub_cifti(ismember({sub_cifti.ses},sub_ses{j}));
            try
                run_indiv(sub_ses_cifti, inputs);
            catch
                fprintf('\tERROR');
            end
        end
    else
        try
            run_indiv(sub_cifti, inputs);
        catch
            fprintf('\tERROR');
        end
    end
    fprintf('\n');
end

end

function out = run_indiv(cifti, inputs)
% run
out = individual_pipeline(cifti,'icaAroma',inputs.icaAroma,'confounds',inputs.confounds,'bandpass',inputs.bandpass,...
            'scrubThresh',inputs.scrubThresh,'scrubBefore',inputs.scrubBefore,'scrubAfter',inputs.scrubAfter,...
            'parc',inputs.parc,'res',inputs.res,'gradient',inputs.gradient,'increment',inputs.increment);
% save
if ~isempty(inputs.outDir) && inputs.saveIndiv
    if isfield(cifti,'ses')
        filename = fullfile(inputs.outDir, sprintf('sub-%s_ses-%s%s.mat',cifti(1).sub,cifti(1).ses,inputs.outSuffix));
    else
        filename = fullfile(inputs.outDir, sprintf('sub-%s%s.mat',cifti(1).sub,inputs.outSuffix));
    end
    if ~inputs.saveConn
        out.conn = [];
    end
    save(filename,'out');
end
end
