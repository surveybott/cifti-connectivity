% function to load HCP-stype cifti files and use BrainSpace toolbox to
% compute gradients - https://brainspace.readthedocs.io/en/latest/
%
% also uses https://github.com/Washington-University/cifti-matlab
function out = individual_pipeline(cifti,varargin)
% handle inputs
p = inputParser;
p.addParameter('parc','schaefer');
p.addParameter('res',1000);
p.addParameter('gradient',true);
p.addParameter('plot',false);
p.addParameter('scrubThresh',0.3,@isnumeric);
p.addParameter('scrubBefore',1,@isnumeric);
p.addParameter('scrubAfter',1,@isnumeric);
p.addParameter('increment',[],@isnumeric);
p.addParameter('confounds',{});
p.addParameter('icaAroma',true);
p.addParameter('bandpass',[0.008 0.1]);
p.addParameter('tr',[]);
p.parse(varargin{:});

inputs = p.Results;
out = [];

warning('off','all')

% load conte69 template (HCP cortex mapping)
[conte.lh, conte.rh] = load_conte69();

% load parcellation
if ~isempty(inputs.res) && inputs.res
    parc = load_parcellation(inputs.parc,inputs.res);
    parc = parc.(sprintf('%s_%d',inputs.parc,inputs.res));
end

% load ciftis, split off surface (with mask where data exist ~29/32k), denoise, concatenate 
gii = [];
mask = [];
scrubTs = [];
for i=1:numel(cifti)
   img = cifti_read(cifti(i).file);
   % cortical surface
   [tmp_lh, tmp_mask_lh] = cifti_struct_dense_extract_surface_data(img,'CORTEX_LEFT');
   [tmp_rh, tmp_mask_rh] = cifti_struct_dense_extract_surface_data(img,'CORTEX_RIGHT');
   tmp = [tmp_lh; tmp_rh];
   img.cdata = [];
   clear tmp_lh tmp_rh
   tmp_mask = [tmp_mask_lh; tmp_mask_rh];
   clear tmp_mask_lh tmp_mask_rh
   % setup denoising
   if inputs.icaAroma
       melodic = load(cifti(i).melodic); % tpts x components
       badComps = load(cifti(i).aroma);
       % add constant
       melodic = [melodic ones(size(melodic,1),1)];
   end
   if exist(cifti(i).confounds,'file')
       confounds = readtable(cifti(i).confounds,'FileType','delimitedtext');
       if ~isempty(inputs.confounds)
           noise = confounds{:,ismember(confounds.Properties.VariableNames,inputs.confounds)};
           noise(isnan(noise)) = 0;
           % add constant
           noise = [noise ones(size(noise,1),1)];
       end
   else
       error("Counfounds file doesn't exist");
   end
   % denoise (each vertex!)
   if inputs.icaAroma || ~isempty(inputs.confounds)
       for v=find(tmp_mask)'
           y = tmp(v,:)';
           if inputs.icaAroma
               b = regress(y,melodic);
               y = y - melodic(:,badComps)*b(badComps);
           end
           if ~isempty(inputs.confounds)
               [~,~,y] = regress(y,noise);
           end
           tmp(v,:) = y;
       end
   end
   % bandpass freq filter
   if ~isempty(inputs.bandpass)
       if ~isempty(inputs.tr)
           hz = 1 / inputs.tr;
       elseif strcmpi(img.diminfo{2}.seriesUnit,'SECOND') || strcmpi(img.diminfo{2}.seriesUnit,'SECONDS')
           hz = 1 / img.diminfo{2}.seriesStep;
       else
           fprintf('BANDPASS FILTERING NOT APPLIED, COULDN''T FIND TR');
           hz = 0;
       end
       if hz
           tmp(tmp_mask,:) = bandpass(tmp(tmp_mask,:)',inputs.bandpass,hz)';
       end
   end
   % scrubbing
   if inputs.scrubThresh
       fd = confounds.framewise_displacement;
       % fd is a 1D array, make row vector, add '0' to make full length
       if size(fd,1) > size(fd,2)
           fd = fd';
       end
       fd(isnan(fd)) = 0;
       rm = false(size(fd));
       % scrub bad volumes
       rm = rm | fd > inputs.scrubThresh;

       % scrub volumes around bad volumes (scrub before and scrub after)
       before = false(size(rm));
       for j=1:inputs.scrubBefore
           before = before | [rm(j+1:end) false(1,j)];
       end

       after = false(size(rm));
       for j=1:inputs.scrubAfter
           after = after | [false(1,j) rm(1:end-j)];
       end
       rm = rm | before | after;
       scrubTs = [scrubTs rm];
   end
   gii = [gii tmp];
   mask = [mask tmp_mask];

   % volumetric structures
   %[tmp, ~, tmp_mask] = cifti_struct_dense_extract_volume_all_data(img);
end
% do scrubbing
if ~isempty(scrubTs)
   gii = gii(:,~scrubTs);
end
clear tmp*

mask = all(mask,2);
gii(~mask,:) = NaN;

% setup increments (multiple connectomes / gradients of with increasing amts of data)
if isempty(inputs.increment)
    vols = size(gii,2);
else
    vols = inputs.increment:inputs.increment:size(gii,2);
    if vols(end) < size(gii,2)
        vols(end+1) = size(gii,2);
    end
end

% correlation matrix
g = {};
gradients = [];
for i=1:numel(vols)
    if ~isempty(inputs.res) && inputs.res
        parcellated = full2parcel(gii(:,1:vols(i))',parc)';
        conn(:,:,i) = corr(parcellated);
    else
        conn(:,:,i) = corr(gii(:,1:vols(i))');
    end
    idx(:,i) = ~all(isnan(conn(:,:,i))); % deal with missing parcels
    % do gradient mapping
    if inputs.gradient
        g{i} = GradientMaps('kernel','cosine','approach','dm');
        try
            g{i} = g{i}.fit(conn(idx(:,i),idx(:,i),i));

            % fix missing parcels and plot
            gradients(:,:,i) = nan(numel(idx(:,i)),size(g{i}.gradients{1},2));
            gradients(idx(:,i),:,i) = g{i}.gradients{1};
            if inputs.plot
                scree_plot(g.lambda{1});
                plot_hemispheres(gradients(:,1:3,i),{conte.lh conte.rh},'parcellation',parc); % first 2
                gradient_in_euclidean(gradients(:,1:3,i),{conte.lh conte.rh},parc);
            end
        catch
            fprintf('Error computing gradient');
        end
    end
end
% outputs
out.inputs = inputs;
out.cifti = cifti;
out.scrub = scrubTs;
out.fd = fd
out.vols = vols;
out.conn = conn;
out.idx = idx;
out.g = g;
out.gradients = gradients;
end
