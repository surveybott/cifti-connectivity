# cifti-connectivity
Denoising, parcellation, connectivity matrices and gradient mapping on CIFTI data (currently cortex only) from fmriprep

## Setup
The follwing dependencies (properly installed and added to the MATLAB path) are required:
- cifti-matlab: [GitHub](https://github.com/Washington-University/cifti-matlab)
- BrainSpace: [docs](https://brainspace.readthedocs.io/en/latest/pages/install.html#matlab-installation), [download](https://github.com/MICA-MNI/BrainSpace/releases)
### Parcellations
Parcellation is handled by BrainSpace. Additional parcellation definitions (e.g. the [Schaefer 1000](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal)) are included here in the [parcellations](parcellations/) folder and need to be copied into `<BRAINSPACE_DIRECTORY>/shared/parcellations` prior to use

**Warning**: The Schaefer 1000 parcellation actually has **998** definied parcels (no 533 or 903). Check outputs and indices carefully, especially when mapping back onto the brain (fsLR/conte69 32k)

## [group_pipeline.m](group_pipeline.m)
finds functional cifti and confound data (assumed fmriprep derivatives) and runs [`individual_gradient.m`](individual_gradient.m) on each subject (in parallel)

### Denoising options 
  - `'confounds'` list corresponds to header of `*_confounds.tsv` fmriprep file. Voxel-wise regression will be performed with these confounds 
  - `'icaAroma'=true` will automatically add noise components to the above regression
  - `'bandpass'` filtering
  - `'scrubThresh'` and `'scrubBefore'`/`'scrubAfter'` to excise high motion volumes based on framewise displacement (fd)

### Parcellation and Connectome
After denoising, a parcellation (`'parc'` and `'res'` parameters) is applied, mean timecourses and covariance (Pearson's) computed

### Example Usage
The below overloads many default parameters (*which exist in the code for **example** purposes*) and performs a simple regression of fd before computing connectivity

`gradient_pipeline('dataDir', '~/bids/derivatives', 'confounds', {'framewise_displacement'}, 'icaAroma', false, 'bandpass', [], 'scrubThresh', 0)`


