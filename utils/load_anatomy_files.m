load('Q:\tactstim_resting\New folder\brainstorm_db\tactstim_resting\anat\@default_subject\tess_cortex_pial_low')
destSurfFile = '@default_subject/tess_cortex_pial_low.mat';
HeadModelType = 'surface';
HeadModelFunction = 'lcmv'; 
atlas_name = 'Mindboggle';
atlas_number = 6;
% ROIs from this atlas that are part of our dataset
ROIs = [3,4,5,6,11,12,15,16,19,20,21,22,23,24,25,26,27,28,29,30,...
    33,34,35,36,37,38,41,42,45,46,47,48,51,52,53,54,55,56,57,58,59,60]; %42 regions

database_anat = 'Q:\tactstim_resting\New folder\brainstorm_db\tactstim_resting\';
surface_file = destSurfFile;
surface_file_data = load([database_anat 'anat' filesep surface_file]);