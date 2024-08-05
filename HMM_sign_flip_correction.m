
dataset_path = 'Q:\tactstim_resting\code\tactstim_rest_data\shuffled_DATA_PNAI_trial_preStim_rest_250_Hz_created_on_04_Feb_2022_13_01_09.mat';
load(dataset_path)
dataset_name = strsplit(dataset_path,'\');
dataset_name = dataset_name{end};
time_path = 'Q:\tactstim_resting\code\tactstim_rest_data\shuffled_Time_PNAI_trial_preStim_rest_250_Hz_created_on_04_Feb_2022_13_01_09.mat';
load(time_path)
data = DATA;
clear DATA

%% Sign flip correction
disp('Running sign flip correction on data...')
disp('You can pause and stop execution if you dont want to flip')
% pause(5)
options.maxlag = 10;
% options.noruns = 10;
[flips,scorepath,covmats_unflipped] = findflip(data,T,options);
data = flipdata(data,T,flips);
save(['sign_flip_corrected_' dataset_name],'data','-v7.3')
save(['flips_and_scorepath_and_covmats_' dataset_name],'flips','scorepath','covmats_unflipped','-v7.3')
clear options