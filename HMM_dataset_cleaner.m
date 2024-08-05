%% 
% Open script HMM_source_modelling
% Within the HMM_source_modelling script run first three sections
% Final section to run is titled: Cleaning up repetitive data blocks saved in brainstorm

clearvars -except original_sFiles
original_sFiles = original_sFiles';

%% Load dataset to clean
% NEVER LOAD SHUFFLED DATASETS and VARIABLES 
% LOAD ORIGINAL NON-SHUFFLED DATA
dataset_path = 'Q:\tactstim_resting\code\tactstim_rest_data\';
dataset_name = 'data_subject_PNAI_trial_preStim_rest_250_Hz_created_on_02_Jan_2022_12_52_01.mat';
time_name = 'Time_PNAI_trial_preStim_rest_250_Hz_created_on_02_Jan_2022_12_52_01.mat';
behavior_name = 'Behavior_PNAI_trial_preStim_rest_250_Hz_created_on_02_Jan_2022_12_52_01.mat';


load([dataset_path, dataset_name])
load([dataset_path, time_name])
load([dataset_path, behavior_name])

BEHAVIOR_cpy = BEHAVIOR(:,2);
%%

all_paths = cellfun(@(x) cell_bhvr(x), BEHAVIOR_cpy,'UniformOutput', false);
valid_block_names = cellfun(@(x) cell_path_split(x), original_sFiles, 'UniformOutput',false);
valid_trial_wrapper = @(x) cell_valid_trials(x,valid_block_names);
valid_trials = cellfun (valid_trial_wrapper, all_paths,'UniformOutput',false);

%% Now remove all extra data

isempty_wrapper = @(x) cellfun(@isempty,x,'UniformOutput',false);
valid_trials_index = cellfun(isempty_wrapper, valid_trials,'UniformOutput',false);

corrected_data = cellfun(@(x,y) cell_rmv_repeat_trials(x,y), ...
    data_subject', valid_trials_index, 'UniformOutput',false);
corrected_time = cellfun(@(x,y) cell_rmv_repeat_trials(x,y), ...
    T', valid_trials_index, 'UniformOutput',false);
corrected_behavior = cellfun(@(x,y) cell_rmv_repeat_trials(x,y), ...
    BEHAVIOR_cpy, valid_trials_index, 'UniformOutput',false);

%% Update
data_subject = corrected_data';
T = corrected_time';
BEHAVIOR_cpy = corrected_behavior;
BEHAVIOR(:,2) = BEHAVIOR_cpy;

%% Save data
clck = char(datetime);
clck = strrep(clck,':','_');
clck = strrep(clck,' ','_');
clck = strrep(clck,'-','_');

source_model = 'PNAI';
trial_type = 'preStim_rest';
downsample_rate = 250;

save([dataset_path 'Time_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'T','-v7.3')
save([dataset_path 'Behavior_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'BEHAVIOR','-v7.3')
save([dataset_path 'data_subject_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'data_subject','-v7.3')



