%%
% Load dataset (Should contain both time, data and sampling frequency
% information)
dataset_path = 'Q:\tactstim_resting\code\tactstim_rest_data\';
model_path = 'Q:\tactstim_resting\code\tactstim_rest_models\';
model_name = '_04_Feb_2022_17_04_27_HMM_model_preStim_rest_250Hz.mat';
single_trial_hmm = 'Q:\tactstim_resting\code\tactstim_rest_data\subject_trial_specific_hmms\';
mdl_nm = strsplit(model_name,'.');
base_model_name = mdl_nm{1, 1};

%% Load 

load(fullfile(model_path, model_name))
load(HMM_model.time_path)
load(HMM_model.freq_path)
load(HMM_model.behavior_path)


%% 
% Load subject and trial shuffled specific data
% Load flip loading 
% flip data dependent upon subject number
% 
load([dataset_path ...
    'shuffled_data_subject_PNAI_trial_preStim_rest_250_Hz_created_on_04_Feb_2022_13_01_09.mat']);
load([dataset_path ...
    'flips_and_scorepath_and_covmats_shuffled_DATA_PNAI_trial_preStim_rest_250_Hz_created_on_04_Feb_2022_13_01_09.mat']);

%%
cnt_resp1 = 1;
cnt_resp2 = 1;
gamma_resp1 = [];
gamma_resp2 = [];

sampling_freq = downsample_rate;
options_mt = struct('Fs',sampling_freq); % Sampling rate
options_mt.fpass = [1 45];  % band of frequency you're interested in
options_mt.tapers = [4 7]; % taper specification - leave it with default values
% options_mt.win = 2 * sampling_freq; % multitaper window
options_mt.to_do = [1 0]; % turn off pdc
options_mt.embeddedlags = HMM_model.hmm.train.embeddedlags;

size_wrapper = @(x) size(x,1);
total_trials = cell2mat(cellfun(size_wrapper, data_subject, 'UniformOutput', 0));
total_trials = sum(total_trials);
trial_number = 1;
h = waitbar(0, 'Single level spectral props');
for subj = 1:length(data_subject)
    f = flips(subj,:);
    for trials = 1:length(data_subject{1,subj})
        
        
        
        data_subj_trial = data_subject{1,subj}{trials,1};
        T_subj_trial = T{1,subj}(1,trials);
        d = flipdata(data_subj_trial,T_subj_trial,f);
%         data_subj_trial = {d};
%         T_subj_trial = {T_subj_trial};
        response_subj_trial = BEHAVIOR{subj, 2}(trials).Response;
        
%         clear d
        
        % Load the estimated hmm for this specific trial and subject
        load([single_trial_hmm base_model_name '\' 'single_trial_hmms\' 'hmm_subj_' num2str(subj) '_trial_'...
            num2str(trials) '_' response_subj_trial])
        
        gamma = hmm.Gamma_subj_trial;
        if length(gamma) == 1962
        
        switch response_subj_trial

            case 'Response1'
              resp1_string = ['hmm_subj_' num2str(subj) '_trial_'...
                    num2str(trials) '_' response_subj_trial];
              disp(resp1_string)
              
              disp_resp1_string = strsplit(resp1_string,'_');
              disp_resp1_string = strjoin(disp_resp1_string,' ');
                
              waitbar(trial_number/total_trials,h,disp_resp1_string)
                
                resp1_data{cnt_resp1,1} = d;
                resp1_T{cnt_resp1,1} = T_subj_trial;
                gamma_resp1 = [gamma_resp1; gamma];
                
                single_trial_spectral_analysis(resp1_data{cnt_resp1,1}, ...
                    resp1_T{cnt_resp1,1},...
                    gamma, options_mt,...
                    [single_trial_hmm base_model_name '\' 'single_trial_spectra\' resp1_string '_spectra'])
                
                cnt_resp1 = cnt_resp1 + 1;    
                trial_number = trial_number + 1;
                
            case 'Response2'
                
                resp2_string = ['hmm_subj_' num2str(subj) '_trial_'...
                    num2str(trials) '_' response_subj_trial];

                disp(resp2_string)
                disp_resp2_string = strsplit(resp2_string,'_');
                disp_resp2_string = strjoin(disp_resp2_string,' ');
                
                waitbar(trial_number/total_trials,h,disp_resp2_string)
                
                resp2_data{cnt_resp2,1} = d;
                resp2_T{cnt_resp2,1} = T_subj_trial;
                gamma_resp2 = [gamma_resp2; gamma];
                
                 single_trial_spectral_analysis(resp2_data{cnt_resp2,1}, ...
                    resp2_T{cnt_resp2,1},...
                    gamma, options_mt,...
                    [single_trial_hmm base_model_name '\' 'single_trial_spectra\' resp2_string '_spectra'])
                
                cnt_resp2 = cnt_resp2 + 1;    
                trial_number = trial_number + 1;
                                        
        end

        end

    end
end


%% 

% Group level spectral characterisation of states separated by response
% types

clearvars -except resp1_data resp1_T gamma_resp1 resp2_data resp2_T gamma_resp2
sampling_freq = downsample_rate;
options_mt = struct('Fs',sampling_freq); % Sampling rate
options_mt.fpass = [1 45];  % band of frequency you're interested in
options_mt.tapers = [4 7]; % taper specification - leave it with default values
% options_mt.win = 2 * sampling_freq; % multitaper window
options_mt.to_do = [1 0]; % turn off pdc
options_mt.embeddedlags = HMM_model.hmm.train.embeddedlags;

[group_spectra_resp1] = ...
    trial_level_spectral_analysis(resp1_data,resp1_T,gamma_resp1,options_mt);

[group_spectra_resp2] = ...
    trial_level_spectral_analysis(resp2_data,resp2_T,gamma_resp2,options_mt);

% HMM_trial_spectra.fitmt_trial_resp1 = trial_spectra_resp1;
% HMM_trial_spectra.fitmt_trial_resp2 = trial_spectra_resp2;
HMM_trial_spectra.fitmt_group_resp1 = group_spectra_resp1;
clear group_spectra_resp1
HMM_trial_spectra.fitmt_group_resp2 = group_spectra_resp2;
clear group_spectra_resp2

%%

% Load nnmf factors into which to project the above data

% load([model_path base_model_name '_spectral_analysis\nnmf_profiles_and_group_projections_3' ])
% clear fitmt_subj_fact_4b fitmt_group_fact_4b

spectral_dir = dir([single_trial_hmm base_model_name '\' 'single_trial_spectra\']);
h = waitbar(0, 'Single level spectral props projection to 4 bands');
trial_number = 1;
total_trials = length(spectral_dir) - 2;

for spdr = 1:length(spectral_dir)

    if ~spectral_dir(spdr).isdir
    
        disp_resp_string = strsplit(spectral_dir(spdr).name,'_');
        disp_resp_string = strjoin(disp_resp_string,' ');    
        waitbar(trial_number/total_trials,h, disp_resp_string)

        load(fullfile(spectral_dir(spdr).folder, spectral_dir(spdr).name))

        [state] = spectral_projection_response...
            (trial_spectra, supply_profiles);

%         trial_spectra_4b.psd = psd_proj_response;
%         trial_spectra_4b.coh = coh_proj_response;

        trial_spectra_4b.state = state;

        save(fullfile(single_trial_hmm, base_model_name, 'single_trial_projection_spectra', spectral_dir(spdr).name),"trial_spectra_4b")
        
        clear trial_spectra_4b
        clear trial_spectra

        trial_number = trial_number + 1;

    end

end

%% Test the projections






