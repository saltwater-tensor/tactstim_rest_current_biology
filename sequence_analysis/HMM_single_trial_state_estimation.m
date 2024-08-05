%% Single trial level states
% Subject specific states

% Load model 
% Load dataset (Should contain both time, data and sampling frequency
% information)
model_path = 'Q:\tactstim_resting\code\tactstim_rest_models\';
model_name = '_04_Feb_2022_17_04_27_HMM_model_preStim_rest_250Hz.mat';
dataset_path = 'Q:\tactstim_resting\code\tactstim_rest_data\';
mdl_nm = strsplit(model_name,'.');
load([model_path model_name]);
load(HMM_model.dataset_path)
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

output_path = 'Q:\tactstim_resting\code\tactstim_rest_data\subject_trial_specific_hmms\';

base_model_name = mdl_nm{1, 1};
mkdir([output_path base_model_name]);
output_path = [output_path base_model_name '\'];

for subj = 1:length(data_subject)
    f = flips(subj,:);
    for trials = 1:length(data_subject{1,subj})
        data_subj_trial = data_subject{1,subj}{trials,1};
        T_subj_trial = T{1,subj}(1,trials);
        d = flipdata(data_subj_trial,T_subj_trial,f);
        data_subj_trial = {d};
        T_subj_trial = {T_subj_trial};
        response_subj_trial = BEHAVIOR{subj, 2}(trials).Response;
        
        clear d
        
        % All the data standardisation steps are done befor the dual
        % estimation. Nothing is done within this step
        [hmm_subj_trial, Gamma_subj_trial, vpath_subj_trial] = hmmdual(data_subj_trial,T_subj_trial,HMM_model.hmm);
        hmm.hmm_subj_trial = hmm_subj_trial;
        hmm.Gamma_subj_trial = Gamma_subj_trial;
        hmm.vpath = vpath_subj_trial;
        hmm.response = response_subj_trial;
        
        save([output_path 'hmm_subj_' num2str(subj) '_trial_'...
            num2str(trials) '_' response_subj_trial],...
            "hmm")
        
        
    end
end

%% Get covariance matrix 
output_path_covs = [output_path 'covariance_matrices\'];
mkdir(output_path_covs)
hmm_files = dir(output_path);

for dr = 1:length(dir(output_path))
    if ~hmm_files(dr).isdir
        load(fullfile(hmm_files(dr).folder, hmm_files(dr).name))
        covmat = NaN(hmm.hmm_subj_trial.train.ndim, hmm.hmm_subj_trial.train.ndim,...
            hmm.hmm_subj_trial.K);
    
        for k = 1:hmm.hmm_subj_trial.K
            covmat(:,:,k) = getFuncConn(hmm.hmm_subj_trial,k);
        end
    
        save([output_path_covs 'covmat_' hmm_files(dr).name],...
        "covmat")
    end
end

%% Calculate state wise and response wise distance 
cov_files = dir(output_path_covs);

cov_resp1 = [];
cov_resp2 = [];
state = 5;
for dr = 1:length(dir(output_path_covs))
    dr
    if ~cov_files(dr).isdir

        load(fullfile(cov_files(dr).folder, cov_files(dr).name))
        nm = strsplit(cov_files(dr).name,'_');
        nm = nm{7};
        nm = strsplit(nm,'.');
        nm = nm{1};
        switch nm
            case 'Response1'
%                 cov_resp1 = cat(3,cov_resp1,covmat(:,:,state));
                cov_resp1 = cat(4,cov_resp1,covmat(:,:,:));
            case 'Response2'
%                 cov_resp2 = cat(3,cov_resp2,covmat(:,:,:));
                  cov_resp2 = cat(4,cov_resp2,covmat(:,:,:));
        end
    
    end
end

%% Riemannian statewise and response wise
for state = 1:6

    C1(:,:,state) = mean_covariances(squeeze(cov_resp1(:,:,state,1:end)),'riemann');
    C2(:,:,state) = mean_covariances(squeeze(cov_resp2(:,:,state,1:end)),'riemann');

end

% Response1
total_states = 6;
C_all_1 = zeros(total_states,total_states);
for k = 1:1:total_states

    for k2 = 1:1:total_states

        a = distance_riemann(C1(:,:,k), C1(:,:,k2));
        C_all_1(k,k2) = a;
    end   
end
figure(1)
heatmap(C_all_1,'Colormap',jet)
% Response2
total_states = 6;
C_all_2 = zeros(total_states,total_states);
for k = 1:1:total_states

    for k2 = 1:1:total_states

        a = distance_riemann(C2(:,:,k), C2(:,:,k2));
        C_all_2(k,k2) = a;
    end   
end
figure(2)
heatmap(C_all_2,'Colormap',jet)

% Response 1 vs Response 2
total_states = 6;
C_all_3 = zeros(total_states,total_states);
for k = 1:1:total_states

    for k2 = 1:1:total_states

        a = distance_riemann(C1(:,:,k), C2(:,:,k2));
        C_all_3(k,k2) = a;
    end   
end
figure(3)
heatmap(C_all_3,'Colormap',jet)

%% Response1 distances (within)
cnt1 = 1;
resp1_distances = [];
for resp1a = 1:size(cov_resp1,3)
    for resp1b = 1:size(cov_resp1,3)
        if ~(resp1a==resp1b)
            a = distance_riemann(squeeze(cov_resp1(:,:,resp1a)),squeeze(cov_resp1(:,:,resp1b)));
            resp1_distances(cnt1) = a;
            cnt1 = cnt1 + 1;
        end

    end

end
resp1_distances(isnan(resp1_distances)) = [];

% Response2 distances (within)
cnt2 = 1;
resp2_distances = [];
for resp2a = 1:size(cov_resp2,3)
    for resp2b = 1:size(cov_resp2,3)
        if ~(resp2a==resp2b)
            a = distance_riemann(squeeze(cov_resp2(:,:,resp2a)),squeeze(cov_resp2(:,:,resp2b)));
            resp2_distances(cnt2) = a;
            cnt2 = cnt2 + 1;
        end

    end

end
resp2_distances(isnan(resp2_distances)) = [];

histogram(resp1_distances,100,"Normalization","pdf")
hold on
histogram(resp2_distances,100,"Normalization","pdf")
% Response1 to Response 2 distances (across)

cnt3 = 1;
resp1_2_distances = [];
for resp1a = 1:size(cov_resp1,3)
    for resp2b = 1:size(cov_resp2,3)
        if ~(resp1a==resp2b)
            a = distance_riemann(squeeze(cov_resp1(:,:,resp1a)),squeeze(cov_resp2(:,:,resp2b)));
            resp1_2_distances(cnt3) = a;
            cnt3 = cnt3 + 1;
        end

    end

end
resp1_2_distances(isnan(resp1_2_distances)) = [];
histogram(resp1_2_distances,100,"Normalization","pdf")
