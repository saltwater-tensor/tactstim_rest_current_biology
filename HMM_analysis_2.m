%% Viterbi sequence analysis

% Reshape vpath output of the model for each trial
% Remember that vpath is only defined for embedded time series !
VPATH_SUBJECT = {};
GAMMA_SUBJECT = {};
Gamma = HMM_model.Gamma;
T_embedded = cellfun(@(x) x-(length(HMM_model.hmm.train.embeddedlags)-1),T,...
    'UniformOutput',false);

subject_T = T_embedded{1,1};
subject_total_time = sum(subject_T);
subject_start = 1;
subject_end = subject_total_time;
subject_vpath = vpath(subject_start:subject_end);
subject_gamma = Gamma(subject_start:subject_end,:);
subj_count = 1;

for subj = 2:length(T)+1 
    

        trial_start = 1;
        trial_end = (subject_T(1,1));
        trial_count = 1;

        for trial = 2:length(subject_T)+1

            vpath_trial = subject_vpath(trial_start:trial_end,1);
            gamma_trial = subject_gamma(trial_start:trial_end,:);
            VPATH_SUBJECT{subj_count,1}{trial_count,1} = vpath_trial;
            GAMMA_SUBJECT{subj_count,1}{trial_count,1} = gamma_trial;
            
            vpath_trial = [];
            gamma_trial = [];
            
            trial_start = trial_end+1;
            trial_count = trial_count + 1;
            
            if trial_start > subject_total_time
                continue
            else
                trial_end = trial_end + subject_T(1,trial);
            end
        end

        if subj > length(T)
            if subject_end == length(vpath)
                msgbox('Viterbi path and gamma reshaped','RESHAPE')
            end
            continue
        else
            subject_T = T_embedded{1,subj};
            subject_total_time = sum(subject_T);
            subject_start = subject_end+1;
            subject_end = subject_end + subject_total_time;
            subject_vpath = vpath(subject_start:subject_end);
            subject_gamma = Gamma(subject_start:subject_end,:);
            subj_count = subj_count + 1;
        end
    
end


%% VITERBI PSTH
VPATH_PSTH_RESPONSE1 = {};
VPATH_PSTH_RESPONSE2 = {};
GAMMA_PSTH_RESPONSE1 = {};
GAMMA_PSTH_RESPONSE2 = {};

for subj_bh = 1:size(BEHAVIOR,1)
    bhvr = BEHAVIOR{subj_bh,2};
    response1 = find(ismember({bhvr.Response}, {'Response1'}));
    response2 = find(ismember({bhvr.Response}, {'Response2'}));
    response1_SOA = [bhvr(response1).SOA_ms].';
    response2_SOA = [bhvr(response2).SOA_ms].';
    
    response1_vpaths = VPATH_SUBJECT{subj_bh,1}(response1,1);
    response2_vpaths = VPATH_SUBJECT{subj_bh,1}(response2,1);
    
    response1_gammas = GAMMA_SUBJECT{subj_bh,1}(response1,1);
    response2_gammas = GAMMA_SUBJECT{subj_bh,1}(response2,1);
    
    time_length_max = cellfun(@(x) length(x),response1_vpaths);
    time_length_max = max(time_length_max);
    
    
    
    
    
    % PSTH response 1
    r1_matrix = NaN(time_length_max,HMM_model.hmm.K);
    r1_count = zeros(time_length_max,HMM_model.hmm.K);
    g1_matrix = NaN(time_length_max,HMM_model.hmm.K);
    g1_count = zeros(time_length_max,HMM_model.hmm.K);
    
    for resp_length = 1:length(response1)
        
        r1 = response1_vpaths{resp_length};
        r1_matrix_tmp = NaN(time_length_max,HMM_model.hmm.K);
        r1_matrix_tmp(find(r1==1),1) = 1;
        r1_matrix_tmp(find(r1==2),2) = 1;
        r1_matrix_tmp(find(r1==3),3) = 1;
        r1_matrix_tmp(find(r1==4),4) = 1;
        r1_matrix_tmp(find(r1==5),5) = 1;
        r1_matrix_tmp(find(r1==6),6) = 1;
        % Indices for which entries are found and will be updated
        r1_count_tmp = double(~isnan(r1_matrix_tmp));
%         r2_matrix = r2_matrix + r2_matrix_tmp;
        r1_matrix = cat(3,r1_matrix, r1_matrix_tmp);
        r1_matrix = sum(r1_matrix,3,'omitnan');
        r1_count = r1_count + r1_count_tmp;
        
        g1 = response1_gammas{resp_length};
        g1_matrix_tmp = NaN(time_length_max,HMM_model.hmm.K);
        g1_count_tmp = zeros(time_length_max,HMM_model.hmm.K);
        
        g1_matrix_tmp(1:size(g1,1),:) = g1;
        g1_count_tmp = double(~isnan(g1_matrix_tmp));
%         g2_matrix = g2_matrix + g2_matrix_tmp;
        g1_matrix = cat(3,g1_matrix, g1_matrix_tmp);
        g1_matrix = sum (g1_matrix,3,'omitnan');
        g1_count = g1_count + g1_count_tmp;
        
        
    end
    
    VPATH_PSTH_RESPONSE1{subj_bh,1}= r1_matrix;
    VPATH_PSTH_RESPONSE1{subj_bh,2}= r1_count;
    
    % Averaging across all responses 1 for a given subject
%     g1_matrix_test1 = g1_matrix./length(response1);
%     g1_matrix_test2 = g1_matrix./g1_count;
    g1_matrix = g1_matrix./g1_count;
    GAMMA_PSTH_RESPONSE1{subj_bh,1}= g1_matrix;
    
    
    
    % PSTH response 2
    r2_matrix = NaN(time_length_max,HMM_model.hmm.K);
    r2_count = zeros(time_length_max,HMM_model.hmm.K);
    g2_matrix = NaN(time_length_max,HMM_model.hmm.K);
    g2_count = zeros(time_length_max,HMM_model.hmm.K);
    
    for resp_length = 1:length(response2)
        
        r2 = response2_vpaths{resp_length};
        r2_matrix_tmp = NaN(time_length_max,HMM_model.hmm.K);
        r2_count_tmp = zeros(time_length_max,HMM_model.hmm.K);
        
        
        r2_matrix_tmp(find(r2==1),1) = 1;
        r2_matrix_tmp(find(r2==2),2) = 1;
        r2_matrix_tmp(find(r2==3),3) = 1;
        r2_matrix_tmp(find(r2==4),4) = 1;
        r2_matrix_tmp(find(r2==5),5) = 1;
        r2_matrix_tmp(find(r2==6),6) = 1;
        % Indices for which entries are found and will be updated
        r2_count_tmp = double(~isnan(r2_matrix_tmp));
%         r2_matrix = r2_matrix + r2_matrix_tmp;
        r2_matrix = cat(3,r2_matrix, r2_matrix_tmp);
        r2_matrix = sum(r2_matrix,3,'omitnan');
        r2_count = r2_count + r2_count_tmp;
        
        
        % Averaging across all responses for a given subject
        g2 = response2_gammas{resp_length};
        g2_matrix_tmp = NaN(time_length_max,HMM_model.hmm.K);
        g2_count_tmp = zeros(time_length_max,HMM_model.hmm.K);
        
        g2_matrix_tmp(1:size(g2,1),:) = g2;
        g2_count_tmp = double(~isnan(g2_matrix_tmp));
%         g2_matrix = g2_matrix + g2_matrix_tmp;
        g2_matrix = cat(3,g2_matrix, g2_matrix_tmp);
        g2_matrix = sum (g2_matrix,3,'omitnan');
        g2_count = g2_count + g2_count_tmp;
        
    end
    
    VPATH_PSTH_RESPONSE2{subj_bh,1}= r2_matrix;
    VPATH_PSTH_RESPONSE2{subj_bh,2}= r2_count;
    
    % Averaging across all responses for a given subject
%     g2_matrix = g2_matrix./length(response2);
    g2_matrix = g2_matrix./g2_count;
    GAMMA_PSTH_RESPONSE2{subj_bh,1}= g2_matrix;
    
end

%% Gamma analysis

response1_viterbi_path_SUBJ_AVRG = NaN(time_length_max,HMM_model.hmm.K);
response2_viterbi_path_SUBJ_AVRG = NaN(time_length_max,HMM_model.hmm.K);

response1_gamma_SUBJ_AVRG = NaN(time_length_max,HMM_model.hmm.K);
response2_gamma_SUBJ_AVRG = NaN(time_length_max,HMM_model.hmm.K);

for gm_subj = 1:length(GAMMA_PSTH_RESPONSE1)
    
    response1_viterbi_path_SUBJ_AVRG = cat(3, response1_viterbi_path_SUBJ_AVRG, ...
        VPATH_PSTH_RESPONSE1{gm_subj,1});
    
    response2_viterbi_path_SUBJ_AVRG = cat(3, response2_viterbi_path_SUBJ_AVRG,...
        VPATH_PSTH_RESPONSE2{gm_subj,1});
    
    response1_gamma_SUBJ_AVRG = cat(3, response1_gamma_SUBJ_AVRG, ...
        GAMMA_PSTH_RESPONSE1{gm_subj});
    
    response2_gamma_SUBJ_AVRG = cat(3, response2_gamma_SUBJ_AVRG,...
        GAMMA_PSTH_RESPONSE2{gm_subj});
    
end
response1_gamma_SUBJ_AVRG = sum(response1_gamma_SUBJ_AVRG,3,'omitnan');
response2_gamma_SUBJ_AVRG = sum(response2_gamma_SUBJ_AVRG,3,'omitnan');
response1_viterbi_path_SUBJ_AVRG = sum(response1_viterbi_path_SUBJ_AVRG,3,'omitnan');
response2_viterbi_path_SUBJ_AVRG = sum(response2_viterbi_path_SUBJ_AVRG,3,'omitnan');

response1_gamma_SUBJ_AVRG = response1_gamma_SUBJ_AVRG./length(GAMMA_PSTH_RESPONSE1);
response2_gamma_SUBJ_AVRG = response2_gamma_SUBJ_AVRG./length(GAMMA_PSTH_RESPONSE2);
response1_viterbi_path_SUBJ_AVRG = response1_viterbi_path_SUBJ_AVRG./length(GAMMA_PSTH_RESPONSE1);
response2_viterbi_path_SUBJ_AVRG = response2_viterbi_path_SUBJ_AVRG./length(GAMMA_PSTH_RESPONSE2);
