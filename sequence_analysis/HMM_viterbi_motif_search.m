%%
% Load dataset (Should contain both time, data and sampling frequency
% information)
dataset_path = 'Q:\tactstim_resting\code\tactstim_rest_data\';
model_path = 'Q:\tactstim_resting\code\tactstim_rest_models\';
model_name = '_04_Feb_2022_17_04_27_HMM_model_preStim_rest_250Hz.mat';
single_trial_hmm = 'Q:\tactstim_resting\code\tactstim_rest_data\subject_trial_specific_hmms\';
mdl_nm = strsplit(model_name,'.');
base_model_name = mdl_nm{1, 1};

%%
hmm_files = dir(fullfile(single_trial_hmm, base_model_name, 'single_trial_hmms'));
v = [1,4,5,6];
p = perms(v);

for dr = 1:length(hmm_files)
    dr
    if ~hmm_files(dr).isdir
            
            load(fullfile(hmm_files(dr).folder, hmm_files(dr).name))
            file_name = strsplit(hmm_files(dr).name,'.');
            file_name = file_name{1};
            vpath = hmm.vpath;
            [unique_viterbi_transitions] = unique_transitions(vpath);
            unique_viterbi_transitions = strrep(unique_viterbi_transitions,'  ','');
%             k1 = strfind(unique_viterbi_transitions,'4615');
%             k2 = strfind(unique_viterbi_transitions,'6415');
%             k3 = strfind(unique_viterbi_transitions,'1546');
            hmm.unique_vpath = unique_viterbi_transitions;
            save(fullfile(single_trial_hmm, base_model_name, 'single_trial_unique_viterbis', [file_name '_unq_vtrb']), "unique_viterbi_transitions")
            
    end

end

%%

% We will test all possible permutations for the 4 states that were found
% to be interesting

hmm_files_unq_vtrb = dir(fullfile(single_trial_hmm, base_model_name, 'single_trial_unique_viterbis'));
v = [1,4,5,6];
p = perms(v);
motif_numbers = zeros(2, size(p,1));
motifs_all_trials = zeros(length(hmm_files_unq_vtrb)-2,size(p,1),2);

for dr = 1:length(hmm_files_unq_vtrb)
    dr
    if ~hmm_files_unq_vtrb(dr).isdir
            
            load(fullfile(hmm_files_unq_vtrb(dr).folder, hmm_files_unq_vtrb(dr).name))
            file_name = strsplit(hmm_files_unq_vtrb(dr).name,'.');
            response_name = strsplit(file_name{1},'_');
            response_name = response_name{6};
            file_name = file_name{1};
            if exist("unique_viterbi_transitions","var")
                hmm_unq_vtrb.unique_sequence = unique_viterbi_transitions;
            elseif exist("hmm_unq_vtrb","var")
                unique_viterbi_transitions = hmm_unq_vtrb.unique_sequence;
            end
            
            motif_numbers_trialwise = zeros(size(p,1),1);
            
            % Finding all motifs
            for prm = 1:size(p,1)
                current_permutation = p(prm,:);
                current_permutation = num2str(current_permutation);
                current_permutation = convertCharsToStrings(strrep(current_permutation,'  ',''));
                k1 = strfind(unique_viterbi_transitions,current_permutation);
                motif_numbers_trialwise(prm) = length(k1);
                % Accumulating motif search across all trials
                switch response_name
                    
                    case 'Response1'
                        motif_numbers(1,prm) = motif_numbers(1,prm) + length(k1);
                        motifs_all_trials(dr-2,prm,1) = length(k1);

                    case 'Response2'
                        motif_numbers(2,prm) = motif_numbers(2,prm) + length(k1);
                        motifs_all_trials(dr-2,prm,2) = length(k1);
                end
                
            end


            hmm_unq_vtrb.motif_numbers_trialwise = motif_numbers_trialwise;
            
            % Rewriting
%             save(fullfile(single_trial_hmm, base_model_name, 'single_trial_unique_viterbis', [file_name]), "hmm_unq_vtrb")
            
            clear unique_viterbi_transitions

    end

end

%% Testing sequencing

% Even if you test for a specific permutation correct for multiple
% comparisions by taking into account all possible permutations

% Response1 vs Response2
motif_number = 24;
motif_specific_data = squeeze(motifs_all_trials(:,motif_number,:));
response1_motifs = motif_specific_data(:,1);
response1_motifs(response1_motifs == 0) = [];
response2_motifs = motif_specific_data(:,2);
response2_motifs(response2_motifs == 0) = [];

original_difference = sum(response1_motifs) - sum(response2_motifs);
original_id_vector = [ones(length(response1_motifs),1); 2*ones(length(response2_motifs),1)];
original_motif_vector = [response1_motifs; response2_motifs];

for iters = 1:10000
    % Keep the total number of trials for which the motif was found
    % constant
    % Randomly assign the label to these trials now and calculate the
    % difference 
    shuffled_id_vector = randsample(2,length(response1_motifs)+length(response2_motifs),true,[0.5 0.5]);
    response1_motifs_shuffled = original_motif_vector(shuffled_id_vector == 1);
    response2_motifs_shuffled = original_motif_vector(shuffled_id_vector == 2);

    shuffled_difference(iters) = sum(response1_motifs_shuffled) - sum(response2_motifs_shuffled);

end

p_val_uncorrected_greater(motif_number) = length((find(original_difference>shuffled_difference)==1));
p_val_uncorrected_greater(motif_number) = eps + (iters-p_val_uncorrected_greater(motif_number))/iters;

%% Response vs random null

% We will keep the total number of trials for a given response type
% constant and also the length for each trial and run the whole unique
% transitions pipeiline and then find motifs

total_response1_trials = 4938;
total_response2_trials = 8478-total_response1_trials;
trial_length = 1962;

% non markov
motifs_all_trials_random_response1 = null_sequence_motif_search(total_response1_trials,p);
motif_all_trial_random_response2 = null_sequence_motif_search(total_response2_trials,p);


%%  Markov testing
% Load group level model
load(fullfile(model_path,model_name))

% Choose Dirichlet priors from this model
dir_init_params = HMM_model.hmm.Dir_alpha;
dir_transition_params = HMM_model.hmm.Dir2d_alpha;
clear HMM_model

motifs_all_trials_random_response1 = null_sequence_motif_search(total_response1_trials,p,...
    true,dir_init_params,dir_transition_params);

motif_all_trial_random_response2 = null_sequence_motif_search(total_response2_trials,p,...
    true,dir_init_params,dir_transition_params);

% for timing purposes only
%  f = @() null_sequence_motif_search(total_response1_trials,p,...
%     true,dir_init_params,dir_transition_params);
% Tt = timeit(f)