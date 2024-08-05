motifs_all_trials_random_response1 = null_sequence_motif_search(total_response1_trials,p,...
    true,dir_init_params,dir_transition_params);

%% for timing purposes only
 f = @() null_sequence_motif_search(total_response1_trials,p,...
    true,dir_init_params,dir_transition_params);
Tt = timeit(f)