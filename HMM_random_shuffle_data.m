source_model = 'PNAI';
trial_type = 'preStim_rest';
downsample_rate = 250;

clck = char(datetime);
clck = strrep(clck,':','_');
clck = strrep(clck,' ','_');
clck = strrep(clck,'-','_');



%% REMOVE TRIALS WITH SPECIFIC SOA CONDITIONS
soa_removal = 180;

for sub = 1:size(data_subject,2)
    
    current_sub_data = data_subject{1,sub};
    current_bhvr = BEHAVIOR{sub,2};
    current_time = T{1,sub};
    
    all_SOAs =  [current_bhvr(1:length(current_bhvr)).SOA_ms].';
    correct_SOAs = find(all_SOAs < soa_removal);
    
    current_sub_data_2 = current_sub_data(1,correct_SOAs);
    current_bhvr_2 = current_bhvr(1,correct_SOAs);
    current_time_2 = current_time(1,correct_SOAs);
    
    data_subject{1,sub} = current_sub_data_2;
    BEHAVIOR{sub,2} = current_bhvr_2;
    T{1,sub} = current_time_2;
    
    clear current_sub_data_2 current_bhvr_2 current_time_2 current_time current_bhvr current_sub_data
    
end




%% RANDOM SHUFFLE DATASET

for shuffle_sub = 1:size(data_subject,2)
    
    % Shuffle response1 and response2 trials within a subject so that the HMM is not
    % stuck in some local minima for a specific response type
    s = RandStream('dsfmt19937');
    trial_shuffle = randperm(s,size(BEHAVIOR{shuffle_sub,2},2));
    shuffled_behavior = BEHAVIOR{shuffle_sub,2};
    shuffled_behavior = shuffled_behavior(trial_shuffle);
    BEHAVIOR{shuffle_sub,2} = shuffled_behavior;
    shuffled_data = data_subject{1,shuffle_sub};
    shuffled_data = shuffled_data(trial_shuffle);
    data_subject{1,shuffle_sub} = shuffled_data';
    shuffled_T =  T{1,shuffle_sub};
    shuffled_T = shuffled_T(trial_shuffle);
    T{1,shuffle_sub} = shuffled_T;
    
    DATA{shuffle_sub,1} = cell2mat(data_subject{1,shuffle_sub});
    
end

%%

save(['shuffled_DATA_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'DATA','-v7.3')
save(['shuffled_Time_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'T','-v7.3')
save(['shuffled_Behavior_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'BEHAVIOR','-v7.3')
save(['shuffled_data_subject_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'data_subject','-v7.3')
save(['shuffled_sampling_frequency_' source_model '_trial_' trial_type '_' num2str(downsample_rate) '_Hz_' 'created_on_' clck],'downsample_rate','-v7.3')
