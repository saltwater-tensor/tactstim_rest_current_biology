function [data_embed_SUBJECT] = HMM_trial_wise_rearrange(X_embed, T_embed)

%% The script rearranges embedded data into trial wise epoch

data_embed_SUBJECT = {};
subject_T = T_embed{1,1};
subject_total_time = sum(subject_T);
subject_start = 1;
subject_end = subject_total_time;
subject_DATA = X_embed{1,1};
subj_count = 1;

for subj = 2:length(T_embed)+1 
    

        trial_start = 1;
        trial_end = (subject_T(1,1));
        trial_count = 1;

        for trial = 2:length(subject_T)+1

            data_trial = subject_DATA(trial_start:trial_end,:);
            data_embed_SUBJECT{subj_count,1}{trial_count,1} = data_trial;
            data_trial = [];
            trial_start = trial_end+1;
            trial_count = trial_count + 1;
            
            if trial_start > subject_total_time
                continue
            else
                trial_end = trial_end + subject_T(1,trial);
            end
        end

        if subj > length(T_embed)
            
%                 msgbox('Embedded data has been reshaped','RESHAPE')
            
            continue
        else
            subject_T = T_embed{1,subj};
            subject_total_time = sum(subject_T);
            subject_start = 1;
            subject_end = subject_total_time;
            subject_DATA = X_embed{subj,1};
            subj_count = subj_count + 1;
        end
    
end

end