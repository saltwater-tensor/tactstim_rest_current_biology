d = dir('Q:\tactstim_resting\code\tactstim_rest_data\subject_trial_specific_hmms\_04_Feb_2022_17_04_27_HMM_model_preStim_rest_250Hz\single_trial_hmms');

for i = 3:length(d)
    i
    load(fullfile(d(i).folder, d(i).name))
    vpath = hmm.vpath;
    vpath = double(vpath);
    hmm.vpath = vpath;
    save(( d(i).name),"hmm")
end