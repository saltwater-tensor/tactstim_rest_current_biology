%% Load data and time
dataset_path = 'Q:\tactstim_resting\code\tactstim_rest_data\sign_flip_corrected_shuffled_DATA_PNAI_trial_preStim_rest_250_Hz_created_on_04_Feb_2022_13_01_09';
load(dataset_path)
time_path = 'Q:\tactstim_resting\code\tactstim_rest_data\shuffled_Time_PNAI_trial_preStim_rest_250_Hz_created_on_04_Feb_2022_13_01_09';
load(time_path)
behavior_path = 'Q:\tactstim_resting\code\tactstim_rest_data\shuffled_Behavior_PNAI_trial_preStim_rest_250_Hz_created_on_04_Feb_2022_13_01_09';
freq_path = 'Q:\tactstim_resting\code\tactstim_rest_data\shuffled_sampling_frequency_PNAI_trial_preStim_rest_250_Hz_created_on_04_Feb_2022_13_01_09.mat';
load(freq_path)
OutputPath = ['Q:\tactstim_resting\code\tactstim_rest_models' filesep];

%%
tic
options = [];

%% Hmm options
if isempty(options)
    
    if ~isempty(options)
        K_states = options.K;
    else
        K_states = input('SUPPLY THE APPROPRIATE NUMBER OF STATES!!');
    end

    model_name = input('Enter a name for your model for eg, standard_trials_100Hz');
    % filename = dataset_name;
    
    
    run_mar = 0;
    options.K = K_states ;%input('enter desired num of states');
    try
        options.Fs = sampling_freq; %input('enter sampling freq in Hz');
    catch
        options.Fs = downsample_rate;
    end
    if run_mar
        run_mar_tag = 'MAR';
        options.order = 11; % (sampling_freq/oreder = smallest_detectable_frequency) set to detect ateast beta
    else
        run_mar_tag = 'NO_MAR';
        options.order = 0; % (sampling_freq/oreder = smallest_detectable_frequency) set to detect ateast beta
    end

    options.timelag = 2;
    options.orderoffset = 1;
    options.covtype = 'full';
    options.zeromean = 1;
    options.DirichletDiag = 10; %default
    options.pca = 84; 
    options.standardise = 1;
    options.onpower = 0;
    options.filter = [];
    options.detrend = 0;
    options.downsample = 0;
    options.symmetricprior = 1; 
    options.dropstates = 1;
    options.tol = 1e-5;
    options.cyc = 1000;
    options.inittype = 'hmmmar';
    options.initrep = 5;
    options.initcyc = 100;
    options.repetitions = 1;
    options.updateGamma = 1;
    options.decodeGamma = 1;
    options.useParallel = 0;
    options.embeddedlags = -7:7; %[-10 8 6 4 2 0 2 4 6 8 10];
    options.embedd_lag_tags = 'embed_lags';
    
    % Stochastic gradient paramters
    options.BIGNbatch = 6; 
    options.BIGNinitbatch = options.BIGNbatch ; options.BIGtol = 1e-7; options.BIGcyc = 100; % or 
    options.BIGundertol_tostop = 5; options.BIGdelay = 1; 
    options.BIGforgetrate = 0.7; options.BIGbase_weights = 0.9;
    
    
end


%% 
% Saving options
clck = char(datetime);
clck = strrep(clck,':','_');
clck = strrep(clck,' ','_');
clck = strrep(clck,'-','_');
save(strcat(OutputPath,'_',clck,'_',model_name,'_','options_file'),'options')

addpath(genpath('Q:\software\HMM_MAR_2019'));

[hmm, Gamma, Xi, vpath, GammaInit, residuals, fehist, feterms, rho] = hmmmar (data,T,options);

fe = hmmfe(data,T,hmm,Gamma,Xi);
HMM_model.hmm = hmm;
HMM_model.Gamma = Gamma;
HMM_model.Xi = Xi;
HMM_model.vpath = vpath;
HMM_model.GammaInit = GammaInit;
HMM_model.residuals = residuals;
HMM_model.fehist = fehist;
HMM_model.feterms = feterms;
HMM_model.rho = rho;
sign_flag = 1;
HMM_model.sign_corrected = sign_flag;
HMM_model.dataset_name = dataset_path;
HMM_model.dataset_location = pwd;
HMM_model.free_energy = fe;
HMM_model.options = options;
HMM_model.dataset_path = dataset_path;
HMM_model.model_path = OutputPath;
HMM_model.time_path = time_path;
HMM_model.behavior_path = behavior_path;
HMM_model.freq_path = freq_path;

save(strcat(OutputPath,'_',clck,'_','HMM_model_',model_name),...
    'HMM_model','-v7.3')

end_time = toc;
save(strcat(OutputPath,'_',clck,'_',model_name,'_','total_run_time'),'end_time')

rmpath(genpath('Q:\software\HMM_MAR_2019'));