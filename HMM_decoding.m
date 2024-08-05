%% HMM decoding analysis
dataset_path = 'Q:\tactstim_resting\code\tactstim_rest_data\';
dataset_path = [dataset_path ...
    'sign_flip_corrected_shuffled_DATA_PNAI_trial_preStim_rest_250_Hz_created_on_03_Jan_2022_10_33_02.mat'];
XamObj = matfile(dataset_path);
XamObj.Properties.Writable = true;
Xamsize = size(XamObj.data,1);
clck = char(datetime);
clck = strrep(clck,':','_');
clck = strrep(clck,' ','_');
clck = strrep(clck,'-','_');
Output_path = 'Q:\tactstim_resting\code\tactstim_rest_data\decoding_data\';
Output_path = [Output_path clck '\'];
mkdir([Output_path]);
addpath('Q:\tactstim_resting\code\tactstim_rest_code\')



%% load model data and time
model_directory = 'Q:\tactstim_resting\code\tactstim_rest_models\';
model_name = '_06_Jan_2022_10_00_37_HMM_model_preStim_rest_250Hz.mat';
mdl_nm = strsplit(model_name,'.');
load([model_directory model_name])
% load(HMM_model.dataset_path)
load(HMM_model.time_path)
load(HMM_model.behavior_path)
T_full = T;
BEHAVIOR_full = BEHAVIOR;

%% 
for xms = 1:Xamsize
    
    data = XamObj.data(xms,:);
    T = T_full(1, xms);
    BEHAVIOR = BEHAVIOR_full(xms,:);
    %% Obtain the high dim PCA matrix, covariance matrix, and time delay embedded dataset
   
    
    options = HMM_model.options;
    options = format_options(options);
    
    [X_embed,T_embed,A_high_dim,A,B,C,e] = highdim_pca_custom(data,T,options.pca,...
        options.embeddedlags,options.standardise,...
        options.onpower,options.varimax,options.detrend,...
        options.filter,options.leakagecorr,options.Fs,options.As);
    options.A = A;
    options.embeddedlags = [];
    
    % data has been standardised i.e. mean removed and divided by standard
    % deviation on a trial by trial basis for each subject 
    % We do not need to do it again for subsequent projections and matrix
    % calculations
    options.standardise = 0;
    clear data T
    
    %% Get high dimensional covariance matrices of HMM states from your loaded model
    
    total_states = HMM_model.hmm.K;
    Cov_mats = NaN([total_states, HMM_model.hmm.train.pca, HMM_model.hmm.train.pca]);
    for k = 1:1:total_states
        % Takes the low dimensional pcs_dims by pc_dims matrix 
        % projects it back to the high dimensional num_lags X contacts space
        % using PC loading matrix A 
        Cov_mats(k,:,:) = getAutoCovMat(HMM_model.hmm,k);
    end
    
    %% Removing short trials and get back response variable
    
    for xem = 1:length(X_embed)
        [tuda_data_updated, tuda_time, tuda_Y] = data_time_behavior_update(X_embed(xem,1), ...
            T_embed(1, xem), BEHAVIOR(xem,:));
        X_embed_update{xem,1} = tuda_data_updated;
        T_embed_update{1,xem} = tuda_time';
        Y_embed_update{xem,1} = tuda_Y;
        X_embed{xem,1} = [];
        T_embed{1,xem} = [];
    end
    
    clear X_embed T_embed
     
    %% Projecting to embedded PC space and calculating activations
    %  Calcualte activation vector
    % At a time point t activation vector is calculated as follows:
    % Bt = Ndimensional vector at time t which is the delay embedded pc
    % projected vector
    % Cm = Ndimensional X Ndimensional cross covariance matrix that defines an
    % HMM state
    % Avt = Activation vector = Cm*Bt
    
    for xx = 1:length(X_embed_update)
        
        % Projecting to embedded PC space
        % Each PC projection is then standardised on trial by trial basis
        [X_embed_pca_file, T_embed_pca_file] = loadfile_custom(X_embed_update{xx,1}, T_embed_update{1,xx}, options);
    %     X_embed_pca{xx,1} = X_embed_pca_file;
    %     T_embed_pca{1,xx} = T_embed_pca_file;
        X_embed_pca_activations_file = [];
        for cvm = 1:total_states
            % Calculation of activation vector
            X_embed_pca_activations_file(:,:,cvm) = (squeeze(Cov_mats(cvm,:,:))*X_embed_pca_file')';
        end
        X_embed_pca_activations{xx,1} = X_embed_pca_activations_file;
        T_embed_pca_activations{1,xx} = T_embed_pca_file;
    %     clear X_embed_pca_activations_file
        
    end
    
    
    %% Trial wise rearrangement of embedded dataset
    
    data_activations_statewise_subject = cell(length(T_embed_pca_activations),total_states);
    
    for sbj = 1:size(T_embed_pca_activations,2)
        for cvm = 1:total_states
        
          data_embed_pca_SUBJECT_state =  HMM_trial_wise_rearrange({squeeze(X_embed_pca_activations{sbj, 1}(:,:,cvm))}, ...
              T_embed_pca_activations(1,sbj));
          data_activations_statewise_subject{sbj, cvm} = data_embed_pca_SUBJECT_state;
        end
    %     [data_embed_pca_SUBJECT] = HMM_trial_wise_rearrange(X_embed_pca, T_embed_pca);
    end
    % clear X_embed_pca X_embed data
    save([Output_path 'Decode_Subject_' num2str(xms)], ...
        'data_activations_statewise_subject','Y_embed_update', ...
        'T_embed_pca_activations','-v7.3');
    clear data_activations_statewise_subject Y_embed_update T_embed_pca_activations
end