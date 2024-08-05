%% HMM analysis
%%
% Add necessary code libs
addpath(genpath('Q:\tactstim_resting\code\tacstim_rest_code'));


%%
% Load model 
% Load dataset (Should contain both time, data and sampling frequency
% information)
model_path = 'Q:\tactstim_resting\code\tactstim_rest_models\';
model_name = '_04_Feb_2022_17_04_27_HMM_model_preStim_rest_250Hz.mat';
mdl_nm = strsplit(model_name,'.');
load(fullfile(model_path, model_name))
load(HMM_model.dataset_path)
load(HMM_model.time_path)
load(HMM_model.freq_path)
load(HMM_model.behavior_path)

%% Basic HMM properties

[vpath] = HMM_model.vpath;%hmmdecode(data,T,HMM_model.hmm,1);
LifeTimes = getStateLifeTimes (vpath,T,HMM_model.options,5);
Intervals = getStateIntervalTimes (vpath,T,HMM_model.options,5);
FO_2 = getFractionalOccupancy (HMM_model.Gamma,T,HMM_model.options,2);
FO_1_start_aligned = getFractionalOccupancy (HMM_model.Gamma,T,HMM_model.options,1,'start');
maxFO = getMaxFractionalOccupancy(HMM_model.Gamma,T,HMM_model.options);
switchingRate =  getSwitchingRate(HMM_model.Gamma,T,HMM_model.options);

HMM_props.LifeTimes = LifeTimes;
HMM_props.Intervals = Intervals;
HMM_props.FO_2 = FO_2;
HMM_props.FO_1_start_aligned = FO_1_start_aligned;
HMM_props.maxFO = maxFO;
HMM_props.switchingRate = switchingRate;
save([mdl_nm{1} ,'_','HMM_properties'],'HMM_props')



%% Distance between states
total_states = HMM_model.hmm.K;
C_all = zeros(total_states,total_states);
for k = 1:1:total_states
    % Takes the low dimensional pcs_dims by pc_dims matrix 
    % projects it back to the high dimensional num_lags X contacts space
    % using PC loading matrix A 
    C1 = getAutoCovMat(HMM_model.hmm,k);
    C_all_covs(:,:,k) = C1;
    for k2 = 1:1:total_states
          C2 = getAutoCovMat(HMM_model.hmm,k2);
          [E1] = eig(C1,C2);
          logE1 = log(E1);
          sqr = logE1 .* logE1;
          sumsqr = sum(sqr);
          sqrtsumsqr = sqrt(sumsqr);
          abssq = abs(sqrtsumsqr);
          dist1 = abssq;
          if k == k2
              C_all(k,k2) = 0;
          else
              C_all(k,k2) = dist1;
          end 
    end   
end

figure(26)
heatmap(C_all,'Colormap',jet)
% save([mdl_nm{1} ,'_','covariance_distances'],'C_all_covs')

%% Spectral properties

% For data driven spectral analysis please remove awkward trials or short trials
% The short trials were also removed for training the HMM but are present
% in the dataset for completeness sake
Gamma = HMM_model.Gamma;
acc = 0 ;
d = length(HMM_model.hmm.train.embeddedlags) - 1;
gamma_new = [];
for n = 1: size(data,1)
    
    gamma = Gamma(acc + (1:(sum(T{n})-length(T{n})*d)),:);
    acc = acc + size(gamma,1);
    gamma_update = gamma_adjust(gamma,T(n),d);

    [tuda_data_updated, tuda_time, tuda_Y] = ...
        data_time_behavior_update(data(n), T(n), BEHAVIOR(n,:));

    data(n) = {tuda_data_updated};
    T(n) = {tuda_time};
    gamma_new{n,1} = gamma_update;
end
gamma_new = cell2mat(gamma_new);

% A safety check
diff_wrapper = @(x)  sum(x)-length(x)*d;
T_check = cellfun(diff_wrapper, T,'UniformOutput',false);
if ~size(gamma_new,1) == sum(cell2mat(T_check))
    error('Please check previous steps!!!')
end

%% Running multitaper analysis
sampling_freq = downsample_rate;
options_mt = struct('Fs',sampling_freq); % Sampling rate
options_mt.fpass = [1 45];  % band of frequency you're interested in
options_mt.tapers = [4 7]; % taper specification - leave it with default values
% options_mt.win = 2 * sampling_freq; % multitaper window
options_mt.to_do = [1 0]; % turn off pdc
options_mt.embeddedlags = HMM_model.hmm.train.embeddedlags;

% Calculating on a single subject basis
N = size(data,1);
fitmt_subj = cell(N,1);
d = length(options_mt.embeddedlags) - 1;
acc = 0 ;
N = size(data,1);
if exist('gamma_new','var')
    Gamma = gamma_new;
else
    Gamma = HMM_model.Gamma;
end

for n = 1:N
    X = data{n};
    gamma = Gamma(acc + (1:(sum(T{n})-length(T{n})*d)),:);
    acc = acc + size(gamma,1);
    fitmt_subj{n} = hmmspectramt(X,T{n},gamma,options_mt);
%     fitmt_subj{n}.state = rmfield(fitmt_subj{n}.state,'ipsd');
%     fitmt_subj{n}.state = rmfield(fitmt_subj{n}.state,'pcoh');
%     fitmt_subj{n}.state = rmfield(fitmt_subj{n}.state,'phase');
    disp(['Subject ' num2str(n)])
end
HMM_spectra.fitmt_subj = fitmt_subj;

% Calculating on the group level
T = cell2mat(T');
% T = T';
data = cell2mat(data);
fitmt_group = hmmspectramt(data,T,HMM_model.Gamma,options_mt);
HMM_spectra.fitmt_group = fitmt_group;

save([mdl_nm{1} ,'_','HMM_spectra'],'HMM_spectra')

%% Standard frequency band analysis
% We could also perform data driven spectral analysis
% https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#-factorising-the-states-spectra-into-frequency-bands

% Straightforward spectral analysis based on classical frequency band
% definition
bands = [1 7; 8 12; 12 30; 40 Inf];
sp_fit_group_bands = spectbands(HMM_spectra.fitmt_group,bands);

for sb_sp = 1:length(HMM_spectra.fitmt_subj)
    subj = HMM_spectra.fitmt_subj{sb_sp,1};
    sp_fit_subj_bands{sb_sp,1} = spectbands(subj,bands);
    
end

%% Data driven frequency band analysis
for sbj = 1:length(HMM_spectra.fitmt_subj)
    ln(sbj) = length(HMM_spectra.fitmt_subj{sbj, 1}.state(1).f);
end

options_fact = struct();
options_fact.Ncomp = 4; 
options_fact.Base = 'coh';
[fitmt_subj_fact_4b,fitmt_group_fact_4b,supply_profiles] = spectdecompose(HMM_spectra.fitmt_subj,options_fact);

%% You can also save the matrix to factorise

chan = [1:42];
chan1 = [];
chan2 = [];
supply_profiles = [];
[X] = ...
    spectdecompose_custom_updated(fitmt_subj,options_fact,...
    chan,chan1,chan2,supply_profiles);

%%
figure(9022)
[supply_profiles,H] = nnmf(X,4,'algorithm','als');
plot(supply_profiles,'DisplayName','supply_profiles')

%%

[X,fitmt_subj_fact_4b,fitmt_group_fact_4b,sp_profiles_4b] = ...
    spectdecompose_custom_updated(HMM_spectra.fitmt_subj,options_fact,chan,chan1,chan2,supply_profiles);
  
  
%% Visualisation

% Load default anatomy by changing the path
load('Q:\tactstim_resting\New folder\brainstorm_db\tactstim_resting\anat\@default_subject\tess_cortex_pial_low')
output_path_spectral_figs = ...
    'Q:\tactstim_resting\code\tactstim_rest_models\_04_Feb_2022_17_04_27_HMM_model_preStim_rest_250Hz_spectral_analysis\nnmf_profiles_and_group_projections_3\0.001\';
if ~exist(output_path_spectral_figs,"dir")
    mkdir(output_path_spectral_figs)
end
% Load nnmf profiles and projections
load('nnmf_profiles_and_group_projections_3.mat')
% Call script for visualisation
% Save figures that have been generated
HMM_visualisation

%% For viterbi and gamma sequence analysis on the group level 
% HMM_analysis_2

%% For permutation testing of spectral matrices obtained
% HMM_analysis_3
