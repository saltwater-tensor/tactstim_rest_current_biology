
clearvars -except T BEHAVIOR GAMMA_SUBJECT VPATH_SUBJECT HMM_model
X = [];
Y = [];
model_name = '_06_Jan_2022_10_00_37_HMM_model_preStim_rest_250Hz.mat';
mdl_nm = strsplit(model_name,'.');
load(model_name)
%% We will not use short trials
se_time = cellfun(@transpose,T,'UniformOutput',0);
se_time = cell2mat(se_time');
max_time = max(se_time);
max_time_embed = max_time - 14;
%%
for subj_bh = 1:size(BEHAVIOR,1)
    
    bhvr = BEHAVIOR{subj_bh,2}';
    t = T{subj_bh}';
    
    response1 = find(ismember({bhvr.Response}, {'Response1'}));
    response1_times = t(response1);
    response1_gammas = GAMMA_SUBJECT{subj_bh,1}(response1,1);
    response1_shorts = find(response1_times < max_time);
    response1_times(response1_shorts) = [];
    response1_gammas(response1_shorts) = [];
    
    
    response2 = find(ismember({bhvr.Response}, {'Response2'}));
    response2_times = t(response2);
    response2_gammas = GAMMA_SUBJECT{subj_bh,1}(response2,1);
    response2_shorts = find(response2_times < max_time);
    response2_times(response2_shorts) = [];
    response2_gammas(response2_shorts) = [];
    
    
    
    A1 = cell2mat(response1_gammas);
    B1 = reshape(A1,[size(response1_gammas,1),max_time_embed*HMM_model.hmm.K]);
    Y1 = ones(size(B1,1),1);

    A2 = cell2mat(response2_gammas);
    B2 = reshape(A2,[size(response2_gammas,1),max_time_embed*HMM_model.hmm.K]);
    Y2 = 2*ones(size(B2,1),1);

    X = [X;B1;B2];
    Y = [Y;Y1;Y2];

end

%%
prd_nm = 1;
for  prd_nm_2 = 1:HMM_model.hmm.K
    for prd_nm_1 = 1:1962
        Predictor_Name{1,prd_nm} = ['State_' num2str(prd_nm_2) '_T_' num2str(prd_nm_1)];
        prd_nm = prd_nm + 1;
    end   
end 
%% Directly fit a model
X_zscore = zscore(X);

Mdl_cv_1 = fitctree(X_zscore,Y,'PredictorSelection','curvature',...
    'CrossVal','on','NumVariablesToSample','all','PredictorName',Predictor_Name,...
    'MaxNumSplits',500,'MinLeafSize',2);

 view(Mdl_cv_1.Trained{1},'Mode','graph')
 classErrorDefault = kfoldLoss(Mdl_cv_1);
 imp_1 = predictorImportance(Mdl_cv_1.Trained{1});

 Mdl_2 = fitctree(X_zscore,Y,'PredictorSelection','curvature',...
   'PredictorName',Predictor_Name,...
    'MaxNumSplits',500,'MinLeafSize',2);
 imp = predictorImportance(Mdl_2);

 imp = reshape(imp,[1962,6]);
 figure
 hold on
 bar(imp(:,1))
 bar(imp(:,2))
 bar(imp(:,3))
 bar(imp(:,4))
 bar(imp(:,5))
 bar(imp(:,end))
 ylabel('Estimates');
 xlabel('Predictors');
 legend on

%% Hyperparameter optimisation single decision tree
clearvars -except X X_zscore Y Predictor_Name
rng(1)

% Parameters to optimize for the decision tree
params = hyperparameters('fitctree',X_zscore,Y);
% Max Splits for the decision tree
params(2).Optimize = true;

% Not optimising min leaf size
params(1).Optimize = true;
% params(1).Range = [1 10];

% Hyperparameter optimization options
c = cvpartition(length(Y),'KFold',10) ;
optimizeParams_options = struct(...
    'Optimizer','bayesopt',...
    'MaxObjectiveEvaluations',30,...
    'CVpartition',c...
    );

tree_Mdl_hyperoptimize = fitctree(X_zscore,Y,'PredictorNames',...
    Predictor_Name,'OptimizeHyperparameters', params,...
    'HyperparameterOptimizationOptions',optimizeParams_options);

%% Hyperparameter optimisation for an ensemble

% Single tree does not optimize for the number of predictor variables to
% use
clearvars -except X Y Predictor_Name

% Parameters to optimize for the ensemble
params = hyperparameters('fitcensemble',X,Y,'Tree');

params(1).Optimize = true;
% Not optimising min leaf size
params(4).Optimize = false;
% Max Splits for the decision tree
params(5).Range = [1 20];
params(5).Optimize = false;
% Number of variables to select
% Max Splits for the decision tree
params(7).Range = [1 11772];
params(7).Optimize = false;

optimizeParams_options = struct(...
    'Optimizer','bayesopt',...
    'MaxObjectiveEvaluations',70....,
    );

ensemble_Mdl_hyperoptimize = fitcensemble(X,Y,'PredictorNames',...
    Predictor_Name,'OptimizeHyperparameters', params,'Learners','tree',...
    'HyperparameterOptimizationOptions',optimizeParams_options);


%  Mdl = fitcensemble(X,Y,"Learners","tree");
