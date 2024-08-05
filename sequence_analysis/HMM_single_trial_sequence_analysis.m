%%

%% paths
model_path = 'Q:\tactstim_resting\code\tactstim_rest_models\';
model_name = '_04_Feb_2022_17_04_27_HMM_model_preStim_rest_250Hz.mat';
dataset_path = 'Q:\tactstim_resting\code\tactstim_rest_data\';
mdl_nm = strsplit(model_name,'.');
load([model_path model_name]);
load(HMM_model.behavior_path)
output_path = 'Q:\tactstim_resting\code\tactstim_rest_data\subject_trial_specific_hmms\';
folder_type = 'single_trial_hmms';
base_model_name = mdl_nm{1, 1};
output_path = fullfile(output_path, base_model_name, folder_type);
output_path = [output_path '\'];

p_value_to_test = 0.05;
%% 
hmm_files = dir(output_path);
p_resp1 = [];
p_resp2 = [];

cnt_resp1 = 1;
cnt_resp2 = 1;

for dr = 1:length(dir(output_path))
    dr
    if ~hmm_files(dr).isdir
        load(fullfile(hmm_files(dr).folder, hmm_files(dr).name))
        vpath = hmm.vpath;
        T = length(vpath);
        gamma = hmm.Gamma_subj_trial;

        %
        subj_trial_string = strsplit(hmm_files(dr).name,'_');
        subj = str2num(subj_trial_string{3});
%         SOA = BEHAVIOR{subj, 2}(trials).Response;

        LifeTimes = getStateLifeTimes (vpath,T,hmm.hmm_subj_trial.train);
        LifeTimes = cellfun(@transpose, LifeTimes, 'UniformOutput', false)';
        
        Intervals = getStateIntervalTimes (vpath,T,hmm.hmm_subj_trial.train);
        Intervals = cellfun(@transpose,Intervals, 'UniformOutput', false)';
        
        switch_rate = getSwitchingRate(gamma,T,hmm.hmm_subj_trial.train);
        FO = getFractionalOccupancy(gamma,T,hmm.hmm_subj_trial.train);
               
        nm = strsplit(hmm_files(dr).name,'_');
        nm = nm{6};
        nm = strsplit(nm,'.');
        nm = nm{1}; 
        
        switch nm
            case 'Response1'
                p_resp1 = cat(3,p_resp1,hmm.hmm_subj_trial.P);
                LifeTimes_resp1{cnt_resp1,:} = LifeTimes;
                Intervals_resp1{cnt_resp1,:} = Intervals;
                FO_resp1(cnt_resp1,:) = FO;
                switch_rate_resp1(cnt_resp1,:) = switch_rate;
                cnt_resp1 = cnt_resp1 + 1;
                
            case 'Response2'
                p_resp2 = cat(3,p_resp2,hmm.hmm_subj_trial.P);
                LifeTimes_resp2{cnt_resp2,:} = LifeTimes;
                Intervals_resp2{cnt_resp2,:} = Intervals;
                FO_resp2(cnt_resp2,:) = FO;
                switch_rate_resp2(cnt_resp2,:) = switch_rate;
                cnt_resp2 = cnt_resp2 + 1;
                
        end
    end
end

%% Exploratory histrograms for transition probabilities
% Permutation testing RESPONSE 1 VS RESPONSE 2
% No within response testing takes place
close all
figure(1)
plot_num = 1;
hold on
from_state = 1;
to_state = 5;
% vec1 = squeeze(p_resp1(from_state,to_state,:));
% vec2 = squeeze(p_resp2(from_state,to_state,:));
% [] = ttest2(vec1, vec2, "Vartype","unequal")
for from_state = 1:6
    for to_state = 1:6
        subplot(6,6,plot_num)
        histogram(p_resp1(from_state,to_state,:),100,"Normalization","pdf")
        hold on
        histogram(p_resp2(from_state,to_state,:),100,"Normalization","pdf")
        vec1 = squeeze(p_resp1(from_state,to_state,:));
        vec2 = squeeze(p_resp2(from_state,to_state,:));
        [h_t, p_t] = ttest2(vec1, vec2, "Vartype","unequal","Alpha",0.05);
        [p_perm, observeddifference, effectsize] = permutationTest(vec1, vec2, 5000);
        p_tbl_ttest(from_state,to_state) = p_t;
        p_tbl_permtest(from_state,to_state) = p_perm;
        plot_num = plot_num + 1;
        
        % Sided test
        [p_perm_larger, observeddifference, effectsize] = permutationTest(vec1, vec2, 5000,'sidedness','larger');
        [p_perm_smaller, observeddifference, effectsize] = permutationTest(vec1, vec2, 5000,'sidedness','smaller');
        p_tbl_permtest_larger(from_state,to_state) = p_perm_larger;
        p_tbl_permtest_smaller(from_state,to_state) = p_perm_smaller;
    end
end
% Correcting for multiple comparisons two sided tests
p_ttest_mafdr = mafdr(p_tbl_ttest(:),'BHFDR',true);
p_ttest_mafdr = reshape(p_ttest_mafdr,[6,6]);
p_permtest_mafdr = mafdr(p_tbl_permtest(:),'BHFDR',true);
p_permtest_mafdr = reshape(p_permtest_mafdr,[6,6]);

p_disp_ttest = p_tbl_ttest;
p_disp_permtest = p_tbl_permtest;
p_disp_ttest(p_disp_ttest > p_value_to_test) = NaN;
p_disp_permtest(p_disp_permtest > p_value_to_test) = NaN;

figure(2)
heatmap(p_disp_ttest,'ColorMap',jet)
figure(3)
heatmap(p_disp_permtest,'ColorMap',jet)

p_disp_ttest_mafdr = p_ttest_mafdr;
p_disp_permtest_mafdr = p_permtest_mafdr;
p_disp_ttest_mafdr(p_disp_ttest_mafdr > p_value_to_test) = NaN;
p_disp_permtest_mafdr(p_disp_permtest_mafdr > p_value_to_test) = NaN;

figure('Name','T test two sided P matrix','NumberTitle','off') 
heatmap(p_disp_ttest_mafdr,'ColorMap',jet)

figure('Name','Permutation test two sided P matrix','NumberTitle','off') 
heatmap(p_disp_permtest_mafdr,'ColorMap',jet)

% Correcting for multiple comparisons one sided tests 
p_permtest_mafdr_larger = mafdr(p_tbl_permtest_larger(:),'BHFDR',true);
p_permtest_mafdr_larger = reshape(p_permtest_mafdr_larger,[6,6]);

p_permtest_mafdr_smaller = mafdr(p_tbl_permtest_smaller(:),'BHFDR',true);
p_permtest_mafdr_smaller = reshape(p_permtest_mafdr_smaller,[6,6]);

p_disp_permtest_mafdr_larger = p_permtest_mafdr_larger;
p_disp_permtest_mafdr_larger(p_disp_permtest_mafdr_larger > p_value_to_test) = NaN;
p_disp_permtest_mafdr_smaller = p_permtest_mafdr_smaller;
p_disp_permtest_mafdr_smaller(p_disp_permtest_mafdr_smaller > p_value_to_test) = NaN;


figure('Name','trans prob R1 > R2 larger perm test','NumberTitle','off') 
heatmap(p_disp_permtest_mafdr_larger,'ColorMap',jet)
xlabel('to state')
ylabel('from state')
figure('Name','trans prob R1 < R2 smaller perm test','NumberTitle','off') 
heatmap(p_disp_permtest_mafdr_smaller,'ColorMap',jet)
xlabel('to state')
ylabel('from state')

%% lifetime 
 
 lftimes_mat_resp1 = rearrange_temporal_cell(LifeTimes_resp1);
 lftimes_mat_resp2 = rearrange_temporal_cell(LifeTimes_resp2);
 p_lftimes = temporal_test(lftimes_mat_resp1, lftimes_mat_resp2,p_value_to_test);   
 
%% intervals
 intervals_mat_resp1 = rearrange_temporal_cell(Intervals_resp1);
 intervals_mat_resp2 = rearrange_temporal_cell(Intervals_resp2);
 p_intrvl_times = temporal_test(intervals_mat_resp1, intervals_mat_resp2,p_value_to_test);
 
%% FO
p_FO = temporal_test(FO_resp1, FO_resp2,p_value_to_test);

%% switching rate
% Switching rate per minute.
switch_rate_resp1 = (switch_rate_resp1/8)*60;
switch_rate_resp2 = (switch_rate_resp2/8)*60;

switch_rate_resp1 = (switch_rate_resp1/60)*8;
switch_rate_resp2 = (switch_rate_resp2/60)*8;

%
switch_rate_resp1 = switch_rate_resp1/1947;
switch_rate_resp2 = switch_rate_resp2/1947;

switch_rate_resp1 = switch_rate_resp1/4;
switch_rate_resp2 = switch_rate_resp2/4;


switch_rate_resp1 = switch_rate_resp1*1947;
switch_rate_resp2 = switch_rate_resp2*1947;



[p_perm, observeddifference, effectsize] = permutationTest(switch_rate_resp1,switch_rate_resp2, 5000);
[p_perm_smaller_swtch, observeddifference, effectsize] = permutationTest(switch_rate_resp1,switch_rate_resp2, ...
    5000,'sidedness','smaller');
fprintf('%0.11f\n',p_perm_smaller_swtch)


%% ANOVA TESTING

%% Intervals 
clear Y_resp1 Y_state_resp1 Y_response_resp1
clear Y_resp2 Y_state_resp2 Y_response_resp2

[Y_resp1_intervals,Y_state_resp1_intervals,Y_response_resp1_intervals] = creat_anovan_variables(intervals_mat_resp1,'resp1');

[Y_resp2_intervals,Y_state_resp2_intervals,Y_response_resp2_intervals] = creat_anovan_variables(intervals_mat_resp2,'resp2');

Y_intervals = [Y_resp1_intervals;Y_resp2_intervals];
Y_States_intervals = [Y_state_resp1_intervals;Y_state_resp2_intervals];
Y_Response_intervals = [Y_response_resp1_intervals;Y_response_resp2_intervals];
[p,tbl,stats_intervals,terms] = anovan(Y_intervals,{Y_States_intervals Y_Response_intervals},'model','full','varnames',{'States','Response'});

figure(2)
[results_all,m_intervals,~,gnames_all] = multcompare(stats_intervals,'Dimension',[1 2]);

figure(3)
[results_states,~,~,gnames_states] = multcompare(stats_intervals,'Dimension',[1]);

figure(4)
[results_response,~,~,gnames_response] = multcompare(stats_intervals,'Dimension',[2]);

tbl = array2table(results_all,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A")=gnames_all(tbl.("Group A"));
tbl.("Group B")=gnames_all(tbl.("Group B"));

%% Lifetimes
clear Y_resp1 Y_state_resp1 Y_response_resp1
clear Y_resp2 Y_state_resp2 Y_response_resp2

[Y_resp1_lifetimes,Y_state_resp1_lifetimes,Y_response_resp1_lifetimes] = creat_anovan_variables(lftimes_mat_resp1,'resp1');

[Y_resp2_lifetimes,Y_state_resp2_lifetimes,Y_response_resp2_lifetimes] = creat_anovan_variables(lftimes_mat_resp2,'resp2');

Y_lifetimes = [Y_resp1_lifetimes;Y_resp2_lifetimes];
Y_States_lifetimes = [Y_state_resp1_lifetimes;Y_state_resp2_lifetimes];
Y_Response_lifetimes = [Y_response_resp1_lifetimes;Y_response_resp2_lifetimes];
[p,tbl,stats_lifetimes,terms] = anovan(Y_lifetimes,{Y_States_lifetimes Y_Response_lifetimes},'model','full','varnames',{'States','Response'});

figure(2)
[results_all,m_lifetimes,~,gnames_all] = multcompare(stats_lifetimes,'Dimension',[1 2]);

figure(3)
[results_states,~,~,gnames_states] = multcompare(stats_lifetimes,'Dimension',[1]);

figure(4)
[results_response,~,~,gnames_response] = multcompare(stats_lifetimes,'Dimension',[2]);

tbl = array2table(results_all,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A")=gnames_all(tbl.("Group A"));
tbl.("Group B")=gnames_all(tbl.("Group B"));

%% FO
clear Y_resp1 Y_state_resp1 Y_response_resp1
clear Y_resp2 Y_state_resp2 Y_response_resp2

[Y_resp1_FO,Y_state_resp1_FO,Y_response_resp1_FO] = creat_anovan_variables(FO_resp1,'resp1');

[Y_resp2_FO,Y_state_resp2_FO,Y_response_resp2_FO] = creat_anovan_variables(FO_resp2,'resp2');

Y_FO = [Y_resp1_FO;Y_resp2_FO];
Y_States_FO = [Y_state_resp1_FO;Y_state_resp2_FO];
Y_Response_FO = [Y_response_resp1_FO;Y_response_resp2_FO];
[p,tbl,stats_FOs,terms] = anovan(Y_FO,{Y_States_FO Y_Response_FO},'model','full','varnames',{'States','Response'});

figure(2)
[results_all,m_FO,~,gnames_all] = multcompare(stats_FOs,'Dimension',[1 2]);

figure(3)
[results_states,~,~,gnames_states] = multcompare(stats_FOs,'Dimension',[1]);

figure(4)
[results_response,~,~,gnames_response] = multcompare(stats_FOs,'Dimension',[2]);


tbl = array2table(results_all,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A")=gnames_all(tbl.("Group A"));
tbl.("Group B")=gnames_all(tbl.("Group B"));

%% Figures

% We will plot statewise due to a large difference in scale. Same scale
% plots obscures visuals
m_var = (4*m_intervals([2,8],:));
plot_temporal_props(m_var)

% plot for switching rate
switch_rate_resp1 = switch_rate_resp1*10000;
switch_rate_resp2 = switch_rate_resp2*10000;
m_var(1,1) = mean(switch_rate_resp1);
m_var(1,2) = var(switch_rate_resp1)/sqrt(length(switch_rate_resp1));
m_var(2,1) = mean(switch_rate_resp2);
m_var(2,2) = var(switch_rate_resp2)/sqrt(length(switch_rate_resp2));
plot_temporal_props(m_var)

