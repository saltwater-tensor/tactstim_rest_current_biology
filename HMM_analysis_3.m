%% Permutation testing

%% Within state testing
within_state_relative_tests = specttest(sp_fit_subj_bands,100,1,1);
within_state_absolute_tests = specttest(sp_fit_subj_bands,100,1,0);
within_state_relative_tests_significance = spectsignificance(within_state_relative_tests,0.05);
within_state_absolute_tests_significance = spectsignificance(within_state_absolute_tests,0.05);

HMM_permutation_tests.within_tests(1) = within_state_relative_tests;
HMM_permutation_tests.within_tests(2) = within_state_absolute_tests;

%% Absolute state testing
across_state_relative_tests = specttest(sp_fit_subj_bands,500,0,1);
across_state_absolute_tests = specttest(sp_fit_subj_bands,500,0,0);
across_state_relative_tests_significance = spectsignificance(across_state_relative_tests,0.05);
across_state_absolute_tests_significance = spectsignificance(across_state_absolute_tests,0.05);

HMM_permutation_tests.across_tests(1) = across_state_relative_tests;
HMM_permutation_tests.across_tests(2) = across_state_absolute_tests;

%%
save([mdl_nm{1} ,'_','HMM_permutation_tests'],'HMM_permutation_tests')