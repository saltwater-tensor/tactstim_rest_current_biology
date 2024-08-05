%% Anatomical parcellation

ROIs = [3,4,5,6,11,12,15,16,19,20,21,22,23,24,25,26,27,28,29,30,...
    33,34,35,36,37,38,41,42,45,46,47,48,51,52,53,54,55,56,57,58,59,60];
myLabel = cell(length(ROIs));

% for i =1:1:length(ROIs)
%   myLabel{i} = Atlas(atlas_number).Scouts(ROIs(i)).Label;
% end

total_parcellations = 12;

% Original ROI indices

frontal = [7,8,11,12,21,22,23,24,25,26,33,34,35,36];
frontal_right = [8,12,22,24,26,34,36];
frontal_left = flip([7,11,21,23,25,33,35]);

medial_PFC = [1,2,15,16];
medial_PFC_right = [2,16];
medial_PFC_left = flip([1,15]);

temporal = [17,18,39,40];
temporal_right = [18,40];
temporal_left = flip([17,39]);

sensory_motor = [27,28,29,30];
sensory_motor_right = [28,30];
sensory_motor_left = flip([27,29]);

parietal = [5,6,19,20,31,32,37,38,41,42];
parietal_right = [6,20,32,38,42];
parietal_left = flip([5,19,31,37,41]);


visual = [3,4,9,10,13,14];
visual_right = [4,10,14];
visual_left = flip([3,9,13]);

re_ROIs = [visual_right,parietal_right,sensory_motor_right,temporal_right,...
    medial_PFC_right,frontal_right,frontal_left,medial_PFC_left,...
    temporal_left,sensory_motor_left,parietal_left,visual_left];

regions = { [frontal],[medial_PFC], [temporal],...
   [sensory_motor], [parietal], [visual] }';

%% Properly load all data
% Load dataset (Should contain both time, data and sampling frequency
% information)
dataset_path = 'Q:\tactstim_resting\code\tactstim_rest_data\';
model_path = 'Q:\tactstim_resting\code\tactstim_rest_models\';
model_name = '_04_Feb_2022_17_04_27_HMM_model_preStim_rest_250Hz.mat';
single_trial_hmm = 'Q:\tactstim_resting\code\tactstim_rest_data\subject_trial_specific_hmms\';
mdl_nm = strsplit(model_name,'.');
base_model_name = mdl_nm{1, 1};
spectral_dir = dir([single_trial_hmm base_model_name '\' 'single_trial_projection_spectra\']);
cnt_resp1 = 1;
cnt_resp2 = 1;

for spdr = 1:length(spectral_dir)

    if ~spectral_dir(spdr).isdir
        spdr
        load(fullfile(spectral_dir(spdr).folder, spectral_dir(spdr).name))
        
        nm = strsplit(spectral_dir(spdr).name,'_');
        nm = nm{6};
        nm = strsplit(nm,'.');
        nm = nm{1};
        
        switch nm
            
            case 'Response1'
                
        	    response1{cnt_resp1,1}.state = trial_spectra_4b.state;
                cnt_resp1 = cnt_resp1 + 1;
                
            case 'Response2'
                
                response2{cnt_resp2,1}.state = trial_spectra_4b.state;
                cnt_resp2 = cnt_resp2 + 1;
                
        end
    end
end



%% COHERENCE ANALYSIS
%% Collecting all data for coherence comparison

mask = ones(size(response1{1, 1}.state(1).coh,2),size(response1{1, 1}.state(1).coh,2));
response1_coh = [];
response2_coh = [];

coh_resp1 = [];
coh_resp2 = [];

for factor_number = 1:4
    
    for region1 = 1:1:length(regions)
        for region2 = 1:1:length(regions)
            
            mask = ones(size(response1{1, 1}.state(1).coh,2),...
                size(response1{1, 1}.state(1).coh,2));
            
            %---------Select specific connections you want to run tests on
            rows = regions{region1,1};
            cols = regions{region2,1};
            mask(rows,cols) = 0;
            
            
            %Just for safety
            mask(cols,rows) = 0;
            mask = ~mask;
            
            mask = triu(mask,1);
            mask = double(mask);
            mask(find(mask == 0)) = NaN;
            
            % Response 1
            for state = 1:1:6
            
                for resp_1_num = 1:1:size(response1,1)
                    
                    data_resp1 = mask.* squeeze(response1...
                        {resp_1_num, 1}.state(state).coh(factor_number,:,:)); 
                    data_resp1 = data_resp1(:);
                    data_resp1 = data_resp1((~isnan(data_resp1)));
                    coh_resp1 = [data_resp1;coh_resp1];
                   
                    
                end
%                  mean_coh_resp1(state) = mean(coh_resp1);
%                  err_coh_resp1(state) = std(coh_resp1) / sqrt( length(coh_resp1));
                 
                 
                 coh_raw_region_resp1{state,factor_number,region1,region2} = (coh_resp1);
                 
                 coh_resp1 = [];
                 
            end


            for state = 1:1:6

                for resp_2_num = 1:1:size(response2,1)
                    
                    data_resp2 = mask.* squeeze(response2...
                        {resp_2_num, 1}.state(state).coh(factor_number,:,:));
                    data_resp2 = data_resp2(:);
                    data_resp2 = data_resp2((~isnan(data_resp2)));
                    coh_resp2 = [data_resp2;coh_resp2];
                    
                    
                end
                %                  mean_coh_resp1(state) = mean(coh_resp1);
                %                  err_coh_resp1(state) = std(coh_resp1) / sqrt( length(coh_resp1));
                
                
                coh_raw_region_resp2{state,factor_number,region1,region2} = (coh_resp2);
                
                coh_resp2 = [];
                
            end   
        end
    end
end

%% Testing
load('HMM_trial_spectra_proj.mat')

%%
% We are only comparing same state and same factor response1 vs response2 
% regions = { [frontal] (1),[medial_PFC](2), [temporal](3),...
%    [sensory_motor](4), [parietal](5), [visual](6) }';

  for state_number = 1:6
      state_number
  for factor_number = 1:4
      factor_number
      for region1 = 1:6
          for region2 = 1:6

              response1_coh_cmpr = coh_raw_region_resp1{state_number,...
                  factor_number,region1,region2};

              response2_coh_cmpr = coh_raw_region_resp2{state_number,...
                  factor_number,region1,region2};
              
              % response1 > response2
              resp1_grtr_resp2 = permutationTest(response1_coh_cmpr,response2_coh_cmpr,...
                  1000,'sidedness','larger','showprogress', 10);
              
              % response1 < response2
              resp1_lssr_resp2 = permutationTest(response1_coh_cmpr,response2_coh_cmpr,...
                  1000,'sidedness','smaller','showprogress', 10);
              
              p_val_resp1_grtr_resp2(state_number,factor_number,region1,region2) = resp1_grtr_resp2;
              
              p_val_resp1_lssr_resp2(state_number,factor_number,region1,region2) = resp1_lssr_resp2;
              
          end
      end
  end
  end

%% Visualise 

% Select a p-value matrix, correct for multiple comparisons and send to
% visualise
% regions = { [frontal] (1),[medial_PFC](2), [temporal](3),...
%    [sensory_motor](4), [parietal](5), [visual](6) }';
% for nnmf profiles 3 : 1 delta, 2 high beta/gamma, 3 beta , 4 alpha


for state_number = 1:6   
    for factor_number = 1:4
        p_val_vis_greater = squeeze(p_val_resp1_grtr_resp2(state_number,factor_number,:,:));
        p_val_vis_greater = mafdr(p_val_vis_greater(:),'BHFDR',true);
        p_val_vis_greater_corr(state_number,factor_number,:,:) = reshape(p_val_vis_greater,[6,6]);

        p_val_vis_lesser = squeeze(p_val_resp1_lssr_resp2(state_number,factor_number,:,:));
        p_val_vis_lesser = mafdr(p_val_vis_lesser(:),'BHFDR',true);
        p_val_vis_lesser_corr(state_number,factor_number,:,:) = reshape(p_val_vis_lesser,[6,6]);      
    end
end

% Visualise
p_thresh = 0.01;
factor_names = {'Delta','HighBetaGam', 'Beta', 'Alpha'};
for state_number = 1:6
    for factor_number = 1:4
        p_grtr = squeeze(p_val_vis_greater_corr(state_number,factor_number,:,:));
        p_grtr_vis = NaN(6,6);
        p_grtr_vis(p_grtr > p_thresh) = 0;
        p_grtr_vis(p_grtr < p_thresh) = 1;
        f = figure();
        f.Position = [1496 375 165 150];
        imagesc(p_grtr_vis);
        figurename = fullfile(pwd, ['grtr_State_' num2str(state_number) '_' factor_names{factor_number}]);
        savefig(f,figurename)
        close(f)
        p_lssr = squeeze(p_val_vis_lesser_corr(state_number,factor_number,:,:));
        p_lssr_vis = NaN(6,6);
        p_lssr_vis(p_lssr > p_thresh) = 0;
        p_lssr_vis(p_lssr < p_thresh) = 1;
        f = figure();
        f.Position = [1496 375 165 150];
        imagesc(p_lssr_vis);
        figurename = fullfile(pwd, ['lssr_State_' num2str(state_number) '_' factor_names{factor_number}]);
        savefig(f,figurename)
        close(f)

    end
end

%%
%%POWER SPECTRA ANALYSIS
%%
pow_resp1 = [];
pow_resp2 = [];

for factor_number = 1:4
    
    for region = 1:1:size(response1{1, 1}.state(1).psd,2)
    
            % Response 1
            for state = 1:1:6
            
                for resp_1_num = 1:1:size(response1,1)
                    
                    data_resp1 = squeeze(response1...
                        {resp_1_num, 1}.state(state).psd(factor_number,region,region)); 
                    pow_resp1 = [data_resp1;pow_resp1];
                   
                    
                end
                                 
                 pow_raw_region_resp1{state,factor_number,region} = (pow_resp1);
                 
                 pow_resp1 = [];
                 
            end


            for state = 1:1:6

                for resp_2_num = 1:1:size(response2,1)
                    
                    data_resp2 = squeeze(response2...
                        {resp_2_num, 1}.state(state).psd(factor_number,region,region));
                    data_resp2 = data_resp2(:);                    
                    pow_resp2 = [data_resp2;pow_resp2];
                    
                    
                end
                %                  mean_coh_resp1(state) = mean(coh_resp1);
                %                  err_coh_resp1(state) = std(coh_resp1) / sqrt( length(coh_resp1));
                
                
                pow_raw_region_resp2{state,factor_number,region} = (pow_resp2);
                
                pow_resp2 = [];
                
            end   
       
    end
end

%% Test the psd
 for state_number = 1:6
      state_number
  for factor_number = 1:4
      factor_number
      for region = 1:size(response1{1, 1}.state(1).coh,2)
          

              response1_psd_cmpr = pow_raw_region_resp1{state_number,...
                  factor_number,region};

              response2_psd_cmpr = pow_raw_region_resp2{state_number,...
                  factor_number,region};
              
              % response1 > response2
              resp1_grtr_resp2 = permutationTest(response1_psd_cmpr,response2_psd_cmpr,...
                  5000,'sidedness','larger','showprogress', 10);
              
              % response1 < response2
              resp1_lssr_resp2 = permutationTest(response1_psd_cmpr,response2_psd_cmpr,...
                  5000,'sidedness','smaller','showprogress', 10);
              
              pow_p_val_resp1_grtr_resp2(state_number,factor_number,region) = resp1_grtr_resp2;
              
              pow_p_val_resp1_lssr_resp2(state_number,factor_number,region) = resp1_lssr_resp2;
              
          
      end
  end
 end

%% vis
% regions = { [frontal] (1),[medial_PFC](2), [temporal](3),...
%    [sensory_motor](4), [parietal](5), [visual](6) }';
% for nnmf profiles 3 : 1 delta, 2 high beta/gamma, 3 beta , 4 alpha
factor_number = 4;
state_number = 4;
% region = 4;
for state_number = 1:6
    for factor_number = 1:4


        p_val_vis_greater = squeeze(pow_p_val_resp1_grtr_resp2(state_number,factor_number,:));
        p_val_vis_greater_corr(state_number,factor_number,:) = mafdr(p_val_vis_greater,'BHFDR',true);
%         
%         h = figure();
%         h.Position = [1708 300 150 400];
%         heatmap(p_val_vis_greater_corr,'Colormap',jet)

        p_val_vis_lesser = squeeze(pow_p_val_resp1_lssr_resp2(state_number,factor_number,:));
        p_val_vis_lesser_corr(state_number,factor_number,:) = mafdr(p_val_vis_lesser,'BHFDR',true);
        
%         f = figure();
%         f.Position = [1708 300 150 400];
%         heatmap(p_val_vis_lesser_corr,'Colormap',jet)

    end
end

% [p_state, p_factor, p_region] = ind2sub(size(pow_p_val_resp1_lssr_resp2),find(pow_p_val_resp1_lssr_resp2<0.01));
% 
% [p_state, p_factor, p_region] = ind2sub(size(pow_p_val_resp1_grtr_resp2),find(pow_p_val_resp1_grtr_resp2<0.01));

%% Produce brain maps

% A brain map should be a matrix 15002 X 2 where 2 has higher power regions
% and 1 has lower power regions
load('Q:\tactstim_resting\New folder\brainstorm_db\tactstim_resting\anat\@default_subject\tess_cortex_pial_low')
destSurfFile = '@default_subject/tess_cortex_pial_low.mat';
HeadModelType = 'surface';
HeadModelFunction = 'lcmv'; 
atlas_name = 'Mindboggle';
atlas_number = 6;
% ROIs from this atlas that are part of our dataset
ROIs = [3,4,5,6,11,12,15,16,19,20,21,22,23,24,25,26,27,28,29,30,...
    33,34,35,36,37,38,41,42,45,46,47,48,51,52,53,54,55,56,57,58,59,60]; %42 regions

database_anat = 'Q:\tactstim_resting\New folder\brainstorm_db\tactstim_resting\';
surface_file = destSurfFile;
surface_file_data = load([database_anat 'anat' filesep surface_file]);
brain_map = zeros(15002,2);
all_brain_maps = zeros(15002,6,4);

map_output_folder = 'Group_analysis\PSD_state_factor_maps';
p_thresh = 0.05;
p_thresh_string = '05';
factor_names = {'Delta','HighBetaGam', 'Beta', 'Alpha'};

% P values to use
pow_p_val_lssr = pow_p_val_resp1_lssr_resp2;
pow_p_val_grtr = pow_p_val_resp1_grtr_resp2;

% pow_p_val_lssr = p_val_vis_lesser_corr;
% pow_p_val_grtr = p_val_vis_greater_corr;


for state_number = 1:6
    for factor_number = 1:4

       brain_map = zeros(15002,3);
       % LESSER
       current_p_vals_lssr = pow_p_val_lssr(state_number, factor_number,:);
       roi_indices_lssr = ind2sub(size(current_p_vals_lssr),find(current_p_vals_lssr<p_thresh));
       ROIs_lssr = ROIs(roi_indices_lssr);
       vertices_lssr = {surface_file_data.Atlas(atlas_number).Scouts(ROIs_lssr).Vertices}.';
       vertices_lssr = cell2mat(cellfun(@transpose, vertices_lssr,'UniformOutput',false));
       
       % putting 10 as value in the lesser dimension 
       % Doing this as our emphasis is more on response2 > response1
       % or response1 < response2
       brain_map(vertices_lssr,1) = 10*ones(length(vertices_lssr),1);
        
       % GREATER
       current_p_vals_grtr = pow_p_val_grtr(state_number, factor_number,:);
       roi_indices_grtr = ind2sub(size(current_p_vals_grtr),find(current_p_vals_grtr<p_thresh));
       ROIs_grtr = ROIs(roi_indices_grtr);
       vertices_grtr = {surface_file_data.Atlas(atlas_number).Scouts(ROIs_grtr).Vertices}.';
       vertices_grtr = cell2mat(cellfun(@transpose, vertices_grtr,'UniformOutput',false));
       
       % putting 5 as value in the higher dimension
       % This is response2 < response1
       % OR response1 > response2
       brain_map(vertices_grtr,2) = 5*ones(length(vertices_grtr),1);
       
       % creating a combined map
       brain_map([vertices_lssr;vertices_grtr],3) = [10*ones(length(vertices_lssr),1);...
           5*ones(length(vertices_grtr),1)];
       
       all_brain_maps(:,state_number, factor_number) = brain_map(:,3);
        
       % Save and create brainstorm output
       kernelMat.ImageGridAmp = brain_map;
       kernelMat.Comment      = ['State_' num2str(state_number) '_factor_' factor_names{factor_number} '_p_' p_thresh_string ];
       kernelMat.sRate      = 1;
       kernelMat.ImageGridTime = 1:size(brain_map,2);
       kernelMat.DataFile = [];
       kernelMat.Time         = 1:size(brain_map,2);
       kernelMat.SurfaceFile = destSurfFile;
       kernelMat.HeadModelType = HeadModelType;
       kernelMat.Function =  HeadModelFunction;
       kernelMat.GoodChannel =    [];
       kernelMat.GridAtlas =    [];
       kernelMat.GridLoc =  [];
       kernelMat.GridOrient =    [];
       % Output filename
       OPTIONS.ResultFile = fullfile([database_anat '\data\' map_output_folder], ...
           ['results_' 'State_' num2str(state_number) '_factor_' factor_names{factor_number} '_p_' p_thresh_string ] );
       % Save file
       save(OPTIONS.ResultFile, '-struct', 'kernelMat');
    
    end
end

% figure_3d_save(str,output_dir)

%% Use all_brain_maps to make a combined figure
% The time axis should be states and should be made for a given factor

alpha_state_rhythm = all_brain_maps(:,[4,6,5,1],4);
write_brainstorm_output(alpha_state_rhythm, 'Alpha_state_4651',...
    map_output_folder, database_anat, destSurfFile)

beta_state_rhythm = all_brain_maps(:,[4,6,5,1],3);
write_brainstorm_output(beta_state_rhythm, 'Beta_state_4651',...
    map_output_folder, database_anat, destSurfFile)

highbeta_state_rhythm = all_brain_maps(:,[4,6,5,1],2);
write_brainstorm_output(highbeta_state_rhythm, 'High_Beta_state_4651',...
    map_output_folder, database_anat, destSurfFile)

%% Within state 
dataset_path = 'Q:\tactstim_resting\code\tactstim_rest_data\';
model_path = 'Q:\tactstim_resting\code\tactstim_rest_models\';
model_name = '_04_Feb_2022_17_04_27_HMM_model_preStim_rest_250Hz.mat';
single_trial_hmm = 'Q:\tactstim_resting\code\tactstim_rest_data\subject_trial_specific_hmms\';
mdl_nm = strsplit(model_name,'.');
base_model_name = mdl_nm{1, 1};
spectral_dir = dir([single_trial_hmm base_model_name '\' 'single_trial_projection_spectra\']);

N = length(spectral_dir)-2;
% fitmt_subj_fact_1d = cell(N,1);
fitmt_subj_fact_4d = cell(N,1);
n = 1;
% factor_number = 2;
K = 6;

% for factor_number = 1:4
    for spdr = 1:length(spectral_dir)
        spdr
        if ~spectral_dir(spdr).isdir

            load(fullfile(spectral_dir(spdr).folder, spectral_dir(spdr).name))

            fitmt_subj_fact_1d{n} = struct();
            fitmt_subj_fact_1d{n}.state = struct();
            for k = 1:K % we don't care about the second component
                fitmt_subj_fact_4d{n}.state(k).psd = trial_spectra_4b.state(k).psd(:,:,:);

            end
            n = n + 1;
        end
    end
%     tests{factor_number} = specttest(fitmt_subj_fact_1d,5000,1,1);
tests = specttest(fitmt_subj_fact_4d,5000,1,1);
% tests_across_states = specttest(fitmt_subj_fact_4d,5000,0,1);
%     significant{factor_number} = spectsignificance(tests{factor_number},0.05);
% end
%%
% load('HMM_trial_spectra_proj_PSD_permutation_test_within_state_results.mat')

% load('HMM_trial_spectra_proj_PSD_permutation_test_across_state_results.mat')

load('HMM_trial_spectra_proj_PSD_permutation_test_within_state_results_2.mat')

p_thresh = 0.01;
p_thresh_string = '01';
factor_names = {'Delta','HighBetaGam', 'Beta', 'Alpha'};
load_anatomy_files
map_output_folder = 'Group_analysis\PSD_state_factor_maps';

all_brain_maps = zeros(15002,6,4);
 
for factor_number = 1:4
    
%     significant = spectsignificance(tests{factor_number},p_thresh);
%     significant = spectsignificance(tests_across_states,p_thresh);

    significant = spectsignificance(tests,p_thresh);

    hgh_corr = significant.higher_corr;
    lwr_corr = significant.lower_corr;
    for k = 1:6
    % for within state testing or results with a single factor
%         state_hghr = squeeze(hgh_corr.state(k).psd);
%         state_lwr = squeeze(lwr_corr.state(k).psd);
        
    % for across state testing or within state results with all factor
    % numbers
        state_hghr = squeeze(hgh_corr.state(k).psd(factor_number,:,:));
        state_lwr = squeeze(lwr_corr.state(k).psd(factor_number,:,:));
        
        
        
        brain_map = zeros(15002,3);
        r1_hghr = [];
        r2_hghr = [];
        r1_lwr = [];
        r2_lwr = [];
        
        %         state_hghr_ROI_numbers =
        [r1_hghr, r2_hghr] = ind2sub(size(state_hghr),find(state_hghr == 1));
        
        
        [r1_lwr, r2_lwr] = ind2sub(size(state_lwr),find(state_lwr == 1));
        
        ROIs_lssr = ROIs(r1_lwr);
        ROIs_grtr = ROIs(r1_hghr);
        
        vertices_lssr = {surface_file_data.Atlas(atlas_number).Scouts(ROIs_lssr).Vertices}.';
        vertices_lssr = cell2mat(cellfun(@transpose, vertices_lssr,'UniformOutput',false));
        vertices_grtr = {surface_file_data.Atlas(atlas_number).Scouts(ROIs_grtr).Vertices}.';
        vertices_grtr = cell2mat(cellfun(@transpose, vertices_grtr,'UniformOutput',false));
        
        brain_map(vertices_lssr,1) =  10*ones(length(vertices_lssr),1);
        brain_map(vertices_grtr,2) =  20*ones(length(vertices_grtr),1);
        brain_map([vertices_lssr;vertices_grtr],3) = [10*ones(length(vertices_lssr),1);...
           20*ones(length(vertices_grtr),1)];
        
        
        all_brain_maps(:, k, factor_number) = brain_map(:,3);
        
        
%         write_brainstorm_output(brain_map, ...
%             ['within_state_' num2str(k) '_factor_' factor_names{factor_number} '_' p_thresh_string '__two'], ...
%             map_output_folder, database_anat, destSurfFile)
        
%         write_brainstorm_output(brain_map, ...
%             ['across_state_' num2str(k) '_factor_' factor_names{factor_number} '_' p_thresh_string], ...
%             map_output_folder, database_anat, destSurfFile)

    end
end

%% Use all_brain_maps to make a combined figure
% The time axis should be states and should be made for a given factor

alpha_state_rhythm = all_brain_maps(:,[4,6,5,1],4);
write_brainstorm_output(alpha_state_rhythm, 'Alpha_state_4651_within_state_p_01',...
    map_output_folder, database_anat, destSurfFile)

beta_state_rhythm = all_brain_maps(:,[4,6,5,1],3);
write_brainstorm_output(beta_state_rhythm, 'Beta_state_4651_within_state_p_01',...
    map_output_folder, database_anat, destSurfFile)

highbeta_state_rhythm = all_brain_maps(:,[4,6,5,1],2);
write_brainstorm_output(highbeta_state_rhythm, 'High_Beta_state_4651_within_state_p_01',...
    map_output_folder, database_anat, destSurfFile)
