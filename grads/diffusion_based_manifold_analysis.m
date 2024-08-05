
model_path = 'Q:\tactstim_resting\code\tactstim_rest_models\';
model_name = '_04_Feb_2022_17_04_27_HMM_model_preStim_rest_250Hz.mat';
dataset_path = 'Q:\tactstim_resting\code\tactstim_rest_data\';
mdl_nm = strsplit(model_name,'.');
load([model_path model_name]);
load(HMM_model.dataset_path)
load(HMM_model.time_path)

%%
T = cellfun(@transpose, T,'Uniformoutput',false);
Gamma_2 = padGamma(HMM_model.Gamma,T,HMM_model.options);
sz_wrapper = @(x) size(x,1);
dt_sz = cell2mat(cellfun(sz_wrapper, data, 'Uniformoutput', false));
assert(sum(dt_sz) == size(Gamma_2,1))
data = cell2mat(data);
gamma_mask = Gamma_2*100;
gamma_mask(gamma_mask < 80) = NaN;

%%
for state = 1:6

    state_cell_data = [];
    % Save a copy for finding the transitions
    gamma_mask2 = gamma_mask;
    gamma_mask2(~isnan(gamma_mask2)) = 1;
    gamma_mask2(isnan(gamma_mask2)) = 0;
    % Detect indices where the state visits start and end for this state
    [unique_gamma_transitions, transition_indices] = unique_transitions(gamma_mask2(:,state));
    transitions = reshape(transition_indices,[2,length(transition_indices)/2]);

    % gamma state is used for extracting data from the time series
    gamma_state = repmat(gamma_mask(:,state),1,42);
    gamma_state(~isnan(gamma_state)) = 1;

    % masking the data
    data_masked = data.*gamma_state;

    % Atlas based ROIs
    ROIs = [3,4,5,6,11,12,15,16,19,20,21,22,23,24,25,26,27,28,29,30,...
        33,34,35,36,37,38,41,42,45,46,47,48,51,52,53,54,55,56,57,58,59,60];
    % region specific data extraction
    sensory_motor_indices = [27,28,29,30];

    sensory_motor_data = data_masked(:,sensory_motor_indices);

    % Rearranding data as cells
    cell_count = 1;
    for trans = 1:size(transitions,2)

        init = transitions(1, trans) + 1;
        stop = transitions(2, trans);
        if (stop-init) > 10
            state_cell_data{cell_count,1} = data_masked(init:stop, :);
            cell_count = cell_count + 1;
        end



    end

    cvr(:,:,state) = cov(data_masked(:,:),'omitrows');
    %     [U,S,V] = svd(cvr(:,:,state));
%     cvr2(:,:,state) = cov(zscore(cell2mat(state_cell_data)));

    data_masked_nonan = data_masked(~isnan(data_masked));
    data_masked_nonan = reshape(data_masked_nonan,...
        length(data_masked_nonan)/size(data_masked,2),size(data_masked,2));

    cvr3(:,:,state) = cov(data_masked_nonan(:,:));

    % If bootstrap testing is required to be performed
%         [cvr_bootrstrap(:,:,state,:), cvr_bootrstrap_cell{state}] = ...
%             covar_bootstrap(data_masked_nonan,100);


end


%     gm = GradientMaps('n_components',4,'approach','dm','kernel','spearman','alignment','pa');
%     gm = gm.fit({cvr(:,:,1) cvr(:,:,2) cvr(:,:,3) cvr(:,:,4) cvr(:,:,5) cvr(:,:,6)});
%     scree_plot(gm.lambda{1});
%     hold on
%      [h,C] = gradient_in_euclidean(gm.gradients{1}(:,1:3));
%     
%      str_cell = {'motor1', 'motor2', 'motor3', 'motor4'};
%      text(gm.gradients{1}(sensory_motor_indices,1), gm.gradients{1}(sensory_motor_indices,2),...
%          gm.gradients{1}(sensory_motor_indices,3), str_cell)

%% Make plots for the gradient
database_anat = 'Q:\tactstim_resting\New folder\brainstorm_db\tactstim_resting\';
destSurfFile = '@default_subject/tess_cortex_pial_low.mat';
surface_file = destSurfFile;
surface_file_data = load([database_anat 'anat' filesep surface_file]);
atlas_name = 'Mindboggle';
atlas_number = 6;
% Atlas based ROIs
% These are the Scout numbers from the Mindboggle atlas. Directly from the
% atlas!
ROIs = [3,4,5,6,11,12,15,16,19,20,21,22,23,24,25,26,27,28,29,30,...
    33,34,35,36,37,38,41,42,45,46,47,48,51,52,53,54,55,56,57,58,59,60];

% Store the indices of the above ROIs vector and partition them into
% different regions so ROIs(sensory_motor_indices) would produce 41 42 45
% 46 and if you surface_file_data.Atlas(atname).ScoutsROIs(41 42 45
% 46).Label should produce post central L and R and then pre central L and
% R
% OR YOU COULD DO
% surface_file_data.Atlas(atname).ScoutsROIs(ROIs(sensory_motor_indices).Label


% ROI colors
colors = lines(7);
color_matrix = zeros(42,3);

visual_indices = [3,4,9,10,13,14];
visual_color = colors(7,:);
color_matrix(visual_indices,:) = repmat(visual_color,...
    [length(visual_indices),1]);

parietal_indices = [5,6,19,20,31,32,37,38,41,42];
pareital_color = colors(6,:);
color_matrix(parietal_indices,:) = repmat(pareital_color,...
    [length(parietal_indices),1]);

sensory_motor_indices = [27,28,29,30];
sensory_motor_color = colors(5,:);
color_matrix(sensory_motor_indices,:) = repmat(sensory_motor_color,...
    [length(sensory_motor_indices),1]);

temporal_indices = [17,18,39,40];
temporal_color = colors(4,:);
color_matrix(temporal_indices,:) = repmat(temporal_color,...
    [length(temporal_indices),1]);

mpfc_indices = [1,2,15,16];
mpfc_color = colors(3,:);
color_matrix(mpfc_indices,:) = repmat(mpfc_color,...
    [length(mpfc_indices),1]);

frontal_indices = [7,8,11,12,21,22,23,24,25,26,33,34,35,36];
frontal_color = colors(2,:);
color_matrix(frontal_indices,:) = repmat(frontal_color,...
    [length(frontal_indices),1]);

for atname = 1:1:length(surface_file_data.Atlas)
    if strcmpi(surface_file_data.Atlas(atname).Name,atlas_name)
        Scouts = surface_file_data.Atlas(atname).Scouts;
        for sc = 1:1:length(ROIs)
            V_label{sc} = Scouts(ROIs(sc)).Label;
        end
    end
end

%%
addpath(genpath('Q:\software\BrainSpace'));

% The sparsity threshold should be 90% this controls the connections in the
% graph
% diffusion alpha should be 0.5
% diffusion time should be 0 i.e. one step
% Check for these parameters directly in the Gradient maps script

% gm1 = GradientMaps('n_components',4,'approach','dm','kernel','spearman','alignment','pa');
% gm1 = gm1.fit({cvr(:,:,1) cvr(:,:,2) cvr(:,:,3) cvr(:,:,4) cvr(:,:,5) cvr(:,:,6)});

% gm2 = GradientMaps('n_components',4,'approach','dm','kernel','spearman','alignment','pa');
% gm2 = gm2.fit({cvr2(:,:,1) cvr2(:,:,2) cvr2(:,:,3) cvr2(:,:,4) cvr2(:,:,5) cvr2(:,:,6)});

% max you can go is 5
gm3 = GradientMaps('n_components',4,'approach','dm','kernel','spearman','alignment','pa');
gm3 = gm3.fit({cvr3(:,:,1) cvr3(:,:,2) cvr3(:,:,3) cvr3(:,:,4) cvr3(:,:,5) cvr3(:,:,6)});

% % If bootstrap testing is performed
% gm_bootstrap = GradientMaps('n_components',4,'approach',...
%     'dm','kernel','spearman','alignment','pa');
% for st = 1:6
%     gm_bootstrap_gradients{st} = gm_bootstrap.fit(cvr_bootrstrap_cell{st},...
%         'reference',gm3.gradients{st});
% end

scree_plot(gm3.lambda{1});
% hold on

%%
map_output_folder = 'Group_analysis\gradient_state_maps';
gm = gm3;
for st = 1:6
    brain_map = zeros(15002,3);
    %     [h,C] = gradient_in_euclidean(gm.gradients{st}(:,1:3));
    state_gradient = gm.gradients{st}(:,1:3);
    state_gradient = zscore(state_gradient,0,1);
    
    state_gradient_aligned = gm.aligned{st}(:,1:3);
    state_gradient_aligned = zscore(state_gradient_aligned,0,1);
    %  str_cell = {'motor1', 'motor2', 'motor3', 'motor4'};
    str_cell = V_label([sensory_motor_indices,visual_indices,parietal_indices]);
    
    
    %     h.subplot = subplot(3,2,st);
    
    % 3D plot
    h = figure(st);
    h.Position = [1379 88 380 271];
    ax = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
    scatter3(ax,state_gradient(:,1),state_gradient(:,2),...
        state_gradient(:,3),100,color_matrix,'filled','Marker','o');
    text(state_gradient([sensory_motor_indices,visual_indices,parietal_indices],1),...
        state_gradient([sensory_motor_indices,visual_indices,parietal_indices],2)...
        ,state_gradient([sensory_motor_indices,visual_indices,parietal_indices],3),str_cell)
    
    % All the gradients are aligned to the first gradient
%     h2 = figure(st*100);
%     h2.Position = [1379 88 380 271];
%     ax2 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
%     scatter3(ax2,state_gradient_aligned(:,1),state_gradient_aligned(:,2),...
%         state_gradient_aligned(:,3),100,color_matrix,'filled','Marker','o');
%     text(state_gradient_aligned(sensory_motor_indices,1),state_gradient_aligned(sensory_motor_indices,2)...
%         ,state_gradient_aligned(sensory_motor_indices,3),str_cell)
    
    % 2D plot
    %     h = figure(st);
    %     h.Position = [1379 88 380 271];
    %     ax = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
    %     %     ax = gca;
    %     scatter(ax,state_gradient(:,1),state_gradient(:,2)...
    %         ,100,color_matrix,'filled','Marker','o');
    %     text(state_gradient(sensory_motor_indices,1),state_gradient(sensory_motor_indices,2)...
    %         ,state_gradient(sensory_motor_indices,3),str_cell)
    
    
        for roi_nums = 1:42
    
            current_ROI = ROIs(roi_nums);
            current_gradient_val_1 = state_gradient(roi_nums,1);
            current_gradient_val_2 = state_gradient(roi_nums,2);
            current_gradient_val_3 = state_gradient(roi_nums,3);
            vertices_roi = {surface_file_data.Atlas(atlas_number).Scouts(current_ROI).Vertices}.';
            vertices_roi = cell2mat(cellfun(@transpose, vertices_roi,'UniformOutput',false));
    
            % Coloring the vertices
            brain_map(vertices_roi,1) = current_gradient_val_1...
                *ones(length(vertices_roi),1);
            brain_map(vertices_roi,2) = current_gradient_val_2...
                *ones(length(vertices_roi),1);
            brain_map(vertices_roi,3) = current_gradient_val_3...
                *ones(length(vertices_roi),1);
    
    
    
        end

%         str_gradient = ['Gradient_state_NEW' num2str(st) ];
%         write_brainstorm_output(brain_map, str_gradient,...
%             map_output_folder, database_anat, destSurfFile)

end

%%

output_dir = 'Q:\tactstim_resting\code\tactstim_rest_models\_04_Feb_2022_17_04_27_HMM_model_preStim_rest_250Hz_spectral_analysis\nnmf_profiles_and_group_projections_3\gradient_maps';
str = [];

figure_3d_save(str,output_dir)


%% Bootstrap based test stats
gm_full = gm3;
num_of_bootstraps = 100;
bootstrap_corr = [];
bootstrap_corr_align = [];

% Change the eigenvector or the gradient to test here
gradient_number = 3;

for stats_state = 1:6
    counter = 1;
    state_type = {};
    for cmpr_state = 1:6
        for bstr = 1:num_of_bootstraps
            vec_base = gm_bootstrap_gradients{1,stats_state}.gradients{1,bstr}(:,gradient_number);
            vec_base_align = gm_bootstrap_gradients{1,stats_state}.aligned{1,bstr}(:,gradient_number);
            % Compare to ground truth
            vec_cmpr =  gm_full.gradients{cmpr_state}(:,gradient_number);
            
            
            bootstrap_corr(stats_state,cmpr_state,bstr) =  corr(vec_base,vec_cmpr);
            bootstrap_corr_align(stats_state,cmpr_state,bstr) =  corr(vec_base_align,vec_cmpr);
            state_type(counter,1) = {['State' num2str(cmpr_state)]};
            counter = counter + 1;
        end
        
        
    end
    
%     figure(1)
%     boxplot(abs(squeeze(bootstrap_corr(stats_state,:,:))'))
%     
%     figure(2)
%     boxplot(squeeze(bootstrap_corr_align(stats_state,:,:))')
% 
%     figure(3)
%     boxplot(abs(squeeze(bootstrap_corr_align(stats_state,:,:))'))
    
    close all

    % anova testing
    anova_data = abs(squeeze(bootstrap_corr_align(stats_state,:,:))');
    anova_data = anova_data(:);
    [p,tbl,stats,terms] = anovan(anova_data,{state_type},'model','full','varnames',{'States'});
    figure(4)
    [results_all,mean_val,~,gnames_all] = multcompare(stats);

    % Now it will go to the next state
end

%% 
dist_mat_3D = [];
dist_mat_2D = [];

for st = 1:6

    euclidean_grads = gm3.gradients{st}(:,1:3);
    for nodes = 1:42
        
        % 3D distance i.e. using all 3 gradients
        current_point = euclidean_grads(nodes,:);
        current_point = repmat(current_point,[42,1]);
        remaining_points = euclidean_grads(1:42,:);
        difference = (current_point - remaining_points).*(current_point - remaining_points);
        euclidean_dist_3D = sqrt(sum(difference,2));
        dist_mat_3D(st,nodes,:) = euclidean_dist_3D;
        
        % 2D distance using top two principal gradients
        current_point = euclidean_grads(nodes,1:2);
        current_point = repmat(current_point,[42,1]);
        remaining_points = euclidean_grads(1:42,1:2);
        difference = (current_point - remaining_points).*(current_point - remaining_points);
        euclidean_dist_2D = sqrt(sum(difference,2));
        dist_mat_2D(st,nodes,:) = euclidean_dist_2D;
        
    end
end

%% Stats for distance

% Testing for significance of a single connection
K = 6;
ndim = 42;
ind2 = triu(true(ndim),1);
p = ndim * (ndim-1) / 2; 
X = zeros(K,p); 

% for k = 1:K
%     d = dist_mat(k,ind2)';
%     [p_val, p_val_corr] = permtest(d,10000);
%     X(k,:) = d(:);
% 
% end

% Cluster based testing
% We can use the same cluster based statistics test
P_VALUE_LARGER = 0.0000000001;
P_VALUE_SMALLER = 1e-12;
dist_mat_testing = dist_mat_3D; %dist_mat_3D

for k = 1:6
    
%     % Significantly larger distances
%     Graph = squeeze(dist_mat_testing(k,:,:));
% 
%     % Correct for average distance in the state
%     mean_dist = sum(sum(Graph))/(42*42);
%     Graph = Graph - mean_dist;
% 
%     tmp = squash(tril(Graph));
%     inds2 = find(tmp>1e-10);
%     data = tmp(inds2);
%     S2 = [];
%     S2.data = data;
%     S2.do_fischer_xform = false;
%     S2.do_plots = 1;
%     S2.pvalue_th = P_VALUE_LARGER/length(S2.data);
%     graph_ggm = teh_graph_gmm_fit(S2);
%     th = graph_ggm.normalised_th;
%     Graph = graph_ggm.data';
%     Graph(Graph<th) = NaN;
%     graphmat = zeros(ndim, ndim);
%     graphmat(inds2) = Graph;
%     graph_significant_larger(k,:,:) = graphmat;
        
    
    
    
    % Significantly shorter distances
    Graph = squeeze(dist_mat_testing(k,:,:));
    
    % Correct for average distance in the state
    mean_dist = sum(sum(Graph))/(42*42);
    Graph = Graph - (mean_dist);
    Graph = -1*Graph;
    
    tmp = squash(tril(Graph,-1));
%     inds2 = find(tmp<1e-10 & tmp~=0);
    inds2 = find(tmp>1e-10);
    data = tmp(inds2);
    S2 = [];
    S2.data = data;
    S2.do_fischer_xform = false;
    S2.do_plots = 1;
    
    % Bonferroni correction
    S2.pvalue_th = P_VALUE_LARGER/861;  %length(S2.data);
    
    % Non corrected
%     S2.pvalue_th = P_VALUE_SMALLER;
    
    
    graph_ggm = teh_graph_gmm_fit(S2);
    th = graph_ggm.normalised_th;
    Graph = graph_ggm.data';
    Graph(Graph<th) = NaN;
    graphmat = zeros(ndim, ndim);
    graphmat(inds2) = Graph;
    graph_significant_smaller(k,:,:) = graphmat;
    
    
end

%% Plotting results

% h = figure(100);
for k = 1:6
    state_gradient = gm3.gradients{k}(:,1:3);
    state_gradient = zscore(state_gradient,0,1);
    h = figure(k*100);
    h.Position = [690 101 535 438]; %  Position: [690 101 521 438] [1379 88 380 271]
    ax = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
    scatter3(ax,state_gradient(:,1),state_gradient(:,2),...
        state_gradient(:,3),20,color_matrix,'filled','Marker','o');
    xlabel('1st gradient')
    ylabel('2nd gradient')
    zlabel('3rd gradient')
    ax.CameraPosition = [8.3220 -30.2429 7.8337];
    ax.CameraPositionMode = 'auto';
    ax.CameraTarget = [0 0.3359 0];
    ax.CameraTargetMode = 'auto';
    ax.CameraUpVector = [0 0 1];
    ax.CameraUpVectorMode = 'auto';
    ax.CameraViewAngle = 7.7384;
    ax.CameraViewAngleMode = 'auto';
    ax.XLim = [-3 3];
    ax.YLim = [-3 3];
    ax.ZLim = [-3 3];
    hold on
%    
    
%     str_cell = V_label([sensory_motor_indices,visual_indices,parietal_indices]);
%         str_cell = V_label([sensory_motor_indices]);
%     text(state_gradient([sensory_motor_indices,visual_indices,parietal_indices],1),...
%         state_gradient([sensory_motor_indices,visual_indices,parietal_indices],2)...
%         ,state_gradient([sensory_motor_indices,visual_indices,parietal_indices],3),str_cell)
    
% Arrow based display 
% arrowmat1 = [-2,-2,-2,-2;2,2,2,2;2,1,0,-1]';
% arrowmat2 = [state_gradient([sensory_motor_indices],1),state_gradient([sensory_motor_indices],2),...
%     state_gradient([sensory_motor_indices],3)];
% daspect([1 1 1])
% arrow3(arrowmat1,arrowmat2,'r');
% 
% text([-2 -2 -2 -2],...
%         [2 2 2 2]...
%         ,[2 1 0 -1],str_cell,'Color','Black','FontSize',12)

% all labels
str_cell = V_label([sensory_motor_indices,visual_indices,frontal_indices,mpfc_indices,parietal_indices,...
    temporal_indices]);
text(state_gradient([sensory_motor_indices,visual_indices,frontal_indices,mpfc_indices,parietal_indices,...
    temporal_indices],1),...
        state_gradient([sensory_motor_indices,visual_indices,frontal_indices,mpfc_indices,parietal_indices,...
    temporal_indices],2)...
        ,state_gradient([sensory_motor_indices,visual_indices,frontal_indices,mpfc_indices,parietal_indices,...
    temporal_indices],3),str_cell,...
        'Color','Black','FontSize',12,'FontWeight','bold','HorizontalAlignment','center')


%     text(state_gradient([sensory_motor_indices,visual_indices,parietal_indices],1),...
%         state_gradient([sensory_motor_indices,visual_indices,parietal_indices],2)...
%         ,state_gradient([sensory_motor_indices,visual_indices,parietal_indices],3),str_cell)

    
    %     gsig_larger = squeeze(graph_significant_larger(k,:,:));
    gsig_smaller = squeeze(graph_significant_smaller(k,:,:));

    for node1 =  1:42  %  [sensory_motor_indices]
        for node2 = 1:42  %  [sensory_motor_indices]

            %             gsig_entry_larger = gsig_larger(node1,node2);
            %             if ~(isnan(gsig_entry_larger) || gsig_entry_larger == 0)
            %
            %                 point1_coords = state_gradient(node1,:);
            %                 point2_coords = state_gradient(node2,:);
            %                 X = [point1_coords(1);point2_coords(1)];
            %                 Y = [point1_coords(2);point2_coords(2)];
            %                 Z = [point1_coords(3);point2_coords(3)];
            %                 plot3(ax,X,Y,Z,'b','LineWidth',1)
            %
            %
            %
            %
            %             end

            gsig_entry_smaller = gsig_smaller(node1,node2);
            if ~(isnan(gsig_entry_smaller) || gsig_entry_smaller == 0)

                point1_coords = state_gradient(node1,:);
                point2_coords = state_gradient(node2,:);
                PC_sig = [point1_coords;point2_coords];
                X = [point1_coords(1);point2_coords(1)];
                Y = [point1_coords(2);point2_coords(2)];
                Z = [point1_coords(3);point2_coords(3)];
                plot3(ax,X,Y,Z,'k','LineWidth',2)

                scatter3(PC_sig(:,1),PC_sig(:,2),...
                    PC_sig(:,3),150,color_matrix([node1,node2],:),'filled','Marker','o','MarkerEdgeColor',[0 0 0]);

            end
        end
    end
    clear ax
%     scatter3(ax,state_gradient(:,1),state_gradient(:,2),...
%             state_gradient(:,3),20,color_matrix,'filled','Marker','o');
       
        % TO PRINT NAMES OF THE REGIONS
%         str_cell = V_label([sensory_motor_indices]);
% 
%         text(state_gradient([sensory_motor_indices],1),...
%         state_gradient([sensory_motor_indices],2)...
%         ,state_gradient([sensory_motor_indices],3),str_cell)

    % PRINCIPAL GRADIENT PLOT
    hold off
    
    h2 = figure(k*1000);
    h2.Position = [1379 88 180 180];
    ax2 = axes('Position',[0.1300 0.1100 0.7750 0.8150]);
    [sorted, idx_sorted] = sort(state_gradient(:,1));
    colors_sorted = color_matrix(idx_sorted,:);
    scatter(ax2,1:42,sorted,100,colors_sorted,'filled','Marker','o');
    V_label_2([1:42]) = {""};
    V_label_2([sensory_motor_indices]) = V_label([sensory_motor_indices]);
    sorted_text_cells = V_label_2(idx_sorted);
    
%     text(1:42,sorted,sorted_text_cells)
    hold on
    for node1 =  1:42  %  [sensory_motor_indices]
        for node2 = 1:42  %  [sensory_motor_indices]

            %             gsig_entry_larger = gsig_larger(node1,node2);
            %             if ~(isnan(gsig_entry_larger) || gsig_entry_larger == 0)
            %
            %                 point1_coords = state_gradient(node1,:);
            %                 point2_coords = state_gradient(node2,:);
            %                 X = [point1_coords(1);point2_coords(1)];
            %                 Y = [point1_coords(2);point2_coords(2)];
            %                 Z = [point1_coords(3);point2_coords(3)];
            %                 plot3(ax,X,Y,Z,'b','LineWidth',1)
            %
            %
            %
            %
            %             end

            gsig_entry_smaller = gsig_smaller(node1,node2);
            if ~(isnan(gsig_entry_smaller) || gsig_entry_smaller == 0)
                sorted_node1 = find(node1 == idx_sorted);
                sorted_node2 = find(node2 == idx_sorted);
                node1_val = sorted(sorted_node1);
                node2_val = sorted(sorted_node2);
                scatter(ax2,[sorted_node1,sorted_node2],[node1_val,node2_val],100,...
                    colors_sorted([sorted_node1,sorted_node2],:),...
                    'filled','Marker','o');
                
                
                
                
            end
        end
    end
    

end

%% Number of significant connections

% Although we see that the number of connections (long distance or short
% distance) are more or less in some states but when making the statement
% that one is a more global state and another a more local state we have to
% have some concrete evidence for this.


% Swapping pairs

msgbox("PERMUTATION TESTING SWAPPING DISTANCE PAIRS")

% pairs_larger_distance_matrix_p_val = NaN(6,6);
pairs_smaller_distance_matrix_p_val = NaN(6,6);
% pairs_larger_distance_effect_direction = NaN(6,6);
pairs_smaller_distance_effect_direction = NaN(6,6);

dist_mat_testing = dist_mat_3D; %dist_mat_3D
for k1 = 1:6
    distance_mat_one = squeeze(dist_mat_testing(k1,:,:));

    for k2 = setdiff(1:6,k1)
    
        distance_mat_two = squeeze(dist_mat_testing(k2,:,:));
        [P_VAL] =  perm_global_local_test(distance_mat_one, distance_mat_two);

%         pairs_larger_distance_matrix_p_val(k1,k2) = P_VAL.p_value_larger_dist;
%         pairs_larger_distance_effect_direction(k1,k2) = P_VAL.true_larger_dist_diff;
        pairs_smaller_distance_matrix_p_val(k1,k2) = P_VAL.p_value_smaller_dist;
        pairs_smaller_distance_effect_direction = P_VAL.true_smaller_dist_diff;
    
    end
end

% Multiple hypothesis correction

% Hypothesis 1
% Our hypothesis is that a state is global or local with respect to other
% remaining states. We don't wish to stay that there is one global state
% only. Hence the multiple correction has to be done state wise and only
% for that specific row and not across all comparisons performed.

% Hypothesis 2 : Global or local state across all comparisons made
% If we need to correct on a global level then only use the lower or upper
% triangle portion of the p_val matrix

% Hypothesis 1
for k = 1:6
    
%     current_p_val_vector = pairs_larger_distance_matrix_p_val(k,:);
    % mafdr takes into account nan values
%     current_p_val_vector_corr = mafdr(current_p_val_vector,'BHFDR',true);
%     pairs_larger_distance_matrix_p_val_hypothesis1_corr(k,:) = current_p_val_vector_corr;
    
    
    current_p_val_vector = pairs_smaller_distance_matrix_p_val(k,:);
    % mafdr takes into account nan values
    current_p_val_vector_corr = mafdr(current_p_val_vector,'BHFDR',true);
    pairs_smaller_distance_matrix_p_val_hypothesis1_corr(k,:) = current_p_val_vector_corr;

end

% Hypothesis 2
tmp = tril(ones(6,6),-1);
inds2 = find(tmp~=0);

% pairs_larger_distance_matrix_p_val_hypothesis2 = NaN(6,6);
% pairs_larger_distance_matrix_p_val_hypothesis2(inds2) = pairs_larger_distance_matrix_p_val(inds2);

pairs_smaller_distance_matrix_p_val_hypothesis2 = NaN(6,6);
pairs_smaller_distance_matrix_p_val_hypothesis2(inds2) = pairs_smaller_distance_matrix_p_val(inds2);

% pairs_larger_distance_matrix_p_val_hypothesis2_corr = NaN(6,6);
% pairs_larger_distance_matrix_p_val_hypothesis2_corr(inds2) = mafdr(pairs_larger_distance_matrix_p_val_hypothesis2(inds2),...
%     'BHFDR',true);

pairs_smaller_distance_matrix_p_val_hypothesis2_corr = NaN(6,6);
pairs_smaller_distance_matrix_p_val_hypothesis2_corr(inds2) = mafdr(pairs_smaller_distance_matrix_p_val_hypothesis2(inds2),...
    'BHFDR',true);


%% Swapping single points

msgbox("PERMUTATION TESTING SWAPPING SINGLE POINTS")

% singlepnt_larger_distance_matrix_p_val = NaN(6,6);
singlepnt_smaller_distance_matrix_p_val = NaN(6,6);
% singlepnt_larger_distance_effect_direction = NaN(6,6);
singlepnt_smaller_distance_effect_direction = NaN(6,6);

for k1 = 1:6
    state_gradients_one = gm3.gradients{k1}(:,1:3);

    for k2 = setdiff(1:6,k1)

        state_gradients_two = gm3.gradients{k2}(:,1:3);
        [P_VAL] = perm_global_local_point_subs(state_gradients_one,...
            state_gradients_two);
        
%         singlepnt_larger_distance_matrix_p_val(k1,k2) = P_VAL.p_value_larger_dist;
%         singlepnt_larger_distance_effect_direction(k1,k2) = P_VAL.true_larger_dist_diff;
        singlepnt_smaller_distance_matrix_p_val(k1,k2) = P_VAL.p_value_smaller_dist;
        singlepnt_smaller_distance_effect_direction = P_VAL.true_smaller_dist_diff;
        

    end
end

% Multiple hypothesis correction same hypothesis as above
% Hypothesis 1
for k = 1:6
    
%     current_p_val_vector = singlepnt_larger_distance_matrix_p_val(k,:);
    % mafdr takes into account nan values
%     current_p_val_vector_corr = mafdr(current_p_val_vector,'BHFDR',true);
%     singlepnt_larger_distance_matrix_p_val_hypothesis1_corr(k,:) = current_p_val_vector_corr;
    
    
    current_p_val_vector = singlepnt_smaller_distance_matrix_p_val(k,:);
    % mafdr takes into account nan values
    current_p_val_vector_corr = mafdr(current_p_val_vector,'BHFDR',true);
    singlepnt_smaller_distance_matrix_p_val_hypothesis1_corr(k,:) = current_p_val_vector_corr;
    
end

% Hypothesis 2
tmp = tril(ones(6,6),-1);
inds2 = find(tmp~=0);

% singlepnt_larger_distance_matrix_p_val_hypothesis2 = NaN(6,6);
% singlepnt_larger_distance_matrix_p_val_hypothesis2(inds2) = singlepnt_larger_distance_matrix_p_val(inds2);

singlepnt_smaller_distance_matrix_p_val_hypothesis2 = NaN(6,6);
singlepnt_smaller_distance_matrix_p_val_hypothesis2(inds2) = singlepnt_smaller_distance_matrix_p_val(inds2);

% singlepnt_larger_distance_matrix_p_val_hypothesis2_corr = NaN(6,6);
% singlepnt_larger_distance_matrix_p_val_hypothesis2_corr(inds2) = mafdr(singlepnt_larger_distance_matrix_p_val_hypothesis2(inds2),...
%     'BHFDR',true);

singlepnt_smaller_distance_matrix_p_val_hypothesis2_corr = NaN(6,6);
singlepnt_smaller_distance_matrix_p_val_hypothesis2_corr(inds2) = mafdr(singlepnt_smaller_distance_matrix_p_val_hypothesis2(inds2),...
    'BHFDR',true);


%% Coherence and distance analysis
% Hypothesis: Higher coherence means smaller distances 
% Result : No relationship, p-values only governed by outliers
spectral_profile_name = 'nnmf_profiles_and_group_projections_3.mat';
load([model_path mdl_nm{1} '_spectral_analysis' '\' spectral_profile_name]);

%%
tmp_mat = ones(42,42);
tmp_mat = tril(tmp_mat,-1);
inds2 = find(tmp_mat == 1);
for st_spec = 1:6
    % Get distance matrix for this state 
    dist_mat_state = squeeze(dist_mat_3D(st_spec,:,:));
    B = dist_mat_state(inds2);
    C = squeeze(graph_significant_larger(st_spec,:,:));
    inds_larger = find(C~=0 & ~isnan(C));
    D = squeeze(graph_significant_smaller(st_spec,:,:));
    inds_smaller = find(D~=0 & ~isnan(D));
    
    coh_mat_all_factors = fitmt_group_fact_4b.state(st_spec).coh; 
    
    for fct = 2:4
        coh_mat_factor_specific = tril(squeeze(coh_mat_all_factors(fct,:,:)),-1);
        
        % Coherence between all pairs
        A = coh_mat_factor_specific(inds2);
        mdl = fitlm(A,B);
        plot(mdl);
        mdl
        
        % Coherence betweeen pairs which have significantly shorter or
        % longer distances within each state
        % Result: No relationship!
        mdl_larger = fitlm(C(inds_larger),coh_mat_factor_specific(inds_larger));
        plot(mdl_larger);
        mdl_larger
        
        mdl_smaller = fitlm(D(inds_smaller),coh_mat_factor_specific(inds_smaller));
        plot(mdl_smaller);
        mdl_smaller
        
        
    end

end

%% 


%%

% Xe = [channel 1,   lag 1
%       channel 1,   lag 2
%       channel 1,   lag ...
%       channel 1,   lag N
%       channel 2,   lag 1
%       channel 2,   lag 2
%       channel 2,   lag ...
%       channel 2,   lag N
%       ....................
%       ....................
%       ....................
%       ....................
%       channel M,   lag 1
%       channel M,   lag 2
%       channel M,   lag ...
%       channel M,   lag N


% Embedding data
  options = HMM_model.options;
  options = format_options(options);
  [X_embed,T_embed,A_high_dim,A,B,C,e] = highdim_pca_custom(data,T,options.pca,...
        options.embeddedlags,options.standardise,...
        options.onpower,options.varimax,options.detrend,...
        options.filter,options.leakagecorr,options.Fs,options.As);

% Creating a selection vector based on the 42 brain regions
brain_regions = 1:42;

% Repeating the brain region vector for number of delays
brain_regions = repmat(brain_regions,[15,1]);
% Vectorising the matrix
brain_regions = brain_regions(:);

% Creating vector for delays and across brain regions
options = HMM_model.options;
E = options.embeddedlags;
E = repmat(E,[1,42])';

% Creating a brain mask
brain_regions_mask = NaN(630,1);

% Selecting specific brain regions to study and masking the brain regions
brain_regions_mask(find(brain_regions == sensory_motor_indices(3))) = 1;
brain_regions_mask = repmat(brain_regions_mask,[1,size(X_embed{1,1},1)])';



%% Experimental section one subject only
% Only using one sensory motor index for phase space analysis
% Remember we are trying to reduce the phase space dimensions so we will
% preserve time.

X_embed_sn = X_embed{1,1}.*brain_regions_mask;
X_emblo = X_embed_sn(~isnan(X_embed_sn));
X_emblo = reshape(X_emblo,[size(X_embed{1,1},1)],15);

% Masking and selecting based only on state
G = HMM_model.Gamma(1:size(X_embed{1,1},1),:)*100;
G(G < 75) = NaN;
G_state = G(:,4); % specify the state here
G_state(~isnan(G_state)) = 1;
X_emblo = X_emblo.*repmat(G_state,[1,15]);
X_emblo = double(X_emblo);

% Can use any sort of diffusion map, if it produces warning or error then
% do not sparsify the connectivity graph or use a kernel
gm3 = GradientMaps('n_components',3,'approach','dm','kernel','none','alignment','none');

% Manually create the variable for which to calculate the covariance
gm3 = gm3.fit({cov(X_emblo12')});
gradient_in_euclidean(gm3.gradients{1,1}(:,1:3))
tm_C = jet(size(gm3.gradients{1,1},1));
figure(2)
scatter3(gm3.gradients{1,1}(:,1),gm3.gradients{1,1}(:,2),gm3.gradients{1,1}(:,3),150,tm_C,'filled','Marker','o');
hold on
scatter3(gm3.gradients{1,1}(1,1),gm3.gradients{1,1}(1,2),gm3.gradients{1,1}(1,3),200,'Marker','+');
scatter3(gm3.gradients{1,1}(end,1),gm3.gradients{1,1}(end,2),gm3.gradients{1,1}(end,3),200,'Marker','*');
comet3(gm3.gradients{1,1}(:,1),gm3.gradients{1,1}(:,2),gm3.gradients{1,1}(:,3),0.99)

%Result:
% The delay embedding that we use creates highly auto-correlated time series.
% For a single region it is not correct to asses attractor dyanmics or manifold using this time series.
% It will produce smooth time organised manifolds for even random segments of
% data.



