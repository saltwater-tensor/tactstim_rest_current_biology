function scatter_color = ring_visualisation_psd_control_2(psd_value,tests,k,factor_number)

p_thresh = 0.01;
p_thresh_string = '01';
factor_names = {'Delta','HighBetaGam', 'Beta', 'Alpha'};

significant = spectsignificance(tests,p_thresh);
hgh_corr = significant.higher_corr;
lwr_corr = significant.lower_corr;

state_hghr = squeeze(hgh_corr.state(k).psd(factor_number,:,:));
state_lwr = squeeze(lwr_corr.state(k).psd(factor_number,:,:));

scatter_color = zeros(size(state_lwr,1),3);
r1_hghr = [];
r2_hghr = [];
r1_lwr = [];
r2_lwr = [];

[r1_hghr, r2_hghr] = ind2sub(size(state_hghr),find(state_hghr == 1));


[r1_lwr, r2_lwr] = ind2sub(size(state_lwr),find(state_lwr == 1));

nodes_lssr = (r1_lwr);
nodes_grtr = (r1_hghr);

% Significant connections

% PSD scaled colors
% psd_grtr = psd_value(nodes_grtr);
% max_grtr = max(psd_grtr);
% min_grtr = min(psd_grtr);
% grtr_linspace = linspace(min_grtr,max_grtr,100);
% grtr_clrs = autumn(100);
% 
% % This part can be more efficient
% for g = 1:length(nodes_grtr)
%     idx = nearest_osl(grtr_linspace, psd_value(nodes_grtr(g)));
%     scatter_color(nodes_grtr(g),:) = grtr_clrs(idx,:);
% end
% 
% psd_lssr = psd_value(nodes_lssr);
% max_lssr = max(psd_lssr);
% min_lssr = min(psd_lssr);
% lssr_linspace = linspace(min_lssr,max_lssr,100);
% lssr_clrs = winter(100);
% 
% for l = 1:length(nodes_lssr)
%     idx = nearest_osl(grtr_linspace, psd_value(nodes_lssr(l)));
%     scatter_color(nodes_lssr(l),:) = lssr_clrs(idx,:);
% end


% Binary color scheme blue-lower red-higher black-non significant
scatter_color(nodes_grtr,:) = repmat([1 0 0],length(nodes_grtr),1);
scatter_color(nodes_lssr,:) = repmat([0 0 1],length(nodes_lssr),1);


end