% - kernel (default: normalized angle)
%       - 'p','pearson'
%       - 'sm','spearman'
%       - 'g','gaussian'
%       - 'na','normalized angle'
%       - 'cs','cosine similarity'
%       - '','none'
%       - a function handle

%   - approach (default: diffusion embedding)
%       - 'dm','diffusion embedding'
%       - 'le','laplacian eigenmap'
%       - 'pca','principal component analysis'

gm = GradientMaps('n_components',4,'approach','dm','kernel','spearman','alignment','pa');
gm = gm.fit(cvr2);
scree_plot(gm.lambda{1});
gradient_in_euclidean(gm.gradients{1}(:,1:3));


[mappedA, mapping] = compute_mapping(state_cell_data{1,1}', 'DiffusionMaps', 3);
