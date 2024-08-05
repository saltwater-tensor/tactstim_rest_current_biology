function [P_VAL] = perm_global_local_point_subs(original_state_gradients_one,...
    original_state_gradients_two)
Nperms = 5000;

% The logic here is to swap single points randomly across the two
% manifolds, calculate distance and then look at the difference in the
% number of significant long distance and short distance connections

%% -------- This portion can be refactored---------------------------------

% Calculate euclidean distance
[dist_mat_3D, dist_mat_2D] = manifold_euclidean_distance(original_state_gradients_one);

% Find larger and smaller significant connections for the first manifold
% [gsig_larger_one] = graph_significant_connections(dist_mat_3D,1,...
%     0.01,0.001);
[gsig_smaller_one] = graph_significant_connections(dist_mat_3D,0,...
    0.01,0.01);

% tmp_larger_one = squash(tril(gsig_larger_one));
% inds2_larger_one = find(tmp_larger_one>1e-10);
tmp_smaller_one = squash(tril(gsig_smaller_one,-1));
inds2_smaller_one = find(tmp_smaller_one>1e-10);

% Calculate euclidean distance
[dist_mat_3D, dist_mat_2D] = manifold_euclidean_distance(original_state_gradients_two);

% Find larger and smaller significant connections for the second manifold
% [gsig_larger_two] = graph_significant_connections(dist_mat_3D,1,...
%     0.01,0.05);
[gsig_smaller_two] = graph_significant_connections(dist_mat_3D,0,...
    0.01,0.01);

% tmp_larger_two = squash(tril(gsig_larger_two));
% inds2_larger_two = find(tmp_larger_two>1e-10);
tmp_smaller_two = squash(tril(gsig_smaller_two,-1));
inds2_smaller_two = find(tmp_smaller_two>1e-10);


% True differences
% true_diff_larger = length(inds2_larger_one) - length(inds2_larger_two);
true_diff_smaller = length(inds2_smaller_one) - length(inds2_smaller_two);

% ------------------------------------------------------------------------

%%
clc
p = size(dist_mat_3D,1);
for nperms = 1:Nperms
    nperms
    random_vector = randi(2,1,p*2)';
    random_vector = reshape(random_vector,[p,2]);

    % Generate shuffled manifold by randomly permuting points from the two
    % original manfiolds
    for P = 1:p
        if random_vector(P,1) == 1

            shuffled_state_gradients_one(P,:) = original_state_gradients_one(P,:);

        elseif random_vector(P,1) == 2

            shuffled_state_gradients_one(P,:) = original_state_gradients_two(P,:);

        end

        if random_vector(P,2) == 1

            shuffled_state_gradients_two(P,:) = original_state_gradients_one(P,:);

        elseif random_vector(P,2) == 2

            shuffled_state_gradients_two(P,:) = original_state_gradients_two(P,:);

        end

    end

    % Calculate euclidean distance
    [dist_mat_3D_shuffled, dist_mat_2D_shuffled] = manifold_euclidean_distance(shuffled_state_gradients_one);
    % Find larger and smaller significant connections for the first shuffled manifold
%     [shuffled_gsig_larger_one] = graph_significant_connections(dist_mat_3D_shuffled,1,...
%         0.01,0.05);
    [shuffled_gsig_smaller_one] = graph_significant_connections(dist_mat_3D_shuffled,0,...
        0.01,0.01);

%     tmp_larger_one = squash(tril(shuffled_gsig_larger_one));
%     inds2_larger_one = find(tmp_larger_one>1e-10);
    tmp_smaller_one = squash(tril(shuffled_gsig_smaller_one));
    inds2_smaller_one = find(tmp_smaller_one>1e-10);

    % Calculate euclidean distance
    [dist_mat_3D_shuffled, dist_mat_2D_shuffled] = manifold_euclidean_distance(shuffled_state_gradients_two);
    % Find larger and smaller significant connections for the first shuffled manifold
%     [shuffled_gsig_larger_two] = graph_significant_connections(dist_mat_3D_shuffled,1,...
%         0.01,0.05);
    [shuffled_gsig_smaller_two] = graph_significant_connections(dist_mat_3D_shuffled,0,...
        0.01,0.01);

%     tmp_larger_two = squash(tril(shuffled_gsig_larger_two));
%     inds2_larger_two = find(tmp_larger_two>1e-10);
    tmp_smaller_two = squash(tril(shuffled_gsig_smaller_two));
    inds2_smaller_two = find(tmp_smaller_two>1e-10);

    % Shuffled differences
%     shuffle_diff_larger(nperms) = length(inds2_larger_one) - length(inds2_larger_two);
    shuffle_diff_smaller(nperms) = length(inds2_smaller_one) - length(inds2_smaller_two);

    clc

end



% if abs(true_diff_larger) > 0
%     P_VAL.p_value_larger_dist = sum(abs(true_diff_larger) < abs(shuffle_diff_larger))/Nperms;
%     if  P_VAL.p_value_larger_dist == 0
%         P_VAL.p_value_larger_dist = eps;
%     end
%     P_VAL.true_larger_dist_diff = true_diff_larger;
% else
%     P_VAL.p_value_larger_dist = 1;
%     P_VAL.true_larger_dist_diff = 0;
% end

if abs(true_diff_smaller) > 0
    P_VAL.p_value_smaller_dist = sum(abs(true_diff_smaller) < abs(shuffle_diff_smaller))/Nperms;
    if P_VAL.p_value_smaller_dist == 0
        P_VAL.p_value_smaller_dist = eps;
    end
    P_VAL.true_smaller_dist_diff = true_diff_smaller;
else
    P_VAL.p_value_smaller_dist = 1;
    P_VAL.true_smaller_dist_diff = 0;
end



end