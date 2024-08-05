function [dist_mat_3D, dist_mat_2D] = manifold_euclidean_distance(euclidean_grads)


for nodes = 1:size(euclidean_grads,1)

    % 3D distance i.e. using all 3 gradients
    current_point = euclidean_grads(nodes,:);
    current_point = repmat(current_point,[42,1]);
    remaining_points = euclidean_grads(1:42,:);
    difference = (current_point - remaining_points).*(current_point - remaining_points);
    euclidean_dist_3D = sqrt(sum(difference,2));
    dist_mat_3D(nodes,:) = euclidean_dist_3D;
    
    % 2D distance using top two principal gradients
    current_point = euclidean_grads(nodes,1:2);
    current_point = repmat(current_point,[42,1]);
    remaining_points = euclidean_grads(1:42,1:2);
    difference = (current_point - remaining_points).*(current_point - remaining_points);
    euclidean_dist_2D = sqrt(sum(difference,2));
    dist_mat_2D(nodes,:) = euclidean_dist_2D;

end


end