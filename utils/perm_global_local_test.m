function [P_VAL] = perm_global_local_test(distance_mat_one, distance_mat_two)

Nperms = 1000;
P_VALUE_SMALLER = 1e-9;
P_VALUE_LARGER = 0.01;
% graph_significant_connections(dist_mat_testing,side_flag,...
%     P_VALUE_LARGER,P_VALUE_SMALLER)
% Original matrices and results
% [gsig_larger_one] = graph_significant_connections(distance_mat_one,1,0.01,0.05);
% [gsig_larger_two] = graph_significant_connections(distance_mat_two,1,0.01,0.05);
[gsig_smaller_one] = graph_significant_connections(distance_mat_one,0,P_VALUE_LARGER,P_VALUE_SMALLER);
[gsig_smaller_two] = graph_significant_connections(distance_mat_two,0,P_VALUE_LARGER,P_VALUE_SMALLER);

% tmp_larger_one = squash(tril(gsig_larger_one));
% inds2_larger_one = find(tmp_larger_one>1e-10);

% tmp_larger_two = squash(tril(gsig_larger_two));
% inds2_larger_two = find(tmp_larger_two>1e-10);

% true_diff_larger = length(inds2_larger_one) - length(inds2_larger_two);

tmp_smaller_one = squash(tril(gsig_smaller_one));
inds2_smaller_one = find(tmp_smaller_one>1e-10);

tmp_smaller_two = squash(tril(gsig_smaller_two));
inds2_smaller_two = find(tmp_smaller_two>1e-10);

true_diff_smaller = length(inds2_smaller_one) - length(inds2_smaller_two);

% 
ind2 = tril(true(size(distance_mat_one,1)),-1);
ind2 = find(ind2==1);
p = size(distance_mat_one,1) * (size(distance_mat_one,1)-1) / 2; 

% Original state distance vectors
state_one = distance_mat_one(ind2);
state_two = distance_mat_two(ind2);
original_distance_vector = [state_one,state_two];

for nperms = 1:Nperms
    nperms    
    random_vector = randi(2,1,p*2)';
    random_vector = reshape(random_vector,[p,2]);

    for P = 1:p
        state_shuffled_one(P,1) = original_distance_vector(P,random_vector(P,1));
        state_shuffled_two(P,1) = original_distance_vector(P,random_vector(P,2));
        
    end

    shuffled_distance_one = zeros(size(distance_mat_one));
    shuffled_distance_one(ind2) = state_shuffled_one;
    shuffled_distance_one = shuffled_distance_one + shuffled_distance_one';
    shuffled_distance_two = zeros(size(distance_mat_two));
    shuffled_distance_two(ind2) = state_shuffled_two;
    shuffled_distance_two = shuffled_distance_two + shuffled_distance_two';

%     gsig_larger_one = [];
%     gsig_larger_two = [];
    gsig_smaller_one = [];
    gsig_smaller_two = [];

%     if abs(true_diff_larger) > 0
% 
%         [gsig_larger_one] = graph_significant_connections(shuffled_distance_one,1,P_VALUE_LARGER,P_VALUE_SMALLER);
%         [gsig_larger_two] = graph_significant_connections(shuffled_distance_two,1,P_VALUE_LARGER,P_VALUE_SMALLER);
%         tmp_larger_one = squash(tril(gsig_larger_one));
%         inds2_larger_one = find(tmp_larger_one>1e-10);
% 
%         tmp_larger_two = squash(tril(gsig_larger_two));
%         inds2_larger_two = find(tmp_larger_two>1e-10);
% 
%         shuffle_diff_larger(nperms) = length(inds2_larger_one) - length(inds2_larger_two);
% 
%     end
    

    if abs(true_diff_smaller) > 0

        [gsig_smaller_one] = graph_significant_connections(shuffled_distance_one,0,P_VALUE_LARGER,P_VALUE_SMALLER);
        [gsig_smaller_two] = graph_significant_connections(shuffled_distance_two,0,P_VALUE_LARGER,P_VALUE_SMALLER);

        tmp_smaller_one = squash(tril(gsig_smaller_one));
        inds2_smaller_one = find(tmp_smaller_one>1e-10);

        tmp_smaller_two = squash(tril(gsig_smaller_two));
        inds2_smaller_two = find(tmp_smaller_two>1e-10);

        shuffle_diff_smaller(nperms) = length(inds2_smaller_one) - length(inds2_smaller_two);

    end

    


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