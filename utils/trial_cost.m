function [mincost_vec] = trial_cost(gamma,cost_matrix)
    

    for idx = 1:size(gamma,1)-1
        start_simplex = gamma(idx,:);
        end_simplex = gamma(idx+1,:);
        [mincost,P] = solveSBP(start_simplex',end_simplex',cost_matrix);
        mincost_vec(idx) = mincost;
%         heatmap(P)
    end
    
    
    
end