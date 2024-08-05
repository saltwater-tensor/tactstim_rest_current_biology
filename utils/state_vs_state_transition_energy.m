function energy_matrix = state_vs_state_transition_energy(cost_matrix)

% Logic for the analysis:
% If you plot or analyse the probability values for the posterior
% probability distirbution of states from the hmm output (stored in the
% variable gamma) you would notice that (by the very nature of the HMM
% algorithm) that at a given time point the probability of a single state
% will be close to 1 while that of others might be close to zero or very
% small. 
% This would be valid for every time point.
% Treat the probability distribution at time point 't' as the initial
% probability distribution then the distribution at 't+1' as the final
% probability distribution. 
% Now under an ideal scenario to make a guaranteed or 'optimal'switch from 
% state 1 at t to state 2 at 't+1' the transition probability matrix has to
% be setup such that p(1,2) = 1 and zero or close to zero for other entries
% But that is not the case
% For every trial we have a transition probability matrix where this is not
% the case.
% So in order to quantify this we can measure how difficult it is to
% transition from state 1 to state 2 under the 'non-optimal' conditions.
% This quantification can be done by solving the Schrodinger's bridge
% problem (SBP)
% Solve the Schrodinger Bridges problem to transition from probability
% distribution at 't' to some other probability distribution at 't+1' where
% the cost matrix is the transition probability matrix obtained from the
% hmm fit to the data which will output the a new 'optimal' transition
% matrix
% Then simply compute the KL divergence between this new optimal matrix and
% your actual cost matrix 

% Or you can directly compute the KL distance between your cost matrix and
% the 'optimal' transition matrix since you already know the solution
% Solving the SBP is a more rigorous way to do so and is also valid when
% you don't know the 'optimal' transition matrix

% Remember if we solve the optimal transport problem via Sinkhorn to get an
% optimal transport plan between the 'optimal' transition matrix and the
% real hmm fit transition matrix, that would answer the transport plan
% between these two matrices not really the cost issue.

    energy_matrix = NaN(6,6);
    simplexes = ones(6,1);
    simplexes = diag(simplexes);
    for i = 1:6
        for j = 1:6
            start_simplex = simplexes(i,:);
            end_simplex = simplexes(j,:);
             [mincost,P] = solveSBP(start_simplex',end_simplex',cost_matrix);
             energy_matrix(i,j) = mincost;
    
        end
    end

%     heatmap(energy_matrix,'Colormap',jet)   
%     heatmap(avrg_transition_r2-avrg_transition_r1,'Colormap',jet) 
%     heatmap(round(avrg_transition_r2-avrg_transition_r1),'Colormap',jet)

end