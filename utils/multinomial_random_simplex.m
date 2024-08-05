function [transition_matrix] = multinomial_random_simplex(Dirichlet_concentration_params)

% Dir2d_alpha = hmm.hmm_subj_trial.Dir2d_alpha;
% Dir2d_alpha = hmm.hmm_subj_trial.Dir_alpha;

for alp = 1:size(Dirichlet_concentration_params,1)
    % Generate y_i (s) according to gamma distribution 
    y_i_state = gamrnd(Dirichlet_concentration_params(alp,:),1);
    normalizing_constant = sum(y_i_state);
    x_i_state = y_i_state/normalizing_constant;

    % if size(Dir2d_alpha,1) == 1 then it would be the initial state
    % probability matrix
    transition_matrix(alp,:) = x_i_state;
end
prob = round(sum(sum(transition_matrix,2)));
prob = int16(prob);
sz = size(Dirichlet_concentration_params,1);
sz = int16(sz);

assert(prob == sz)
