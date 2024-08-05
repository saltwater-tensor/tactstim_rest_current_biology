function [random_sequence] = markov_random_sequence(trial_length,dir_init_params, dir_transition_params)

% Generate a random initial state probability vector
initial_multinomial_distr = multinomial_random_simplex(dir_init_params);
pd_init = makedist('Multinomial','Probabilities',initial_multinomial_distr);

% Generate a random transition probability matrix
[transition_matrix] = multinomial_random_simplex(dir_transition_params);
for st = 1:size(transition_matrix,1)
    pd_transition = makedist('Multinomial','Probabilities',transition_matrix(st,:));
    transition_matrix_prob(st) = pd_transition; 
end

% Generate a state sequence
current_state = random(pd_init);
random_sequence(1,1) = current_state;
for st = 2:trial_length
    
    % Select the transition probability row according to the current state
    current_transition_pd = transition_matrix_prob(current_state);
    % Generate next state according to the transitions above
    pd_next = random(current_transition_pd);
    random_sequence(st,1) = pd_next;
    current_state = pd_next;
    

end

end