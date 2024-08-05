function random_motif_matrix = null_sequence_motif_search(varargin)


if nargin < 3
    total_response1_trials = varargin{1};
    p = varargin{2};
end

if nargin > 2

total_response1_trials = varargin{1};
    p = varargin{2};
    markov = varargin{3};
    dir_init_params = varargin{4};
    dir_transition_params = varargin{5};

end
% h = waitbar(0,'Processing trials');

% 0.5 has been experimentally determined
total_time = 0.6*total_response1_trials; % in seconds

for tr = 1:total_response1_trials

    time_rem = 0.05*(total_response1_trials-tr);
    time_rem = round(time_rem/60) %minutes
%     waitbar(tr/total_response1_trials,h,...
%         ['Processing trial ' num2str(tr) ' time remaining ' num2str(time_rem) ' minutes'])

    if markov
        R1 = markov_random_sequence(1962,dir_init_params,dir_transition_params);
    else  
        R1 = randsample(6,1962,true);
    end
    
    [unique_R1] = unique_transitions(R1);
    unique_R1 = strrep(unique_R1,'  ','');

    for prm = 1:size(p,1)

        current_permutation = p(prm,:);
        current_permutation = num2str(current_permutation);
        current_permutation = convertCharsToStrings(strrep(current_permutation,'  ',''));
        k1 = strfind(unique_R1,current_permutation);
        random_motif_matrix(tr,prm) = length(k1);
    end


end

end
