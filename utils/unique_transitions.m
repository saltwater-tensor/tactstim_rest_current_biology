function [unique_viterbi_transitions, transition_index] = unique_transitions(viterbi_trial_sequence)
    
    a = viterbi_trial_sequence;
    a_f = flip(a);
    % original
    chain_length = 1;
    for trm = 2:length(a)
        prv_letter = a(trm-1);
        current_letter = a(trm);
        if prv_letter ~= current_letter
            b(chain_length) = prv_letter;
            transition_index(chain_length) = trm-1;

            chain_length = chain_length + 1;
        end
    
        if trm == (length(a)-1)
            b(chain_length) = current_letter;
        end
    
    end
    
    %flipped
    chain_length = 1;
    for trm = 2:length(a_f)
        prv_letter = a_f(trm-1);
        current_letter = a_f(trm);
        if prv_letter ~= current_letter
            b_f(chain_length) = prv_letter;
            chain_length = chain_length + 1;
        end
    
        if trm == (length(a)-1)
            b_f(chain_length) = current_letter;
        end
    
    end
    b_f_corrected = flip(b_f);
    
    
    if length(b_f) == length(b)
        final_transition_list = b;
    elseif length(b_f) > length(b)
        final_transition_list = b_f_corrected;
    elseif length(b) > length(b_f)
        final_transition_list = b;
    end

    final_transition_list = num2str(final_transition_list);
    unique_viterbi_transitions = convertCharsToStrings(final_transition_list);
    
end
