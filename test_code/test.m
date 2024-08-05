%% Testing code
% This test code was written to to find 
% a) Presence of a specific sequence or a substring within a large string
% b) To find unique transitions from a viterbi sequence i.e. remove self
% transitions 

%% part a

clear
s = RandStream('mlfg6331_64');
R1 = randsample(s,'123456',1967,true,[0.3 0.2 0.1 0.2 0.1 0.1]);
R2 = randsample(s,'123456',1967,true);
k1 = strfind(R1,'1264');
k2 = strfind(R2,'1264');

R1 = randsample('123456',1967,true,[0.3 0.2 0.1 0.2 0.1 0.1]);
R2 = randsample('123456',1967,true);
k1 = strfind(R1,'1264');
k2 = strfind(R2,'1264');

convertCharsToStrings(ans)
k1 = strfind(ans,'1264');

%% part b

a = randsample(2,200,true,[0.1 0.9]);
a_f = flip(a);

% original
chain_length = 1;
for trm = 2:length(a)
    prv_letter = a(trm-1);
    current_letter = a(trm);
    if prv_letter ~= current_letter
        b(chain_length) = prv_letter;
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

convertCharsToStrings(final_transition_list)