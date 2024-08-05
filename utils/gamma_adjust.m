function gamma_update = gamma_adjust(gamma,T,d)

tuda_time = cellfun(@transpose,T,'UniformOutput',0);
tuda_time = cell2mat(tuda_time');
tuda_time_mask = [];
%Excluding shorter trials
max_time = max(tuda_time);
short_trials = find(tuda_time<max_time);

for tdtime = 1:length(tuda_time)
    time_length = tuda_time(tdtime);
    time_length = time_length - d;
    if ismember(short_trials,tdtime)

        time_vector = ones(time_length,1);

    else
        time_vector = NaN(time_length,1);
    end
    
    tuda_time_mask = [tuda_time_mask;time_vector];

end

tuda_gamma_updated = gamma.*tuda_time_mask;
NaN_mat = isnan(tuda_gamma_updated);
[r,c] = find(NaN_mat~=1);
gamma_update = gamma;
gamma_update(r,:) = [];
T_update = T;
T_update{1,1}(short_trials) = [];

if size(gamma_update,1) ~= (sum(T_update{1})-length(T_update{1})*d)
    error('gamma and time do not match!!!')
end

end