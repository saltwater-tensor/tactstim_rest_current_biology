function lftimes_mat_resp1 = rearrange_temporal_cell(LifeTimes_resp1)

for st = 1:6
    for trial = 1:length(LifeTimes_resp1)
        try
            lftimes_mat_resp1(trial,st) = {LifeTimes_resp1{trial, 1}{st, 1}};
        catch
            lftimes_mat_resp1(trial,st) = {NaN};
        end
    end    
end