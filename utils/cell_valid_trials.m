function [valid_trials] = cell_valid_trials(x,y)

    valid_wrapper = @(x) find(ismember(x , y));
    valid_trials = cellfun(valid_wrapper,x,'UniformOutput',false);

end