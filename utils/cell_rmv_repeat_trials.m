function [X] = cell_rmv_repeat_trials(x,y)

   Y = cell2mat(y);
   X = x;
   X (Y == 1) =[];
    
end