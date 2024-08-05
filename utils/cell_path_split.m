function [taskblock] = cell_path_split(x)
    taskblock = strsplit(x,'\');
    if length(taskblock) == 1
        taskblock = strsplit(x,'/');
    end

    if length(taskblock) > 6
        taskblock = taskblock{8};
    end
    
    if length(taskblock) < 4
        taskblock = taskblock{2};
        taskblock = strsplit(taskblock,'@raw');
        taskblock = taskblock{2};
    end
end