function [pathsplit] = cell_bhvr(x)

    original_file_location = {x.original_file_location}.';
    pathsplit = cellfun(@(x) cell_path_split(x), original_file_location,'UniformOutput',false);
    
end