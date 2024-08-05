function [] = write_brainstorm_output(brain_map, str, opdir, database_anat, destSurfFile)

kernelMat.ImageGridAmp = brain_map;
kernelMat.Comment      = str;
kernelMat.sRate      = 1;
kernelMat.ImageGridTime = 1:size(brain_map,2);
kernelMat.DataFile = [];
kernelMat.Time         = 1:size(brain_map,2);
kernelMat.SurfaceFile = destSurfFile;
kernelMat.HeadModelType = 'surface';
kernelMat.Function =  'lcmv';
kernelMat.GoodChannel =    [];
kernelMat.GridAtlas =    [];
kernelMat.GridLoc =  [];
kernelMat.GridOrient =    [];
% Output filename
OPTIONS.ResultFile = fullfile([database_anat '\data\' opdir], ...
    ['results_' str] );
% Save file
save(OPTIONS.ResultFile, '-struct', 'kernelMat');




end