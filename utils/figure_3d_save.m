function [] = figure_3d_save( str, output_dir)

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');

for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = get(FigHandle, 'Name');
    FigName = strrep(FigName,'/','_');
    FigName = strrep(FigName,': ','_');
    %   FigHandle.Renderer = 'painters';
    if sum(ismember(FigName,['3D_@']))
        
        if isempty(str)
            saveas(FigHandle, [output_dir '\' FigName '.png']);
        elseif ~isempty(str)
            FigName = [str num2str(iFig)];
            saveas(FigHandle, [output_dir '\' FigName '.png']);
        end
        
        
    end
end

end