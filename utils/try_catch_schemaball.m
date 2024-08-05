
% This script is for try catch error only
% This script should be never run individually

if ndim == 48
    graph2 = graph + graph';
    graph2 = graph2(re_ROIs_48_index_complete,re_ROIs_48_index_complete);
    graph2 = tril(graph2);
    graph_schemaball = normalize(graph2,1,'range',[0 1]);

    myLabel_2_schemaball = myLabel_2(:,1);

    %Numbering based labels
    for LBL = 1:1:24
        myLabel_2_schemaball{LBL,1} = [num2str(LBL) 'R'];
    end
    llb = 24;
    for LBL = 25:1:48
        myLabel_2_schemaball{LBL,1} = [num2str(llb) 'L'];
        llb = llb -1;
    end

    C = figure((k*100 + band));
    set(gcf,'NumberTitle','off') %don't show the figure number
    set(gcf,'Name',['State ' num2str(k) ' band ' num2str(band)]) %select the name you want
    schemaball(graph_schemaball,myLabel_2_schemaball,[],colr_mat_2_schm,schemaball_size,colr_mat_2_orig,...
        ['State ' num2str(k) ' band ' num2str(band)]);
    clear graph graphmat th grap_ggm

elseif ndim == 42
    graph2 = graph + graph';
    graph2 = nan([42,42]);
    graph2 = graph2(re_ROIs_42_index_complete,re_ROIs_42_index_complete);
    graph2 = tril(graph2);
    graph_schemaball = normalize(graph2,1,'range',[0 1]);
    myLabel_2_schemaball = myLabel_2(4:45,1);
    if size(colr_mat_2_orig,1)~=42
        colr_mat_2_orig = colr_mat_2_orig(4:45,:);
    end
    if size(colr_mat_2_schm,1)~=42
        colr_mat_2_schm = colr_mat_2_schm(4:45,:);
    end

    psd_scatter_color_rearranged = psd_scatter_color(re_ROIs_42_index_complete,:);
    %Numbering based labels
    %             for LBL = 1:1:21
    %                 myLabel_2_schemaball{LBL,1} = [num2str(LBL) 'R'];
    %             end
    %             llb = 24;
    %             for LBL = 22:1:42
    %                 myLabel_2_schemaball{LBL,1} = [num2str(llb) 'L'];
    %                 llb = llb -1;
    %             end
    C = figure((k*100 + band));
    colr_vals = jet(100);
    Cc = colr_vals([1,100],:,:);
    set(gcf,'NumberTitle','off') %don't show the figure number
    set(gcf,'Name',['State ' num2str(k) ' band ' num2str(band)]) %select the name you want
    schemaball(graph_schemaball,myLabel_2_schemaball,Cc,colr_mat_2_schm,schemaball_size,colr_mat_2_orig,...
        ['State ' num2str(k) ' band ' num2str(band)],psd_scatter_color_rearranged, output_path_spectral_figs);

    % To create glass brain figures for the coherence rings
    %                     glassbrain_tactstim
    clear graph graphmat th grap_ggm
    close all





end