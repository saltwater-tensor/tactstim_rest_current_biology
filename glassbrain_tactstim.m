database ='Q:\tactstim_resting\New folder\brainstorm_db\tactstim_resting\';

ROIs_new = [3,4,5,6,11,12,15,16,19,20,21,22,23,24,25,26,27,28,29,30,...
    33,34,35,36,37,38,41,42,45,46,47,48,51,52,53,54,55,56,57,58,59,60];

surface_file = '@default_subject\tess_cortex_pial_low.mat';

surface_file_data = load([database 'anat' filesep surface_file]);
a= surface_file_data;
Atlas_name = surface_file_data.Atlas(6).Name;
Scouts = surface_file_data.Atlas(6).Scouts;
% load('D:\Abhinav_Sharma\RS_peri_MEGLFP\brainstorm_db\MEG_LFP_peri\anat\@default_subject\tess_cortex_pial_low.mat')

% a= surface_file_data;
% load('Q:\MEG_lfp_peri_analysis\hmm\Mindboggle_analysis\OFF\Whole_brain_stn_lfp_medication_OFF_06_Jan_2020_18_39_55_HMM_model_pca_NO_MAR_Motor_cortex_LFP_all_embed_lags\supply_profiles_4b_type4_filtered1to45\keyfunc_dot.mat')
% close all

% view_surface_matrix(Vertices, Faces, SurfAlpha=0, SurfColor=[], hFig=[], isFem=0, SurfaceFile=[])
% view_surface_matrix(surface_file_data.Vertices  , surface_file_data.Faces...
%     , 0, [], [], 0,surface_file)
%%
SF = 'Q:\tactstim_resting\New folder\brainstorm_db\tactstim_resting\anat\@default_subject\tess_cortex_pial_low.mat';

%%
r = [];
c = [];
adj_mat = [];
X = [];
Y = [];
Z = [];
connsx = [];
connsy = [];
connsz = [];
unique_label_points = [];
text_coordinate = [];
coordinate_names = {};

ROIs_new = re_ROIs_2;
% ROIs_new = re_ROIs_48_index_complete;
adj_mat = graph_schemaball;

[r,c] = find(adj_mat > 0 );

vector1 = linspace(0,1,100);
vector2 = floor(linspace(3,20,100));
vector3 = floor(linspace(1,10,100));
color_mat = hot(length(vector1));

if ~isempty(r)
    
for r1 = 1:1:length(r)
    
    
    point1 = r(r1);
    ROI_actual_1 = ROIs_new(point1);
    ROI_actual_1_seed_vertex = Scouts(ROI_actual_1).Seed;
    point2 = c(r1);
    ROI_actual_2 = ROIs_new(point2);
    ROI_actual_2_seed_vertex = Scouts(ROI_actual_2).Seed;
    
    X = [a.Vertices(ROI_actual_1_seed_vertex,1);...
        a.Vertices(ROI_actual_2_seed_vertex,1)];
    
    Y = [a.Vertices(ROI_actual_1_seed_vertex,2);...
        a.Vertices(ROI_actual_2_seed_vertex,2)];
    
    Z = [a.Vertices(ROI_actual_1_seed_vertex,3);...
        a.Vertices(ROI_actual_2_seed_vertex,3)];
    
    connsx{r1,1} = X;
    connsy{r1,1} = Y;
    connsz{r1,1} = Z;
    
    connection_strength = adj_mat(point1,point2);
    vec = connection_strength-vector1;
    vec_min = min(abs(vec));
    [p1,p2] = find(abs(vec) == vec_min);
    
    markerSize_strength(r1) = vector2(p2);
    lineWidth_strength(r1) = vector3(p2);
    linecolor_strength(r1,:) = color_mat(p2,:);
    
end
unique_label_points = [r;c];
unique_label_points = unique(unique_label_points);

[C,ia,ic] = unique(r);
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];

for u = 1:1:length(unique_label_points)
    
    point = unique_label_points(u);
    ROI_actual = ROIs_new(point);
    ROI_actual_seed_vertex = Scouts(ROI_actual).Seed;
    ROI_actual_seed_label = Scouts(ROI_actual).Label;
    
    text_coordinate(u,1) = a.Vertices(ROI_actual_seed_vertex,1);
    text_coordinate(u,2) = a.Vertices(ROI_actual_seed_vertex,2);
    text_coordinate(u,3) = a.Vertices(ROI_actual_seed_vertex,3);
    coordinate_names{u,1} = ROI_actual_seed_label;
    
end

end
%%
[hFig, iDS, iFig] = view_surface(SF,0.7,[],[]);
if ~isempty(r)

    ax = hFig.CurrentAxes;
    hold(ax)
    for rr = 1:1:length(connsx)

         X = connsx{rr};
         Y = connsy{rr};
         Z = connsz{rr};
         markerSize =  markerSize_strength(rr);
         lineWidth = lineWidth_strength(rr);
         color_val = linecolor_strength(rr,:);
%          plot3(X,Y,Z,'Marker','o','MarkerSize',markerSize,'MarkerFaceColor','r',...
%              'LineWidth',lineWidth,'Color',color_val,'MarkerEdgeColor','none')
         
         
        current_row_val = r(rr);
        [num_conns] =  find(ismember(C,current_row_val));
        number_of_connections = a_counts(num_conns);
        
        if number_of_connections > 1
            plot3(X,Y,Z,'Marker','o','MarkerSize',markerSize,'MarkerFaceColor','c',...
             'LineWidth',lineWidth,'Color',color_val,'MarkerEdgeColor','none')
        else
%             plot3(X,Y,Z,'Marker','o','MarkerSize',markerSize,'MarkerFaceColor','r',...
%              'LineWidth',lineWidth,'Color',color_val,'MarkerEdgeColor','none')
            continue
        end
        
    end
% clear X Y Z connsx connsy connsz

%     Display text
    for ut = 1:1:length(coordinate_names)

        text(text_coordinate(ut,1),text_coordinate(ut,2),text_coordinate(ut,3),...
            coordinate_names{ut},'Color','white','FontSize',10,'FontWeight','normal')


    end
%     waitfor(msgbox('Set transparency to 60%'));
%     figure_3d('FigureKeyPressedCallback',   gcf, keyEvent)
%     
%     figure_3d_save(['State ' num2str(k)],[' band ' num2str(band) ' '])
% waitfor(msgbox('Save all figures'));
% close all
end


%%
function [] = figure_3d_save( state,band )

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');

for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  FigName = strrep(FigName,'/','_');
  FigName = strrep(FigName,': ','_');
%   FigHandle.Renderer = 'painters';
  if sum(ismember(FigName,['3D_@']))
    FigName = [state band FigName num2str(iFig)];
    saveas(FigHandle, [pwd '\' FigName '.fig']);
    break
  end
end

end