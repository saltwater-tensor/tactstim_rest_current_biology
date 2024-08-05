function [graph_significant] = graph_significant_connections(dist_mat_testing,side_flag,...
    P_VALUE_LARGER,P_VALUE_SMALLER)

if side_flag
    % Significantly larger distances
    Graph = dist_mat_testing;
    tmp = squash(tril(Graph,-1));
    inds2 = find(tmp>1e-10);
    data = tmp(inds2);
    S2 = [];
    S2.data = data;
    S2.do_fischer_xform = false;
    S2.do_plots = 0;
    S2.pvalue_th = P_VALUE_LARGER/length(S2.data);
    graph_ggm = teh_graph_gmm_fit(S2);
    th = graph_ggm.normalised_th;
    Graph = graph_ggm.data';
    Graph(Graph<th) = NaN;
    graphmat = zeros(size(dist_mat_testing));
    graphmat(inds2) = Graph;
    graph_significant = graphmat;
    close all

else
    % Significantly shorter distances
    Graph = dist_mat_testing;

    % Correct for average distance in the state
%     mean_dist = sum(sum(Graph))/(42*42);
%     Graph = Graph - abs(mean_dist);

    % Correct for average distance in the state
    mean_dist = sum(sum(Graph))/(42*42);
    Graph = Graph - (mean_dist);
    Graph = -1*Graph;


    tmp = squash(tril(Graph,-1));

    %     inds2 = find(tmp<1e-10 & tmp~=0);
    inds2 = find(tmp>1e-10);

    data = tmp(inds2);
    S2 = [];
    S2.data = data;
    S2.do_fischer_xform = false;
    S2.do_plots = 0;
    S2.pvalue_th = P_VALUE_SMALLER/861;%length(S2.data);
%     S2.pvalue_th = P_VALUE_SMALLER;
    graph_ggm = teh_graph_gmm_fit(S2);
    th = graph_ggm.normalised_th;
    Graph = graph_ggm.data';
    Graph(Graph<th) = NaN;
    graphmat = zeros(size(dist_mat_testing));
    graphmat(inds2) = Graph;
    graph_significant = graphmat;
    close all

end



end