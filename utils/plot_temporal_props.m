function [] = plot_temporal_props(m_var)

h = figure();
h.Position = [1526 295 215 200];
pl = subplot(1,1,1);
errorbar(m_var([1,2],1),m_var([1,2],2),'s','MarkerFaceColor','black','MarkerSize',10,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.79,'CapSize',8)
pl.XLim = [0.5 2.5];
pl.LineWidth = 1.5;

% Y values for optimal display
% max_val = max(m_var([1],1)+ m_var([1],2), m_var([2],1)+ m_var([2],2));
% max_val = (max_val) + 0.01;
% min_val = min(m_var([1],1)- m_var([1],2), m_var([2],1)- m_var([2],2));
% min_val = (min_val) - 0.01;
% numticks = 4;
% ytickval = linspace(min_val,max_val,numticks);
% % ytickval = round(ytickval);
% pl.YTick = (ytickval);
% pl.YLim = [floor(min_val) ceil(max_val)];

% X options
pl.XTick = [1, 2];
pl.XTickLabel = {'Resp1', 'Resp2'};
% pl.FontWeight = 'bold';
pl.FontSize = 10;
pl.TickLength = [0.03,0.025];
pl.Box = 'off';


end