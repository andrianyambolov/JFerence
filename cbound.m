function [pplot, fil] = cbound(Y, varargin)
% Plots error bands.

defaultAlpha = 0.2;
defaultLineWidth = 4;
co = get(gca,'ColorOrder');
defaultColor = co(get(gca,'ColorOrderIndex'),:);
defaultFaceColor = defaultColor;
defaultLineStyle = '-';
defaultXAxis = 1:size(Y,1);
defaultMarker = 'none';
% defaultAxisLabels = ;

p = inputParser;
% validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p,'Y', @ismatrix);
addOptional(p,'alpha',defaultAlpha);
addOptional(p,'LineWidth',defaultLineWidth);
addOptional(p,'Color', defaultColor);
addOptional(p,'FaceColor', defaultColor);
addOptional(p,'LineStyle', defaultLineStyle);
addOptional(p,'Marker', defaultMarker);
addOptional(p,'XAxis', defaultXAxis);
parse(p, Y, varargin{:});

T = size(p.Results.Y,1); 

x_axis = p.Results.XAxis;
x_plot =[x_axis, fliplr(x_axis)];
hold on
if size(Y,2) == 3
    y_plot=[Y(:,1)', flipud(Y(:,3))'];
    fil = fill(x_plot, y_plot, 1, 'facecolor', p.Results.FaceColor, 'edgecolor',  p.Results.FaceColor, 'facealpha', p.Results.alpha);
    pplot = plot(x_axis, Y(:,2),'LineStyle', p.Results.LineStyle, 'Marker', p.Results.Marker, 'Color', p.Results.Color, 'linewidth',  p.Results.LineWidth);
elseif size(Y,2) == 1
    pplot = plot(x_axis, Y(:,1),'LineStyle', p.Results.LineStyle, 'Marker', p.Results.Marker, 'Color', p.Results.Color, 'linewidth',  p.Results.LineWidth);
    fil = [];
else
    y_plot=[Y(:,1)', flipud(Y(:,2))'];
    fil = fill(x_plot, y_plot, 1, 'facecolor', p.Results.FaceColor, 'edgecolor', p.Results.FaceColor, 'facealpha', p.Results.alpha);
    pplot = [];
end
axis tight
end