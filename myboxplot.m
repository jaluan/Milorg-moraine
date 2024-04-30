function myboxplot(quantiles,y,dy,col,dir)

% Jane Lund Andersen, March 2023

switch dir
    case 'horizontal'
        line([quantiles(1) quantiles(5)],[y y],'color',col,'Linewidth',1.5,'handlevisibility','off') %whisker
        patch([quantiles(2) quantiles(4) quantiles(4) quantiles(2)],...
            [y-dy y-dy y+dy y+dy],col,'handlevisibility','off') %box
        line([quantiles(3) quantiles(3)],[y-dy y+dy],'color','k','Linewidth',1.5,'handlevisibility','off') %median

    case 'vertical'
        line([y y],[quantiles(1) quantiles(5)],'color',col,'Linewidth',1.5,'handlevisibility','off') %whisker
        patch([y-dy y-dy y+dy y+dy],...
            [quantiles(2) quantiles(4) quantiles(4) quantiles(2)],col,...
            'handlevisibility','off') %box
        line([y-dy y+dy],[quantiles(3) quantiles(3)],'color','k','Linewidth',1.5,'handlevisibility','off') %median
end