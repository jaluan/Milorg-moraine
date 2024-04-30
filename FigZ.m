% Produce manuscript figure Z

close all;
set(groot','defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

% addpath('Functions','Functions/export_fig')

%% load MC results
snr=[1,5];
str = num2str(snr(1));
if length(snr) > 1 %if multiple samples were modelled together
    for i=2:length(snr)
        str = [str,['-',num2str(snr(i))]];
    end
end

mname = ['models/MC_rockfall_sample_',str,'_v2_temp1_RFmax18ka_max5m_LSDn_sig2.mat'];
load(mname,'model');

%% Find accepted model values
uval = []; restot = []; Be = []; Al = []; C=[];
for nw = 1:model.Nwalk
    I = find(model.walker{nw}.status == 1);
    if (~isempty(I))
        uval = [uval(:,:);(model.walker{nw}.up(:,I)')];
        restot = [restot,(model.walker{nw}.restot(I)')];
        Be = [Be(:,:);model.walker{nw}.gm(1:3:6,I)'];
        Al = [Al(:,:);model.walker{nw}.gm(2:3:6,I)'];
        C = [C(:,:);model.walker{nw}.gm(3:3:6,I)'];
    end
end

%% Colors and styles
cols=[.1,.5,.3; .4,.3,.2; .5,.3,.7; .3,.3,.8; .7,.7,.2];
style={'-','-.',':','--','-'};

%% Summary figure
figure()
tiledlayout(6,8,"TileSpacing","compact",'Padding','compact')

nexttile(1,[3,6])
% scatter(uval(:,1),uval(:,2)*1e3,6,restot,'Marker','.')
scatter(uval(:,1),uval(:,2)*1e3,'Marker','.','MarkerEdgeColor',0.5*[1 1 1],'DisplayName','17MFM-01 + 05')
hold on
scatter(uval(restot<=2,1),uval(restot<=2,2)*1e3,'Marker','.',...
    'MarkerEdgeColor',0.2*[1 1 1],'handlevisibility','off')
% [X,Y] = meshgrid(uval(restot<=2,1),uval(restot<=2,2)*1e3);
% contour(X,Y,N,40);
prct=[0.05 0.25 0.5 0.75 0.95]; quant=quantile(uval(:,1),prct);
myboxplot(quant,3.5,0.2,0.7*[1 1 1],'horizontal')
prct=[0.05 0.25 0.5 0.75 0.95]; quant=quantile(uval(:,2)*1e3,prct);
myboxplot(quant,0.05,0.01,0.7*[1 1 1],'vertical')
xlabel('Duration of rockwall exposure (Myr)')
ylabel('Time of rockfall event (ka)')
set(gca,'FontSize',14), box on
xlim([0 2]), ylim([2 18])
text(0.01,17,'A','FontSize',20)
% colorbar

nexttile(25,[3,3]), hold on
for i=1:model.Nsnr
    [f,xi] = ksdensity(uval(:,3+i)); %compute kernel density
    h(snr(i))=line(xi,f,'color',cols(snr(i),:),'LineWidth',1.5,'LineStyle',style{snr(i)},...
        'DisplayName',char(model.data{i}.name));
end
% legend()
xlabel('Depth in rockwall (m)') %$D_{RW}$
ylabel('Frequency'), set(gca,'YTick',[],'FontSize',14), box on
xlim([0 2])
text(0.05,2.3,'B','FontSize',20)

nexttile(28,[3,3])
[f,xi] = ksdensity(uval(:,3)); %compute kernel density
line(xi,f,'color',0.5*[1 1 1],'LineWidth',1.5,'LineStyle',style{1},'DisplayName','17MFM-01 + 05');
xlabel('Elevation of rock fall (m a.s.l.)') %$E_{RF}$
ylabel('Frequency'), set(gca,'YTick',[],'FontSize',14), box on
% text(1420,0.0028,'c','FontSize',20)

map = parula; %colormap;
map(1,:) = [1 1 1]; %[1,1,1]
colormap(map);
t = linspace(0,2*pi); %radians

for i=1:2
    nexttile(7+(snr(i)-1)*8);
    
    Beint = linspace(min(Be(:,i)),max(Be(:,i)),40)';
    Alint = linspace(min(Al(:,i)),max(Al(:,i)),40)';
    Cint = linspace(min(C(:,i)),max(C(:,i)),40)';
    N = hist3([Al(:,i),Be(:,i)],{Alint Beint}); N = N/sum(N(:));
    [X,Y] = meshgrid(Beint,Alint);
    contourf(X,Y,N,40,'lineColor','none'); hold on
    %Plot measured error on top - ellipses:
    cx = model.data{i}.N10; cy = model.data{i}.N26; %centerpoint
    dcx = model.data{i}.dN10; dcy = model.data{i}.dN26; %1sd
    a1 = dcx ; b1 = dcy; x1 = a1*cos(t); y1 = b1*sin(t); %1sd
    a2 = 2*dcx ; b2 = 2*dcy; x2 = a2*cos(t); y2 = b2*sin(t); %2sd
    plot(x1+cx,y1+cy,'k','LineWidth',2), hold on %plot 1sd
    plot(x2+cx,y2+cy,'k','LineWidth',2) %plot 2sd
    % yl1=ylim; xl1=xlim;
    if snr(i)==3
        text(min(Be(:,i)),mean(Al(:,i)),'$^{26}$Al','FontSize',20,'Rotation',90,'HorizontalAlignment','right') 
    end
    if snr(i)==5
        text(mean(Be(:,i)),min(Al(:,i)-1.5e6),'$^{10}$Be','FontSize',20,'HorizontalAlignment','center')
    end
    set(gca,'XTick',[],'YTick',[],'FontSize',14,'visible','off'), box on
    
    nexttile(8+(snr(i)-1)*8);
    N = hist3([C(:,i),Be(:,i)],{Cint Beint}); N = N/sum(N(:));
    [X,Y] = meshgrid(Beint,Cint);
    contourf(X,Y,N,40,'lineColor','none'); hold on
    %Plot measured error on top - ellipses:
    cy = model.data{i}.N14; %centerpoint
    dcy = model.data{i}.dN14; %1sd
    b1 = dcy; y1 = b1*sin(t); %1sd
    b2 = 2*dcy; y2 = b2*sin(t); %2sd
    plot(x1+cx,y1+cy,'k','LineWidth',2), hold on %plot 1sd
    plot(x2+cx,y2+cy,'k','LineWidth',2) %plot 2sd
    if snr(i)==5
        text(mean(Be(:,i)),min(C(:,i))-4e4,'$^{10}$Be','FontSize',20,'HorizontalAlignment','center')
    end
    set(gca,'XTick',[],'YTick',[],'FontSize',14,'visible','off'), box on
    text(min(Be(:,i))+1e5,max(C(:,i)),char(model.data{i}.name),'HorizontalAlignment','right','FontSize',14)
end

%% Load second inversion
snr=2:4;
str = num2str(snr(1));
if length(snr) > 1 %if multiple samples were modelled together
    for i=2:length(snr)
        str = [str,['-',num2str(snr(i))]];
    end
end
mname = ['models/MC_rockfall_sample_',str,'_v2_temp1_RFmax18ka_max5m_LSDn_sig2.mat'];
load(mname,'model');

%% Find accepted model values
uval = []; restot = []; Be = []; Al = []; C=[];
for nw = 1:model.Nwalk
    I = find(model.walker{nw}.status == 1);
    if (~isempty(I))
        uval = [uval(:,:);(model.walker{nw}.up(:,I)')];
        restot = [restot,(model.walker{nw}.restot(I)')];
        Be = [Be(:,:);model.walker{nw}.gm(1:3:9,I)'];
        Al = [Al(:,:);model.walker{nw}.gm(2:3:9,I)'];
        C = [C(:,:);model.walker{nw}.gm(3:3:9,I)'];
    end
end

%% Add to figure
nexttile(1,[3,6])
scatter(uval(:,1),uval(:,2)*1e3,'Marker','.','MarkerEdgeColor',0.5*[1 .5 1],...
    'DisplayName','17MFM-02 + 03 + 04','handlevisibility','off')
hold on
scatter(uval(restot<=2,1),uval(restot<=2,2)*1e3,'Marker','.',...
    'MarkerEdgeColor',0.2*[1 .5 1],'handlevisibility','off')
prct=[0.05 0.25 0.5 0.75 0.95]; quant=quantile(uval(:,1),prct);
myboxplot(quant,4.5,0.2,brighten(cols(3,:),0.4),'horizontal')
prct=[0.05 0.25 0.5 0.75 0.95]; quant=quantile(uval(:,2)*1e3,prct);
myboxplot(quant,0.05,0.01,brighten(cols(3,:),0.4),'vertical')
hdl(1) = plot(nan, nan,'.','Color',0.5*[1 1 1],'MarkerSize',20,'DisplayName','17MFM-01 + 05');
hdl(2) = plot(nan, nan,'.','Color',brighten(0.5*[1 .5 1],0.2),'MarkerSize',20,'DisplayName','17MFM-02 + 03 + 04');
legend(hdl)

nexttile(25,[3,3]), hold on
for i=1:model.Nsnr
    [f,xi] = ksdensity(uval(:,3+i)); %compute kernel density
    h(snr(i))=line(xi,f,'color',cols(snr(i),:),'LineWidth',1.5,'LineStyle',style{snr(i)},...
        'DisplayName',char(model.data{i}.name));
end
legend([h(1) h(2) h(3) h(4) h(5)])

nexttile(28,[3,3])
[f,xi] = ksdensity(uval(:,3)); %compute kernel density
line(xi,f,'color',cols(3,:),'LineWidth',1.5,'LineStyle',style{2},'DisplayName','17MFM-02 + 03 + 04');
legend()
xlim([1450 1950])
text(min(xlim)+15,max(ylim)-0.00005,'C','FontSize',20,'Horiz','left', 'Vert','top')

for i=1:3
    nexttile(7+(snr(i)-1)*8);
    
    Beint = linspace(min(Be(:,i)),max(Be(:,i)),40)';
    Alint = linspace(min(Al(:,i)),max(Al(:,i)),40)';
    Cint = linspace(min(C(:,i)),max(C(:,i)),40)';
    N = hist3([Al(:,i),Be(:,i)],{Alint Beint}); N = N/sum(N(:));
    [X,Y] = meshgrid(Beint,Alint);
    contourf(X,Y,N,40,'lineColor','none'); hold on
    %Plot measured error on top - ellipses:
    cx = model.data{i}.N10; cy = model.data{i}.N26; %centerpoint
    dcx = model.data{i}.dN10; dcy = model.data{i}.dN26; %1sd
    a1 = dcx ; b1 = dcy; x1 = a1*cos(t); y1 = b1*sin(t); %1sd
    a2 = 2*dcx ; b2 = 2*dcy; x2 = a2*cos(t); y2 = b2*sin(t); %2sd
    plot(x1+cx,y1+cy,'k','LineWidth',2), hold on %plot 1sd
    plot(x2+cx,y2+cy,'k','LineWidth',2) %plot 2sd
    % yl1=ylim; xl1=xlim;
    if snr(i)==3
        text(min(Be(:,i)),mean(Al(:,i)),'$^{26}$Al','FontSize',20,'Rotation',90,'HorizontalAlignment','right') 
    end
    if snr(i)==5
        text(mean(Be(:,i)),min(Al(:,i)-1.5e6),'$^{10}$Be','FontSize',20,'HorizontalAlignment','center')
    end
    set(gca,'XTick',[],'YTick',[],'FontSize',14,'visible','off'), box on
    
    nexttile(8+(snr(i)-1)*8);
    N = hist3([C(:,i),Be(:,i)],{Cint Beint}); N = N/sum(N(:));
    [X,Y] = meshgrid(Beint,Cint);
    contourf(X,Y,N,40,'lineColor','none'); hold on
    %Plot measured error on top - ellipses:
    cy = model.data{i}.N14; %centerpoint
    dcy = model.data{i}.dN14; %1sd
    b1 = dcy; y1 = b1*sin(t); %1sd
    b2 = 2*dcy; y2 = b2*sin(t); %2sd
    plot(x1+cx,y1+cy,'k','LineWidth',2), hold on %plot 1sd
    plot(x2+cx,y2+cy,'k','LineWidth',2) %plot 2sd
    if snr(i) == 3
        text(max(Be(:,i)),mean(C(:,i)),'$^{14}$C','FontSize',20,'Rotation',90,'HorizontalAlignment','right')
    end
    set(gca,'XTick',[],'YTick',[],'FontSize',14,'visible','off'), box on
    text(min(Be(:,i))+1e5,max(C(:,i)),char(model.data{i}.name),'HorizontalAlignment','right','FontSize',14)
end

nexttile(47);
text(0,1,'D','FontSize',20)
set(gca,'Xtick','','Ytick','','visible','off')

set(gcf,'units','normalized','position',[.1,.3,.75,.6]);

%% Set interpreter to default
set(groot','defaulttextinterpreter','default');
set(groot, 'defaultAxesTickLabelInterpreter','default'); 
set(groot, 'defaultLegendInterpreter','default');