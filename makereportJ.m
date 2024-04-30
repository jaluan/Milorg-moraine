function makereportJ(snr)

% function makereportJ(snr,CCurve,emergence,ss)

% This function plots: 
% a. Measured vs. modelled Be-Al, 
% b. histogram of d18O-threshold parameter and c. last deglaciation pm,
% d. 2D-histogram of exhumation pathways of sample to surface through time,
% (e-h). If the model contains 14C, the measured vs modelled
% concentrations compared to 10Be and 26Al will be shown,
% i. Distribution of all model parameters with separate coloured curves for 
% each walker.
% j. parameter values vs model-time through burn-in (gray) and post-burn in
% (colours, rejected models are grey). Separate colours for different walkers.

% Inputs: 
%  - snr: scalar or vector containing sample numbers of samples in input data
%       file generated from compile_mn_data.m.

% Note that the inversion must have been generated first by use of 
% rockfallMC-code

% Written by David Lundbek Egholm, Aarhus University
% Modified by Jane Lund Andersen to also i) include 14C, ii) combine with 
% rockfall forward model

close all;
set(groot','defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

addpath('Functions','Functions/export_fig')
CNprop = getCNprop(); %Cosmogenic properties

%% load MC results
str = num2str(snr(1));
if length(snr) > 1 %if multiple samples were modelled together
    for i=2:length(snr)
        str = [str,['-',num2str(snr(i))]];
    end
end

% mname = ['models/MC_rockfall_sample_',str,'_ntemp10_RFmax18ka_max5m.mat'];
% mname = ['models/MC_rockfall_sample_',str,'_v2_temp10_RFmax18ka_max5m.mat'];
mname = ['models/MC_rockfall_sample_',str,'_v2_temp1_RFmax18ka_max5m_LSDn_sig2.mat'];
load(mname,'model');

%% Name of output file
fname = ['models/reports/Report_sample_',str,'_v2_temp1_RFmax18ka_max5m_LSDn_sig2.pdf'];

%% colors
col1 = 0.9*[1,1,1];
col2 = 0.6*[1,1,1];
col3 = 0.3*[1,1,1];

map = colormap;
nc = length(map);
xc = linspace(0,1,nc);
rc = map(:,1);
gc = map(:,2);
bc = map(:,3);
wcol=zeros(model.Nwalk,3);
for i=1:model.Nwalk
    ri = interp1(xc,rc,(i-1)/model.Nwalk);
    gi = interp1(xc,gc,(i-1)/model.Nwalk);
    bi = interp1(xc,bc,(i-1)/model.Nwalk);
    wcol(i,:) = [ri,gi,bi];
end

%% figures showing walkers *************

set(gcf,'units','normalized','position',[.1,.3,.4,.6]);
set(gca,'position',[0,0,1,1],'visible','off');
set(gca,'xlim',[0,1],'ylim',[0,1]);
set(gcf,'Name','Walker information');
% ax0 = gca;

[np,~] = numSubplots(model.Nmp);

%plot models parameters as a function of model time
for i = 1:model.Nmp
    subplot(np(1),np(2),i); hold on; box on; grid on;
    xlabel('model nr.');
    ylabel(model.mp{i}.name);
    set(gca,'xlim',[0,length(model.walker{1}.status)],'ylim',[model.mp{i}.vmin,model.mp{i}.vmax]);
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == -2); %rejected during burnin
        plot(I,model.walker{nw}.up(i,I),'.','color',col1);
        I = find(model.walker{nw}.status == 0); %rejected after burnin
        plot(I,model.walker{nw}.up(i,I),'.','color',col2);
        I = find(model.walker{nw}.status == -1); %accepted during burnin
        plot(I,model.walker{nw}.up(i,I),'.','color',col3);
        I = find(model.walker{nw}.status == 1); %accepted after burnin
        plot(I,model.walker{nw}.up(i,I),'.','color',wcol(nw,:));
    end
end

print('temp3.pdf','-dpdf','-fillpage'); %print temporary pdf, appended below

%% Parameter distributions ************
figure()
set(gcf,'units','normalized','position',[.3,.2,.4,.6]);
set(gcf,'Name','Parameter distributions');

[np,~] = numSubplots(model.Nmp);

%loop parameters
for i=1:model.Nmp
    
    subplot(np(1),np(2),i); 
    hold on; box on; grid on;
    ylabel('Frequency');
    xlabel(model.mp{i}.name);
    set(gca,'xlim',[model.mp{i}.vmin, ...
                    model.mp{i}.vmax]);
    uval = [];
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        if (~isempty(I))
            [f,xi] = ksdensity(model.walker{nw}.up(i,I));
            line(xi,f,'color',wcol(nw,:));
            uval = [uval(:);model.walker{nw}.up(i,I)'];
        end
    end
    if (~isempty(uval))
        [f,xi] = ksdensity(uval); %compute kernel density
    end
    line(xi,f,'color','k','linewidth',2);

end

print('temp2.pdf','-dpdf','-fillpage'); %print temporary pdf, appended below

%% set plot margins
lpm = 0.1;
ddm = 0.02;

%% 14C model-data fit figures 
%14C figures
if (any(ismember(model.data{1,1}.nuclides,3))) %if 14C in model
    % initiate figure
    figure()
    set(gcf,'units','normalized','position',[.3,.2,.4,.6]);
    set(gcf,'Name','14C model-data fit');
    ax0 = gca;
    set(ax0,'position',[0,0,1,1]);
    set(gca,'visible','off');
    dpx = (1-2*lpm-(model.Nsnr-1)*ddm)/model.Nsnr;
    
    for ns = 1:model.Nsnr %loop samples

        axes('position',[lpm+(ns-1)*(ddm+dpx),0.7,dpx,0.25]);
        hold on; box on; grid on;
        title(model.data{ns}.name);

        n0 = (ns-1)*model.Nds(ns);
        
        % retrieve modelled nuclide distributions
        Be = [];
        Al = [];
        C = [];
        for nw = 1:model.Nwalk
            I = find(model.walker{nw}.status == 1);
            Be = [Be(:);model.walker{nw}.gm(n0+1,I)'];
            Al = [Al(:);model.walker{nw}.gm(n0+2,I)'];
            C = [C(:);model.walker{nw}.gm(n0+3,I)'];
        end
        %calculate ratios
        CBe = C./Be;
        CAl = C./Al;
        
        Beint = linspace(min(Be),max(Be),40)';
        CBeint = linspace(0.5*min(CBe),1.5*max(CBe),40)';
        Alint = linspace(min(Al),max(Al),40)';
        CAlint = linspace(0.5*min(CAl),1.5*max(CAl),40)';

        %C-Be banana
        xlabel(['Normalized NBe (atoms/g), sample ',num2str(ns)]);
        if (ns == 1), ylabel('N14C/N10Be');
        else, set(gca,'yticklabel',[]);
        end

        N = hist3([CBe,Be],{CBeint Beint});
        N = N/sum(N(:));
        
        %Normalize to SLHL using surface production rate (spallation+muons)
        nfac = CNprop.PBe0/(model.data{ns}.production.P10spal +...
            model.data{ns}.production.P10_m1 + model.data{ns}.production.P10_m2);

        [X,Y] = meshgrid(Beint*nfac,CBeint);
        contour(X,Y,N,40);
    
        errorbar(model.data{ns}.N10*nfac,model.data{ns}.r1410,model.data{ns}.dN10*nfac,'horizontal','.k');
        errorbar(model.data{ns}.N10*nfac,model.data{ns}.r1410,model.data{ns}.dr1410,'vertical','.k');
        
        %C-Al banana
        axes('position',[lpm+(ns-1)*(ddm+dpx),0.075,dpx,0.25]);
        hold on; box on; grid on;
        
        xlabel(['Normalized NAl (atoms/g), sample ',num2str(ns)]);
        if (ns == 1), ylabel('N14C/N26Al');
        else, set(gca,'yticklabel',[]);
        end

        N = hist3([CAl,Al],{CAlint Alint});
        N = N/sum(N(:));
        
        %Normalize to SLHL using surface production rate (spallation+muons)
        nfac = CNprop.PAl0/(model.data{ns}.production.P26spal + ...
            model.data{ns}.production.P26_m1 + model.data{ns}.production.P26_m2);
        
        [X,Y] = meshgrid(Alint*nfac,CAlint);
        contour(X,Y,N,40);
        
        model.data{ns}.r1426 = model.data{ns}.N14./model.data{ns}.N26;
        model.data{ns}.dr1426 = model.data{ns}.r1426.*sqrt((model.data{ns}.dN26./...
        model.data{ns}.N26).^2+(model.data{ns}.dN14./model.data{ns}.N14).^2);
    
        errorbar(model.data{ns}.N26*nfac,model.data{ns}.r1426,model.data{ns}.dN26*nfac,'horizontal','.k');
        errorbar(model.data{ns}.N26*nfac,model.data{ns}.r1426,model.data{ns}.dr1426,'vertical','.k');
       
        % Histogram/kernel density with modelled concentrations
        axes('position',[lpm+(ns-1)*(ddm+dpx),0.4,dpx,0.2]);
        hold on; box on; grid on;
        xlabel('N14 (atoms/g)'); ylabel ('Probability')
        histogram(C,'Normalization','probability')

        % Measured data on top
        errorbar(model.data{ns}.N14,.015,model.data{ns}.dN14,...
        'horizontal','.k','Linewidth',1.5,'Color',[.7 .2 .2]);
    
    end
    print('temp1.pdf','-dpdf','-fillpage'); %print temporary pdf, appended below
end

%% ******* main report figure *************
% Initiate figure
figure;
set(gcf,'papertype','a4');
set(gcf,'units','centimeters','position',[5,5,21,29.7]);
set(gcf,'Name','Report');
ax0 = gca;
set(ax0,'position',[0,0,1,1]);
set(gca,'visible','off');
text(0.05,0.97,['File: ',mname],'HorizontalAlignment','left','fontsize',14);

dpx = (1-2*lpm-(model.Nsnr-1)*ddm)/model.Nsnr;

% Al-Be banana
for ns = 1:model.Nsnr %loop samples

    axes('position',[lpm+(ns-1)*(ddm+dpx),0.675,dpx,0.25]);
    hold on; box on; grid on;
    set(gca,'ylim',[3,7.5]);
    
    %title
%     Sname=['$',model.data{1,ns}.name{1},'$'];
%     title(Sname);
    title(model.data{ns}.name);
    
    n0 = (ns-1)*model.Nds(ns);
    
    Be = [];
    Al = [];
    
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        Be = [Be(:);model.walker{nw}.gm(n0+1,I)'];
        Al = [Al(:);model.walker{nw}.gm(n0+2,I)'];
    end
    AlBe = Al./Be;
    
    Beint = linspace(min(Be),max(Be),40)';
    %Alint = linspace(min(Al),max(Al),40)';
    AlBeint = linspace(0.5*min(AlBe),1.5*max(AlBe),40)';
    
    xlabel(['Normalized NBe (atoms/g), sample ',num2str(ns)]);
    if (ns == 1), ylabel('NAl/NBe');
    else, set(gca,'yticklabel',[]);
    end
    
    N = hist3([AlBe,Be],{AlBeint Beint});
    N = N/sum(N(:));

    nfac = CNprop.PBe0/(model.data{ns}.production.P10spal +...
        model.data{ns}.production.P10_m1 + model.data{ns}.production.P10_m2);
    %nfac2 = CNprop.PAl0/(model.data{ns}.production.P26spal +...
    %    model.data{ns}.production.P26_m1 + model.data{ns}.production.P26_m2);
    [X,Y] = meshgrid(Beint*nfac,AlBeint);
    contour(X,Y,N,40);
    
    errorbar(model.data{ns}.N10*nfac,model.data{ns}.r2610,model.data{ns}.dN10*nfac,'horizontal','.k');
    errorbar(model.data{ns}.N10*nfac,model.data{ns}.r2610,model.data{ns}.dr2610,'vertical','.k');
    
%     % Histogram/kernel density with modelled N26 concentrations
%     axes('position',[0.65,0.8,0.2,0.1]);
%     hold on; box on; grid on;
%     xlabel('N26 (atoms/g)'); %ylabel ('Probability')
%     histogram(Al,'Normalization','probability')
% 
%     % Measured data on top
%     errorbar(model.data{ns}.N26,.015,model.data{ns}.dN26,...
%     'horizontal','.k','Linewidth',1.5,'Color',[.7 .2 .2]);

end

% glaciation history parameter distibutions

dpx = (1-2*lpm-ddm)/2;

%loop parameters
for i=1:2
    
    axes('position',[lpm+(i-1)*(ddm+dpx),0.37,dpx,0.25]);
    hold on; box on; grid on;
    
    if (i == 1), ylabel('Frequency'); end
    xlabel(model.mp{i}.name);
    set(gca,'xlim',[model.mp{i}.vmin, ...
                    model.mp{i}.vmax]);
    set(gca,'yticklabel',[]);
    uval = [];
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        if (~isempty(I))
            [f,xi] = ksdensity(model.walker{nw}.up(i,I));
            line(xi,f,'color',wcol(nw,:));
            uval = [uval(:);model.walker{nw}.up(i,I)'];
        end
    end
    if (~isempty(uval))
        [f,xi] = ksdensity(uval); %compute kernel density
    end
    line(xi,f,'color','k','linewidth',2);
end

% Crossplot
uval = []; res = [];
for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == 1); %accepted after burnin
        uval = [uval model.walker{nw}.up(:,I)]; %u or up here? u repeats accepted model and corresponds to restot
        res = [res; model.walker{nw}.restot(I)];
end

axes('position',[0.1,0.05,0.39,0.25]); %half-width panel=0.39 wide
hold on; box on;
plot(uval(1,:),uval(2,:),'.')
plot(uval(1,res<=2),uval(2,res<=2),'.')
xlabel(model.mp{1}.name);  
ylabel(model.mp{2}.name);
legend('Accepted','$<2 \sigma$','location','southeast')

axes('position',[0.51,0.05,0.39,0.25]); %half-width=0.39 wide,full=0.8
hold on; box on;
plot(uval(1,:),uval(3,:),'.')   
plot(uval(1,res<=2),uval(3,res<=2),'.')
xlabel(model.mp{1}.name);  
ylabel(model.mp{3}.name);
set(gca,'YAxisLocation','right')


%% save figures to one pdf and delete temporary files
print(fname,'-dpdf','-fillpage');

if (any(ismember(model.data{1,1}.nuclides,3))) %if 14C in model
    append_pdfs(fname, 'temp1.pdf')
    delete('temp1.pdf')
end

append_pdfs(fname, 'temp2.pdf', 'temp3.pdf')
delete('temp2.pdf','temp3.pdf')
%% Set interpreter to default
set(groot','defaulttextinterpreter','default');
set(groot, 'defaultAxesTickLabelInterpreter','default'); 
set(groot, 'defaultLegendInterpreter','default');