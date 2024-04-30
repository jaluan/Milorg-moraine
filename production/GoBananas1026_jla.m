
function out = GoBananas1026_jla(sm,bur_lines)

% 'sm' is a string for the scaling model to use - 'St' for Lal/Stone, 'LS' for LSDn
% 'bur_lines' is a string - 'y' or 'n' to include or not
% Version 1/23 - NL

file = input('Please enter the name of the sample data file: ','s');
FID = fopen(file);
data = textscan(FID,'%s %n %n %n %s %n %n %n %n %n %n %s %n %n %s');
fclose(FID);
dstring='';

%Make the sample structure.

all_sample_name = data{1};
all_lat = data{2};
all_long = data{3};
all_elv = data{4};
all_pressure = data{4};
all_aa = data{5};
all_thick = data{6};
all_rho = data{7};
all_othercorr = data{8};
all_E = data{9};
all_N10 = data{10};
all_delN10 = data{11};
all_be_std_name = data{12};
all_N26 = data{13};
all_delN26 = data{14}; 
all_al_std_name = data{15};

%Load data and constants
% make_consts_LSD;
load consts_LSDn;

num_samples = length(all_lat);

P10_LS = zeros(size(all_lat));
P26_LS = zeros(size(all_lat));

%Calculate total SLHL P rates for 26Al and 10Be
if strcmp(sm,'St')
    P26 = consts.P26_ref_St;
    P10 = consts.P10_ref_St;
elseif strcmp(sm,'LS')
    P26 = consts.P26_ref_LS;
    P10 = consts.P10_ref_LS;
end

for a = 1:num_samples	
	if all_N26(a) ~= 0
		all_isN26(a) = 1;
    else
		all_isN26(a) = 0; 
    end

	if all_delN26(a) ~= 0
		all_isdelN26(a) = 1; 
    else
		all_isdelN26(a) = 0; 
    end
	
	if all_N10(a) ~= 0
		all_isN10(a) = 1;
    else
		all_isN10(a) = 0;
    end

	if all_delN10(a) ~= 0
		all_isdelN10(a) = 1; 
    else
		all_isdelN10(a) = 0; 
    end
	
	% catch mismatches;
	
	if (~all_isN26(a) && ~all_isN10(a))
   		error(['Need either Al-26 or Be-10 concentration - line ' int2str(a)]);
	elseif (all_isN26(a) && ~all_isdelN26(a)) || (~all_isN26(a) && all_isdelN26(a))
    		error(['Need both Al-26 concentration and uncertainty - line ' int2str(a)]);
	elseif (all_isN10(a) && ~all_isdelN10(a)) || (~all_isN10(a) && all_isdelN10(a))
    		error(['Need both Be-10 concentration and uncertainty - line ' int2str(a)]);
    end
end

for a = 1:num_samples
    sample.sample_name{a} = all_sample_name{a};
    sample.lat = all_lat(a);
    sample.long = all_long(a);
    sample.aa = all_aa{a};
    if strcmp(all_aa{a},'std') || strcmp(all_aa{a},'ant')
			% store the elevation value
			sample.elv = all_elv(a);
	elseif  strcmp(all_aa{a},'pre')
			% store the pressure value
			sample.pressure = all_pressure(a);
    end
    sample.thick = all_thick(a);
    sample.E = all_E(a);
    sample.rho = all_rho(a);
    sample.othercorr = all_othercorr(a);
		sample.N26 = all_N26(a);
		sample.delN26 = all_delN26(a);
    if all_isN26(a)
        if (strcmp(all_al_std_name(a),'KNSTD'))
            sample.N26 = all_N26(a).*consts.al_stds_cfs(1);
            sample.delN26 = all_delN26(a).*consts.al_stds_cfs(1);
        elseif (strcmp(all_al_std_name(a),'ZAL94'))
            sample.N26 = all_N26(a).*consts.al_stds_cfs(2);
            sample.delN26 = all_delN26(a).*consts.al_stds_cfs(2);
        elseif (strcmp(all_al_std_name(a),'AL09'))
            sample.N26 = all_N26(a).*consts.al_stds_cfs(3);
            sample.delN26 = all_delN26(a).*consts.al_stds_cfs(3);
        elseif (strcmp(all_al_std_name(a),'ZAL94N'))
            sample.N26 = all_N26(a).*consts.al_stds_cfs(4);
            sample.delN26 = all_delN26(a).*consts.al_stds_cfs(4);
        elseif (strcmp(all_al_std_name(a),'SMAL11'))
            sample.N26 = all_N26(a).*consts.al_stds_cfs(5);
            sample.delN26 = all_delN26(a).*consts.al_stds_cfs(5);
        elseif (strcmp(all_al_std_name(a),'Z92-0222'))
            sample.N26 = all_N26(a).*consts.al_stds_cfs(6);
            sample.delN26 = all_delN26(a).*consts.al_stds_cfs(6);
        elseif (strcmp(all_al_std_name(a),'0'))
            sample.N26 = all_N26(a).*consts.al_stds_cfs(7);	
            sample.delN26 = all_delN26(a).*consts.al_stds_cfs(7);
        end
     end
    
    if all_isN10(a)
        if (strcmp(all_be_std_name(a),'07KNSTD'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(1);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(1);
        elseif (strcmp(all_be_std_name(a),'KNSTD'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(2);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(2);
        elseif (strcmp(all_be_std_name(a),'NIST_Certified'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(3);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(3);
        elseif (strcmp(all_be_std_name(a),'LLNL31000'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(4);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(4);
        elseif (strcmp(all_be_std_name(a),'LLNL10000'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(5);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(5);
        elseif (strcmp(all_be_std_name(a),'LLNL3000'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(6);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(6);
        elseif (strcmp(all_be_std_name(a),'LLNL1000'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(7);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(7);
        elseif (strcmp(all_be_std_name(a),'LLNL300'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(8);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(8);
        elseif (strcmp(all_be_std_name(a),'NIST_30000'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(9);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(9);
        elseif (strcmp(all_be_std_name(a),'NIST_30200'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(10);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(10);
        elseif (strcmp(all_be_std_name(a),'NIST_30300'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(11);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(11);
        elseif (strcmp(all_be_std_name(a),'NIST_30600'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(12);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(12);
        elseif (strcmp(all_be_std_name(a),'NIST_27900'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(13);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(13);
        elseif (strcmp(all_be_std_name(a),'S555'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(14);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(14);
        elseif (strcmp(all_be_std_name(a),'S2007'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(15);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(15);
        elseif (strcmp(all_be_std_name(a),'BEST433'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(16);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(16);
        elseif (strcmp(all_be_std_name(a),'BEST433N'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(17);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(17);
        elseif (strcmp(all_be_std_name(a),'S555N'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(18);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(18);
        elseif (strcmp(all_be_std_name(a),'S2007N'))
            sample.N10 = all_N10(a).*consts.be_stds_cfs(19);
            sample.delN10 = all_delN10(a).*consts.be_stds_cfs(19);

        end
    end
%     % 1a. Thickness scaling factor. 
% 
%     if sample.thick > 0;
%         sample.thickSF = thickness(sample.thick,consts.Lsp,sample.rho);
%     else 
%         sample.thickSF = 1;
%     end;
% 
% 1b. If no pressure entered yet, create it from the elevation
    % using the appropriate atmosphere
    % Change in version 2: This looks up the surface pressure from 
    % NCAR map to use in the standard atmosphere equation.
    % Obviously, it is always better to estimate pressure from local
    % station data. 

    if (~isfield(sample,'pressure'))
        if (strcmp(sample.aa,'std'))
            % Old code
            % sample.pressure = stdatm(sample.elv);
            % New code
            sample.pressure = ERA40atm(sample.lat,sample.long,sample.elv);
        elseif (strcmp(sample.aa,'ant'))
            sample.pressure = antatm(sample.elv);
        end
    end

% Catch confusion with pressure submission. If sample.pressure is already 
% set, it should have a submitted value. If zero, something is wrong. 
% This should never happen in online use. 

if sample.pressure == 0
    error(['Sample.pressure = 0 on sample ' sample.sample_name]);
end
% 
% 
sample.shielding = 1;

disp(sample.sample_name(a));

    scaling26 = ScalingLSDLal(sample,consts,26);
    scaling10 = ScalingLSDLal(sample,consts,10);

    al_results = get_age_LSDLal(sample,consts,26,scaling26,a);
    be_results = get_age_LSDLal(sample,consts,10,scaling10,a);

    if strcmp(sm,'St')  
        P10_St(a) = scaling10.SF_St(1).*P10;
        P26_St(a) = scaling26.SF_St(1).*P26;
    elseif strcmp(sm,'LS')  
        t_app = scaling10.tv(1:find(scaling10.tv>be_results.t_LS,1,'first'));
        weights = exp(-t_app.*consts.l10);
%         weights = exp(-t_app.*consts.l10).*t_app;
        P10_LS(a) = P10.*sum(scaling10.SF_LS(1:find(scaling10.tv>be_results.t_LS,1,'first')).*weights)./sum(weights);
        t_app = scaling26.tv(1:find(scaling26.tv>al_results.t_LS,1,'first'));
        weights = exp(-t_app.*consts.l26);
%         weights = exp(-t_app.*consts.l26).*t_app;
        P26_LS(a) = P26.*sum(scaling26.SF_LS(1:find(scaling26.tv>al_results.t_LS,1,'first')).*weights)./sum(weights);
    end

    % JLA added the following three lines
    sample.production.Lspal=150; %spallation attenuation length
    sf_spal = exp(-sample.thick.*sample.rho./2./sample.production.Lspal); %thickness shielding factor
    sample.production.P26spal=P26_LS(a)*sf_spal*sample.othercorr; %production corrected for sample thickness and topographic shielding

    if strcmp(sm,'LS') 
        GDRc = [-448.004 1189.18 -1152.15 522.061 -103.241 6.89901 0];
        RcEst = polyval(GDRc,cos(d2r(sample.lat)));
        
        if (RcEst < 0) RcEst = 0; end
    end
        
    %Get muon production rates for 26Al and 10Be
    mconsts.Natoms = consts.NatomsQtzO;
    mconsts.sigma0.LS = consts.sigma010.LS;
    mconsts.sigma0.Lm = consts.sigma010.Lm;
    mconsts.sigma0.St = consts.sigma010.St;
    mconsts.fstar.LS = consts.fstar10.LS;
    mconsts.fstar.Lm = consts.fstar10.Lm;
    mconsts.fstar.St = consts.fstar10.St;
%     mconsts.delsigma190 = consts.delsigma190_10; % not used
    mconsts.k_negpartial = consts.k_negpartial_10;
%     mconsts.delk_neg = consts.delk_neg10; % not used
    mconsts.mfluxRef = consts.mfluxRef;
    
    if strcmp(sm,'St')  
        muSt = P_mu_total_alpha1((sample.thick.*sample.rho./2),sample.pressure,mconsts,'yes');
        P10_mu(a) = muSt.P_fast_St + muSt.P_neg_St;
    elseif strcmp(sm,'LS')  
        muLS = P_mu_totalLSD((sample.thick.*sample.rho./2),sample.pressure,RcEst,consts.SPhiInf,mconsts,'yes');
%         muLS = P_mu_total_alpha1((sample.thick.*sample.rho./2),sample.pressure,mconsts,'yes');
        P10_mu(a) = muLS.P_fast_LS + muLS.P_neg_LS;
%         P10_mu(a) = muLS.P_fast_St + muLS.P_neg_St;
    end
    %

    mconsts.Natoms = consts.NatomsQtzSi; %Changed from original NatomsQtzO - 26Al produced from Si, not O
    mconsts.sigma0.LS = consts.sigma026.LS;
    mconsts.sigma0.Lm = consts.sigma026.Lm;
    mconsts.sigma0.St = consts.sigma026.St;
    mconsts.fstar.LS = consts.fstar26.LS;
    mconsts.fstar.Lm = consts.fstar26.Lm;
    mconsts.fstar.St = consts.fstar26.St;
%     mconsts.delsigma0 = consts.delsigma190_26; % not used
    mconsts.k_negpartial = consts.k_negpartial_26;
%     mconsts.delk_neg = consts.delk_neg26; % not used
    mconsts.mfluxRef = consts.mfluxRef;

    if strcmp(sm,'St')  
        muSt = P_mu_total_alpha1((sample.thick.*sample.rho./2),sample.pressure,mconsts,'yes');
        P26_mu(a) = muSt.P_fast_St + muSt.P_neg_St;
    elseif strcmp(sm,'LS')  
        muLS = P_mu_totalLSD((sample.thick.*sample.rho./2),sample.pressure,RcEst,consts.SPhiInf,mconsts,'yes');
%         muLS = P_mu_total_alpha1((sample.thick.*sample.rho./2),sample.pressure,mconsts,'yes');
        P26_mu(a) = muLS.P_fast_LS + muLS.P_neg_LS;
%         P26_mu(a) = muLS.P_fast_St + muLS.P_neg_St;

  % JLA added the following lines
        maxZ = 1200; %1200 g/cm2=4.5m with rho ~2.65-2.7
        p26_muons = fit_Pmu_with_exp_LSD_jla(sample.pressure,0,maxZ,2,mconsts,RcEst,consts.SPhiInf,0);
            sample.production.P26_Lm1=p26_muons.L(1); %attenuation, first exponential term
            sample.production.P26_Lm2=p26_muons.L(2); %attenuation, second exponential term
            sample.production.P26_m1 = p26_muons.P(1); %production, first exponential term
            sample.production.P26_m2 = p26_muons.P(2); %production, second exponential term
            shield_fac26_m1 = exp(-sample.thick.*sample.rho./2./p26_muons.L(1)); %sample thickness correction
            shield_fac26_m2 = exp(-sample.thick.*sample.rho./2./p26_muons.L(2)); %sample thickness correction
            sample.production.P26_m1 = sample.production.P26_m1*shield_fac26_m1*sample.othercorr; % production, first exponential term, corrected for sample thickness and topographic shielding
            sample.production.P26_m2 = sample.production.P26_m2*shield_fac26_m2*sample.othercorr; %production, second exponential term, corrected for sample thickness and topographic shielding

    end
    samples{a}=sample; %jla
end

%Calculate total P rates for 26Al and 10Be
if strcmp(sm,'St')  
    P26_tot = (P26_St) + P26_mu;
    P10_tot = (P10_St) + P10_mu;
    P10_spall_frac = mean(P10_St./P10_tot);
    P26_spall_frac = mean(P26_St./P26_tot);
elseif strcmp(sm,'LS')
    P26_tot = (P26_LS)' + P26_mu;
    P10_tot = (P10_LS)' + P10_mu;
    P10_spall_frac = mean(P10_LS'./P10_tot);
    P26_spall_frac = mean(P26_LS'./P26_tot);
end  

% create data for the simple-exposure line and the simple-erosion line. 
tempt = logspace(2,7,100);
% tempt = (0:1e4:1e7);
be1 = (P10_spall_frac./consts.l10).*(1-exp(-consts.l10.*tempt))+((1-P10_spall_frac)./consts.l10).*(1-exp(-consts.l10.*tempt));
a1 = (P26_spall_frac./consts.l26).*(1-exp(-consts.l26.*tempt))+((1-P26_spall_frac)./consts.l26).*(1-exp(-consts.l26.*tempt));
bound1a = [be1' (a1./be1)']';
tempe = logspace(-9,1,100);
be2 = P10_spall_frac./(consts.l10 + mean(all_rho).*tempe./160)+(1-P10_spall_frac)./(consts.l10 + mean(all_rho).*tempe./1500);
a2 = P26_spall_frac./(consts.l26 + mean(all_rho).*tempe./160)+(1-P26_spall_frac)./(consts.l26 + mean(all_rho).*tempe./1500);
bound2a = [be2' (a2./be2)']';

N10_norm = all_N10./P10_tot';
delN10_norm = N10_norm.*(all_delN10./all_N10);
N26_norm = all_N26./P26_tot';
delN26_norm = N26_norm.*(all_delN26./all_N26);

%Calculate Al/Be Ratio
% r2610 = all_N26./all_N10;
r2610 = N26_norm./N10_norm;
% drdN26 = 1./all_N10;
% drdN10 = -all_N26./(all_N10.^2);
% delR = sqrt( (all_delN10.*drdN10).^2 + (all_delN26.*drdN26).^2 );
% delR_norm = (delR./r2610).*r2610_norm;

drdN26 = 1./N10_norm;
drdN10 = -N26_norm./(N10_norm.^2);
delR = sqrt( (delN10_norm.*drdN10).^2 + (delN26_norm.*drdN26).^2 );

%Do some calculations for burial lines
% ind = [2:20:(length(a1)/2) (length(a1)/2)+20:100:length(a1)];
% ind = [2 4 11 26 51 121 201 401 2001];
ind = 5:5:length(a1);
num_a =  a1(ind);
num_be = be1(ind);
%Create burial time vector
t_bur = (0:5e5:6e6);

decayed_a = zeros(length(t_bur),length(num_a));
decayed_be =  decayed_a;
for i = 1:length(num_a)
    decayed_a(:,i) = num_a(i).*exp(-consts.l26.*t_bur);
    decayed_be(:,i) = num_be(i).*exp(-consts.l10.*t_bur);
end
decayed_r = decayed_a./decayed_be;

%% Banana Plot
f1 = figure('PaperSize',[10 10],'Position',[500 500 600 600]);
steadyExp = semilogx(bound1a(1,:),bound1a(2,:),'k-','LineWidth',1);
hold on;
steadyEro = semilogx(bound2a(1,:),bound2a(2,:),'-','LineWidth',1,'Color','#0072BD');
hold on;
semilogx(decayed_be,decayed_r,'-.','Color',[0.3 0.3 0.3]);
for a = 1:length(all_lat)
    gca;
%     ellipse(all_N10(a),all_delN10(a),all_N26(a),all_delN26(a),1);
    ellipse(N10_norm(a),delN10_norm(a),N26_norm(a),delN26_norm(a),1,'r');
    hold on;
end

semilogx(N10_norm,r2610,'o','MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','MarkerSize',8)
hold on;
if strcmp(bur_lines,'y')
    for i = 2:size(decayed_r,1)
    semilogx(decayed_be(i,:), decayed_r(i,:),'--','Color',[0.3 0.3 0.3]);
    hold on;
    end
else
semilogx(decayed_be(1,:), decayed_r(1,:),'+k');
end

gca;
% axis square
grid on;
%Do some labeling
xlabel('[10Be*] (yr)'); ylabel('[26Al*]/[10Be*]');
% axis([1e3 1.5*max(N10_norm) 0 1.25.*max(r2610)]);
axis([1e3 1.25*max(be1) 0 1.25.*max(a1./be1)]);
% axis([0 1.5*max(N10_norm) 0 1.25.*max(r2610)]);
text(N10_norm.*(1+(0.05.*max(N10_norm)./N10_norm)),r2610,all_sample_name);
% text(N10_norm.*1.05,r2610,all_sample_name);
ax = gca;
ax.FontSize = 12; ax.FontName = 'Arial'; ax.FontWeight = 'bold';

%% Isochron Plot
f2 = figure('PaperSize',[10 10],'Position',[550 450 600 600]);
plot(be1,a1,'-k','LineWidth',1)
hold on;
plot(be2,a2,'-','LineWidth',1,'Color','#0072BD')
hold on;
errorbar(N10_norm,N26_norm,delN26_norm,delN26_norm,delN10_norm,delN10_norm,'o',...
    'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','MarkerSize',8, 'LineWidth',1,...
    'Color','#A2142F');

hold on;
% plot(be1,a1,'-k','LineWidth',1);
grid on;
% axis([0 max(be1) 0 max(a1)]);
axis([0 1.25.*max(N10_norm) 0 1.25.*max(N26_norm)]);
% axis([0 2.2e6 0 1.2e6]);
xlabel('[10Be*] (yr)');
ylabel('[26Al*] (yr)');

plot(decayed_be,decayed_a,'-.','Color',[0.3 0.3 0.3])

% semilogx(decayed_be,decayed_a,'k-.');
if strcmp(bur_lines,'y')
    %Do some calculations for burial lines
    num_a =  a1;
    num_be = be1;
    %Create burial time vector
    t_bur = (0:5e5:6e6);

    decayed_a = zeros(length(t_bur),length(num_a));
    decayed_be =  decayed_a;
    for i = 1:length(num_a)
        decayed_a(:,i) = num_a(i).*exp(-consts.l26.*t_bur);
        decayed_be(:,i) = num_be(i).*exp(-consts.l10.*t_bur);
    end

    for i = 2:length(t_bur)
        hold on;
        plot(decayed_be(i,:),decayed_a(i,:),'--','Color',[0.3 0.3 0.3])
    end
else
end

% multiples = [2 3];
% for i=1:length(multiples)
%     hold on;
%     plot(be1,a1.*multiples(i),':','LineWidth',1,'Color',[0.3 0.3 0.3]);
% end
text(N10_norm.*(1+(0.05.*max(N10_norm)./N10_norm)),N26_norm,all_sample_name);
ax = gca;
ax.FontSize = 12; ax.FontName = 'Arial'; ax.FontWeight = 'bold';

%% Dimensionless Banana

% bound1aDim = [(be1'.*consts.l10) ((a1.*consts.l26)./(be1.*consts.l10))']';
% bound2aDim = [(be2'.*consts.l10) ((a2.*consts.l26)./(be2.*consts.l10))']';
% 
% N10_norm = (all_N10.*consts.l10)./P10_tot';
% delN10_norm = N10_norm.*(all_delN10./all_N10);
% N26_norm = (all_N26.*consts.l26)./P26_tot';
% delN26_norm = N26_norm.*(all_delN26./all_N26);
% 
% %Calculate Al/Be Ratio
% % r2610 = all_N26./all_N10;
% r2610 = (N26_norm)./(N10_norm);
% % drdN26 = 1./all_N10;
% % drdN10 = -all_N26./(all_N10.^2);
% % delR = sqrt( (all_delN10.*drdN10).^2 + (all_delN26.*drdN26).^2 );
% % delR_norm = (delR./r2610).*r2610_norm;
% 
% drdN26 = 1./N10_norm;
% drdN10 = -N26_norm./(N10_norm.^2);
% delR = sqrt( (delN10_norm.*drdN10).^2 + (delN26_norm.*drdN26).^2 );
% 
% %Do some calculations for burial lines
% % ind = [2:20:(length(a1)/2) (length(a1)/2)+20:100:length(a1)];
% % ind = [2 4 11 26 51 121 201 401 2001];
% ind = 5:5:length(a1);
% num_a =  a1(ind).*consts.l26;
% num_be = be1(ind).*consts.l10;
% %Create burial time vector
% t_bur = (0:5e5:6e6);
% 
% decayed_a = zeros(length(t_bur),length(num_a));
% decayed_be =  decayed_a;
% for i = 1:length(num_a)
%     decayed_a(:,i) = num_a(i).*exp(-consts.l26.*t_bur);
%     decayed_be(:,i) = num_be(i).*exp(-consts.l10.*t_bur);
% end
% decayed_r = decayed_a./decayed_be;
% 
% f3 = figure('PaperSize',[10 10],'Position',[600 400 600 600]);
% steadyExp = plot(bound1aDim(1,:),bound1aDim(2,:),'k-','LineWidth',1);
% hold on;
% steadyEro = plot(bound2aDim(1,:),bound2aDim(2,:),'-','LineWidth',1,'Color','#0072BD');
% hold on;
% plot(decayed_be,decayed_r,'-.','Color',[0.3 0.3 0.3]);
% for a = 1:length(all_lat)
%     gca;
% %     ellipse(all_N10(a),all_delN10(a),all_N26(a),all_delN26(a),1);
%     ellipse(N10_norm(a),delN10_norm(a),N26_norm(a),delN26_norm(a),1,'r');
%     hold on;
% end
% 
% plot(N10_norm,r2610,'o','MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','MarkerSize',8)
% hold on;
% if strcmp(bur_lines,'y')
%     for i = 2:size(decayed_r,1)
%     plot(decayed_be(i,:), decayed_r(i,:),'--','Color',[0.3 0.3 0.3]);
%     hold on;
%     end
% else
% plot(decayed_be(1,:), decayed_r(1,:),'+k');
% end
% 
% gca;
% % axis square
% grid on;
% %Do some labeling
% xlabel('[10Be*]'); ylabel('[26Al*]/[10Be*]');
% % axis([1e3 1.5*max(N10_norm) 0 1.25.*max(r2610)]);
% % axis([1e3 1.25*max(be1) 0 1.25.*max(a1./be1)]);
% % axis([0 1.5*max(N10_norm) 0 1.25.*max(r2610)]);
% text(N10_norm.*(1+(0.05.*max(N10_norm)./N10_norm)),r2610,all_sample_name);
% % text(N10_norm.*1.05,r2610,all_sample_name);
% ax = gca;
% ax.FontSize = 12; ax.FontName = 'Arial'; ax.FontWeight = 'bold';





out.ID = all_sample_name;
out.allN10 = all_N10;
out.allN26 = all_N26;
out.r2610 = r2610;
out.delR = delR;
out.P26 = P26;
out.P10 = P10;
out.P14_tot = P26_tot;
out.P10_tot = P10_tot;
out.P26_mu = P26_mu;
out.P10_mu = P10_mu;
out.N10_norm = N10_norm;
out.N26_norm = N26_norm;
out.delN10_norm = delN10_norm;
out.delN26_norm = delN26_norm;
out.bound1a = bound1a;
out.bound2a = bound2a;
out.be = be1;
out.a = a1;
out.t = tempt;
out.decayed_a = decayed_a;
out.decayed_be = decayed_be;
out.decayed_r = decayed_r;
