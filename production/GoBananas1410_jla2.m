
function samples = GoBananas1410_jla2()

% 'sm' is a string for the scaling model to use - 'St' for Lal/Stone, 'LS' for LSDn
% 'bur_lines' is a string - 'y' or 'n' to include or not
% Version 1/23 - NL
% _jla: set up to fit muon production profile with exponentials
% _jla2: tweaked to get production parameters as function of rockfall
% elevation

% file = input('Please enter the name of the sample data file: ','s');
% FID = fopen(file);
% data = textscan(FID,'%s %n %n %n %s %n %n %n %n %n %n %n %n %s');
% fclose(FID);
% dstring='';

%Make the sample structure.

% all_sample_name = data{1};
% all_lat = data{2}; 
% all_long = data{3};
% all_elv = data{4};
% all_pressure = data{4};
% all_aa = data{5};
% all_thick = data{6};
% all_rho = data{7};
% all_othercorr = data{8};
% all_E = data{9};
% all_N14 = data{10};
% all_delN14 = data{11}; 
% all_N10 = data{12};
% all_delN10 = data{13};
% all_be_std_name = data{14};

% all_curve = data{15};

%Load data and constants
% make_consts_LSD;
load consts_LSDn;

num_samples = 1; %length(all_lat);

P10_LS = 0; %zeros(size(all_lat));
P14_LS = 0; %zeros(size(all_lat));

%Calculate total SLHL P rates for 14C and 10Be
% if strcmp(sm,'St')
%     P14 = consts.P14_ref_St;
%     P10 = consts.P10_ref_St;
% elseif strcmp(sm,'LS')
    P14 = consts.P14_ref_LS;
    P10 = consts.P10_ref_LS;
% end

% % Check for nuclides
% for a = 1:num_samples	
% 	if all_N14(a) ~= 0
% 		all_isN14(a) = 1;
%     else
% 		all_isN14(a) = 0; 
%     end
% 
% 	if all_delN14(a) ~= 0
% 		all_isdelN14(a) = 1; 
%     else
% 		all_isdelN14(a) = 0; 
%     end
% 	
% 	if all_N10(a) ~= 0
% 		all_isN10(a) = 1;
%     else
% 		all_isN10(a) = 0;
%     end
% 
% 	if all_delN10(a) ~= 0
% 		all_isdelN10(a) = 1; 
%     else
% 		all_isdelN10(a) = 0; 
%     end
% 	
% 	% catch mismatches;
% 	
% 	if (~all_isN14(a) && ~all_isN10(a))
%    		error(['Need either C-14 or Be-10 concentration - line ' int2str(a)]);
% 	elseif (all_isN14(a) && ~all_isdelN14(a)) || (~all_isN14(a) && all_isdelN14(a))
%     		error(['Need both C-14 concentration and uncertainty - line ' int2str(a)]);
% 	elseif (all_isN10(a) && ~all_isdelN10(a)) || (~all_isN10(a) && all_isdelN10(a))
%     		error(['Need both Be-10 concentration and uncertainty - line ' int2str(a)]);
%     end
% end

for a = 1:num_samples
    sample.sample_name{a} = 'Rockfall'; %all_sample_name{a};
    sample.lat = -74.3; %all_lat(a);
    sample.long = -9.86; %all_long(a);
    sample.aa = 'std'; %all_aa{a};
%     if strcmp(all_aa{a},'std') || strcmp(all_aa{a},'ant')
			% store the elevation value
% 			sample.elv = 1450; %all_elv(a);
% 	elseif  strcmp(all_aa{a},'pre')
% 			% store the pressure value
% 			sample.pressure = all_pressure(a);
%     end
    sample.thick = 0; %all_thick(a);
    sample.E = 0; %all_E(a);
    sample.rho = 2.65; %all_rho(a);
    sample.othercorr = 1; %all_othercorr(a); %topographic shielding
	sample.N14 = 4e5; %all_N14(a);
	sample.delN14 = 4e3; %all_delN14(a);
    
%     if all_isN10(a)
%         if (strcmp(all_be_std_name(a),'07KNSTD'))
            sample.N10 = 1.3e6.*consts.be_stds_cfs(1); %all_N10(a).*consts.be_stds_cfs(1);
            sample.delN10 = 2.4e4.*consts.be_stds_cfs(1); %all_delN10(a).*consts.be_stds_cfs(1);
%         elseif (strcmp(all_be_std_name(a),'KNSTD'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(2);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(2);
%         elseif (strcmp(all_be_std_name(a),'NIST_Certified'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(3);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(3);
%         elseif (strcmp(all_be_std_name(a),'LLNL31000'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(4);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(4);
%         elseif (strcmp(all_be_std_name(a),'LLNL10000'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(5);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(5);
%         elseif (strcmp(all_be_std_name(a),'LLNL3000'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(6);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(6);
%         elseif (strcmp(all_be_std_name(a),'LLNL1000'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(7);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(7);
%         elseif (strcmp(all_be_std_name(a),'LLNL300'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(8);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(8);
%         elseif (strcmp(all_be_std_name(a),'NIST_30000'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(9);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(9);
%         elseif (strcmp(all_be_std_name(a),'NIST_30200'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(10);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(10);
%         elseif (strcmp(all_be_std_name(a),'NIST_30300'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(11);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(11);
%         elseif (strcmp(all_be_std_name(a),'NIST_30600'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(12);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(12);
%         elseif (strcmp(all_be_std_name(a),'NIST_27900'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(13);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(13);
%         elseif (strcmp(all_be_std_name(a),'S555'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(14);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(14);
%         elseif (strcmp(all_be_std_name(a),'S2007'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(15);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(15);
%         elseif (strcmp(all_be_std_name(a),'BEST433'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(16);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(16);
%         elseif (strcmp(all_be_std_name(a),'BEST433N'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(17);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(17);
%         elseif (strcmp(all_be_std_name(a),'S555N'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(18);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(18);
%         elseif (strcmp(all_be_std_name(a),'S2007N'))
%             sample.N10 = all_N10(a).*consts.be_stds_cfs(19);
%             sample.delN10 = all_delN10(a).*consts.be_stds_cfs(19);
%         end
%     end
%     % 1a. Thickness scaling factor. 

    if sample.thick > 0
        sample.thickSF = thickness(sample.thick,consts.Lsp,sample.rho);
    else 
        sample.thickSF = 1;
    end

% 1b. If no pressure entered yet, create it from the elevation
    % using the appropriate atmosphere
    % Change in version 2: This looks up the surface pressure from 
    % NCAR map to use in the standard atmosphere equation.
    % Obviously, it is always better to estimate pressure from local
    % station data. 

%     if (~isfield(sample,'pressure'))
%         if (strcmp(sample.aa,'std'))
            % Old code
            % sample.pressure = stdatm(sample.elv);
            % New code
            elevations = 1450:5:1950;
            sample.production.elevs=elevations;
            for i=1:length(elevations)
                sample.elv=elevations(i);
                sample.pressure = ERA40atm(sample.lat,sample.long,sample.elv);
%         elseif (strcmp(sample.aa,'ant'))
%             sample.pressure = antatm(sample.elv);
%         end
%     end

% Catch confusion with pressure submission. If sample.pressure is already 
% set, it should have a submitted value. If zero, something is wrong. 
% This should never happen in online use. 

    if sample.pressure == 0
        error(['Sample.pressure = 0 on sample ' sample.sample_name]);
    end
    % 
    sample.shielding = 1;
    
%     disp(sample.sample_name(a));

    scaling14 = ScalingLSDLal(sample,consts,14);
    scaling10 = ScalingLSDLal(sample,consts,10);

    c_results = get_age_LSDLal(sample,consts,14,scaling14,a);
    be_results = get_age_LSDLal(sample,consts,10,scaling10,a);

%     if strcmp(sm,'St')  
%         P10_St(a) = scaling10.SF_St(1).*P10;
%         P14_St(a) = scaling14.SF_St(1).*P14;
%     elseif strcmp(sm,'LS')  
        t_app = scaling10.tv(1:find(scaling10.tv>be_results.t_LS,1,'first'));
        weights = exp(-t_app.*consts.l10);
        P10_LS(a) = P10.*sum(scaling10.SF_LS(1:find(scaling10.tv>be_results.t_LS,1,'first')).*weights)./sum(weights);
        t_app = scaling14.tv(1:find(scaling14.tv>c_results.t_LS,1,'first'));
        weights = exp(-t_app.*consts.l14);
        P14_LS(a) = P14.*sum(scaling14.SF_LS(1:find(scaling14.tv>c_results.t_LS,1,'first')).*weights)./sum(weights); %time averaged and decay weighted production
%     end
    
    % JLA added the following four lines
    sample.production.Lspal=150; %spallation attenuation length
    sf_spal = exp(-sample.thick.*sample.rho./2./sample.production.Lspal); %thickness shielding factor
    sample.production.P14spal(i)=P14_LS(a)*sf_spal*sample.othercorr; %production corrected for sample thickness and topographic shielding
    sample.production.P10spal(i)=P10_LS(a)*sf_spal*sample.othercorr; %production corrected for sample thickness and topographic shielding

%     if strcmp(sm,'LS') 
        GDRc = [-448.004 1189.18 -1152.15 522.061 -103.241 6.89901 0];
        RcEst = polyval(GDRc,cos(d2r(sample.lat)));
        
        if (RcEst < 0) RcEst = 0; end
%     end
    
    %Get muon production rates for 14C and 10Be
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
    
%     if strcmp(sm,'St')  
%         muSt = P_mu_total_alpha1((sample.thick.*sample.rho./2),sample.pressure,mconsts,'yes');
%         P10_mu(a) = muSt.P_fast_St + muSt.P_neg_St;
%     elseif strcmp(sm,'LS')  
%         muLS = P_mu_totalLSD((sample.thick.*sample.rho./2),sample.pressure,RcEst,consts.SPhiInf,mconsts,'yes');
%         muLS = P_mu_total_alpha1((sample.thick.*sample.rho./2),sample.pressure,mconsts,'yes');
%         P10_mu(a) = muLS.P_fast_LS + muLS.P_neg_LS;
%         P10_mu(a) = muSt.P_fast_St + muSt.P_neg_St;

        % JLA added the following lines
        maxZ = 1200; %1200 g/cm2=4.5m with rho ~2.65-2.7
        p10_muons = fit_Pmu_with_exp_LSD_jla(sample.pressure,0,maxZ,2,mconsts,RcEst,consts.SPhiInf,0);
            sample.production.P10_Lm1(i)=p10_muons.L(1); %attenuation, first exponential term
            sample.production.P10_Lm2(i)=p10_muons.L(2); %attenuation, second exponential term
            sample.production.P10_m1(i) = p10_muons.P(1); %production, first exponential term
            sample.production.P10_m2(i) = p10_muons.P(2); %production, second exponential term
            shield_fac10_m1 = exp(-sample.thick.*sample.rho./2./p10_muons.L(1)); %sample thickness correction
            shield_fac10_m2 = exp(-sample.thick.*sample.rho./2./p10_muons.L(2)); %sample thickness correction
            sample.production.P10_m1(i) = sample.production.P10_m1(i)*shield_fac10_m1*sample.othercorr; % production, first exponential term, corrected for sample thickness and topographic shielding
            sample.production.P10_m2(i) = sample.production.P10_m2(i)*shield_fac10_m2*sample.othercorr; %production, second exponential term, corrected for sample thickness and topographic shielding
         
%    end
    %

    mconsts.Natoms = consts.NatomsQtzO;
    mconsts.sigma0.LS = consts.sigma014.LS;
    mconsts.sigma0.Lm = consts.sigma014.Lm;
    mconsts.sigma0.St = consts.sigma014.St;
    mconsts.fstar.LS = consts.fstar14.LS;
    mconsts.fstar.Lm = consts.fstar14.Lm;
    mconsts.fstar.St = consts.fstar14.St;
%     mconsts.delsigma0 = consts.delsigma190_14; % not used
    mconsts.k_negpartial = consts.k_negpartial_14;
%     mconsts.delk_neg = consts.delk_neg14; % not used
    mconsts.mfluxRef = consts.mfluxRef;

%     if strcmp(sm,'St')  
%         muSt = P_mu_total_alpha1((sample.thick.*sample.rho./2),sample.pressure,mconsts,'yes');
%         P14_mu(a) = muSt.P_fast_St + muSt.P_neg_St;
%     elseif strcmp(sm,'LS')  
%         muLS = P_mu_totalLSD((sample.thick.*sample.rho./2),sample.pressure,RcEst,consts.SPhiInf,mconsts,'yes');
% %         muLS = P_mu_total_alpha1((sample.thick.*sample.rho./2),sample.pressure,mconsts,'yes');
%         P14_mu(a) = muLS.P_fast_LS + muLS.P_neg_LS;
% %         P14_mu(a) = muLS.P_fast_St + muLS.P_neg_St;

        % JLA added the following lines
        maxZ = 1200; %1200 g/cm2=4.5m with rho ~2.65-2.7
        p14_muons = fit_Pmu_with_exp_LSD_jla(sample.pressure,0,maxZ,2,mconsts,RcEst,consts.SPhiInf,0);
            sample.production.P14_Lm1(i)=p14_muons.L(1); %attenuation, first exponential term
            sample.production.P14_Lm2(i)=p14_muons.L(2); %attenuation, second exponential term
            sample.production.P14_m1(i) = p14_muons.P(1); %production, first exponential term
            sample.production.P14_m2(i) = p14_muons.P(2); %production, second exponential term
            shield_fac14_m1 = exp(-sample.thick.*sample.rho./2./p14_muons.L(1)); %sample thickness correction
            shield_fac14_m2 = exp(-sample.thick.*sample.rho./2./p14_muons.L(2)); %sample thickness correction
            sample.production.P14_m1(i) = sample.production.P14_m1(i)*shield_fac14_m1*sample.othercorr; % production, first exponential term, corrected for sample thickness and topographic shielding
            sample.production.P14_m2(i) = sample.production.P14_m2(i)*shield_fac14_m2*sample.othercorr; %production, second exponential term, corrected for sample thickness and topographic shielding

%     end
    samples{a}=sample; %jla
            end
end

% %Calculate total P rates for 14C and 10Be
% if strcmp(sm,'St')  
%     P14_tot = (P14_St) + P14_mu;
%     P10_tot = (P10_St) + P10_mu;
%     P10_spall_frac = mean(P10_St./P10_tot);
%     P14_spall_frac = mean(P14_St./P14_tot);
% elseif strcmp(sm,'LS')
%     P14_tot = (P14_LS)' + P14_mu;
%     P10_tot = (P10_LS)' + P10_mu;
%     P10_spall_frac = mean(P10_LS'./P10_tot);
%     P14_spall_frac = mean(P14_LS'./P14_tot);
% end  
% 
% % create data for the simple-exposure line and the simple-erosion line. 
% tempt = logspace(2,5.3,100);
% be1 = (P10_spall_frac/consts.l10)*(1-exp(-consts.l10*tempt))+((1-P10_spall_frac)/consts.l10)*(1-exp(-consts.l10*tempt));
% c1 = (P14_spall_frac/consts.l14)*(1-exp(-consts.l14*tempt))+((1-P14_spall_frac)/consts.l14)*(1-exp(-consts.l14*tempt));
% bound1a = [be1' (c1./be1)']';
% tempe = logspace(-9,1,100);
% be2 = P10_spall_frac./(consts.l10 + tempe./160)+(1-P10_spall_frac)./(consts.l10 + tempe./1500);
% c2 = P14_spall_frac./(consts.l14 + tempe./160)+(1-P14_spall_frac)./(consts.l14 + tempe./1500);
% bound2a = [be2' (c2./be2)']';
% 
% N10_norm = all_N10./P10_tot';
% delN10_norm = N10_norm.*(all_delN10./all_N10);
% N14_norm = all_N14./P14_tot';
% delN14_norm = N14_norm.*(all_delN14./all_N14);
% 
% %Calculate C/Be Ratio
% % r1410 = all_N14./all_N10;
% r1410 = N14_norm./N10_norm;
% drdN14 = 1./N10_norm;
% drdN10 = -N14_norm./(N10_norm.^2);
% delR = sqrt( (delN10_norm.*drdN10).^2 + (delN14_norm.*drdN14).^2 );
% % delR_norm = (delR./r1410).*r1410_norm;
% 
% %Do some calculations for burial lines
% % ind = [2:20:(length(c1)/2) (length(c1)/2)+20:100:length(c1)];
% ind = 5:5:length(c1);
% num_c =  c1(ind);
% num_be = be1(ind);
% %Create burial time vector
% t_bur = (0:2e3:2e4);
% 
% decayed_c = zeros(length(t_bur),length(num_c));
% decayed_be = decayed_c;
% for i = 1:length(num_c)
%     decayed_c(:,i) = num_c(i).*exp(-consts.l14.*t_bur);
%     decayed_be(:,i) = num_be(i).*exp(-consts.l10.*t_bur);
% end
% decayed_r = decayed_c./decayed_be;
% 
% %% Banana Plot
% figure('PaperSize',[10 10]);
% steadyExp = semilogx(bound1a(1,:),bound1a(2,:),'k-','LineWidth',1);
% hold on;
% steadyEro = semilogx(bound2a(1,:),bound2a(2,:),'-','LineWidth',1,'Color','#0072BD');
% hold on;
% semilogx(decayed_be,decayed_r,'-.','Color',[0.3 0.3 0.3]);
% 
% for a = 1:length(all_lat)
%     gca;
% %     ellipse(all_N10(a),all_delN10(a),all_N14(a),all_delN14(a),1);
%     ellipse(N10_norm(a),delN10_norm(a),N14_norm(a),delN14_norm(a),1,'r');
%     hold on;
% end
% semilogx(N10_norm,r1410,'o','MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','MarkerSize',8)
% hold on;
% if strcmp(bur_lines,'y')
%     for i = 2:size(decayed_r,1)
%         semilogx(decayed_be(i,:), decayed_r(i,:),'--','Color',[0.3 0.3 0.3]);
%         hold on;
%     end
% else
%     semilogx(decayed_be(1,:), decayed_r(1,:),'+k');
%     hold on;
% end
% 
% gca;
% % axis square
% grid on;
% %Do some labeling
% xlabel('[10Be*] (yr)'); 
% ylabel('[14C*]/[10Be*]');
% % axis([1e3 1.25*max(be1) 0 1.25.*max(c1./be1)]);
% axis([1e3 1.25*max(be1) 0 1.25.*max(c1./be1)]);
% 
% text(N10_norm.*(1+(0.05.*max(N10_norm)./N10_norm)),r1410,all_sample_name);
% hold on;
% ax = gca;
% ax.FontSize = 12; ax.FontName = 'Arial'; ax.FontWeight = 'bold';
% 
% %% Isochron Plot
% figure('PaperSize',[10 10]);
% hold on;
% plot(be1,c1,'-k','LineWidth',1)
% hold on;
% plot(be2,c2,'-','LineWidth',1,'Color','#0072BD')
% hold on;
% errorbar(N10_norm,N14_norm,delN14_norm,delN14_norm,delN10_norm,delN10_norm,'o',...
%     'MarkerFaceColor','#A2142F','MarkerEdgeColor','#A2142F','MarkerSize',8, 'LineWidth',1,...
%     'Color','#A2142F');
% hold on;
% grid on;
% 
% plot(decayed_be,decayed_c,'-.','Color',[0.3 0.3 0.3])
% hold on;
% 
% axis([0 1.25*max(N10_norm) 0 1.25.*max(N14_norm)]);
% xlabel('[10Be*] (yr)');
% ylabel('[14C*] (yr)');
% % semilogx(decayed_be,decayed_c,'k-.');
% hold on;
% if strcmp(bur_lines,'y')
%     %Do some calculations for burial lines
%     num_c =  c1;
%     num_be = be1;
%     %Create burial time vector
% %     t_bur = (0:2e3:2e4);
% 
%     decayed_c = zeros(length(t_bur),length(num_c));
%     decayed_be =  decayed_c;
%     for i = 1:length(num_c)
%         decayed_c(:,i) = num_c(i).*exp(-consts.l14.*t_bur);
%         decayed_be(:,i) = num_be(i).*exp(-consts.l10.*t_bur);
%     end
% 
%     for i = 2:length(t_bur)
%         hold on;
%         plot(decayed_be(i,:),decayed_c(i,:),'--','Color',[0.3 0.3 0.3])
%     end
% else
% end
% 
% text(N10_norm.*(1+(0.05.*max(N10_norm)./N10_norm)),N14_norm,all_sample_name);
% 
% % multiples = [2 3];
% % for i=1:length(multiples)
% %     hold on;
% %     plot(be1,c1.*multiples(i),':','LineWidth',1.5,'Color',[0.3 0.3 0.3]);
% % end
% ax = gca;
% ax.FontSize = 12; ax.FontName = 'Arial'; ax.FontWeight = 'bold';
% 
% out.ID = all_sample_name;
% out.allN10 = all_N10;
% out.allN14 = all_N14;
% out.r1410 = r1410;
% out.delR = delR;
% out.P14 = P14;
% out.P10 = P10;
% out.P14_tot = P14_tot;
% out.P10_tot = P10_tot;
% out.P14_mu = P14_mu;
% out.P10_mu = P10_mu;
% out.N10_norm = N10_norm;
% out.N14_norm = N14_norm;
% out.delN10_norm = delN10_norm;
% out.delN14_norm = delN14_norm;
% out.bound1a = bound1a;
% out.bound2a = bound2a;
% out.be = be1;
% out.c = c1;
% out.t = tempt;
% out.decayed_c = decayed_c;
% out.decayed_be = decayed_be;
% out.decayed_r = decayed_r;
