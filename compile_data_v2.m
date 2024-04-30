function compile_data_v2()

% This code reads sample data from an excelfile, calculates production
% parameters using 'st' scaling, and saves a data structure for use in
% bedrockMC codes.

% _v2: spallation: LSDn implementation with codes from Nat Lifton
% 
% Input excelfile format:
% Each sample is presented with a line for each nuclide measured. 
% The first line contains sample site information (Sample name and number,
% lat (dd), lon (dd), elevation (m), topographic shielding factor (0-1), 
% sample thickness (cm) and density (g/cm3), depth to top of sample (cm), 
% year of sampling (yr), number of nuclides measured (this determines how 
% many lines of excel-file to read before next sample occurs). 
% Each line has a nuclide identificator (1=10Be, 2=26Al, 3=14C), 
% the cosmogenic nuclide concentration in at/g and related uncertainty.

% This script cites the following literature:
% Balco 2017: Production rate calculations for cosmic-ray-muon-produced 
% 10Be and 26Al benchmarked against geological calibration data. Quaternary
% Geochronology 39 (150-173).
% Balco & Shuster 2009: Balco, G. and Shuster, D.L., 2009. Production rate
% of cosmogenic 21Ne in quartz estimated from 10Be, 26Al, and 21Ne 
% concentrations in slowly eroding Antarctic bedrock surfaces. Earth and 
% Planetary Science Letters, 281(1-2), pp.48-58.
% Borchers et al., 2016: Borchers, B., Marrero, S., Balco, G., Caffee, M., 
% Goehring, B., Lifton, N., Nishiizumi, K., Phillips, F., Schaefer, J. and 
% Stone, J., 2016. Geological calibration of spallation production rates in 
% the CRONUS-Earth project. Quaternary Geochronology, 31, pp.188-198.
% Fernandez-Mosquera 2010: ?
% Marrero, S. M., Phillips, F. M., Borchers, B., Lifton, N., Aumer, R., & 
% Balco, G. (2016). Cosmogenic nuclide systematics and the CRONUScalc 
% program. Quaternary Geochronology, 31, 160-187.

% Codes used for calculation of production rates retrieved from:
% https://bitbucket.org/cronusearth/cronus-calc/src/master/ in july 2021
% Website: https://cronus.cosmogenicnuclides.rocks/2.1/
% See Marrero et al., 2016 for details

% Written by Jane Lund Andersen, Aarhus University, 2019-2022

close all; clear

% Add paths to functions called in this script
addpath('Functions')

%% Read data info from Excelfile
excelfile = 'data/milorgmoraine_input.xlsx';
[num,text,~] = xlsread(excelfile);

ns = max(num(:,1)); %number of samples/models in sheet

ij=1; %row number of first sample
for i=1:ns %loop through all samples
    sample{i}.Nnuc = num(ij,10); %Number of nuclides measured for sample i
    sample{i}.nuclides = num(ij:ij+sample{i}.Nnuc-1,11); %nuclide identifications for sample i
    
    %Sample specific data
    sample{i}.name = text(ij+1,1);
    sample{i}.lat=num(ij,2);
    sample{i}.lon=num(ij,3);
    sample{i}.elev=num(ij,4);
    % Calculate atm pressure (hPa) from elevation
        p = ERA40atm(sample{i}.lat,sample{i}.lon,sample{i}.elev);    
        sample{i}.pressure = p; 
    sample{i}.shield=num(ij,5); %Topographic shielding
    sample{i}.thick=num(ij,6); %Sample thickness (cm)
    sample{i}.density = num(ij,7); %(g/cm3)
    rho = sample{i}.density; %Density of rock sample (g/cm3) for internal use in function
    sample{i}.depth=num(ij,8); %surface samples = 0
    sample{i}.sampleyr = num(ij,9);
    
%     %% Site-specific production parameters
%     
%     % Define depths below surface z/rho [cm/(g/cm3)]=[g/cm^2] for fitting of production profiles
%     D_m = 100; %Max depth (m), changed below
%     z_m = linspace(0,10,100);
%     z_D = D_m*z_m.^3/10*rho; %denser depth-grid near surface
%     
%     % Spallation attenuation length calculated from CronusCalc functions 
%     % based on site cutoff rigidity, not considering terrain shielding [g/cm2]
%     % Lspal=attenuationlength(sample{i}.lat,sample{i}.lon,sample{i}.elev,p);
%     Lspal = 155; %Alternative, constant value [g/cm2]
%     sample{i}.production.Lspal=Lspal; %[g/cm2]
%     
%     % Thickness correction
%     sf_spal = exp(-sample{i}.thick/2*rho/Lspal); %Factor to correct production 
%     % rate for thickness of sample, sets surface production = production midway 
%     % through sample. Make sure this is not already factored in to site-specific 
%     % production rates
%     
%     maxZ = 1200; %maxdepth (g/cm2) for muon-production profile used for fitting
%     % of exponentials below. If this depth is very large, exponential terms will
%     % be dominated by fast muon production, which isn't ideal. 1200 g/cm2=4.5m
%     % with rho ~2.65-2.7, Test effect of this choice
% 
% 
%     nuclide specific data
        for j=1:sample{i}.Nnuc %loop over nuclides measured for sample i
            nuclide=num(ij+j-1,11); %get nuclide identification
            switch nuclide %calculate relevant production parameters only
                case 1 %10Be
                sample{i}.N10 = num(ij+j-1,12);
                sample{i}.dN10 = num(ij+j-1,13);
%                 
%                 % Spallation surface production in atoms / g qtz / year
%                 % 4.01 = Be10-spallation at SLHL from Borchers et al.,
%                 % 2016, Table 7, for the st scaling framework
%                 sample{i}.production.P10spal = 4.01.*stone2000(sample{i}.lat,p,1); 
%                 sample{i}.production.P10spal = sample{i}.production.P10spal*sf_spal*sample{i}.shield; %sample thickness and topographic shielding correction
%     
%                 %Muon production following Balco 2017, 
%                 %Muon production parameters, from BCO fit, Model 1A, alpha=1; f_star*f_C*f_D
%                 mc10.k_neg = 0.00191 .* 0.704 .* 0.1828; %summary cross-section for negative muon capture (at/muon)
%                 mc10.sigma0 = 0.280e-30; %x-section for fast muon production at 1 Gev
%                 mc10.Natoms = 2.006e22; %Oxygen atoms pr gram Quartz
%                 
%                 % Fit muon production profile calculated with P_mu_total_alpha1.m 
%                 % with two exponential terms following Balco 2017, Eq. 7/Fig. 15
%                 p10_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,2,mc10,0);
%                 sample{i}.production.P10_Lm1=p10_muons.L(1); %attenuation, first exponential term
%                 sample{i}.production.P10_Lm2=p10_muons.L(2); %attenuation, second exponential term
%                 sample{i}.production.P10_m1 = p10_muons.P(1); %production, first exponential term
%                 sample{i}.production.P10_m2 = p10_muons.P(2); %production, second exponential term
%                 shield_fac10_m1 = exp(-sample{i}.thick/2*rho/p10_muons.L(1)); %sample thickness correction
%                 shield_fac10_m2 = exp(-sample{i}.thick/2*rho/p10_muons.L(2)); %sample thickness correction
%                 sample{i}.production.P10_m1 = sample{i}.production.P10_m1*shield_fac10_m1*sample{i}.shield; % production, first exponential term, corrected for sample thickness and topographic shielding
%                 sample{i}.production.P10_m2 = sample{i}.production.P10_m2*shield_fac10_m2*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding
%                 
                case 2 %26Al
                sample{i}.N26 = num(ij+j-1,12);
                sample{i}.dN26 = num(ij+j-1,13);
%                 
%                 % Spallation surface production in atoms / g qtz / year
%                 % 27.93 = Al26-spallation at SLHL from Borchers et al.,
%                 % 2016, Table 7, for the st scaling framework
%                 sample{i}.production.P26spal = 27.93.*stone2000(sample{i}.lat,p,1);
%                 sample{i}.production.P26spal = sample{i}.production.P26spal*sf_spal*sample{i}.shield; %sample thickness and topographic shielding correction
%                 
%                 %Muon production following Balco 2017, 
%                 %Muon production parameters, from BCO fit, Model 1A, alpha=1; f_star*f_C*f_D
%                 mc26.k_neg = 0.0133 .* 0.296 .* 0.6559; %summary cross-section for negative muon capture (at/muon)
%                 mc26.sigma0 = 3.89e-30; %x-section for fast muon production at 1 Gev
%                 mc26.Natoms = 1.003e22; %Si atoms pr gram Quartz
%                               
%                 % Fit muon production profile calculated with P_mu_total_alpha1.m 
%                 % with two exponential terms following Balco 2017, Eq. 7/Fig. 15
%                 p26_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,2,mc26,0);
%                 sample{i}.production.P26_Lm1=p26_muons.L(1); %attenuation, first exponential term
%                 sample{i}.production.P26_Lm2=p26_muons.L(2); %attenuation, second exponential term
%                 sample{i}.production.P26_m1 = p26_muons.P(1); %production, first exponential term
%                 sample{i}.production.P26_m2 = p26_muons.P(2); %production, second exponential term
%                 shield_fac26_m1 = exp(-sample{i}.thick/2*rho/p26_muons.L(1)); %sample thickness correction
%                 shield_fac26_m2 = exp(-sample{i}.thick/2*rho/p26_muons.L(2)); %sample thickness correction
%                 sample{i}.production.P26_m1 = sample{i}.production.P26_m1*shield_fac26_m1*sample{i}.shield; %production, first exponential term, corrected for sample thickness and topographic shielding
%                 sample{i}.production.P26_m2 = sample{i}.production.P26_m2*shield_fac26_m2*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding
%                 
                case 3 %14C
                sample{i}.N14 = num(ij+j-1,12);
                sample{i}.dN14 = num(ij+j-1,13);
% 
%                 % Spallation surface production in atoms / g qtz / year
%                 sample{i}.production.P14spal = 15.1*stone2000(sample{i}.lat,p,1);
%                 sample{i}.production.P14spal = sample{i}.production.P14spal*sf_spal*sample{i}.shield; %sample thickness and topographic shielding correction
%     
%                 mc14.Natoms = 2.006e22;
%                 mc14.k_neg = 0.116 .* 0.704 .* 0.1828; % From Leymon High fit
%                 mc14.sigma0 = 0.45e-27./190; % From Heisinger, alpha = 1
%                 
%                 % Fit muon production profile calculated with P_mu_total_alpha1.m 
%                 % with two exponential terms following Balco 2017, Eq. 7/Fig. 15
%                 p14_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,2,mc14,0);
%                 sample{i}.production.P14_Lm1=p14_muons.L(1); %attenuation, first exponential term
%                 sample{i}.production.P14_Lm2=p14_muons.L(2); %attenuation, second exponential term
%                 sample{i}.production.P14_m1 = p14_muons.P(1); %production, first exponential term
%                 sample{i}.production.P14_m2 = p14_muons.P(2); %production, second exponential term
%                 shield_fac14_m1 = exp(-sample{i}.thick/2*rho/p14_muons.L(1)); %sample thickness correction
%                 shield_fac14_m2 = exp(-sample{i}.thick/2*rho/p14_muons.L(2)); %sample thickness correction
%                 sample{i}.production.P14_m1 = sample{i}.production.P14_m1*shield_fac14_m1*sample{i}.shield; %production, first exponential term, corrected for sample thickness and topographic shielding
%                 sample{i}.production.P14_m2 = sample{i}.production.P14_m2*shield_fac14_m2*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding

%                % Add LSDn production created through GoBananas1410_jla.m and GoBananas1026_jla.m
               load LSprod.mat samples
               sample{i}.production=samples{i}.production;

            end
        end
        
        %% Calculate nuclide ratios
        if isfield(sample{i},'N26')
            sample{i}.r2610 = sample{i}.N26./sample{i}.N10;
            sample{i}.dr2610 = sample{i}.r2610.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN26./sample{i}.N26).^2);
        end       
        if isfield(sample{i},'N14')
            sample{i}.r1410 = sample{i}.N14./sample{i}.N10;
            sample{i}.dr1410 = sample{i}.r1410.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN14./sample{i}.N14).^2);
            sample{i}.r1426 = sample{i}.N14./sample{i}.N26;
            sample{i}.dr1426 = sample{i}.r1426.*sqrt((sample{i}.dN26./sample{i}.N26).^2+(sample{i}.dN14./sample{i}.N14).^2);
        end
        
    ij=ij+sample{i}.Nnuc; %Look for next sample in this row
end

save('./data/milorgmoraine_mn_data_v2.mat','sample','excelfile')