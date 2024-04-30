function gm=rockfall_forward_v3(up,model,CNprop,Prockfall)
% inputs: 
% up(1): RWED is Rock Wall Exposure Duration(time in kyr), 
% up(2): Tm is the time of rockfall/time of subaerial exposure on moraine,
% up(3): Drw is Rock Wall Depth (in m or cm),
% model = model setup (generated from rockfallMCvJ1.m), 
% CNprop = structure with halflives and decay constants of radioactive 
%       nuclides (from getCNprop.m)

% output
% gm: cosmogenic nuclide concentrations [Be-10, Al-26, C-14] at end of
% simulation for each sample

% This function will:
% 1. calculate production profile in headwall given exposure time (RWED) and 
% assumed elevation/steepness
% 2. interpolate nuclide concentrations at depth of sample (Drw)
% 3. add nuclides produced since deposition on moraine - assuming stable 
% elevation (Tm)
% Note: Muogenic production is not taken into account at present

% Written by Jane Lund Andersen, Sep/22 - Feb/23
% _v2: adds elevation of rockfall as a parameter
% _v3: add muon production

%% Unpack generic parameters
RWED = up(1)*1e6; %yr
Tm = up(2)*1e6; %a
Lspal = model.data{1}.production.Lspal; %same for all samples

%% Calculate production profile in rockwall assuming elevation/steepness
% lat=-74.3; lon=-9.86; 
elv=round(up(3)/5)*5; %Elevation of rockfall rounded to nearest 5 (m)
shield=0.72; %Topographic shielding of rockwall surface assuming dip of 70 degrees
% p = ERA40atm(lat,lon,elv); % atm pressure (hPa) at assumed pos. of rockfall
zs=0:.1:20; %depth vector from 0 to 20 m depth in rockwall in 0.1 m increments
rho = model.data{1}.density; %density assumed same for all samples

I=find(elv==Prockfall{1}.production.elevs);
P10spal=Prockfall{1}.production.P10spal(I)*shield;
P10_m1=Prockfall{1}.production.P10_m1(I)*shield;
P10_m2=Prockfall{1}.production.P10_m2(I)*shield;
P10_Lm1=Prockfall{1}.production.P10_Lm1(I)*shield;
P10_Lm2=Prockfall{1}.production.P10_Lm2(I)*shield;

P10z=P10spal*exp(-rho.*zs*100/Lspal)+P10_m1*exp(-rho.*zs*100/P10_Lm1)+ ...
    P10_m2*exp(-rho.*zs*100/P10_Lm2); % *100 changes depth in m to cm

P14spal=Prockfall{1}.production.P14spal(I)*shield;
P14_m1=Prockfall{1}.production.P14_m1(I)*shield;
P14_m2=Prockfall{1}.production.P14_m2(I)*shield;
P14_Lm1=Prockfall{1}.production.P14_Lm1(I)*shield;
P14_Lm2=Prockfall{1}.production.P14_Lm2(I)*shield;

P14z=P14spal*exp(-rho.*zs*100/Lspal)+P14_m1*exp(-rho.*zs*100/P14_Lm1)+ ...
    P14_m2*exp(-rho.*zs*100/P14_Lm2); % *100 changes depth in m to cm

P26spal=Prockfall{1}.production.P26spal(I)*shield;
P26_m1=Prockfall{1}.production.P26_m1(I)*shield;
P26_m2=Prockfall{1}.production.P26_m2(I)*shield;
P26_Lm1=Prockfall{1}.production.P26_Lm1(I)*shield;
P26_Lm2=Prockfall{1}.production.P26_Lm2(I)*shield;

P26z=P26spal*exp(-rho.*zs*100/Lspal)+P26_m1*exp(-rho.*zs*100/P26_Lm1)+ ...
    P26_m2*exp(-rho.*zs*100/P26_Lm2); % *100 changes depth in m to cm

% P10rws = 4.01.*stone2000(lat,p,1)*shield; %Be-10
% P10z=P10rws*exp(-rho.*zs*100/Lspal); % *100 changes depth in m to cm
% P26rws = 27.93.*stone2000(lat,p,1)*shield; %Al-26
% P26z=P26rws*exp(-rho.*zs*100/Lspal);
% P14rws = 15.1.*stone2000(lat,p,1)*shield; %C-14 (NB: NL estimate is 13.5 rather than 15.1, for LSDn though)
% % P14rws = 13.5.*stone2000(lat,p,1)*shield;
% P14z=P14rws*exp(-rho.*zs*100/Lspal);

%% Calculate and output nuclide concentrations for each input sample
for i=1:model.Nsnr %loop samples

    n0 = (i-1)*model.Nsmp+model.Mmp; %parameter number start

    % Unpack sample-dependent parameters
    Drw = up(n0+1); %m

    % Initialize
    N10 = 0.0;
    N26 = 0.0;
    N14 = 0.0;

    % Calculate nuclide concentrations at depth of sample (Drw) in
    % rockwall
    P10d=interp1(zs,P10z,Drw); %interpolate production at depth
    N10=P10d/CNprop.lambda_Be*(1.0-exp(-RWED*CNprop.lambda_Be));
    P26d=interp1(zs,P26z,Drw); %interpolate production at depth
    N26=P26d/CNprop.lambda_Al*(1.0-exp(-RWED*CNprop.lambda_Al));
    P14d=interp1(zs,P14z,Drw); %interpolate production at depth
    N14=P14d/CNprop.lambda_C*(1.0-exp(-RWED*CNprop.lambda_C));

    % Add nuclides produced since deposition on moraine (Tm) - assuming 
    % stable elevation since deposition
    P14s = model.data{i}.production.P14spal; %at/g/yr
    P14_m1 = model.data{i}.production.P14_m1; %at/g/yr
    P14_m2 = model.data{i}.production.P14_m2; %at/g/yr

    P10s = model.data{i}.production.P10spal; %at/g/yr
    P10_m1 = model.data{i}.production.P10_m1; %at/g/yr
    P10_m2 = model.data{i}.production.P10_m2; %at/g/yr

    P26s = model.data{i}.production.P26spal; %at/g/yr
    P26_m1 = model.data{i}.production.P26_m1; %at/g/yr
    P26_m2 = model.data{i}.production.P26_m2; %at/g/yr

    N14 = N14*exp(-Tm*CNprop.lambda_C) + ... %decay during timestep
        (P14s/CNprop.lambda_C).*(1-exp(-CNprop.lambda_C.*Tm)) + ... %spallation production
        (P14_m1/CNprop.lambda_C).*(1-exp(-CNprop.lambda_C.*Tm)) + ... %muon pathway #1
        (P14_m2/CNprop.lambda_C).*(1-exp(-CNprop.lambda_C.*Tm)); %muon pathway #2
    N10 = N10*exp(-Tm*CNprop.lambda_Be) + ... %decay during timestep
        (P10s/CNprop.lambda_Be).*(1-exp(-CNprop.lambda_Be.*Tm))+ ... %spallation production
        (P10_m1/CNprop.lambda_Be).*(1-exp(-CNprop.lambda_Be.*Tm)) + ... %muon pathway #1
        (P10_m2/CNprop.lambda_Be).*(1-exp(-CNprop.lambda_Be.*Tm)); %muon pathway #2
    N26 = N26*exp(-Tm*CNprop.lambda_Al) + ... %decay during timestep
        (P26s/CNprop.lambda_Al).*(1-exp(-CNprop.lambda_Al.*Tm))+ ... %spallation production
        (P26_m1/CNprop.lambda_Al).*(1-exp(-CNprop.lambda_Al.*Tm)) + ... %muon pathway #1
        (P26_m2/CNprop.lambda_Al).*(1-exp(-CNprop.lambda_Al.*Tm)); %muon pathway #2

% ouput modelled data    
    for j=1:model.data{i}.Nnuc %loop nuclides
        nuclide=model.data{i}.nuclides(j); %get nuclide identification
            if nuclide == 1 %10Be
                gm((i-1)*model.Nds+j) = N10(1);
            elseif nuclide == 2 %26Al
                gm((i-1)*model.Nds+j) = N26(1);
            elseif nuclide == 3 %14C
                gm((i-1)*model.Nds+j) = N14(1);
            end
    end
end