function results = get_age_LSDLal(sample, consts, nuclide,geomag,i)

% This function calculates the exposure age of a sample and
% packages the results.
%
% syntax : results = get_c_be_age(sample,consts,nuclide);
%
% argument sample is the structure assembled upstream by al_be_age_one
% or al_be_age_many. Many required fields, see hard-copy docs for details.
%
% argument consts is typically the al_be_consts structure derived from
% make_al_be_consts_vxx.m. Many required fields, see hard-copy docs for
% details.
%
% argument nuclide is 14 or 10. Number not string. 
%
% Many dependencies.
% 
% results is a structure with many fields:
%
% Non-scaling-scheme-dependent results:
%
% results.flags: Error messages, mostly about saturation
% results.main_version: version of this function
% results.muon_version: version of P_mu_total called internally
% results.P_mu: surface production rate by muons (atoms/g/yr)
% results.thick_sf: thickness scaling factor (nondimensional)
% results.tv: time vector against which to plot Rc and P
%
% Scaling-scheme-dependent results: five of each of these fields, one for
% each scaling scheme. XX in field names below indicates a two-letter code
% identifying each scaling scheme, as follows:
% St (Stone,2000); Du (Dunai, 2001); De (Desilets, 2006);
% Li (Lifton, 2005); and Lt (paleomagnetically corrected implementation 
% of Lal(1991)/Stone(2000). 
%
% results.P_XX: P(t) at site in scaling scheme XX (atoms/g/yr) (vector)
% results.t_XX: Exposure age WRT scaling scheme XX (yr)
% results.FSF_St: Effective scaling factor WRT scaling scheme XX 
% (the effective scaling factor is the SF which, when put into the simple
% age equation, yields the correct age)
% results.delt_int_XX: Internal uncertainty WRT scaling scheme XX (yr)
% results.delt_ext_XX: External uncertainty WRT scaling scheme XX (yr)
%
% Also:
% results.SF_St_nominal: Stone(2000) scaling factor (historical interest)
% 
% See hard-copy documentation for discussion of the actual method. 
%
% Written by Greg Balco -- UW Cosmogenic Nuclide Lab and Berkeley
% Geochronology Center
% balcs@u.washington.edu
% 
% Modified by Brent Goehring -- Lamont-Doherty Earth Observatory
% goehring@ldeo.columbia.edu
% and Nat Lifton -- Purdue University
% nlifton@purdue.edu

% September, 2010
% 
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math
%
% Copyright 2010, University of Washington, Columbia University and the
% Purdue University
% All rights reserved
% Developed in part with funding from the National Science Foundation.

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).
% What version is this?

ver = '2.0';


% Activate secret production rate calibration feature. Not used in 
% the online calculators. 

calFlag = 0;
if isfield(sample,'truet')
    % It's a calibration sample and we want the time-integrated
    % SF, not the age. 
    calFlag = 1;
end

% 0. Select appropriate values for nuclide of interest

if nuclide == 10
    % Atoms/g measurement
    N = sample.N10; delN = sample.delN10; 
    % Production rates from spallation for six schemes
    P_ref_LS = consts.P10_ref_LS; delP_ref_LS = consts.delP10_ref_LS;    
    P_ref_LD = consts.P10_ref_LD; delP_ref_LD = consts.delP10_ref_LD;    
    P_ref_Lm = consts.P10_ref_Lm; delP_ref_Lm = consts.delP10_ref_Lm;
    P_ref_St = consts.P10_ref_St; delP_ref_St = consts.delP10_ref_St;
    % Decay constant
    l = consts.l10; 
    % Fraction spallation for Stone, 2000 scheme - historical interest
    Fsp = consts.Fsp10;
    % constants structure for muon production rate
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
elseif nuclide == 14
    % Atoms/g measurement
    N = sample.N14; delN = sample.delN14; 
    % Production rates for 6 schemes
    P_ref_LS = consts.P14_ref_LS; delP_ref_LS = consts.delP14_ref_LS;    
    P_ref_LD = consts.P14_ref_LD; delP_ref_LD = consts.delP14_ref_LD;    
    P_ref_Lm = consts.P14_ref_Lm; delP_ref_Lm = consts.delP14_ref_Lm;
    P_ref_St = consts.P14_ref_St; delP_ref_St = consts.delP14_ref_St;
    % Decay constant
    l = consts.l14; dell = consts.dell14; 
    % Fraction spallation for Stone, 2000 scheme - historical interest
    Fsp = consts.Fsp14;
    % constants structure for muon production rate
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
elseif nuclide == 26
    % Atoms/g measurement
    N = sample.N26; delN = sample.delN26; 
    % Production rates for 6 schemes
    P_ref_LS = consts.P26_ref_LS; delP_ref_LS = consts.delP26_ref_LS;    
    P_ref_LD = consts.P26_ref_LD; delP_ref_LD = consts.delP26_ref_LD;    
    P_ref_Lm = consts.P26_ref_Lm; delP_ref_Lm = consts.delP26_ref_Lm;
    P_ref_St = consts.P26_ref_St; delP_ref_St = consts.delP26_ref_St;
    % Decay constant
    l = consts.l26; dell = consts.dell26; 
    % Fraction spallation for Stone, 2000 scheme - historical interest
    Fsp = consts.Fsp26;
    % constants structure for muon production rate
    mconsts.Natoms = consts.NatomsQtzSi;
    mconsts.sigma0.LS = consts.sigma026.LS;
    mconsts.sigma0.Lm = consts.sigma026.Lm;
    mconsts.sigma0.St = consts.sigma026.St;
    mconsts.fstar.LS = consts.fstar26.LS;
    mconsts.fstar.Lm = consts.fstar26.Lm;
    mconsts.fstar.St = consts.fstar26.St;
%     mconsts.delsigma190 = consts.delsigma190_26; % not used
    mconsts.k_negpartial = consts.k_negpartial_26;
%     mconsts.delk_neg = consts.delk_neg26; % not used
    mconsts.mfluxRef = consts.mfluxRef;
elseif nuclide == 3
    % Atoms/g measurement
    N = sample.N3; delN = sample.delN3; 
    % Production rates from spallation for six schemes
    P_ref_LS = consts.P3_ref_LS; delP_ref_LS = consts.delP3_ref_LS;    
    P_ref_LD = consts.P3_ref_LD; delP_ref_LD = consts.delP3_ref_LD;    
    P_ref_Lm = consts.P3_ref_Lm; delP_ref_Lm = consts.delP3_ref_Lm;
    P_ref_St = consts.P3_ref_St; delP_ref_St = consts.delP3_ref_St;
    % Fraction spallation for Stone, 2000 scheme - historical interest
    Fsp = consts.Fsp3;
end

mconsts.SPhiInf = consts.SPhiInf;

% 1. Obtain the thickness correction and the atmospheric pressure.  

% 1a. Thickness scaling factor. 

if sample.thick > 0.0
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

if (~isfield(sample,'pressure'))
    if (strcmp(sample.aa,'std'))
        % Old code
        % sample.pressure = stdatm(sample.elv);
        % New code
%         sample.pressure = NCEPatm_2(sample.lat,sample.long,sample.elv);
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

% Initialize the result flags. 

results.flags = [];

% 2. Make an initial guess at the age. This mainly serves to limit the
% length of the forward calculation and speed things up slightly.

% First, find P according to Stone/Lal SF, get muon production from H2002

if nuclide ~= 3
    
%     P_mu = zeros(length(sample.lat));
%     P_St = P_ref_St * sample.thickSF * sample.shielding * stone2000(sample.lat,sample.pressure,1);
%     SF_St = sample.thickSF * sample.shielding * stone2000(sample.lat,sample.pressure,1);
%     A = l + sample.rho .* sample.E ./consts.Lsp;
%     
%     if sample.E==0;
%         t_simple = N./P_St;
%     else 
%         t_simple = (-1/A)*log(1-(N * A / (P_St))); 
%     end
%     if (N >= ((P_St+P_mu)/A));
%         % if appears to be saturated in simple-age-world, simple age equation
%         % would evaluate to NaN. Avoid this.
%         % set results to -1; this flags it
%         t_simple = -1;
%     else 
%         % Actually do calculation if possible
%         t_simple = (-1/A)*log(1-(N * A / (P_St+P_mu))); 
%     end
    % Actually get the full data from P_mu_total -- needed later
    % Use middle of sample
% RcEst = 14.9.*((cos(abs(d2r(sample.lat)))).^4); %Dipolar estimate for these purposes

    GDRc = [-448.004 1189.18 -1152.15 522.061 -103.241 6.89901 0];
%   Trajectory-traced dipolar estimate for these purposes
    RcEst = polyval(GDRc,cos(d2r(sample.lat)));
    
%     % %Use mean Solar Modulation Parameter (SPhiInf)

    muLS = P_mu_totalLSD((sample.thick.*sample.rho./2),sample.pressure,RcEst,consts.SPhiInf,mconsts,'yes');
    muSt = P_mu_total_alpha1((sample.thick.*sample.rho./2),sample.pressure,mconsts,'yes');
    
    % 
%     P_mu_LS = muLS.P_fast_LS + muLS.P_neg_LS;
    P_mu_LS = muLS.P_fast_St + muLS.P_neg_St;  % NL 4/22 Test
    P_mu_Lm = muSt.P_fast + muSt.P_neg;
    P_mu_St = muSt.P_fast + muSt.P_neg;
    
%     % Don't double-count muons in the following line
    P_St = P_ref_St * sample.thickSF * sample.shielding * stone2000(sample.lat,sample.pressure,1);
    SF_St = sample.thickSF * sample.shielding * stone2000(sample.lat,sample.pressure,1);
    A = l + sample.rho .* sample.E ./consts.Lsp;

%     P_mu = zeros(length(sample.lat));

    P_mu_LS = P_mu_LS';
    P_mu_LD = P_mu_LS';
    P_mu_Lm = P_mu_Lm';
    P_mu_St = P_mu_St';

    if (N >= ((P_St+P_mu_St)/A))
        % if appears to be saturated in simple-age-world, simple age equation
        % would evaluate to NaN. Avoid this.
        % set results to -1; this flags it
        t_simple = -1;
    else 
        % Actually do calculation if possible
        t_simple = (-1/A)*log(1-(N * A / (P_St+P_mu_St))); 
    end
    P_mu_LS = P_mu_LS';
    P_mu_LD = P_mu_LS';
    P_mu_Lm = P_mu_Lm';
    P_mu_St = P_mu_St';
else
    P_mu_LS = zeros(length(sample.lat));
    P_mu_LD = zeros(length(sample.lat));
    P_mu_Lm = P_mu_LS;
    P_mu_St = P_mu_LS;
    
    P_St = P_ref_St * sample.thickSF * sample.shielding * stone2000(sample.lat,sample.pressure,1);
    SF_St = sample.thickSF * sample.shielding * stone2000(sample.lat,sample.pressure,1);
    A = sample.rho .* sample.E ./consts.Lsp;
    
    if sample.E==0
        t_simple = N./P_St;
    else 
        t_simple = (-1/A)*log(1-(N * A / (P_St))); 
    end
end

% catch for negative longitudes before Rc interpolation
if sample.long < 0; sample.long = sample.long + 360;end;

% Make the time domain for the forward age calculation. The time steps 
% follow the paleomagnetic data:
%

% As linear interpolation is used to find the age, using a 1000-yr time step
% step could affect the accuracy at the 10^-3 level in the final timestep
% given the Be-10 and Al-26 half-lives with typical erosion. This has no 
% effect on the eventual age. This would not be acceptable for short-half-life 
% nuclides, i.e. 14-C. 

% Use the non-time-dependent age t_simple to decide on a length for Rc(t). The
% factor of 1.6 is a W.A.G. Note that because the long-term magnetic field is 
% low, old PMC ages will nearly always be less than the simple age. If this goes
% wrong there is a diagnostic flag. 

mt = max(t_simple .* 1.6,5e4); % mt is max time. Changed to take max between t_simple*1.6 and 50ka - NL 01/23

% Make the time vector
tv = geomag.tv;

% clip to limit computations...
if nuclide == 10
    if mt < 0 % Saturated WRT simple age - use full tv
        mt = 1e7;
%     elseif  mt < 12000;
%         mt = 12000; % Don't allow unreasonably short times
    elseif mt > 1e7
        mt = 1e7;
    end
elseif nuclide == 26
    if mt < 0 % Saturated WRT simple age - use full tv
        mt = 5e6;
%     elseif  mt < 12000;
%         mt = 12000; % Don't allow unreasonably short times
    elseif mt > 5e6 
        mt = 5e6;
    end
elseif nuclide == 14
    if mt < 0 % Saturated WRT simple age - use tv to 50 ka
        mt = 50000;
%     elseif  mt < 12000;
%         mt = 12000; % Don't allow unreasonably short times
    elseif mt > 50000 
        mt = 50000;
    end
elseif nuclide == 3
    if mt < 0 % Saturated WRT simple age - use full tv
        mt = 1e7;
%     elseif  mt < 12000;
%         mt = 12000; % Don't allow unreasonably short times
    elseif mt > 1e7 
        mt = 1e7;
    end
end    
% now chop off tv;
clipindex = find(tv < mt, 1, 'last' );
% clipindex = find(tv <= mt, 1, 'last' );
% if not a calibration sample, do the clip
clipindex=clipindex+1; %make sure the array isn't clipped too short
if calFlag == 0
    tv = tv(1:clipindex); % don't actually pad out to what mt is...lose 1 time step
end


% Calculate the unweighted P(t) separately to be sent back in the results.
% This is the surface production rate taking account of thickness. 
% P_St is already calculated
if calFlag == 1

    P_St = geomag.SF_St(i)'.*P_ref_St;% 2/14/22 NL changed from (i,:)
    P_LS = geomag.SF_LS(i,:).*P_ref_LS;
    P_LD = geomag.SF_LD(i,:).*P_ref_LD;
    P_Lm = geomag.SF_Lm(i,:).*P_ref_Lm;
else
    P_St = geomag.SF_St.*P_ref_St;
    P_LS = geomag.SF_LS.*P_ref_LS;
    P_LD = geomag.SF_LD.*P_ref_LD;
    P_Lm = geomag.SF_Lm.*P_ref_Lm;
end


% Also calculate production by muons. 
% This code uses only a highly simplified attenuation-length approximation
% for the depth dependence of production by muons. This is OK here because
% good exposure-dating sites must have low erosion rates. This
% approximation isn't good enough for erosion-rate calculations, and it's
% not used in the erosion rate calculator. 

Lmu = 1500; 


% Special clipping action for calibration measurements.
% Calibration section modified by Nat Lifton 4/09 - returns calculated
% concentrations instead of ages, for minimization by pfit10v2n

if calFlag == 1
    % First, chop off tv
    clipindex = find(tv <= sample.truet, 1, 'last' );
    tv2 = tv(1:clipindex);
    if tv2(end) < sample.truet
        tv2 = [tv2 sample.truet];
    end
    % Now shorten the P's commensurately 
%     P_St2 = interp1(tv,P_St,tv2);    
    P_LS2 = interp1(tv,P_LS,tv2);    
    P_LD2 = interp1(tv,P_LD,tv2);    
    P_Lm2 = interp1(tv,P_Lm,tv2);

    % Give back Pmu
    results.P_mu_LS = P_mu_LS;
    results.P_mu_LD = P_mu_LD;
    results.P_mu_Lm = P_mu_Lm;
    results.P_mu_St = P_mu_St;
    if nuclide == 3
        dcf = 1;
    else
        dcf = exp(-tv2.*l); % decay factor;
    end
        
    % recover the scaling factors
    % Note that sample.thick and sample.shielding should be zero, i.e.
    % already taken out of N someplace upstream
    results.tv = tv2; % this provided for stable isotope calibration

    %Calculated concentrations
    dpfs = exp(-tv2.*sample.E.*sample.rho./consts.Lsp); % spallation depth dependence
    dpfm = exp(-tv2.*sample.E.*sample.rho./Lmu); % muon depth dependence approximation

    N_LS = cumtrapz(tv2,(P_LS2.*dcf.*dpfs + P_mu_LS.*dcf.*dpfm));
    N_LD = cumtrapz(tv2,(P_LD2.*dcf.*dpfs + P_mu_LD.*dcf.*dpfm));
    N_Lm = cumtrapz(tv2,(P_Lm2.*dcf.*dpfs + P_mu_Lm.*dcf.*dpfm));
    N_St = cumtrapz(tv2,(P_St.*dcf.*dpfs + P_mu_St.*dcf.*dpfm));

    results.N_LS = max(N_LS);
    results.N_LD = max(N_LD);
    results.N_Lm = max(N_Lm);
    results.N_St = max(N_St);
    
    % Simplified Error Propagation as below
    
    sfa = ['LS';'LD';'Lm';'St'];
    
    for a = 1:4 % Do everything four times
        % extract t, Pref, delPref for SF
        eval(['tPref = P_ref_' sfa(a,:) ';']);
        eval(['tdelPref = delP_ref_' sfa(a,:) ';']);
        % do most of computation
        if sample.E==0 && nuclide == 3
            FP = N./sample.truet;
            dNdt = -N./(sample.truet.^2);
            delt = sample.deltruet;
        elseif sample.E==0 && nuclide == 14
            FP = N./sample.truet;
            dNdt = -N./(sample.truet.^2);
            delt = sample.deltruet;
        else    
            FP = (N.*A)./(1 - exp(-A.*sample.truet));
            dNdt = FP.*exp(-A.*sample.truet);
            delt = sample.deltruet;
        end
        % make respective delt's
    %    eval(['delN_ext_' sfa(a,:) ' = sqrt(dNdt.^2 * delN.^2 + dNdP.^2 * delFP.^2);']);
        eval(['delN_int_' sfa(a,:) ' = sqrt(dNdt.^2 * delt.^2);']);
       % eval(['FP_' sfa(a,:) ' = FP;']);
        eval(['results.delN_int_' sfa(a,:) ' = delN_int_' sfa(a,:) ';']);
    end
    results.atmdepth = sample.pressure.*1.019716;

    % Done
    return;
end
    
% Calculate N(t) including decay and erosion

if nuclide == 3
    dcf = 1; % decay factor;
else
    dcf = exp(-tv.*l);
end

dpfs = exp(-tv.*sample.E.*sample.rho./consts.Lsp); % spallation depth dependence
dpfm = exp(-tv.*sample.E.*sample.rho./Lmu); % muon depth dependence approximation

if calFlag == 0
      P_LS = P_LS(1:clipindex); % trim production arrays to agree with age array
      P_LD = P_LD(1:clipindex); % trim production arrays to agree with age array
      P_Lm = P_Lm(1:clipindex); % trim production arrays to agree with age array
%       P_St = P_St(1:clipindex); % trim production arrays to agree with age array
end
      
N_LS = cumtrapz(tv,(P_LS.*dcf.*dpfs + P_mu_LS.*dcf.*dpfm));
N_LD = cumtrapz(tv,(P_LD.*dcf.*dpfs + P_mu_LD.*dcf.*dpfm));
N_Lm = cumtrapz(tv,(P_Lm.*dcf.*dpfs + P_mu_Lm.*dcf.*dpfm));
N_St = cumtrapz(tv,(P_St.*dcf.*dpfs + P_mu_St.*dcf.*dpfm));

% Look for saturation with respect to various scaling factors -- 
% If not saturated, get the age by reverse-interpolation.
% Note that this is not necessarily rigorous and doesn't attempt to take
% account of uncertainties in deciding if a sample is saturated. If you
% need rigorous analysis of close-to-saturation measurements, this code is
% not for you. 

if nuclide==10;nstring='Be-10';
elseif nuclide==14;nstring='C-14';
elseif nuclide==26;nstring='Al-26';
elseif nuclide==3;nstring='He-3';
end
    
if N > max(N_LS) 
    flag = ['Sample ' sample.sample_name ' -- ' nstring ' appears to be saturated WRT LSDn Lifton et al 2014 SF.'];
    results.flags = [results.flags '<br>' flag];
    t_LS = 0;
else
    t_LS = interp1(N_LS,tv,N);
end

if N > max(N_LD) 
    flag = ['Sample ' sample.sample_name ' -- ' nstring ' appears to be saturated WRT LSDn Lifton et al 2014 Dipolar SF.'];
    results.flags = [results.flags '<br>' flag];
    t_LD = 0;
else
    t_LD = interp1(N_LD,tv,N);
end

if N > max(N_Lm)
    flag = ['Sample ' sample.sample_name ' -- ' nstring ' appears to be saturated WRT Lal/Stone PMAG SF.'];
    results.flags = [results.flags '<br>' flag];
    t_Lm = 0;
else
    t_Lm = interp1(N_Lm,tv,N);
end

if N > max(N_St) 
    flag = ['Sample ' sample.sample_name ' -- ' nstring ' appears to be saturated WRT Stone(2000) SF.'];
    results.flags = [results.flags '<br>' flag];
    t_St = 0;
else
    t_St = interp1(N_St,tv,N);
end

% Error propagation scheme. 
% This is highly simplified. We approximate the error by figuring out what the
% effective production rate is (disregarding the special depth dependence
% for muons) which gives the right age in the simple age equation. 
% Error in this taken to be linear WRT the reference production rate.
% Then we linearly propagate errors through the S.A.E. 
% We ignore the nominal error in production by muons. This is OK
% because it's small compared to the total error in the reference
% production rate. Future versions will use Monte Carlo error analysis. 

sfa = ['LS';'LD';'Lm';'St'];

if nuclide == 3
    for a = 1:4 % Do everything six times
    % extract t, Pref, delPref for SF
    eval(['tt = t_' sfa(a,:) ';']);
    eval(['tPref = P_ref_' sfa(a,:) ';']);
    eval(['tdelPref = delP_ref_' sfa(a,:) ';']);
    
    % do most of computation
    if sample.E==0
        FP = N./tt;
        delFP = (tdelPref / tPref) .* FP;
        dtdN = 1./FP;
        dtdP = -N./(FP.^2);
    else
        FP = (N.*A)./(1-exp(-A.*tt));
        delFP = (tdelPref / tPref) * FP;
        dtdN = 1./(FP - N.*A);  
        dtdP = -N./(FP.*FP - N.*A.*FP);
    end
    % make respective delt's
    eval(['delt_ext_' sfa(a,:) ' = sqrt( dtdN.^2 * delN.^2 + dtdP.^2 * delFP.^2);']);
    eval(['delt_int_' sfa(a,:) ' = sqrt(dtdN.^2 * delN.^2);']);
    eval(['FP_' sfa(a,:) ' = FP;']);
    end
else
    for a = 1:4 % Do everything three times
    % extract t, Pref, delPref for SF
    eval(['tt = t_' sfa(a,:) ';']);
    eval(['tPref = P_ref_' sfa(a,:) ';']);
    eval(['tdelPref = delP_ref_' sfa(a,:) ';']);
    if tt > 0 % Not saturated, is an age
        % do most of computation
        FP = (N.*A)./(1 - exp(-A.*tt));
        delFP = (tdelPref / tPref) * FP;
        dtdN = 1./(FP - N.*A);  
        dtdP = -N./(FP.*FP - N.*A.*FP);
        % make respective delt's
        eval(['delt_ext_' sfa(a,:) ' = sqrt( dtdN.^2 * delN.^2 + dtdP.^2 * delFP.^2);']);
        eval(['delt_int_' sfa(a,:) ' = sqrt(dtdN.^2 * delN.^2);']);
        eval(['FP_' sfa(a,:) ' = FP;']);
    else % t set to 0, was saturated
        eval(['delt_ext_' sfa(a,:) ' = 0;']);
        eval(['delt_int_' sfa(a,:) ' = 0;']);
        eval(['FP_' sfa(a,:) ' = 0;']);
    end
    end
end

% 5. Results structure assignment

% Thickness scaling factor
results.thick_sf = sample.thickSF;

% Muons

results.P_mu_LS = P_mu_LS;
results.P_mu_Lm = P_mu_Lm;
results.P_mu_St = P_mu_St;

% Time template
results.tv = tv;

% Results x 4 by scaling factor

for a = 1:4
    if a < 4 % No Rc record for non-time-dependent SF
        eval(['results.Rc_' sfa(a,:) ' = geomag.Rc_' sfa(a,:) ';']);
    end
    eval(['results.P_' sfa(a,:) ' = P_' sfa(a,:) ';']); % vector for LS, LD, Lm, scalar for St    
    eval(['results.t_' sfa(a,:) ' = t_' sfa(a,:) ';']); % age
    eval(['results.FSF_' sfa(a,:) ' = FP_' sfa(a,:) './(results.thick_sf.*sample.shielding.*P_ref_' sfa(a,:) ');']); % effective SF
    eval(['results.FP_' sfa(a,:) ' = FP_' sfa(a,:) ';']); % effective PR
    eval(['results.delt_int_' sfa(a,:) ' = delt_int_' sfa(a,:) ';']);
    eval(['results.delt_ext_' sfa(a,:) ' = delt_ext_' sfa(a,:) ';']);
end

% Lal/Stone SF for historical interest
results.SF_St_nominal = stone2000(sample.lat,sample.pressure,Fsp);
results.atmdepth = sample.pressure.*1.019716;


% Version
results.main_version = ver;

