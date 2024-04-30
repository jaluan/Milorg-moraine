function scaling = ScalingLSDLal(sites,consts,nuclide)

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
% of Lal(1991)/Stone(2000). LS = Sato et al. (2008) Analytical FLuxes and
% Bob Reedy Cross-Section Data (9/2010)
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
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).
% What version is this?

sfa = ['LS';'LD';'Lm';'St'];

% If no pressure entered yet, create it from the elevation
% using the appropriate atmosphere
% Change in version 2: This looks up the surface pressure from 
% NCAR map to use in the standard atmosphere equation.
% Obviously, it is always better to estimate pressure from local
% station data. 

num_samples = length(sites.lat);
% Make the time vector

calFlag = 0;
Rcflg = 'n';

% Age Relative to t0=2010
tv = [0:10:50 60:100:3060 3260:200:75060 76060:1000:800060 802060:2000:2000060 logspace(log10(2002060),7,200)];

LDRc = zeros(num_samples,length(tv));
LSRc = zeros(num_samples,length(tv));
LmRc = zeros(num_samples,length(tv));

if isfield(sites,'truet')
    % It's a calibration sample and we want the time-integrated
    % SF, not the age. 
    calFlag = 1;
end

if calFlag == 1
    % First, chop off tv
    clipindex = find(tv <= max(max(sites.truet),50060), 1, 'last' );
    tv2 = tv(1:clipindex);
    if tv2(end) < max(max(sites.truet),50060)
        tv2 = [tv2 max(max(sites.truet),50060)];
    end
    % Now shorten the P's commensurately 
    LDRc2 = zeros(num_samples,length(tv2));    
    LSRc2 = zeros(num_samples,length(tv2));
    LmRc2 = zeros(num_samples,length(tv2));
    %And the age vector
    tv = tv2;
    LDRc = LDRc2;
    LSRc = LSRc2;
    LmRc = LmRc2;
end

% Need solar modulation for Lifton SF's

this_SPhi = zeros(size(tv)) + consts.SPhiInf; % Solar modulation potential for Sato et al. (2008)
this_SPhi(1:length(consts.SPhi)) = consts.SPhi; % Solar modulation potential for Sato et al. (2008)
w = 0.066; % water content for Sato & Niita (2006), gravimetric, about 14% volumetric per Fred Phillips

% interpolate an M for tv > 0 for dipole field. 10/1/15 Pavon et al. SHA.DIF.14k
temp_M1 = interp1(consts.tRc,consts.MM0_Pavon,tv(1:91));
temp_M2 = interp1(consts.tM,consts.GLOPADMM0,tv(92:end));
temp_M = cat(2,temp_M1,temp_M2);

for a = 1:num_samples
    clear sample;
    % Load the input data structure
%     sample.site = sites.site(a);
%     sample.sample_name = sites.site(a);
    sample.lat = sites.lat(a);
    sample.long = sites.long(a);
    if sites.elv(a) < 0
        sample.aa = 'pre';
        sample.pressure = sites.pressure(a);
    elseif strcmp(sites.aa(a),'ant')
        sample.elv = sites.elv(a);
        sample.aa = 'ant';
    else
        sample.elv = sites.elv(a);
        sample.aa = 'std';
    end
    sample.thick = sites.thick(a);
    sample.rho = sites.rho(a);
    sample.shielding = sites.shielding(a);
    sample.E = sites.E(a);
%    sample.age = sites.truet(a);
    
    % 1a. Thickness scaling factor. 

    if sample.thick > 0.0
        sample.thickSF = thickness(sample.thick,consts.Lsp,sample.rho);
    else 
        sample.thickSF = 1;
    end
    
    % 1b. Pressure correction
    
    if (~isfield(sample,'pressure'))
        if (strcmp(sample.aa,'std'))
            % Old code
%             sample.pressure = stdatm(sample.elv);
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
    
    % catch for negative longitudes before Rc interpolation
    if sample.long < 0; sample.long = sample.long + 360;end


%     mt(a) = t_simple .* 1.6; % mt is max time. 

    % Make up the Rc vectors.
    % Start with Lifton et al. 2005
    % First 15 from the data blocks


% TIME-VARYING GAD FIELD (comparable to Lm)
%   Lifton et al. 2014 -  Fit to Trajectory-traced GAD dipole field as f(M/M0), as long-term average. This is the one I'm using now...
    GDRc = [-448.004 1189.18 -1152.15 522.061 -103.241 6.89901 0];

    % For PavonRc - 10/1/15

        [longi,lati,tvi] = meshgrid(sample.long,sample.lat,tv(1:91));
        LSRc(a,1:91) = interp3(consts.lon_Rc,consts.lat_Rc,consts.tRc,consts.PavonRc,longi,lati,tvi);
        LSRc(a,92:end) = temp_M(92:end).*polyval(GDRc,cos(d2r(sample.lat)));
    
       zIndex = LSRc(a,:) < 0;
       LSRc(a,zIndex) = 0;

       LDLat = abs(90-angdist(sample.lat,sample.long,consts.lat_pp_Pavon,consts.lon_pp_Pavon));
       LDRc(a,1:91) = temp_M(1:91).*polyval(GDRc,cos(d2r(LDLat)));
       LDRc(a,92:end) = temp_M(92:end).*polyval(GDRc,cos(d2r(sample.lat)));

       zIndex = LDRc(a,:) < 0;
       LDRc(a,zIndex) = 0;
%        LmRc = LDRc;

       %FINE TO APPLY THIS per Tatsuhiko Sato, Oct 2015
        if Rcflg == 'y'   
        %     %   Calculate effective Rc skymap and from latitude and
        %     %   trajectory-traced Rc if requested
                LSRc(a,:) = interp2(consts.Rcveff,consts.lat,consts.LatRc,LSRc(a,:),sites.lat(a));
                LDRc(a,:) = interp2(consts.Rcveff,consts.lat,consts.LatRc,LDRc(a,:),sites.lat(a));
        end

        % Approximate paleo-pole-positions and field strengths from Pavon for <14 ka - 7/11/14
    %     Original 14.9 GV references 1950. New value (14.31) is referenced to 2010 DGRF.
%         LDRc(a,1:91) = 14.9.*(consts.MM0_Pavon).*((cos(abs(d2r(LDLat)))).^4);              
%         LDRc(a,92:end) = 14.9.*temp_M(92:end).*((cos(abs(d2r(sample.lat)))).^4);

        LmRc(a,1:91) = 14.31.*(consts.MM0_Pavon).*((cos(abs(d2r(LDLat)))).^4);              
        LmRc(a,92:end) = 14.31.*temp_M(92:end).*((cos(abs(d2r(sample.lat)))).^4);

    % Calculate the unweighted P(t) separately to be sent back in the results.
    % This is the surface production rate taking account of thickness. 
    % P_St is already calculated
% 
    LSD = LSDscaling(sample.pressure,LSRc(a,:),this_SPhi,w,consts,nuclide);
    LSDgd = LSDscaling(sample.pressure,LDRc(a,:),this_SPhi,w,consts,nuclide);

    if nuclide == 3
        SF_LS(a,:) = LSD.He.*sample.thickSF.*sample.shielding;
        SF_LD(a,:) = LSDgd.He.*sample.thickSF.*sample.shielding;
    elseif nuclide == 10
        SF_LS(a,:) = LSD.Be.*sample.thickSF.*sample.shielding;
        SF_LD(a,:) = LSDgd.Be.*sample.thickSF.*sample.shielding;
    elseif nuclide == 14
        SF_LS(a,:) = LSD.C.*sample.thickSF.*sample.shielding;
        SF_LD(a,:) = LSDgd.C.*sample.thickSF.*sample.shielding;
    elseif nuclide == 26
        SF_LS(a,:) = LSD.Al.*sample.thickSF.*sample.shielding;
        SF_LD(a,:) = LSDgd.Al.*sample.thickSF.*sample.shielding;
    else    %Total nucleon flux scaling factors as f(Rc)
        SF_LS(a,:) = LSD.sp.*sample.thickSF.*sample.shielding;
        SF_LD(a,:) = LSDgd.sp.*sample.thickSF.*sample.shielding;
    end

    SF_LSDeth(a,:) = LSD.eth.*sample.thickSF.*sample.shielding;
    SF_LSDth(a,:) = LSD.th.*sample.thickSF.*sample.shielding;
    
%     
    SF_Lm(a,:) = stone2000Rcsp(sample.pressure,LmRc(a,:)).*sample.thickSF.*sample.shielding;
    SF_St(a) = sample.thickSF * sample.shielding * stone2000(sample.lat,sample.pressure,1);
    % 5. Results structure assignment

%     Thickness scaling factor
    scaling.thick_sf(a) = sample.thickSF;
        

    % Results x 4 by scaling factor
    eval(['scaling.SF_' sfa(4,:) '(a) = SF_' sfa(4,:) '(a);']); % vector for LS,LD,Lm, scalar for St    

    for b = 1:3
%         if b < 3 % No Rc record for non-time-dependent SF
            eval(['scaling.Rc_' sfa(b,:) '(a,:) = ' sfa(b,:) 'Rc(a,:);']);
            eval(['scaling.SF_' sfa(b,:) '(a,:) = SF_' sfa(b,:) '(a,:);']); % vector for De,Du,Li,LS,Lm, scalar for St    
%         end
       %eval(['scaling.SF_' sfa(a,:) ' = SF_' sfa(a,:) ';']); % scalar for St   
       if b == 3
            scaling.SF_LSDeth(a,:) = SF_LSDeth(a,:);
            scaling.SF_LSDth(a,:) = SF_LSDth(a,:);
        end
    end
end

% Lal/Stone SF for historical interest
%geomag.SF_St_nominal = stone2000(sample.lat,sample.pressure,Fsp);

% Time template
scaling.tv = tv;
%     geomag.tv(a,:) = tv(a,:);

