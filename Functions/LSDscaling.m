function out = LSDscaling(h,Rc,SPhi,w,consts,nuclide) 

% Implements the Lifton Sato et al scaling scheme for spallation.
%
% Syntax: scalingfactor = LiftonSatoSX(h,Rc,SPhi,w,consts);
%
% Where:
%   h = atmospheric pressure (hPa)
%   Rc = cutoff rigidity (GV)
%   SPhi = solar modulation potntial (Phi, see source paper)
%   w = fractional water content of ground (nondimensional)
%   
%
% Vectorized. Send in scalars or vectors of common length. 
%
% Modified by Shasta Marrero (NMT) to include the reactions for Ti & Fe to
% produce chlorine-36. July 2011.
%
% Written by Nat Lifton 2011, Purdue University
% Based on code by Greg Balco -- UW Cosmogenic Nuclide Lab
% balcs@u.washington.edu
% April, 2007
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math
%
%Modified by Allie Koester, 2020 Purdue University
%koestea@purdue.edu
%
% Copyright 2001-2007, University of Washington
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).

% convert pressure to atmospheric depth

X = h.*1.019716;

%Per Tatsuhiko Sato, personal communication, 2013, convert annually averaged Usoskin 
%solar modulation potential to Sato Force Field Potential due to different
%assumed Local Interstellar Spectrum and other factors

% SPhi = 1.1381076.*SPhi - 1.2738468e-4.*SPhi.^2; already included in
% constants file

% Sref = 509.6; %1950-1850 mean
% Sref = 416.518; %Long-term mean
Sref = SPhi(1); %587.4 MV mean 2001-2010 from Usoskin et al 2011. USE THIS!
% Sref = 400; %Idealized solar minimum value, similar to long-term mean (11.4 ka)
% Changed from SPhi(1) 12/8/11 - more consistent with other calculations
% assuming idealized reference states such as SLHL pressure and Rc = 0. BUT
% does not normalize flux at t0 to 1 at SLHL - more like 0.92 - results in
% really high production predictions. Need to use the value from the 2001-2010.
% Go back to SPhi(1)... 12/15/11

% Select reference values for nuclide of interest or flux

if nuclide == 3
    HeRef = consts.P3nRef_q + consts.P3pRef_q;
elseif nuclide == 10
    BeRef = consts.P10nRef_q + consts.P10pRef_q;
elseif nuclide == 14
    CRef = consts.P14nRef_q + consts.P14pRef_q;
elseif nuclide == 26
    AlRef = consts.P26nRef_q + consts.P26pRef_q;
else    
    SpRef = consts.nfluxRef + consts.pfluxRef;
    % Sato et al. (2008) Reference hadron flux integral >1 MeV
end

EthRef = consts.ethfluxRef;
ThRef = consts.thfluxRef;

% Full version with cross sections and includes thermal and
% epithermal fluxes

NSite = Neutrons(h,Rc,SPhi,w,consts,nuclide); 
[ethflux,thflux] = NeutronsLowE(h,Rc,SPhi,w);
PSite = Protons(h,Rc,SPhi,consts,nuclide);

%Nuclide-specific scaling factors as f(Rc)
if nuclide == 3
    Site.He = (NSite.P3n + PSite.P3p)./HeRef;
elseif nuclide == 10
    Site.Be = (NSite.P10n + PSite.P10p)./BeRef;
elseif nuclide == 14
    Site.C = (NSite.P14n + PSite.P14p)./CRef; %add in other fluxes from elements?
elseif nuclide == 26
    Site.Al = (NSite.P26n + PSite.P26p)./AlRef;
else    %Total nucleon flux scaling factors as f(Rc)
    Site.sp = ((NSite.nflux + PSite.pflux))./SpRef; % Sato et al. (2008) Reference hadron flux integral >1 MeV
end

Site.E = NSite.E;%Nucleon flux energy bins
Site.eth = ethflux./EthRef; %Epithermal neutron flux scaling factor as f(Rc)
Site.th = thflux./ThRef;%Thermal neutron flux scaling factor as f(Rc)


% % Full version with cross sections calculates absolute 3He, 10Be, 14C, and 26Al production rates
% 
% [nflux,P3n,P10n,P14n,P26n,P36Can,P36Kn,P36Tin,P36Fen] = NeutronsXS(h,Rc,SPhi,w,consts);
% [ethflux,thflux] = NeutronsLowE(h,Rc,SPhi,w);
% [pflux,P3p,P10p,P14p,P26p,P36Cap,P36Kp,P36Tip,P36Fep] = ProtonsXS(h,Rc,SPhi,consts);
% Site.sp = ((nflux + pflux))./SpRef;
% Site.He = (P3n + P3p);
% Site.Be = (P10n + P10p);
% Site.C = (P14n + P14p);
% Site.Al = (P26n + P26p);
% Site.ClCa = (P36Can + P36Cap)./ClCaRef;
% Site.ClK = (P36Kn + P36Kp)./ClKRef;
% Site.ClTi = (P36Tin + P36Tip)./ClTiRef;
% Site.ClFe = (P36Fen + P36Fep)./ClFeRef;
% Site.eth = ethflux./EthRef;
% Site.th = thflux./ThRef;

out = Site;

