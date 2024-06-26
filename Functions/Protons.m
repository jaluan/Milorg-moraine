function P = Protons(h,Rc,s,consts,nuclide)

% Sato et al. (2008) Neutron Spectrum
% Analytical Function Approximation (PARMA)
% Implemented in MATLAB by Nat Lifton, 2013
% Purdue University, nlifton@purdue.edu

% Copyright 2013, Purdue University
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 3,
% as published by the Free Software Foundation (www.fsf.org).

x = h.*1.019716; % Convert pressure (hPa) to atm depth (g/cm2)
E = consts.E;
% E = logspace(0,5.3010,200);
% E = [1.1295 11.295 112.95 1129.5 11295];

A = 1;
Z = 1;
Ep = 938.27; % Rest mass of a proton
U = (4-1.675).*pi.*A./Z.*1e-7; % Unit conversion factor

% Flatten low rigidities. 

lowRc = find(Rc < 1.0);
Rc(lowRc) = 1.0 + zeros(size(lowRc));

% Primary spectrum

smin = 400; %units of MV
smax = 1200; %units of MV

a1 = 2.1153;
a2 = 4.4511e-1;
a3 = 1.0064e-2;
a4 = 3.9564e-2;
a5 = 2.9236;
a6 = 2.7076;
a7 = 1.2663e4;
a8 = 4.8288e3;
a9 = 3.2822e4;
a10 = 7.4378e3;
a11 = 3.4643;
a12 = 1.6752;
a13 = 1.3691;
a14 = 2.0665;
a15 = 1.0833e2;
a16 = 2.3013e3;

Etoa = E + a1.*x;
Rtoa = 0.001.*sqrt((A.*Etoa).^2 + 2.*A.*Ep.*Etoa)./Z;

Elis = zeros(1,length(E));
Beta = zeros(1,length(E));
Rlis = zeros(1,length(E));
phiTOA = zeros(1,length(E));
phiLIS = zeros(1,length(E));
phiSec = zeros(1,length(E));
phiPtot = zeros(1,length(E));
p10p = zeros(1,length(E));

% Secondary Spectrum

c11 = 1.2560;
c12 = 3.2260e-3;
c13 = -4.8077e-6;
c14 = 2.2825e-9;
c21 = 4.3783e-1;
c22 = -5.5759e-4;
c23 = 7.8388e-7;
c24 = -3.8671e-10;
c31 = 1.8102e-4;
c32 = -5.1754e-7;
c33 = 7.5876e-10;
c34 = -3.8220e-13;
c41 = 1.7065;
c42 = 7.1608e-4;
c43 = -9.3220e-7;
c44 = 5.2665e-10;

b1 = c11 + c12.*x + c13.*x.^2 + c14.*x.^3;
b2 = c21 + c22.*x + c23.*x.^2 + c24.*x.^3;
b3 = c31 + c32.*x + c33.*x.^2 + c34.*x.^3;
b4 = c41 + c42.*x + c43.*x.^2 + c44.*x.^3;

h11min = 2.4354e-3;
h11max = 2.5450e-3;
h12min = -6.0339e-5;
h12max = -7.1807e-5;
h13min= 2.1951e-3;
h13max = 1.4580e-3;
h14min = 6.6767;
h14max = 6.9150;
h15min = 9.3228e-1;
h15max = 9.9366e-1;
h21min = 7.7872e-3;
h21max = 7.6828e-3;
h22min = -9.5771e-6;
h22max = -2.4119e-6;
h23min = 6.2229e-4;
h23max = 6.6411e-4;
h24min = 7.7842;
h24max = 7.7461;
h25min = 1.8502;
h25max = 1.9431;
h31min = 9.6340e-1;
h31max = 9.7353e-1;
h32min = 1.5974e-3;
h32max = 1.0577e-3;
h33min = -7.1179e-2;
h33max = -2.1383e-2;
h34min = 2.2320;
h34max = 3.0058;
h35min = 7.8800e-1;
h35max = 9.1845e-1;
h41min = 7.8132e-3;
h41max = 7.3482e-3;
h42min = 9.7085e-11;
h42max = 2.5598e-5;
h43min = 8.2392e-4;
h43max = 1.2457e-3;
h44min = 8.5138;
h44max = 8.1896;
h45min = 2.3125;
h45max = 2.9368;

h51 = 1.9100e-1;
h52 = 7.0300e-2;
h53 = -6.4500e-1;
h54 = 2.0300;
h55 = 1.3000;
h61 = 5.7100e-4;
h62 = 6.1300e-6;
h63 = 5.4700e-4;
h64 = 1.1100;
h65 = 8.3700e-1;

% Combine primary and secondary spectra

for a = 1:length(Rc)
    Elis = Etoa + s(a).*Z./A;
    Beta = sqrt(1-(Ep./(Ep + Elis.*A)).^2); % Particle speed relative to light
    Rlis = 0.001.*sqrt((A.*Elis).^2 + 2.*A.*Ep.*Elis)./Z;
    C = a7 + a8./(1 + exp((Elis - a9)./a10));

    phiTOA = (C.*(Beta.^a5)./(Rlis.^a6)).*(Rtoa./Rlis).^2;
    phiPri = (U./Beta).*phiTOA.*(a2.*exp(-a3.*x) + (1 - a2).*exp(-a4.*x));
    
    g1min = h11min + h12min.*Rc(a) + h13min./(1 + exp((Rc(a) - h14min)./h15min));
    g1max = h11max + h12max.*Rc(a) + h13max./(1 + exp((Rc(a) - h14max)./h15max));
    g2min = h21min + h22min.*Rc(a) + h23min./(1 + exp((Rc(a) - h24min)./h25min));
    g2max = h21max + h22max.*Rc(a) + h23max./(1 + exp((Rc(a) - h24max)./h25max));
    g3min = h31min + h32min.*Rc(a) + h33min./(1 + exp((Rc(a) - h34min)./h35min));
    g3max = h31max + h32max.*Rc(a) + h33max./(1 + exp((Rc(a) - h34max)./h35max));
    g4min = h41min + h42min.*Rc(a) + h43min./(1 + exp((Rc(a) - h44min)./h45min));
    g4max = h41max + h42max.*Rc(a) + h43max./(1 + exp((Rc(a) - h44max)./h45max));

    phiPmin = g1min.*(exp(-g2min.*x) - g3min.*exp(-g4min.*x)); %length of Rc
    phiPmax = g1max.*(exp(-g2max.*x) - g3max.*exp(-g4max.*x)); %length of Rc

    g5 = h51 + h52.*Rc(a) + h53./(1 + exp((Rc(a) - h54)./h55));
    g6 = h61 + h62.*Rc(a) + h63./(1 + exp((Rc(a) - h64)./h65));

    f3 = g5 + g6.*x;
    f2 = (phiPmin - phiPmax)./(smin.^f3 - smax.^f3);
    f1 = phiPmin - f2.*smin.^f3;

    phiP = f1 + f2.*s(a).^f3;
    
    phiSec = (phiP.*b1.*E.^b2)./(1 + b3.*E.^b4);
    
    Ec = (sqrt((1000.*Rc(a).*Z).^2 + Ep.^2) - Ep)./A;
    Es = a13.*(Ec - a14.*x);
    Es1 = max(a15,Es);
    Es2 = max(a16,Es);
    
    phiPtot = phiPri.*(tanh(a11.*(E./Es1 - 1)) + 1)./2 + ...
        phiSec.*(tanh(a12.*(1 - E./Es2)) + 1)./2;
    
    clipindex = find(E <= 1, 1, 'last' ); %Make sure the clip index is consistent with the definition of E above
    
    if nuclide == 3
%     Quartz
        P.P3p(a) = (trapz(E,phiPtot.*consts.OpxHe3T).*consts.NatomsQtzO + trapz(E,phiPtot.*consts.SipxHe3T).*consts.NatomsQtzSi).*1e-27.*3.1536e7;    
% %     Olivine
% %         Forsterite
%         N.P3p(a) = (trapz(E,phiPtot.*consts.MgpxHe3T).*consts.NatomsOlFoMg...
%             + trapz(E,phiPtot.*consts.SipxHe3T).*consts.NatomsOlFoSi...
%             + trapz(E,phiPtot.*consts.OpxHe3T).*consts.NatomsOlFoO).*1e-27.*3.1536e7;  
% %         Fayalite
%         N.P3p(a) = (trapz(E,phiPtot.*consts.FepxHe3T).*consts.NatomsOlFaFe...
%             + trapz(E,phiPtot.*consts.SipxHe3T).*consts.NatomsOlFaSi...
%             + trapz(E,phiPtot.*consts.OpxHe3T).*consts.NatomsOlFaO).*1e-27.*3.1536e7;  
% %         80% Fo 20% Fa
%         N.P3p(a) = (trapz(E,phiPtot.*consts.MgnxHe3T).*consts.NatomsOlFo80Mg...
%             + trapz(E,phiPtot.*consts.FepxHe3T).*consts.NatomsOlFo80Fe...
%             + trapz(E,phiPtot.*consts.SipxHe3T).*consts.NatomsOlFo80Si...
%             + trapz(E,phiPtot.*consts.OpxHe3T).*consts.NatomsOlFo80O).*1e-27.*3.1536e7; 
% %     Pyroxene
% %         Enstatite
%         N.P3p(a) = (trapz(E,phiPtot.*consts.MgpxHe3T).*consts.NatomsOPxEnMg...
%             + trapz(E,phiPtot.*consts.SipxHe3T).*consts.NatomsOPxEnSi...
%             + trapz(E,phiPtot.*consts.OpxHe3T).*consts.NatomsOPxEnO).*1e-27.*3.1536e7;  
% %         Ferrosilite
%         N.P3p(a) = (trapz(E,phiPtot.*consts.FepxHe3T).*consts.NatomsOPxFsFe...
%             + trapz(E,phiPtot.*consts.SipxHe3T).*consts.NatomsOPxFsSi...
%             + trapz(E,phiPtot.*consts.OpxHe3T).*consts.NatomsOPxFsO).*1e-27.*3.1536e7;  
% %         Wollastonite
%         N.P3p(a) = (trapz(E,phiPtot.*consts.CapxHe3T).*consts.NatomsCPxCa...
%             + trapz(E,phiPtot.*consts.SipxHe3T).*consts.NatomsCPxSi...
%             + trapz(E,phiPtot.*consts.OpxHe3T).*consts.NatomsCPxO).*1e-27.*3.1536e7;  
    elseif nuclide == 10
%     Quartz
        P.P10p(a) = (trapz(E,phiPtot.*consts.O16pxBe10).*consts.NatomsQtzO...
            + trapz(E,phiPtot.*consts.SipxBe10).*consts.NatomsQtzSi).*1e-27.*3.1536e7;
%         P.P10p(a) = (trapz(E,phiPtot.*consts.O16pxBe10)...
%             + trapz(E,phiPtot.*consts.SipxBe10)./2).*consts.NatomsQtzO.*1e-27.*3.1536e7;        
%         P.P10p(a) = (trapz(E,phiPtot.*consts.O16pxBe10)).*consts.NatomsQtzO.*1e-27.*3.1536e7;
%         P.P10p(a) = (trapz(E,phiPtot.*consts.SipxBe10)).*consts.NatomsQtzSi.*1e-27.*3.1536e7;

% %     Magnetite
% %       P.P10p(a) = (trapz(E,phiPtot.*consts.FepxBe10).*consts.NatomsMagFe...
% %           + trapz(E,phiPtot.*consts.O16pxBe10).*consts.NatomsMagO).*1e-27.*3.1536e7;   
%         P.P10p(a) = trapz(E,phiPtot.*consts.FepxBe10).*consts.NatomsMagFe.*1e-27.*3.1536e7;   
%         P.P10p(a) = trapz(E,phiPtot.*consts.O16pxBe10).*consts.NatomsMagO.*1e-27.*3.1536e7;   
    elseif nuclide == 14    
        P.P14p(a) = (trapz(E,phiPtot.*consts.O16pxC14)+ trapz(E,phiPtot.*consts.SipxC14./2)).*consts.NatomsQtzO.*1e-27.*3.1536e7;
    elseif nuclide == 26
        P.P26p(a) = trapz(E,phiPtot.*consts.SipxAl26).*consts.NatomsQtzSi.*1e-27.*3.1536e7; 

    else
        P.pflux(a) = trapz(E(clipindex:end),phiPtot(clipindex:end));
    end
end

P.E = E;

% Plot it

%     figure;clf;
%     loglog(E,phiPtot(1,:));hold on;
%     figure;clf;
%     loglog(E,phiSec(1,:));hold on;



