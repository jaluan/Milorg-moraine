function rockfallMCvJ2(snr,sig)

% function rockfallMCvJ1(snr)
% Inverse MCMC model exploring rockfall timing and moraine residence time 
% for samples with measured cosmogenic nuclides from blue-ice moraine in 
% Milorgfjella. Delineates the ensemble of parameter values with lowest 
% data misfit.

% Inputs: 
% snr: scalar or vector containing sample numbers of samples in input data
% file generated from compile_data.m. If multiple numbers are given, the
% exposure history of the samples on the rockwall and on the moraine 
% will be assumed to be the same, while exhumation depth within rockwall 
% may vary.
% sig: controls uncertainty on measurements: 1=1 sigma, 2=2 sigma etc.

% Written by David Lundbek Egholm, Aarhus University
% Modified by Jane Lund Andersen to also i) include 14C, ii) combine with 
% rockfall forward model

% Versions
% vJ2: also includes elevation of rockfall as a free parameter in model.

close all;
addpath('Functions')

%% **** load data *****
load('./Data/milorgmoraine_mn_data_v2.mat','sample') %_v2: LSDn scaling
load Prockfall.mat Prockfall

%Cosmogenic halflives etc.
CNprop = getCNprop();

%%%%%% Model setup %%%%%%%%%%

%calculate number of samples input
model.Nsnr = length(snr);

%number of sample specific parameters
model.Nsmp = 1; %(Rockwall depth)
    
%initialize
model.Temp = 1.0; %'temperature'; >1 to increase nuclide uncertainties, if=1 no significance, overwritten below
model.Mmp = 3; %number of common model parameters

%Rockwall exposure duration
model.mp{1}.name = 'Rockwall exposure duration (Myr)';
model.mp{1}.vmin = 1e-3; %Myr
model.mp{1}.vmax = 2; %Myr

%Time of rockfall, here equivalent of time on moraine
model.mp{2}.name = 'Time of rockfall (Ma)'; 
model.mp{2}.vmin = 2e-3; %Ma 
model.mp{2}.vmax = 18e-3; %Ma

%Elevation of rockfall in m
model.mp{3}.name = 'Rockfall elevation (m)'; 
model.mp{3}.vmin = 1450; %m 
model.mp{3}.vmax = 1950; %m

%loop samples to set model parameter boundaries for exhumation and get sample info
for i=1:model.Nsnr
    model.mp{model.Mmp+(i-1)*model.Nsmp+1}.name =  ['Z in rockwall (m), 17MFM-0',num2str(snr(i))];
    model.mp{model.Mmp+(i-1)*model.Nsmp+1}.vmin = 0;
    model.mp{model.Mmp+(i-1)*model.Nsmp+1}.vmax = 5; %m
    
    %add sample info to models
    model.data{i} = sample{snr(i)};
    
    %number of data per sample (nuclides)
    model.Nds(i) = sample{snr(i)}.Nnuc;

end

%save some general MCMC paramers
model.Nwalk = 10; %10; % number of walkers, reduce for testing purposes
model.burnin = 5e3; %4e3; % length of burn-in phase
model.Nmod = 50e3; %40e3; % target number of accepted models per walker
model.Nmax = 1e6; %100e4; % maximum number of models, stops inversion from running indefinitely if something is wrong

%number of total model parameters
model.Nmp = model.Mmp + model.Nsnr*model.Nsmp;

%resample model parameters into vectors
for i=1:model.Nmp
    umin(i) = model.mp{i}.vmin;
    umax(i) = model.mp{i}.vmax;
    du0(i) = model.mp{i}.vmax - model.mp{i}.vmin;
end
umin = umin(:);
umax = umax(:);

%data and covariance
for i=1:model.Nsnr %loop samples
    for j=1:sample{snr(i)}.Nnuc %loop nuclides
        nuclide=sample{snr(i)}.nuclides(j); %get nuclide identification
            if nuclide == 1 %10Be
                    dobs((i-1)*model.Nds+j) = model.data{i}.N10;
                    sigd((i-1)*model.Nds+j) = model.data{i}.dN10*sig;
            elseif nuclide == 2 %26Al
                    dobs((i-1)*model.Nds+j) = model.data{i}.N26;
                    sigd((i-1)*model.Nds+j) = model.data{i}.dN26*sig;
            elseif nuclide == 3 %14C
                    dobs((i-1)*model.Nds+j) = model.data{i}.N14;
                    sigd((i-1)*model.Nds+j) = model.data{i}.dN14*sig;
            else
                    warning('Nuclide not implemented')
            end
    end
end

Cobs = model.Temp*diag(sigd.^2);
Cobsinv = inv(Cobs);

%%%%%% Initiate random model generation %%%%%%%

%Initialize random sequence
rng('default');

%set walker starting positions
if model.Nwalk > 1 %multiple walkers
    for i=1:model.Nmp
        wini(i,:) = 0.8*(randperm(model.Nwalk)-1)/(model.Nwalk-1)+0.1;
    end
else %only one walker
    wini = randi([0,4],model.Nmp,model.Nwalk)/4;
end

%loop walkers
for nw = 1:model.Nwalk

    %walker starting position - for initial parameter vector
    for i = 1:model.Nmp
        u(i) = (1-wini(i,nw))*model.mp{i}.vmin + wini(i,nw)*model.mp{i}.vmax;
    end     
 
    %initialize
    res_current = 1e20; %current residual - first model run
    acount = 0; %acceptance count
    bcount = 0; %burn-in count
    rcount = 0; %reject count
    accrat = 0; %acceptance rate
    status = zeros(model.Nmax,1); %model status
    up_rec = zeros(model.Nmp,model.Nmax); %proposed parameters
    u_rec = zeros(model.Nmp,model.Nmax);
    restot_rec = zeros(model.Nmax,1);
    accrat_rec = zeros(model.Nmax,1);
    k_rec = zeros(model.Nmax,1);
    accfac = 1e-2;
    
    %run models
    mi = 0; %model iteration
    k = 0.01; %initial step length

    while ((mi < model.Nmax)&&(acount < model.Nmod))
        
        mi = mi + 1; %update model iteration

%         disp(['nw = ',num2str(nw),'/',num2str(model.Nwalk),' mi = ',num2str(mi),'/',num2str(model.Nmax),' bcount = ',num2str(bcount),' acount = ',num2str(acount),' accrat = ',num2str(accrat),' k = ',num2str(k)]);
    
       %***** step length ******
       
       %set acceptance ratio
        if (mi > 100)
            accrat = (sum(abs(status((mi-100):(mi-1))))+1)/100; %update based on status of last 100 models
        elseif (bcount < model.burnin)
            accrat = 0.1; %first 100 models of burnin
        else 
            accrat = 0.3; %only gets here if model.burnin is <100?
        end
        
        %steplength
        if (bcount < 0.5*model.burnin) %first half of burn-in
            
            k = k*((1-accfac) + accfac*accrat/0.1);
            
            %take larger steps if too few models are accepted intially
            if ((mi > 100)&&(bcount < 2)) k = 1.0; %limited to 0.5 below
            elseif ((mi > 10)&&(bcount < 2)) k = 0.01*mi;
            end
            
        elseif (bcount < model.burnin) %second half of burn-in

            k = k*((1-accfac) + accfac*accrat/0.2);
           
            if ((mi > 100)&&(bcount < 2)) k = 1.0; %Does not enter here unless burn-in is <200?
            elseif ((mi > 10)&&(bcount < 2)) k = 0.01*mi; %and here if burn-in is <20?
            end
            
        
        elseif (acount < model.Nmod) %after burn-in

            k = k*((1-accfac) + accfac*accrat/0.3);    
            
        end  
        
        if (k > 0.5) k = 0.5; end %steplength limited to 0.5
        
        %model 'temperature' gradually reduced from 11 to 1 during burn-in
        if (bcount < model.burnin) model.Temp = 1.0 + 10.0*(model.burnin-bcount)/model.burnin;
        else model.Temp = 1.0;
        end
        
        %********* propose new parameters ************
        
         %random step
        du = 0.5*randn(model.Nmp,1).*du0(:);
        
        %proposed model
        up = u(:) + k*du(:);
     
        
        %retake if outside parameter boundaries
        while (any(up(:) < umin(:))||any(up(:) > umax(:)))
        
            %random step
            du = 0.5*randn(model.Nmp,1).*du0(:);

            %proposed model
            up = u(:) + k*du(:);
            
        end
        
        %********** Forward model *****************
%         [gmp] = rockfall_forward_v2(up,model,CNprop);
        [gmp] = rockfall_forward_v3(up,model,CNprop,Prockfall);
            
        %Acceptance critieria
        %restot = (dobs(:)-gmp(:))'*Cobsinv*(dobs(:)-gmp(:));
        restot = (dobs(:)-gmp(:))'*Cobsinv*(dobs(:)-gmp(:))/model.Temp;
        rfrac = exp(-0.5*restot)/exp(-0.5*res_current);
        alpha = rand(1);
    
    
        %if model is accepted
        if ((alpha < rfrac)||(mi == 1))
        
            u = up;
            gm = gmp;
            res_current = restot;
        
            %accepted after burn-in
            if (bcount > model.burnin)
        
                status(mi) = 1;
                acount = acount + 1;
                
            %if accepted during burn-in    
            else
        
                status(mi) = -1;
                bcount = bcount + 1;
            
            end
        
        %rejected
        else
        
            status(mi) = 0;
            rcount = rcount + 1;                 
    
        end
       
        %save things
        up_rec(:,mi) = up(:);
        u_rec(:,mi) = u(:);
        gm_rec(:,mi) = gm(:);
        restot_rec(mi) = res_current;
        accrat_rec(mi) = accrat;
        k_rec(mi) = k;
                
    end

    %change status for models that are rejected during burn-in
    Imin=find(status == 1,1); %find index of first accepted model after burn-in
    I = find(status(1:Imin) == 0); %find indices of rejected models in burn-in
    status(I) = -2; %set status of these models to -2
    
    %save things
    model.walker{nw}.status = status(1:mi);
    model.walker{nw}.up = up_rec(:,1:mi);
    model.walker{nw}.u = u_rec(:,1:mi);
    model.walker{nw}.gm = gm_rec(:,1:mi);
    model.walker{nw}.restot = restot_rec(1:mi);
    model.walker{nw}.acount = acount;
    model.walker{nw}.bcount = bcount; 
    model.walker{nw}.rcount = rcount;
    model.walker{nw}.accrate = accrat_rec(1:mi);
    model.walker{nw}.kstep = k_rec(1:mi);

end

%%save output
%sample numbers
str = num2str(snr(1));
if length(snr) > 1
    for i=2:length(snr)
        str = [str,['-',num2str(snr(i))]];
    end
end
%filename
savefile = ['models/MC_rockfall_sample_',str,'_v2_temp1_RFmax18ka_max5m_LSDn_sig',num2str(sig),'.mat'];

%save
save(savefile,'model');
