%% PFig5.m — Standard error of all-India rainfall by month (Figure 5)
%
% Reproduces Figure 5 of:
%   Seshadri, A. K. (2018). Statistics of spatial averages and optimal
%   averaging in the presence of missing data. Spatial Statistics, 25, 1-21.
%   doi:10.1016/j.spasta.2018.04.002
%
% Estimates the standard error (SE) of spatially averaged rainfall for
% individual months from April through November, plus multi-month averages
% for ISMR (June-September) and April-November. The SE quantifies
% uncertainty in reports of the all-India average due to measurement noise
% and squared bias from finite spatial sampling:
%
%   SE = { beta^T (D1 + sigma_e^2 I) beta }^{1/2}
%
% where D1 is the bias matrix and sigma_e is measurement noise standard
% deviation. Weights beta are chosen to minimise MSE (alpha = 1).
%
% The 0.25 x 0.25 deg gridded product is treated as the true spatial
% average. Two values of measurement uncertainty are considered
% (sigma_e = 0.5 and 1.0 mm/day). The plot shows that measurement error
% has little effect on SE after spatial averaging; SE comes mainly from
% bias. For ISMR the SE is approximately 0.29 mm/day (mean ~7 mm/day,
% standard deviation ~0.39 mm/day).
%
% For multi-month periods (indices 9 and 10), sigma_e is scaled down by
% sqrt(number of months) to account for temporal averaging. However this
% does not affect the plotted results, since only single-month periods are
% shown in the paper
%
% Data:
%   - 0.25 x 0.25 deg gridded rainfall (meanhrraindat.mat) — true ISMR.
%   - 1 x 1 deg gridded rainfall (indiadat.mat) — 357 point observations.
%   - Lat/lon coordinates (indialatlon.mat).
%
% Requirements: MATLAB Optimization Toolbox (for quadprog).

clear, close all

%% Input periods
monthstartlist = [4 5 6 7 8 9 10 11 6 4]; % each period is of 1 month, last two periods are of 4 and 8 months respectively
monthendlist = [4 5 6 7 8 9 10 11 9 11]; % period duration = monthendlist - monthstartlist + 1

%% Measurement standard deviation
sigmaelist = [0.5 1]; % mm/day

%% probability of available observations
alpha = 1;

%% load & read data

% 0.25 deg data
load ../data/meanhrraindat % meanhrrain
meanhrrain = meanhrrain(1:111,:,:);

% 1 deg rainfall
inddat = load('../data/indiadat.mat');
indlatlon = load('../data/indialatlon.mat');

% initialize
datmat = NaN(111,357,365);

for i = 1901:2011
    dati = inddat.indiarainmodel(i); % year i
    datmati = cell2mat(dati); % 357 x 365 array
    datmat(i-1900,:,:) = single(datmati); % year, loc, day
end


latlon = indlatlon.indialatlon;

lat = single(latlon(:,1)); lon = single(latlon(:,2));

%% Loop for input period

varmat = NaN(numel(monthstartlist),numel(sigmaelist));
semat = NaN(numel(monthstartlist),numel(sigmaelist));
meanrainmat = NaN(numel(monthstartlist),numel(sigmaelist));

for maincount = 1:numel(monthstartlist)
    monthstart = monthstartlist(maincount);
    monthend = monthendlist(maincount);
    
    
    % calculate indices for 1 deg rainfall
    numdaysinmonth = [31 28 31 30 31 30 31 31 30 31 30 31];
    cumsumdays = cumsum(numdaysinmonth);
    cumsumdays = cat(2,0,cumsumdays);
    
    startindex = cumsumdays(monthstart) + 1;
    endindex = cumsumdays(monthend + 1);
    
    % calculate subset of year and time-average
    precip = datmat(:,:,startindex:endindex);
    precipp = squeeze(mean(precip,3)); % {year, loc}
    
    % ismr
    ismr = NaN(111,1);
    for i = 1:size(precipp,1)
        ismr(i) = sum(precipp(i,:)'.*cosd(lat))/sum(cosd(lat));
    end
    
    % calculate indices for 0.25 deg rainfall
    numdaysinmonth25 = [30 31 30 31 31 30 31 30];
    cumsumdays25 = cumsum(numdaysinmonth25);
    cumsumdays25 = cat(2,0,cumsumdays25);
    
    monthstart25 = monthstart - 3;
    monthend25 = monthend - 3;
    
    startindex25 = cumsumdays25(monthstart25) + 1;
    endindex25 = cumsumdays25(monthend25 + 1);
    
    precipphr = meanhrrain(:,startindex25:endindex25);
    
    
    %% ISMR
    ismrhr = mean(precipphr,2);
    
    %% Loop for measurement standard deviation
    
    for sigmaecount = 1:numel(sigmaelist)
        sigmae = sigmaelist(sigmaecount); % mm/day
        
        if maincount == 9
            sigmae = sigmae / sqrt(4); % (scaling stdev for 4 month period); does not affect paper plots, only 1:8 are plotted
        elseif maincount == 10
            sigmae = sigmae / sqrt(8); % (scaling stdev for 8 month period); does not affect paper plots, only 1:8 are plotted
        end
        
        % mean rainfall
        pmean = mean(precipp,1);
        
        % covariance estimation
        Xp = precipp - repmat(pmean,[111 1]);
        
        Sv = 1/110*(Xp'*Xp);
        
        Sr = Sv + sigmae^2*eye(357);
        
        %% Average
        pmean = pmean';
        muv = mean(pmean);
        
        %% Matrices
        Fv = diag((pmean-muv).^2);
        A = diag(pmean-muv);
        
        for i = 1:size(precipp,1)
            dmi = precipp(i,:)'- ismrhr(i); % 
            if i == 1
                Msum = dmi*dmi';
            else
                Msum = Msum + dmi*dmi';
            end
        end
        M = Msum / 111;
        
        %% Variance matrix
        Srd = diag(diag(Sr));
        C = Sr + (1-alpha)/alpha*(Srd+Fv);
        
        %% Optimal sampling
        
        % inputs for quadprog
        fq = zeros(357,1);
        Aeqq = zeros(357,357); Aeqq(1,:) = 1;
        beqq = zeros(357,1); beqq(1) = 1;
        lbq = zeros(357,1);
        
        opts = optimoptions('quadprog','Algorithm','interior-point-convex','Display','iter', 'TolFun',1e-9);
        
        % min. MSE
        Hqm = double(C + M);
        [betam,fvalm,exitflagm,outputm,lambdam] = ...
            quadprog(Hqm,fq,[],[],Aeqq,beqq,lbq,[],[],opts);
        
        % calculate variance
        varm = betam'*C*betam;
        
        % calculate bias squared
        biasm = betam'*M*betam;
        
        % calculate MSE
        msem = varm + biasm;
        
        % MSE from uniform sampling
        u = 1/357*ones(357,1);
        mseunif = u'*(C+M)*u;
        biasunif = u'*M*u;
        varunif = u'*C*u;
        
        %% For each year, calculate optimal average (Assuming no missing data)
        pmeanopt = NaN(size(precipp,1),1);
        
        for count = 1:size(precipp,1)
            precippi = precipp(count,:)';
            pmeanopt(count) = betam'*precippi;
        end
        
        pmeanunif = mean(precipp,2);
        
        %% calculate variance accurately
        varsigma2 = (1/111*sum(precipp.*precipp,1) - alpha*(1/111*sum(precipp,1)).^2)';
        Srnd = Sr;  I = eye(357); idiag = find(I); Srnd(idiag) = 0;
        Evilist = mean(precipp,1)';
        
        varm1 = 1/alpha*sum(betam.^2.*varsigma2);
        varm2 = betam'*Srnd*betam;
        varm3 = 1/alpha*sigmae^2*(betam'*betam);
                
        %% Store outputs
        varmat(maincount,sigmaecount) = varm1 + varm2 + varm3;
        semat(maincount,sigmaecount) = biasm + varm3;
        meanrainmat(maincount,sigmaecount) = mean(ismrhr);
    end
end

%% Plot figure 5

figure(5), semilogy(monthstartlist(1:8), meanrainmat(1:8), 'kx-', 'LineWidth', 1.5), hold on,...
    semilogy(monthstartlist(1:8), sqrt(varmat(1:8,1)), 'b+-', 'LineWidth', 1.0), hold on,...
    semilogy(monthstartlist(1:8), sqrt(varmat(1:8,2)), 'bo--', 'LineWidth', 1.0), hold on,...
    semilogy(monthstartlist(1:8), sqrt(semat(1:8,1)), 'r+-', 'LineWidth', 1.0), hold on,...
    semilogy(monthstartlist(1:8), sqrt(semat(1:8,2)), 'ro--', 'LineWidth', 1.0), hold on,...
    xlabel('month'), ylabel('rainfall (mm/day)'), legend('mean','Stdev: \sigma_e = 0.5 mm/day','Stdev: \sigma_e = 1.0 mm/day','SE: \sigma_e = 0.5 mm/day','SE: \sigma_e = 1.0 mm/day',...
    'Location','North'),...
    ylim([ 0 10]), grid on

