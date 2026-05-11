%% PFig2_3.m — Optimal weights for bias, variance, and MSE (Figures 2-3)
%
% Reproduces Figures 2 and 3 of:
%   Seshadri, A. K. (2018). Statistics of spatial averages and optimal
%   averaging in the presence of missing data. Spatial Statistics, 25, 1-21.
%   doi:10.1016/j.spasta.2018.04.002
%
% Figure 2: Optimal weights for minimising squared bias in the estimate of
%   ISMR (alpha = 1). Panel (a) shows the spatial distribution of weights
%   over India; panel (b) plots weights against the squared deviation of
%   each location's expected rainfall from the spatial average, showing that
%   large weights occur only where this deviation is small (Section 3.1).
%
% Figure 3: Optimal weights for minimising variance (panels a, c) and MSE
%   (panels b, d) for alpha = 1.0 (top) and alpha = 0.8 (bottom). Weights
%   are computed via quadratic programming. Variance minimisation favours
%   locations with small local variance and covariance; MSE minimisation
%   balances both bias and variance considerations.
%
% The script loops over two values of alpha (1.0 and 0.8), solves three
% quadratic programs per alpha (min bias, min variance, min MSE), runs
% 5000-member Monte Carlo ensembles to validate, and stores the optimal
% weight vectors for plotting.
%
% Data:
%   - 0.25 x 0.25 deg gridded rainfall (meanhrraindat.mat) — true ISMR.
%   - 1 x 1 deg gridded rainfall (indiadat.mat) — 357 point observations.
%   - Lat/lon coordinates (indialatlon.mat).
%
% Requirements: MATLAB Optimization Toolbox (for quadprog),
%               Mapping Toolbox (for contourm), coastlines dataset.

clear, close all

%% Coastlines
% load coast % older versions of Matlab
% latcoast = lat; loncoast = long; % older versions of Matlab

load coastlines
latcoast = coastlat; loncoast = coastlon;

%% probability of available observations
minusalphalist = [0 0.20]; %
alphalist = 1 - minusalphalist;

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

%% average over input period

% input period
monthstart = 6;
monthend = 9;

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


%% Covariance matrix

% measurement standard deviation
sigmae = 1; % mm / day

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
    dmi = precipp(i,:)'- ismrhr(i); 
    if i == 1
        Msum = dmi*dmi';
    else
        Msum = Msum + dmi*dmi';
    end
end
M = Msum / 111;

%% Alpha loop

vardat = NaN(numel(alphalist),1);
biasdat = NaN(numel(alphalist),1);
varmodel = NaN(numel(alphalist),1);
varmodelmat = NaN(5,numel(alphalist),1);


biasmodel = NaN(numel(alphalist),1);

varmodelquad = NaN(numel(alphalist),1);
biasmodelquad = NaN(numel(alphalist),1);

betavmat = NaN(357,numel(alphalist));
betabmat = NaN(357,numel(alphalist));
betammat = NaN(357,numel(alphalist));

for alphacount = 1:numel(alphalist)
    alpha = alphalist(alphacount);
    
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
    
    % min. variance
    Hqv = double(C);
    [betav,fvalv,exitflagv,outputv,lambdav] = ...
        quadprog(Hqv,fq,[],[],Aeqq,beqq,lbq,[],[],opts);
    
    % min. bias squared
    Hqb = double(M);
    [betab,fvalb,exitflagb,outputb,lambdab] = ...
        quadprog(Hqb,fq,[],[],Aeqq,beqq,lbq,[],[],opts);
    
    % min. MSE
    Hqm = double(C + M);
    [betam,fvalm,exitflagm,outputm,lambdam] = ...
        quadprog(Hqm,fq,[],[],Aeqq,beqq,lbq,[],[],opts);
    
    % calculate variance
    varv = betav'*C*betav;
    varb = betab'*C*betab;
    varm = betam'*C*betam;
    
    % calculate bias squared
    biasv = betav'*M*betav;
    biasb = betab'*M*betab;
    biasm = betam'*M*betam;
    
    % calculate MSE
    msev = varv + biasv;
    mseb = varb + biasb;
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
    
    
    %% For each year, calculate optimal average (with missing data, with probability 1-alpha)
    
        
    Nensemb = 5000;
    
    varensemblist = NaN(Nensemb,1);
    biasensemblist = NaN(Nensemb,1);
    
    for outcount = 1:Nensemb
        pmeanoptm = NaN(size(precipp,1),1);
        pmeanoptmne = NaN(size(precipp,1),1);
        for count = 1:size(precipp,1)
            scount = binornd(1,alpha,357,1);
            precippi = (precipp(count,:) + randn(1,357)*sigmae)';
            pmeanoptm(count) = sum(betam.*precippi.*scount)/sum(betam.*scount);
            
            precippine = precipp(count,:)';
            pmeanoptmne(count) = sum(betam.*precippine.*scount)/sum(betam.*scount);
        end
        
        if outcount == 1
            pmeanoptensemb = pmeanoptm;
        else
            pmeanoptensemb = pmeanoptensemb + pmeanoptm;
        end
        
        
        varensemblist(outcount,1) = var(pmeanoptm);
        
        
        pmeanoptm(count) = sum(betam.*precippi.*scount)/sum(betam.*scount);
        biasensemblist(outcount,1) = mean((pmeanoptmne - ismrhr).^2);
    end
    
    pmeanoptensemb = pmeanoptensemb / Nensemb;
    
    vardat(alphacount) = mean(varensemblist);
    biasdat(alphacount) = mean((pmeanoptensemb-ismrhr).^2);
    
    % calculate variance more accurately
    varsigma2 = (1/111*sum(precipp.*precipp,1) - alpha*(1/111*sum(precipp,1)).^2)';
    Srnd = Sr;  I = eye(357); idiag = find(I); Srnd(idiag) = 0;
    Evilist = mean(precipp,1)';
    
    varm1 = 1/alpha*sum(betam.^2.*varsigma2);
    varm2 = betam'*Srnd*betam;
    varm3 = 1/alpha*sigmae^2*(betam'*betam);
    varm4 = -2*(1-alpha)/alpha*(betam'*Evilist)*sum(betam.^2.*Evilist);
    varm5 = (1-alpha)/alpha*(betam'*Evilist)^2*(betam'*betam);
    
    varmodelmat(1,alphacount) = varm1;
    varmodelmat(2,alphacount) = varm2;
    varmodelmat(3,alphacount) = varm3;
    varmodelmat(4,alphacount) = varm4;
    varmodelmat(5,alphacount) = varm5;
    
    varmodel(alphacount) = varm1 + varm2 + varm3 + varm4 + varm5;
    varmodelquad(alphacount) = betam'*C*betam;
    
    biasmodelquad(alphacount) = biasm;
    
    % calculate bias
    bias2 = (1-alpha)/alpha*((betam'*Evilist)*(betam'*betam) - sum(betam.^2.*Evilist));
    
    biasmodela = 0;
    for yearcount = 1:111
        precippi = squeeze(precipp(yearcount,:))';
        bias1i = betam'*precippi - ismrhr(yearcount);
        biasmodela = biasmodela + (bias1i + bias2)^2;
    end
    
    biasmodela = biasmodela / 111;
    
    biasmodel(alphacount) = biasmodela;
    
    %% Assign values of beta
    betabmat(:,alphacount) = betab;
    betavmat(:,alphacount) = betav;
    betammat(:,alphacount) = betam;
    
end

%% Figure 2
dmup = (mean(precipp,1) - mean(ismrhr)).^2;

rangebetab = (max(betab) - min(betab)) / 20;

ibetab = betab > rangebetab;

% bias
figure(2), 

subplot(1,2,1),
[latu,lonu,betabplot] = vect2matgeophys(lat,lon,squeeze(betabmat(:,1)));

xlimplot = [min(lonu) max(lonu)];
ylimplot = [min(latu) max(latu)];

aspectr = [1 (max(latu) - min(latu)) / (max(lonu) - min(lonu)) 1];


caxisb = [0 max(max(betabmat))];...
    flagout1 = contourm(double(latu),double(lonu),betabplot'); hold on, plot(loncoast,latcoast), ...
    xlim(xlimplot), ylim(ylimplot), clim(caxisb), colorbar,...
    title(strcat('Min. Bias; \alpha=', num2str(alphalist(1),2))),...
    pbaspect(aspectr), set(gca,'box','on');


subplot(1,2,2),
plot(dmup(ibetab),betab(ibetab),'k+'),...
    xlabel('squared deviation from Ev (mm/day)^{2}'), ylabel('optimal weights \beta_i for minimizing Bias'),...
     pbaspect(aspectr), set(gca,'box','on'), grid on, xlim([-10 500])


% variance
figure(3)
subplot(2,2,1), ...
    [latu,lonu,betavplot] = vect2matgeophys(lat,lon,squeeze(betavmat(:,1)));
caxisv = [0 max(max(betavmat))];
flagout2 = contourm(double(latu),double(lonu),betavplot'); hold on, plot(loncoast,latcoast), ...
    xlim(xlimplot), ylim(ylimplot), clim(caxisv), colorbar,...
    title(strcat('Min. Variance; \alpha=', num2str(alphalist(1),2))),...
    pbaspect(aspectr), set(gca,'box','on');

figure(3),

subplot(2,2,3), ...
    [latu,lonu,betavplot] = vect2matgeophys(lat,lon,squeeze(betavmat(:,2)));
flagout5 = contourm(double(latu),double(lonu),betavplot'); hold on, plot(loncoast,latcoast), ...
    xlim(xlimplot), ylim(ylimplot), clim(caxisv), colorbar,...
    title(strcat('Min. Variance; \alpha=', num2str(alphalist(2)))),...
    pbaspect(aspectr), set(gca,'box','on');

% MSE
figure(3), subplot(2,2,2), ...
    [latu,lonu,betamplot] = vect2matgeophys(lat,lon,squeeze(betammat(:,1)));
caxism = [0 max(max(betammat))];
flagout3 = contourm(double(latu),double(lonu),betamplot'); hold on, plot(loncoast,latcoast), ...
    xlim(xlimplot), ylim(ylimplot), clim(caxism), colorbar,...
    title(strcat('Min. MSE; \alpha=', num2str(alphalist(1),2))),...
    pbaspect(aspectr), set(gca,'box','on');

figure(3), subplot(2,2,4), ...
    [latu,lonu,betamplot] = vect2matgeophys(lat,lon,squeeze(betammat(:,2)));
flagout6 = contourm(double(latu),double(lonu),betamplot'); hold on, plot(loncoast,latcoast), ...
    xlim(xlimplot), ylim(ylimplot), clim(caxism), colorbar,...
    title(strcat('Min. MSE; \alpha=', num2str(alphalist(2)))),...
    pbaspect(aspectr), set(gca,'box','on');