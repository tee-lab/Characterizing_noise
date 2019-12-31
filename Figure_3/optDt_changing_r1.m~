%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by Jitesh Jhawar somewhere in 2017. This code simulates a            %
% stochastic process using the standard algorithm presented by Sir D Gillespie,% 
% 1976. The process simulates the evolution of temporal dynamics of consensus  %
% in a population that contains individuals in two states. These individuals   %
% switch between their states via. different kinds of reactions/interactions.  %
% Following this, we also derive the underying SDE from data and compare the   %
% reconstructed functions with the expected ones.                              % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters 
N = 50; 
s = flip(1./(2*(10:10:200))); 
r2 = 1; 
r3 = 0; 
r4 = 0; 
Tint = 50;
opt_Dt = zeros(length(s),1);
est_tau = zeros(length(s),1);

if r3 == 0
    address = '~/Documents/RSPhilTran/submit/Characterizing_noise/pairwise/optimalDt/varyingN/';
else
    address = '~/Documents/RSPhilTran/submit/Characterizing_noise/ternary/optimalDt/varyingN/';
end

Dt = 1:1:200;
% Dt = [1,Dt];
Drift = cell(length(Dt),1);
Diffusion = Drift; Diffusion_mod = Drift;
dist_drift = zeros(length(Dt),1); dist_diff = dist_drift;

for j = 1:length(s)
    r1 = s;
    [tSample,S] = GS_runner1D(N,r1,r2,r3,r4,Tint);
    
    % calculated drift diffusion for different Dt for same time series and store them in cell array 
    for i = 1:length(Dt)
        T_int = 1;
        [Diffusion_temp,Diffusion_mod_temp,Drift_temp,op] = driftAndDiffusion_const_time(S,T_int,Dt(i));
        Drift{i,1} = Drift_temp;
        Diffusion{i,1} = Diffusion_temp;
        Diffusion_mod{i,1} = Diffusion_mod_temp;
    end
    
    Diffusion_mod{1,1}(Diffusion_mod{1,1}>10) = nan;
    Drift{1,1}(Drift{1,1}>10) = nan;
    exp_drift = -2*r1*op - 2*r4*op + (r3/2)*op.*(1-op.^2); %expected
    exp_diff = (2/sqrt(N))*sqrt(r1+r4+((2*r2+r3-2*r4)*(1-op.^2)/4));
    exp_diff = exp_diff.^2;
    op = op';
    
    %calculate distance between the expected and the extracted
    for i = 1:length(Dt)
        Drift_temp = Drift{i,1};
        dist_drift(i,1) = sqrt(nanmean((Drift_temp - exp_drift').^2))/nanmean(abs(exp_drift));
        Diffusion_temp = Diffusion_mod{i,1};
        idx = find(~isnan(Diffusion_temp));
        dist_diff(i,1) = sqrt(nanmean((Diffusion_temp - exp_diff').^2))/nanmean(exp_diff);
    end
    
    %Find optimal Dt and store
    [M,I] = min(dist_drift);
    opt_Dt(j) = I;
    
    %calculate correlation time for single time series and store
     t_lag = 1000;
     S(isnan(S)) = 0;
     acf = autocorr(S',t_lag);
     t_lag = (1:t_lag+1);
     t_lag = t_lag';
     idx = find(~isnan(acf));
     f = fit(t_lag(idx),acf(idx),'exp1');
     val = coeffvalues(f);
     b = val(1,2);
     a = val(1,1);
     ct = -1*ceil(1/b);
 %     t_acf = tSample(1:length(t_lag));
     est_tau(j) = tSample(ct);
end
save([address 'opt_Dt_changing_N'],'opt_Dt');
save([address 'est_tau_interp_r1'],'est_tau');



