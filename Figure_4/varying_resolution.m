N = 200; r1 = 0.01; r2 = 1; r3 = 0; r4 = 0; res = 1:5:500;
Diffusion_distance = zeros(length(res),1);
Drift_distance = zeros(length(res),1);
if r3 == 0
    address = ['~/Documents/RSPhilTran/submit/Characterizing_noise/pairwise/varying_resolution/N_' num2str(N) '/'];
else
    address = ['~/Documents/RSPhilTran/submit/Characterizing_noise/ternary/varying_resolution/N_' num2str(N) '/'];
end
Dt = 1;
% Dt = [1,Dt];
Drift = cell(length(Dt),1);
Diffusion = Drift; Diffusion_mod = Drift;
dist_drift = zeros(length(Dt),1); dist_diff = dist_drift;
for j = 1:length(res)
    Tint = res(j);
    for i = 1:length(Dt)
        [tSample,S] = GS_runner1D(N,r1,r2,r3,r4,Tint);
        T_int = tSample(end)/(length(tSample));
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
    for i = 1:length(Dt)
        Drift_temp = Drift{i,1};
        dist_drift(i,1) = sqrt(nanmean((Drift_temp - exp_drift').^2))/nanmean(abs(exp_drift));
        Diffusion_temp = Diffusion_mod{i,1};
        idx = find(~isnan(Diffusion_temp));
        dist_diff(i,1) = sqrt(nanmean((Diffusion_temp - exp_diff').^2))/nanmean(exp_diff);
    end
    Diffusion_distance(j) = dist_diff;
    Drift_distance(j) = dist_drift;
end
save([address 'Diffusion_distance'],'Diffusion_distance');
save([address 'Drift_distance'],'Drift_distance');
save([address 'resolution'],'res');
