function[avgDiff,avgDiff_mod,avgDrift,op]  = driftAndDiffusion_const_time(X,t_int,Dt)
%% Drift and Diffusion
%X = (vel_x.^2 + vel_y.^2).^0.5;
inc = 0.01;
op = -1:inc:1.0;
diff = zeros(length(X),1);
drift = zeros(length(X),1);

avgDrift = zeros(length(op),1);
avgDiff = zeros(length(op),1);
avgDiff_mod = zeros(length(op),1);
stdDrift = zeros(length(op),1);
stdDif = zeros(length(op),1);
stdDif_mod = zeros(length(op),1);
% dt = 0.01;
% Dt = 10;
% X(X == 0) = nan;
%     X = filter(ones(1,10)/10,1,X(:,1));
%loop for calculating drift and diffusion coeff
for i = 1:(length(X)-Dt)
    
    drift(i) = (X(i+Dt) - X(i))/(t_int*Dt);
    diff(i) = ((X(i+Dt) - X(i))^2)/(t_int*Dt);
end
%
bin = min(op);

n = 1;
while (bin < max(op))       %binning loop
    c = 0;
    j = 1;
    binDri = zeros(length(X),1);
    binDif = zeros(length(X),1);
    for i = 1:(length(X)-1)
        if (X(i) < (bin + inc) && X(i) >= bin)
            binDif(j) = diff(i);
            binDri(j) = drift(i);
            j = j+1;
        end
    end
    binDri(binDri==0) = nan;
    binDif(binDif==0) = nan;
    avgDrift(n) = nanmean(binDri);
    avgDiff(n) = nanmean(binDif);
    stdDrift(n) = nanstd(binDri);
    stdDif(n) = nanstd(binDif);
    bin = bin + inc;
    n = n + 1;
end

avgDrift(avgDrift==0) = nan;
avgDiff(avgDiff==0) = nan;

bin = min(op);

n = 1;
while (bin < max(op))       %binning loop
    c = 0;
    j = 1;
    binDif_mod = zeros(length(X),1);
    for i = 1:(length(X)-Dt)
        if (X(i) < (bin + inc) && X(i) >= bin)
            binDif_mod(j) = ((X(i+Dt) - X(i)) - (t_int*Dt)*avgDrift(n))^2/(t_int*Dt);
            j = j+1;
        end
    end
    binDif_mod(binDif_mod==0) = nan;
    avgDiff_mod(n) = nanmean(binDif_mod);
    stdDif_mod(n) = nanstd(binDif_mod);
    bin = bin + inc;
    n = n + 1;
end

avgDiff_mod(avgDiff_mod==0) = nan;

% % y = 1-2*op;
% % z = 1.0*(1 - op.^2);
% figure
% scatter(op,avgDrift)
% xlim([min(op),max(op)])
% xlabel('Polarization','FontWeight','bold')
% ylabel('Drift','FontWeight','bold')
% hline = refline([0 0]);
% hline.Color = 'r';
% % % hold on
% % % plot(op,y)
% % % legend('Data derived','f(U) = 1 - 2U')
% figure
% scatter(op,avgDiff)
% xlim([min(op),max(op)])
% xlabel('Polarization','FontWeight','bold')
% ylabel('Diffusion','FontWeight','bold')
% % hold on
% % plot(op,z)
% % legend('Data derived','g^2(U) = (1 - U^2)')


% figure
% scatter(X,drift,'.')
% xlim([min(op),1])
% xlabel('Polarization','FontWeight','bold')
% ylabel('Drift','FontWeight','bold')
% hline = refline([0 0]);
% hline.Color = 'r';
% figure
% plot(X,diff,'.')
% xlim([min(op),1])
% xlabel('Polarization','FontWeight','bold')
% ylabel('Diffusion','FontWeight','bold')

% posStd_Drift = avgDrift + stdDrift;
% negStd_Drift = avgDrift - stdDrift;
% posStd_Dif = avgDiff + stdDif;
% negStd_Dif = avgDiff - stdDif;
%
% drift(isnan(drift)) = 0;
% diff(isnan(diff)) = 0;
% X(isnan(X)) = 0;
% address = '/home/jitesh/Documents/MATLAB/DriftDiffusion/60indi/pooled';
% dlmwrite([address '/raw_drift.txt'],[X,drift],'delimiter','\t');
% dlmwrite([address '/raw_diffusion.txt'],[X,diff],'delimiter','\t');
% dlmwrite([address '/drift.txt'],[op.',avgDrift,posStd_Drift,negStd_Drift],'delimiter','\t');
% dlmwrite([address '/diffusion.txt'],[op.',avgDiff,posStd_Dif,negStd_Dif],'delimiter','\t');

end
