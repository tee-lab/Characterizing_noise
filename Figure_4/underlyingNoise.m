function[noise]  = underlyingNoise(X,Dt,point,inc,Tint)
%% Drift and Diffusion
% inc = 0.01;
op = -1:inc:1.0;
drift = zeros(length(X),1);

avgDrift = zeros(length(op),1);
% dt = 0.01;
% Dt = 10;
% X(X == 0) = nan;
%     X = filter(ones(1,10)/10,1,X(:,1));
%loop for calculating drift and diffusion coeff
for i = 1:(length(X)-Dt)
    
    drift(i) = (X(i+Dt) - X(i))/(Tint*Dt);
end
%
bin = point;

n = 1;
while (bin < point+inc)       %binning loop
    j = 1;
    binDri = zeros(length(X),1);
    for i = 1:(length(X)-1)
        if (X(i) < (bin + inc) && X(i) >= bin)
            binDri(j) = drift(i);
            j = j+1;
        end
    end
    binDri(j:length(binDri)) = nan;
    avgDrift(n) = nanmean(binDri);
    bin = bin + inc;
    n = n + 1;
end

% avgDrift(isnan(avgDrift)) = 0;

bin = point;

n = 1;
while (bin < point+inc)       %binning loop
    c = 0;
    j = 1;
    binDif_mod = zeros(length(X),1);
    for i = 1:(length(X)-1)
        if (X(i) < (bin + inc) && X(i) >= bin)
            binDif_mod(j) = ((X(i+1) - X(i)) - (Tint*Dt)*avgDrift(find(op==point)))/sqrt(Tint);
            j = j+1;
        end
    end
    binDif_mod(j:length(binDif_mod)) = nan;
    idx = find(~isnan(binDif_mod));
    noise = binDif_mod(idx);
    bin = bin + inc;
    n = n + 1;
end
end