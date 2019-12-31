function[tSample,S_res] = GS_runner1D(N,r1,r2,r3,r4,Tint)
%initialize species, reaction propensities
mu  = 8;
C = zeros(mu,1);
A = C;
%reaction rates
 
C(1:2) = r1;
C(3:4) = r2;
C(5:6) = r3;
C(7:8) = r4;

% Tint = 50;%ceil((r2+r3+r4)/((N^0.5)*r1));
Tend = ceil(Tint*500000);
steps = floor(Tend);
rel = 1;    %realizations/repetitions
%Species number storing matrices
S = zeros(steps,rel);
S1 = S; S2 = S1; sum_all = S1; tSample = S1;
for iter = 1:rel
    X1 = ceil(rand(1,1)*N); X2 = N-X1;
    n = 0;
    T = 0;
    Tprint = 0.01;
    while (T < Tend)
        %         loop = loop + 1
        %Reactions
        A(1) = C(1)*X2/N;  %1
        A(2) = C(2)*X1/N;  %2
        A(3) = C(3)*X1*X2/(N^2);   %3
        A(4) = C(4)*X1*X2/(N^2);   %4
        A(5) = C(5)*(X1^2)*X2/(N^3);   %5
        A(6) = C(6)*X1*(X2^2)/(N^3);   %6
        A(7) = C(7)*(X1^2)/(N^2);    %7
        A(8) = C(8)*(X2^2)/(N^2);    %8
        
        A0 = sum(A);
        A0 = A0(1,1);
        
        %Generating random numbers
        R(1) = rand(1,1);
        
        T = T + (log(1/R(1)))/A0;
        if T > Tprint
            n = n  + 1;
            S1(n,iter) = X1; S2(n,iter) = X2; sum_all(n,iter) = (X1+X2)/N;
            S(n,iter) = (X1 - X2)/N;
            Tprint = Tprint + Tint;
            tSample(n,iter) = T;
            
        end
        R(2) = rand(1,1);
        R2A0 = R(2)*A0;
        Sum = 0;
        nu = 1;
        while Sum <= R2A0 %&& nu <= 24
            mu = nu;
            Sum = Sum + A(nu);
            nu = nu + 1;
        end
        if mu == 1 || mu == 4 || mu == 5 || mu == 8
            X2 = X2 - 1;
            X1 = X1 + 1;
        elseif mu == 2 || mu == 3 || mu == 6 || mu == 7
            X1 = X1 - 1;
            X2 = X2 + 1;
        end
    end
%         S(:,iter) = filter(ones(1,10)/10,1,S(:,iter));
%     S(find(tSample == 0)) = nan;
    
end

% S(S==0) = nan;
% tSample(tSample==0) = nan;
tSample = tSample/N;
S = S(find(tSample ~= 0));
tSample = tSample(find(tSample ~= 0));
% xlim([0,1])

S1(S1==0) = nan;
S2(S2==0) = nan;
S_res = interp1(tSample,S,min(tSample):1:max(tSample),'previous');
end
