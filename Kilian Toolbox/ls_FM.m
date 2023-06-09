% 3 equations monetary VAR.M

clear; 

% Quarterly data from FRED
% Real Gross Domestic Product: Billions of Chained 2012 Dollars, SAAR, Q
% GDP Implicit Price Deflator in United States: SA, Q
% Effective federal funds rate, M
% Shadow rate from Altanta Fed
% 1960.II-2019.IV
load GDPC1.csv; drgdp=diff(log(GDPC1(:,2)))*400;
load GDPdeflator.csv; infl=diff(log(GDPdeflator(:,2)))*400;
load FEDFUNDS_SR.csv;

h = 40;                            % Maximum impulse response horizon
p=5;                               % VAR lag order

irate=[];
for i=1:3:length(FEDFUNDS_SR(:,2))
  irate=[irate; mean(FEDFUNDS_SR(i:i+2,2))];
end

end_sample = 192; %2008Q1
% end_sample = length(irate); %2019Q4
%% Same order as in the slides (Fed funds rate ordered third)
y=[drgdp(1:end_sample) infl(1:end_sample) irate(1:end_sample)]; [~,K]=size(y); 

[A,SIGMA,~,V,~]=olsvarc(y,p);   % VAR with intercept	

% Display LS model estimates
disp(V(1:K,1))
disp(A(1:K,1:K))
disp(A(1:K,K+1:2*K))
disp(A(1:K,2*K+1:3*K))
disp(A(1:K,3*K+1:4*K))
disp(SIGMA(1:K,1:K))

% VAR impulse response analysis
horizon=0:h;
B0inv = chol(SIGMA(1:K,1:K),'lower');
[IRF]=irfvar(A,B0inv,p,K,h);
IRF(7,:)=cumsum(IRF(7,:));
IRF(8,:)=cumsum(IRF(8,:));

B0 = inv(B0inv);
psi31_0 = -B0(3,1)/B0(3,3);
psi32_0 = -B0(3,2)/B0(3,3);
disp('Contemporaneous coefficients in the monetary policy rule')
disp('Baseline ordering:')
disp(strcat('$\psi_{31,0}=$)',num2str(psi31_0)))
disp(strcat('$\psi_{32,0}=$)',num2str(psi32_0)))

fig = figure(1);
subplot(1,3,1)
plot(horizon,IRF(7,:),'r-',horizon,zeros(size(horizon)),'linewidth',3);
title('GNP to Monetary Policy Shock','fontsize',18)
ylabel('Percent','fontsize',18)
xlabel('Quarters','fontsize',18)
axis([0 h -4 2])

subplot(1,3,2)
plot(horizon,IRF(8,:),'r-',horizon,zeros(size(horizon)),'linewidth',3);
title('Prices to Monetary Policy Shock','fontsize',18)
ylabel('Percent','fontsize',18)
xlabel('Quarters','fontsize',18)
axis([0 h -4 2])

subplot(1,3,3)
plot(horizon,IRF(9,:),'r-',horizon,zeros(size(horizon)),'linewidth',3);
title('Fed funds to Monetary Policy Shock','fontsize',18)
ylabel('Percent','fontsize',18)
xlabel('Quarters','fontsize',18)
axis([0 h -0.4 2])

%% Alternative order (Fed funds rate ordered first)
y=[irate drgdp infl]; [t,K]=size(y); 

[A,SIGMA,Uhat,V,X]=olsvarc(y,p);   % VAR with intercept	

% VAR impulse response analysis
horizon=0:h;
B0inv_alt = chol(SIGMA(1:K,1:K),'lower');
[IRF]=irfvar(A,B0inv_alt,p,K,h);
IRF(2,:)=cumsum(IRF(2,:));
IRF(3,:)=cumsum(IRF(3,:));

B0_alt = inv(B0inv_alt);
psi12_0 = -B0_alt(1,2)/B0_alt(1,1);
psi13_0 = -B0_alt(1,3)/B0_alt(1,1);

disp('--------------------------------------')
disp('Alternative ordering:')
disp(strcat('$\psi_{12,0}=$)',num2str(psi12_0)))
disp(strcat('$\psi_{13,0}=$)',num2str(psi13_0)))

fig = figure(2);
subplot(1,3,1)
plot(horizon,IRF(2,:),'r-',horizon,zeros(size(horizon)),'linewidth',3);
title('GNP to Monetary Policy Shock','fontsize',18)
ylabel('Percent','fontsize',18)
xlabel('Quarters','fontsize',18)
axis([0 h -2 2])

subplot(1,3,2)
plot(horizon,IRF(3,:),'r-',horizon,zeros(size(horizon)),'linewidth',3);
title('Prices to Monetary Policy Shock','fontsize',18)
ylabel('Percent','fontsize',18)
xlabel('Quarters','fontsize',18)
axis([0 h -2 4])

subplot(1,3,3)
plot(horizon,IRF(1,:),'r-',horizon,zeros(size(horizon)),'linewidth',3);
title('Fed funds to Monetary Policy Shock','fontsize',18)
ylabel('Percent','fontsize',18)
xlabel('Quarters','fontsize',18)
axis([0 h -0.4 2])

