function in = iresponse_proxy(in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'iresponse_proxy' computes the impulse response functions to shock
% identified via instrumental variables 
% Codes are written based on the Mertens Ravn (2013, AER) source codes,
% modified to  
%1) handle different length of VAR and factors
%2) Multiple instruments to explain the same shock

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y         = in.vars(in.p+1:end,:);
[T,n]     = size(Y);
[T_m,n_m] = size(in.proxies);

% number of proxies
k  = 1;
%Assuming proxies start at least p periods later
instrument  = in.proxies(1:end,:);
 
% Identification
%%%%%%%%%%%%%%%%%

% covariance of the reduced form shocks
in.Sigma_m = in.Sigma;

% Instrument on VAR residuals
Phib = [ones(T_m,1) instrument]\in.res(T-T_m-in.T_m_end+1:T-in.T_m_end,:);

% Fitted values of the identified shocks
% here is where thes prior on beta should enter , instead of Phib(:,1)
% m = beta e1 + v
uhat1           =   [ones(T_m,1) instrument]*Phib(:,1); 

% regress the fitted values on the other reduced form shocks
b21ib11_TSLS    =   [ones(T_m,1) uhat1]\in.res(T-T_m-in.T_m_end+1:T-in.T_m_end,2:end);  
b21ib11_TSLS    =   b21ib11_TSLS(2:end,:)';
b21ib11         =   b21ib11_TSLS;

% Identification of b11 and b12 from the covariance matrix of the VAR
Sig11   = in.Sigma_m(1:k,1:k);
Sig21   = in.Sigma_m(k+1:n,1:k);
Sig22   = in.Sigma_m(k+1:n,k+1:n);
ZZp     = b21ib11*Sig11*b21ib11'-(Sig21*b21ib11'+b21ib11*Sig21')+Sig22;
b12b12p = (Sig21- b21ib11*Sig11)'*(ZZp\(Sig21- b21ib11*Sig11));
b11b11p = Sig11-b12b12p;
b11     = sqrt(b11b11p);
in.b1   = [b11; b21ib11*b11];
in.Phib = Phib;

% Impulse Responses
%%%%%%%%%%%%%%%%%%%%
% initial shock: eps(1,1)=1
irs(in.p+1,:) = in.b1(:,1);

for jj = 2:in.irhor%+max(max(VAR.term_spreads_matur),max(VAR.real_rates_init+VAR.real_rates_matur-1))
    lvars = (irs(in.p+jj-1:-1:jj,:))';
    irs(in.p+jj,:) = lvars(:)'*in.Phi(1:in.p * n,:);     
end

in.irs   = irs(in.p+1:in.p+in.irhor,:);
in.uhat1 = uhat1;


if in.compute_F_stat == 1
    % F-test
    %%%%%%%%%%%%%%%%%%%%%%%%%
    XX_m      = [ones(T_m,1) instrument];
    Res_m     = in.res( T-T_m-in.T_m_end+1 : T-in.T_m_end , :)- XX_m * Phib;
    Res_const = in.res( T-T_m-in.T_m_end+1 : T-in.T_m_end,1) - ...
        ones(T_m,1)*(ones(T_m,1)\in.res(T-T_m-in.T_m_end+1 : T-in.T_m_end,1));
    
    SST_m      = Res_const'*Res_const;
    SSE_m      = Res_m(:,1)'*Res_m(:,1);
    in.F_m     = ((SST_m-SSE_m)/n_m)/ ...
        (SSE_m/(length(in.res(T-T_m-in.T_m_end+1 : T-in.T_m_end,1))-(n_m+1)));
    in.R2_m    = (1-SSE_m/(SST_m));
    in.R2adj_m = in.R2_m-...
        (n_m / ((length(in.res(T-T_m-in.T_m_end+1 : T-in.T_m_end,1)) -(n_m+1))))*(1-in.R2_m);
    
    % Calculate robust standard errors
    SS_m = zeros(n_m+1,n_m+1);
    for ii = 1 : T_m
        SS_m = SS_m+1/T_m * XX_m(ii,:)'*XX_m(ii,:)*Res_m(ii,1)^2;
    end
    Avarb_m    = inv(1/T_m * XX_m'*XX_m)*SS_m *inv(1/T_m * XX_m'*XX_m);
    RR_m       = [zeros(n_m,1) eye(n_m)];
    WW_m       = T_m*(RR_m*Phib(:,1))' * inv(RR_m*Avarb_m*RR_m') * (RR_m*Phib(:,1));
    in.F_m_rob = WW_m/n_m;
end