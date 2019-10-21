function olr = integ_fr(Ts, T, tau0, p)

% constants
sigb = 5.67e-8;

%normalized pressure
pn   = p./max(p);

% define tau and get tau difference
tau = tau0.*pn.^2;
dtau = -diff(tau);

% get T and tau at mid points
Tmid = (T(1:end-1) + T(2:end))/2;
taumid = (tau(1:end-1) + tau(2:end))/2;

% compute outgoing longwave radiation
olr = sigb*Ts^4*exp(-tau0)+sum(sigb*Tmid.^4.*exp(-taumid).*dtau);