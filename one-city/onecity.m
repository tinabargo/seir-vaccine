function dydt = onecity(t,y,paramvec)
% right-hand side of the ODE system for the one-city model
% Input parameters:
% - t: time (real)
% - y: [S E I R D Sv Ev Iv Rv Dv W]'
% - paramvec: parameter source [I0 beta r delta kappa gamma alpha]
%   --> remove I0 because of monotonicity issues!
% Output parameter:
% - dydt: right-hand side of the ODE system (vector)

% constants
N = 5000000; % total population
%gamma = 0.05; % recovery rate
I0 = 5000; % initial number of infected individuals
alphav = 0; % death rate for vaccinated individuals

% extract variables from y
S = y(1); E = y(2); I = y(3); R = y(4); D = y(5);
Sv = y(6); Ev = y(7); Iv = y(8); Rv = y(9); Dv = y(10);
W = y(11);

% extract parameters from paramvec
%I0 = paramvec(1); 
beta = paramvec(1); r = paramvec(2); delta = paramvec(3);
kappa = paramvec(4); gamma = paramvec(5); alpha = paramvec(6);

% solve for adjusted transmission rates
thetav = beta/3;
betav = (1-r)*beta;
psiv = betav/3;

% Solve for deltat
Sstar = S - (S/N)*(beta*I + thetav*Iv);
deltat = max(min([delta Sstar]),0);
%disp([t deltat])

Sprime = -(S/N)*(beta*I + thetav*Iv) - deltat;
Eprime = (S/N)*(beta*I + thetav*Iv) - kappa*E;
Iprime = kappa*E - (gamma + alpha)*I;
Rprime = gamma*I;
Dprime = alpha*I;
Svprime = -(Sv/N)*(betav*I + psiv*Iv) + deltat;
Evprime = (Sv/N)*(betav*I + psiv*Iv) - kappa*Ev;
Ivprime = kappa*Ev - (gamma + alphav)*Iv;
Rvprime = gamma*Iv;
Dvprime = alphav*Iv;
W = kappa*(E + Ev);
dydt = [Sprime Eprime Iprime Rprime Dprime Svprime Evprime Ivprime ...
    Rvprime Dvprime W]';