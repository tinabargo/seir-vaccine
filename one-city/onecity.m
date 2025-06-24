function dydt = onecity(t,y,paramvec,constvec)
% right-hand side of the ODE system for the one-city model
% Input parameters:
% - t: time (real)
% - y: [S E I R D Sv Ev Iv Rv Dv W]'
% - paramvec: parameter source [beta r delta kappa gamma alpha]
% - constvec: constant vector [N, I0, alphav]
% Output parameter:
% - dydt: right-hand side of the ODE system (vector)

% extract constants from cvec
N = constvec(1); I0 = constvec(2); alphav = constvec(3);

% extract variables from y
S = y(1); E = y(2); I = y(3); R = y(4); D = y(5);
Sv = y(6); Ev = y(7); Iv = y(8); Rv = y(9); Dv = y(10);
W = y(11);

% extract parameters from paramvec
beta = paramvec(1); r = paramvec(2); delta = paramvec(3);
kappa = paramvec(4); gamma = paramvec(5); alpha = paramvec(6);

% solve for adjusted transmission rates
thetav = beta/3;
betav = (1-r)*beta;
psiv = betav/3;

% Solve for deltat
Stemp = (S/N)*(beta*I + thetav*Iv); Svtemp = (Sv/N)*(betav*I + psiv*Iv);
Sstar = S - Stemp;
% Ensures that deltat >= 0
deltat = max(min([delta Sstar]),0);
% Ensures that S and Sv are nonnegative
if S >= Stemp
    Sout = Stemp;
else
    Sout = S;
end
if Sv >= Svtemp
    Svout = Svtemp;
else
    Svout = Sv;
end

Sprime = -Sout - deltat;
Eprime = Sout - kappa*E;
Iprime = kappa*E - gamma*I;
Rprime = (1-alpha)*gamma*I;
Dprime = alpha*gamma*I;
Svprime = -Svout + deltat;
Evprime = Svout - kappa*Ev;
Ivprime = kappa*Ev - gamma*Iv;
Rvprime = (1-alphav)*gamma*Iv;
Dvprime = alphav*gamma*Iv;
W = kappa*(E + Ev);
dydt = [Sprime Eprime Iprime Rprime Dprime Svprime Evprime Ivprime ...
    Rvprime Dvprime W]';