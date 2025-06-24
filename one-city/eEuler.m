function y = eEuler(t,x,paramvec,constvec,h)
% one-step explicit Euler scheme
% Input parameters:
% - t: time (real)
% - x: [S E I R D Sv Ev Iv Rv Dv W]'
% - paramvec: parameter source [I0 beta r delta kappa gamma alpha]
%   --> remove I0 because of monotonicity issues!
% - h: timestep
% Output parameter:
% - y: values of the state variables at the next timestep (vector)

y = x + h*onecity(t,x,paramvec,constvec);