function s=lhsu(nsample)
% Latin Hypercube Sampling
% Input:
%   nsample : no. of samples
% Output:
%   s       : random sample (nsample,6)
%   parameter vectors are of the form [beta r delta kappa gamma alpha]
%   Budiman (2003)

nvar = 6;
s = zeros(nsample,nvar);

% R0 follows the truncated normal distribution, beta = R0/6.5
% mean = 2.41
% sd = 0.03826531
% range: [0,+\infty)
r0mean = 2.41; r0sd = 0.03826531;
idx = randperm(nsample); ran = rand(nsample,1);
pmin = cdf('Normal',0,r0mean,r0sd);
P = pmin + ((1-pmin)/nsample)*(idx'-ran);
r0vec = icdf('Normal',P,r0mean,r0sd); bvec = r0vec/6.5;
s(:,1) = bvec;

% r follows the uniform distribution:
% range: [0.5,1]
rmin = 0.5; rmax = 1;
idx = randperm(nsample); ran = rand(nsample,1);
P = (idx'-ran)/nsample; rvec = rmin + (rmax-rmin)*P;
s(:,2) = rvec;

% delta follows the beta distribution:
% a = 0.7340199944487569
% b = 112.5201094419424
% scale = 4302834.680007711
% offset = -1.197495370221827e-23
aconst = 0.7340199944487569; bconst = 112.5201094419424;
scale = 4302834.680007711; offset = -1.197495370221827e-23;
idx = randperm(nsample); ran = rand(nsample,1);
P = (idx'-ran)/nsample;
dvec = icdf('Beta',P,aconst,bconst); deltavec = scale*dvec + offset;
s(:,3) = deltavec;

% kappa1 follows the truncated normal distribution, kappa = 1/kappa1
% mean = 5.5
% sd = 0.97
% range: [1,+\infty)
kmean = 5.5; ksd = 0.97;
idx = randperm(nsample); ran = rand(nsample,1);
pmin = cdf('Normal',1,kmean,ksd);
P = pmin + ((1-pmin)/nsample)*(idx'-ran);
kvec = icdf('Normal',P,kmean,ksd); kappavec = 1./kvec;
s(:,4) = kappavec;

% gamma follows the truncated normal distrubution
% mean = 6.5
% sd = 0.77
% range: [4,+\infty)
gmean = 6.5; gsd = 0.77;
idx = randperm(nsample); ran = rand(nsample,1);
pmin = cdf('Normal',4,gmean,gsd);
P = pmin + ((1-pmin)/nsample)*(idx'-ran);
gvec = icdf('Normal',P,gmean,gsd); gammavec = 1./gvec;
s(:,5) = gammavec;

% alpha follows the uniform distribution
% range: [0,0.44]
amin = 0; amax = 0.44;
idx = randperm(nsample); ran = rand(nsample,1);
P = (idx'-ran)/nsample; avec = amin + (amax-amin)*P;
s(:,6) = avec;