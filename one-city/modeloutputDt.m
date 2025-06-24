function ymat = modeloutputDt(matparam)
% computes the model output for each line of parameters in matparam
% Input parameters:
% - matparam: matrix of parameters, each row gives 1 set of parameters
% Output parameter:
% - ymat: matrix of output variables

[rowno,colno] = size(matparam);
ymat = zeros(rowno,20);

% constants
t0 = 0; tf = 140; nstep = 140; dt = (tf-t0)/nstep; tvec = (0:nstep)*dt;
matY = zeros(11,nstep+1); popsize = 5000000; infect0 = 5000; av = 0;
cvec = [popsize infect0 av];

for indrow = 1:rowno
    % extract the parameters [I0 beta r delta kappa gamma alpha]
    pvec = matparam(indrow,:);
    y0 = [(popsize-infect0) 0 infect0 0 0 0 0 0 0 0 infect0]';
    matY(:,1) = y0;
    for ind = 1:nstep
        matY(:,ind+1) = eEuler(tvec(ind),matY(:,ind),pvec,cvec,dt);
    end
    ymat(indrow,:) = matY(5,8:7:end);
end