function ymat = modeloutputWt(matparam)
% computes the model output for each line of parameters in matparam
% Input parameters:
% - matparam: matrix of parameters, each row gives 1 set of parameters
% Output parameter:
% - ymat: matrix of output variables

[rowno,colno] = size(matparam);
ymat = zeros(rowno,20);

% constants
t0 = 0; tf = 140; nstep = 140; dt = (tf-t0)/nstep; tvec = (0:nstep)*dt;
matY = zeros(11,nstep+1); N = 5000000; I0 = 5000;

for indrow = 1:rowno
    % extract the parameters [I0 beta r delta kappa gamma alpha]
    pvec = matparam(indrow,:);
    y0 = [(N-I0) 0 I0 0 0 0 0 0 0 0 I0]';
    matY(:,1) = y0;
    for ind = 1:nstep
        matY(:,ind+1) = eEuler(tvec(ind),matY(:,ind),pvec,dt);
    end
    ymat(indrow,:) = matY(11,8:7:end);
end