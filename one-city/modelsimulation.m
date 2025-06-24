% generate plots of the model

% constants
t0 = 0; tf = 140; nstep = 140; dt = (tf-t0)/nstep; tvec = (0:nstep)*dt;
matY = zeros(11,nstep+1); popsize = 5000000; infect0 = 5000; av = 0;
cvec = [popsize infect0 av];

% values of [beta r delta kappa gamma alpha]
R0 = 2.41; kappa1 = 5.5; gamma1 = 6.5; alph = 0.09;
arr = 0.75; del = 27887;
pvec = [(R0/gamma1) arr del (1/kappa1) (1/gamma1) alph];
yinit = [(popsize - infect0) 0 infect0 0 0 0 0 0 0 0 infect0]';
matY(:,1) = yinit;

for ind = 1:nstep
    matY(:,ind+1) = eEuler(tvec(ind),matY(:,ind),pvec,cvec,dt);
end

totalS = matY(1,:) + matY(6,:);
totalE = matY(2,:) + matY(7,:);
totalI = matY(3,:) + matY(8,:);
totalR = matY(4,:) + matY(9,:);
totalD = matY(5,:) + matY(10,:);
total = totalS + totalE + totalI + totalR + totalD;
% generate plots
figure
hold on
grid on
xlim([0 140])
ylim([0 2.5e6])
plot(tvec,[matY(11,:)' totalD'])
legend("infected count", "death count")
title('Cumulative Infected and Cumulative Death Count Over Time')
xlabel('time')
ylabel('count')
hold off