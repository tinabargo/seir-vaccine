% Write the LHS matrix and output matrices

% generate the LHS matrix
% [beta r delta kappa gamma alpha]
R0 = 2.41; kappa1 = 5.5; gamma1 = 6.5; alph = 0.09; arr = 0.75;
del = 27887;
pbase = [(R0/gamma1) arr del (1/kappa1) (1/gamma1) alph];
nsamp = 700;
matLHS = lhsu(nsamp);
matLHS = [pbase; matLHS]; % add the base values in 1st row
% save LHS matrix for replicability
writematrix(matLHS,'LHSmatrix.csv')

% generate the output matrix
matOW = modeloutputWt(matLHS);
matOD = modeloutputDt(matLHS);

% write output matrices on csv files
writematrix(matOD,'DOutput.csv')
writematrix(matOW, 'WOutput.csv')