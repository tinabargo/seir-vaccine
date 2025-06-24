% Checking monotonicity of input variables with respect to model outputs

nrow = 100; ncol = 6;
temp = ones(nrow+1,1);
tind = 7:7:140;
tlength = length(tind);

% parameter values for [beta r delta kappa gamma alpha]
R0 = 2.41; kappa1 = 5.5; gamma1 = 6.5; alph = 0.09;
arr = 0.75; del = 27887;
pbase = [(R0/gamma1) arr del (1/kappa1) (1/gamma1) alph];

% generate sample values of the parameters using Latin Hypercube Sampling
paramtemp = lhsu(nrow);
paramsource = [pbase; paramtemp];

% input labels [beta r delta kappa gamma alpha]
labelx = ["\beta","r","\delta","\kappa","\gamma","\alpha"];

for ind = 1:ncol
    % generate the matrix of parameters
    Mparam = temp*pbase;

    % replace one column with sample values
    paramvec = paramsource(:,ind);
    Mparam(:,ind) = paramvec;

    modelout = modeloutputWt(Mparam);

    for ts = 1:tlength
        figure
        plot(paramvec,modelout(:,1),'.')
        title('W vs ' + labelx(ind) + ' at t = ' + tind(ts))
        xlabel(labelx(ind))
        ylabel('W')
    end

end