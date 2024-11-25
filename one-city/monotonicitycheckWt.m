% Checking monotonicity of input variables with respect to model outputs

nrow = 100; ncol = 6;
temp = ones(nrow+1,1);
tind = 7:7:140;
tlength = length(tind);

% parameter values for [beta r delta kappa gamma alpha]
%pbase = [1.3 0.5 10000 0.2 0.05 0.05];
%pmin = [0.01 0 100 0.01 0.05 0.01];
%pmax = [1.3 1 100000 0.6 0.25 0.3];

pbase = [1.3 0.5 10000 0.2 0.05 0.5];
pmin = [0.01 0 100 0.01 0.05 0.4];
pmax = [1.3 1 100000 0.6 0.25 0.75];

% input labels [beta r delta kappa gamma alpha]
labelx = ["\beta","r","\delta","\kappa","\gamma","\alpha"];

for ind = 1:ncol
    % generate the matrix of parameters
    Mparam = temp*pbase;

    % replace one column with sample values
    param0 = pmin(ind); paramend = pmax(ind);
    dparam = (paramend - param0)/nrow;
    paramvec = param0 + dparam*(0:nrow)';
    Mparam(:,ind) = paramvec;

    modelout = modeloutputWt(Mparam);

    for ts = 1:tlength
        figure
        plot(paramvec,modelout(:,1))
        title('W vs ' + labelx(ind) + ' at t = ' + tind(ts))
        xlabel(labelx(ind))
        ylabel('W')
    end

end