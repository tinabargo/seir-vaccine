% Sampling via LHS
% Sensitivity Analysis via PRCC

% matrix where we store prcc and p-values
% rows: [beta r delta kappa gamma alpha]
% columns: W(7) to W(140) and D(7) to D(140)
matWPRCC = zeros(6,20); matWPval = zeros(6,20);
matDPRCC = zeros(6,20); matDPval = zeros(6,20);
rn = {'beta','r','delta','kappa','gamma','alpha'};
cn = {'W(7)','W(14)','W(21)','W(28)','W(35)','W(42)','W(49)','W(56)',...
    'W(63)','W(70)','W(77)','W(84)','W(91)','W(98)','W(105)','W(112)',...
    'W(119)','W(126)','W(133)','W(140)','D(7)','D(14)','D(21)','D(28)',...
    'D(35)','D(42)','D(49)','D(56)','D(63)','D(70)','D(77)','D(84)',...
    'D(91)','D(98)','D(105)','D(112)','D(119)','D(126)','D133)','D(140)'};

% reading the stored LHS matrix
matLHS = csvread('LHSmatrix.csv');

% reading the stored output matrices
matOW = csvread('WOutput.csv');
matOD = csvread('DOutput.csv');

% rank-transform the LHS matrix
[LHSr,LHSc] = size(matLHS);
matLHSrank = zeros(LHSr,LHSc);
for ind = 1:LHSc
    matLHSrank(:,ind) = tiedrank(matLHS(:,ind));
end

% rank-transform the output matrices
[matOr,matOc] = size(matOW);
matOWrank = zeros(matOr,matOc); matODrank = zeros(matOr,matOc);
for ind = 1:matOc
    matOWrank(:,ind) = tiedrank(matOW(:,ind));
    matODrank(:,ind) = tiedrank(matOD(:,ind));
end

% input labels [beta r delta kappa gamma alpha]
labelx = ["\beta","r","\delta","\kappa","\gamma","\alpha"];

for indx = 1:LHSc
    % matrix where we store the residuals
    matWRes = zeros(matOr,matOc+1); matDRes = zeros(matOr,matOc+1);

    % prepare matrices for multiple linear regression
    Xtemp = [matLHSrank(:,1:indx-1) matLHSrank(:,indx+1:end)];
    Xind = matLHSrank(:,indx);

    % for computing the model values
    Xnew = [ones(LHSr,1) Xtemp];
    
    % multiple linear regression
    mdlx = fitlm(Xtemp,Xind);
    % compute the residual in Xind
    coeffs = table2array(mdlx.Coefficients);
    yhat = Xnew*coeffs(:,1);
    Xres = Xind - yhat;
    matWRes(:,1) = Xres; matDRes(:,1) = Xres;
    
    for indy = 1:matOc
        YWind = matOWrank(:,indy); YDind = matODrank(:,indy);
        mdlyw = fitlm(Xtemp,YWind); mdlyd = fitlm(Xtemp,YDind);
        % compute the residual in YWind and YDind
        coeffsw = table2array(mdlyw.Coefficients);
        coeffsd = table2array(mdlyd.Coefficients);
        ywhat = Xnew*coeffsw(:,1); ydhat = Xnew*coeffsd(:,1);
        YWres = YWind - ywhat; YDres = YDind - ydhat;
        matWRes(:,indy+1) = YWres; matDRes(:,indy+1) = YDres;
    end

    % linear regression on the residuals
    [ccw,pvw] = corrcoef(matWRes); [ccd,pvd] = corrcoef(matDRes);
    matWPRCC(indx,:) = ccw(1,2:end); matDPRCC(indx,:) = ccd(1,2:end);
    matWPval(indx,:) = pvw(1,2:end); matDPval(indx,:) = pvd(1,2:end);

    % generate plot of pval of W
    figure
    hold on
    plot(7:7:140,pvw(1,2:end),'LineWidth',2.0)
    yline(0.05,'--r')
    xticks(0:7:140)
    grid on
    title('p-values of W vs ' + labelx(indx) + ' over time')
    xlabel('Day')
    ylabel('p-value')
    hold off
    
    % generate plot of PRCC of W
    figure
    plot(7:7:140,ccw(1,2:end),'LineWidth',2.0)
    ylim([-1 1])
    xticks(0:7:140)
    yline(-0.07,'--r')
    yline(0.07,'--r')
    grid on
    title('PRCC of W vs ' + labelx(indx) + ' over time')
    xlabel('Day')
    ylabel('PRCC')

    % generate plot of pval of D
    figure
    hold on
    plot(7:7:140,pvd(1,2:end),'LineWidth',2.0)
    yline(0.05,'--r')
    xticks(0:7:140)
    grid on
    title('p-values of D vs ' + labelx(indx) + ' over time')
    xlabel('Day')
    ylabel('p-value')
    hold off

    % generate plot of PRCC of D
    figure
    % patch([7 140 140 7],[-0.1 -0.1 0.1 0.1],'red')
    plot(7:7:140,ccd(1,2:end),'LineWidth',2.0)
    ylim([-1 1])
    xticks(0:7:140)
    yline(-0.07,'--r')
    yline(0.07,'--r')
    grid on
    title('PRCC of D vs ' + labelx(indx) + ' over time')
    xlabel('Day')
    ylabel('PRCC')
end

Tpval = array2table([matWPval matDPval],'RowNames',rn,'VariableNames',cn);
Tprcc = array2table([matWPRCC matDPRCC],'RowNames',rn,'VariableNames',cn);
writetable(Tpval,'usanalysistime-pval.csv','WriteRowNames',true, ...
    'WriteVariableNames',true,'Delimiter',',')
writetable(Tprcc,'usanalysistime-prcc.csv','WriteRowNames',true, ...
    'WriteVariableNames',true,'Delimiter',',')