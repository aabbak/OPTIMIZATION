function [gbest, bestParams] = gwo(nRun, n, D, maxFes, GEO)
% GWO  Grey Wolf Optimiser (Mirjalili et al., 2014) for hybrid-geoid
% reference-point selection. Same setup as PSO: continuous X in [-5,5]^D,
% binaryTransfer -> objectiveFunction (maximised, alpha = best wolf).

fun  = @(so) objectiveFunction(so, GEO);
dmin = -5;  dmax = 5;
emptyMetrics = struct('RMSE',Inf,'MSE',Inf,'MAE',Inf,'R',0,'R2',0,'m0',Inf,'extRMS',Inf);

traces = struct([]);
objectiveValues=zeros(1,nRun); pointCounts=zeros(1,nRun);
rmseValues=zeros(1,nRun); mseValues=zeros(1,nRun); maeValues=zeros(1,nRun);
corrValues=zeros(1,nRun); r2Values=zeros(1,nRun); m0Values=zeros(1,nRun);
externalRmsValues=zeros(1,nRun); bestSelections=cell(1,nRun); cpuTimes=zeros(1,nRun);

for rx = 1:nRun
    rng('shuffle');  tic;

    X = dmin + (dmax-dmin)*rand(n,D);
    alphaPos=zeros(1,D); alphaScore=-inf; alphaMet=emptyMetrics; alphaNum=D; alphaSel=[];
    betaPos =zeros(1,D); betaScore =-inf;
    deltaPos=zeros(1,D); deltaScore=-inf;
    fes = 0;

    for i = 1:n
        [met,obj,sel] = fun(binaryTransfer(X(i,:)));  fes = fes + 1;
        [alphaPos,alphaScore,alphaMet,alphaNum,alphaSel, betaPos,betaScore, deltaPos,deltaScore] = ...
            updateLeaders(obj,met,sel,X(i,:), alphaPos,alphaScore,alphaMet,alphaNum,alphaSel, ...
                          betaPos,betaScore, deltaPos,deltaScore);
        traces(rx,fes).best=alphaScore; traces(rx,fes).rmse=alphaMet.RMSE; traces(rx,fes).num=alphaNum;
    end

    iter = 1;
    while fes < maxFes
        a = 2 - 2*(fes/maxFes);          % 2 -> 0 over the FE budget
        for i = 1:n
            for j = 1:D
                A1=2*a*rand-a; C1=2*rand; d1=abs(C1*alphaPos(j)-X(i,j)); x1=alphaPos(j)-A1*d1;
                A2=2*a*rand-a; C2=2*rand; d2=abs(C2*betaPos(j) -X(i,j)); x2=betaPos(j) -A2*d2;
                A3=2*a*rand-a; C3=2*rand; d3=abs(C3*deltaPos(j)-X(i,j)); x3=deltaPos(j)-A3*d3;
                X(i,j) = (x1+x2+x3)/3;
            end
        end
        X = max(min(X,dmax),dmin);

        for i = 1:n
            if fes >= maxFes, break; end
            [met,obj,sel] = fun(binaryTransfer(X(i,:)));  fes = fes + 1;
            [alphaPos,alphaScore,alphaMet,alphaNum,alphaSel, betaPos,betaScore, deltaPos,deltaScore] = ...
                updateLeaders(obj,met,sel,X(i,:), alphaPos,alphaScore,alphaMet,alphaNum,alphaSel, ...
                              betaPos,betaScore, deltaPos,deltaScore);
            traces(rx,fes).best=alphaScore; traces(rx,fes).rmse=alphaMet.RMSE; traces(rx,fes).num=alphaNum;
        end
        fprintf('GWO run=%d  FEs=%d  iter=%d  f=%g  RMSE=%.4fcm  points=%d\n', ...
                 rx,fes,iter,alphaScore,alphaMet.RMSE*100,alphaNum);
        iter = iter + 1;
    end

    objectiveValues(rx)=alphaScore; rmseValues(rx)=alphaMet.RMSE; mseValues(rx)=alphaMet.MSE;
    maeValues(rx)=alphaMet.MAE; corrValues(rx)=alphaMet.R; r2Values(rx)=alphaMet.R2;
    m0Values(rx)=alphaMet.m0; externalRmsValues(rx)=alphaMet.extRMS; pointCounts(rx)=alphaNum;
    bestSelections{rx}=alphaSel; cpuTimes(rx)=toc;
end

gbest = max(objectiveValues);  bestParams = alphaPos;

R = packResults('results_gwo','GWO','GWO (hybrid geoid)', maxFes, nRun, ...
    objectiveValues, rmseValues, mseValues, maeValues, corrValues, r2Values, ...
    m0Values, externalRmsValues, pointCounts, cpuTimes, bestSelections, traces);
saveRunResults(R);
end

function [aP,aS,aM,aN,aSel, bP,bS, dP,dS] = updateLeaders(obj,met,sel,pos, ...
         aP,aS,aM,aN,aSel, bP,bS, dP,dS)
% Maintain alpha/beta/delta leaders for a maximised objective.
if obj > aS
    dS=bS; dP=bP;  bS=aS; bP=aP;
    aS=obj; aP=pos; aM=met; aN=numel(sel); aSel=sel;
elseif obj > bS
    dS=bS; dP=bP;  bS=obj; bP=pos;
elseif obj > dS
    dS=obj; dP=pos;
end
end
