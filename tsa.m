function [gbest, bestParams] = tsa(nRun, n, D, maxFes, GEO)
% TSA  Tree-Seed Algorithm (Kiran, 2015) for hybrid-geoid reference-point
% selection. Same setup as PSO: continuous trees in [-5,5]^D,
% binaryTransfer -> objectiveFunction (maximised, B = best tree).

fun  = @(so) objectiveFunction(so, GEO);
dmin = -5;  dmax = 5;
n = max(n,3);
ST   = 0.1;                        % search tendency
lowS = max(1, ceil(0.10*n));       % min seeds per tree
highS= max(lowS, ceil(0.25*n));    % max seeds per tree
emptyMetrics = struct('RMSE',Inf,'MSE',Inf,'MAE',Inf,'R',0,'R2',0,'m0',Inf,'extRMS',Inf);

traces = struct([]);
objectiveValues=zeros(1,nRun); pointCounts=zeros(1,nRun);
rmseValues=zeros(1,nRun); mseValues=zeros(1,nRun); maeValues=zeros(1,nRun);
corrValues=zeros(1,nRun); r2Values=zeros(1,nRun); m0Values=zeros(1,nRun);
externalRmsValues=zeros(1,nRun); bestSelections=cell(1,nRun); cpuTimes=zeros(1,nRun);

for rx = 1:nRun
    rng('shuffle');  tic;

    Trees = dmin + (dmax-dmin)*rand(n,D);
    objTree = -inf(n,1);
    gbest=-inf; gMet=emptyMetrics; gNum=D; gSel=[];
    fes = 0;

    for i = 1:n
        [met,obj,sel]=fun(binaryTransfer(Trees(i,:))); fes=fes+1;
        objTree(i)=obj;
        if obj>gbest, gbest=obj; gMet=met; gNum=numel(sel); gSel=sel; end
        traces(rx,fes).best=gbest; traces(rx,fes).rmse=gMet.RMSE; traces(rx,fes).num=gNum;
    end
    [~,bi]=max(objTree);  B=Trees(bi,:);

    iter = 1;
    while fes < maxFes
        for i = 1:n
            if fes>=maxFes, break; end
            nSeed = round(lowS + rand*(highS-lowS));  nSeed=max(nSeed,1);
            bestSeedObj=-inf; bestSeed=Trees(i,:);
            for s = 1:nSeed
                if fes>=maxFes, break; end
                r=randi(n); while r==i, r=randi(n); end
                seed=Trees(i,:);
                for j=1:D
                    if rand < ST
                        seed(j)=Trees(i,j)+(2*rand-1)*(B(j)       -Trees(r,j));
                    else
                        seed(j)=Trees(i,j)+(2*rand-1)*(Trees(i,j) -Trees(r,j));
                    end
                end
                seed=max(min(seed,dmax),dmin);
                [met,obj,sel]=fun(binaryTransfer(seed)); fes=fes+1;
                if obj>bestSeedObj, bestSeedObj=obj; bestSeed=seed; end
                if obj>gbest, gbest=obj; gMet=met; gNum=numel(sel); gSel=sel; end
                traces(rx,fes).best=gbest; traces(rx,fes).rmse=gMet.RMSE; traces(rx,fes).num=gNum;
            end
            % Greedy selection: replace the tree if the best seed is better
            if bestSeedObj > objTree(i)
                Trees(i,:)=bestSeed; objTree(i)=bestSeedObj;
                if objTree(i) >= max(objTree), B=Trees(i,:); end
            end
        end
        fprintf('TSA run=%d  FEs=%d  iter=%d  f=%g  RMSE=%.4fcm  points=%d\n', ...
                 rx,fes,iter,gbest,gMet.RMSE*100,gNum);
        iter = iter + 1;
    end

    objectiveValues(rx)=gbest; rmseValues(rx)=gMet.RMSE; mseValues(rx)=gMet.MSE;
    maeValues(rx)=gMet.MAE; corrValues(rx)=gMet.R; r2Values(rx)=gMet.R2;
    m0Values(rx)=gMet.m0; externalRmsValues(rx)=gMet.extRMS; pointCounts(rx)=gNum;
    bestSelections{rx}=gSel; cpuTimes(rx)=toc;
end

gbest = max(objectiveValues);  bestParams = B;

R = packResults('results_tsa','TSA','TSA (hybrid geoid)', maxFes, nRun, ...
    objectiveValues, rmseValues, mseValues, maeValues, corrValues, r2Values, ...
    m0Values, externalRmsValues, pointCounts, cpuTimes, bestSelections, traces);
saveRunResults(R);
end
