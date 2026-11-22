function [gbest, bestParams] = abc(nRun, n, D, maxFes, GEO)
% ABC  Artificial Bee Colony (Karaboga, 2005) for hybrid-geoid reference-point
% selection. Same setup as PSO: continuous food sources in [-5,5]^D,
% binaryTransfer -> objectiveFunction (maximised). Selection probabilities use
% fitnessTransform(-objective) so a larger objective maps to a larger, positive
% fitness.

fun  = @(so) objectiveFunction(so, GEO);
dmin = -5;  dmax = 5;
n = max(n,3);                      % at least 3 sources for neighbour selection
limit = round(0.6*n*D);           % abandonment threshold (scout)
emptyMetrics = struct('RMSE',Inf,'MSE',Inf,'MAE',Inf,'R',0,'R2',0,'m0',Inf,'extRMS',Inf);

traces = struct([]);
objectiveValues=zeros(1,nRun); pointCounts=zeros(1,nRun);
rmseValues=zeros(1,nRun); mseValues=zeros(1,nRun); maeValues=zeros(1,nRun);
corrValues=zeros(1,nRun); r2Values=zeros(1,nRun); m0Values=zeros(1,nRun);
externalRmsValues=zeros(1,nRun); bestSelections=cell(1,nRun); cpuTimes=zeros(1,nRun);

for rx = 1:nRun
    rng('shuffle');  tic;

    Foods  = dmin + (dmax-dmin)*rand(n,D);
    fitVal = zeros(n,1);  trial = zeros(n,1);
    gbest=-inf; gMet=emptyMetrics; gNum=D; gSel=[];
    fes = 0;

    for i = 1:n
        [met,obj,sel] = fun(binaryTransfer(Foods(i,:)));
        fitVal(i)=fitnessTransform(-obj);  trial(i)=0;  fes=fes+1;
        if obj>gbest, gbest=obj; gMet=met; gNum=numel(sel); gSel=sel; end
        traces(rx,fes).best=gbest; traces(rx,fes).rmse=gMet.RMSE; traces(rx,fes).num=gNum;
    end

    iter = 1;
    while fes < maxFes
        % ----- Employed bees -----
        for i = 1:n
            if fes>=maxFes, break; end
            k=randi(n); while k==i, k=randi(n); end
            j=randi(D);
            v=Foods(i,:); v(j)=Foods(i,j)+(2*rand-1)*(Foods(i,j)-Foods(k,j));
            v=max(min(v,dmax),dmin);
            [met,obj,sel]=fun(binaryTransfer(v)); fes=fes+1;
            fv=fitnessTransform(-obj);
            if fv>fitVal(i), Foods(i,:)=v; fitVal(i)=fv; trial(i)=0; else, trial(i)=trial(i)+1; end
            if obj>gbest, gbest=obj; gMet=met; gNum=numel(sel); gSel=sel; end
            traces(rx,fes).best=gbest; traces(rx,fes).rmse=gMet.RMSE; traces(rx,fes).num=gNum;
        end

        % ----- Onlooker bees -----
        P = fitVal / sum(fitVal);
        i=1; t=0;
        while t < n
            if fes>=maxFes, break; end
            if rand < P(i)
                t=t+1;
                k=randi(n); while k==i, k=randi(n); end
                j=randi(D);
                v=Foods(i,:); v(j)=Foods(i,j)+(2*rand-1)*(Foods(i,j)-Foods(k,j));
                v=max(min(v,dmax),dmin);
                [met,obj,sel]=fun(binaryTransfer(v)); fes=fes+1;
                fv=fitnessTransform(-obj);
                if fv>fitVal(i), Foods(i,:)=v; fitVal(i)=fv; trial(i)=0; else, trial(i)=trial(i)+1; end
                if obj>gbest, gbest=obj; gMet=met; gNum=numel(sel); gSel=sel; end
                traces(rx,fes).best=gbest; traces(rx,fes).rmse=gMet.RMSE; traces(rx,fes).num=gNum;
            end
            i=i+1; if i>n, i=1; end
        end

        % ----- Scout bee -----
        [mx,si]=max(trial);
        if mx > limit && fes < maxFes
            Foods(si,:)=dmin+(dmax-dmin)*rand(1,D);
            [met,obj,sel]=fun(binaryTransfer(Foods(si,:))); fes=fes+1;
            fitVal(si)=fitnessTransform(-obj); trial(si)=0;
            if obj>gbest, gbest=obj; gMet=met; gNum=numel(sel); gSel=sel; end
            traces(rx,fes).best=gbest; traces(rx,fes).rmse=gMet.RMSE; traces(rx,fes).num=gNum;
        end

        fprintf('ABC run=%d  FEs=%d  iter=%d  f=%g  RMSE=%.4fcm  points=%d\n', ...
                 rx,fes,iter,gbest,gMet.RMSE*100,gNum);
        iter = iter + 1;
    end

    objectiveValues(rx)=gbest; rmseValues(rx)=gMet.RMSE; mseValues(rx)=gMet.MSE;
    maeValues(rx)=gMet.MAE; corrValues(rx)=gMet.R; r2Values(rx)=gMet.R2;
    m0Values(rx)=gMet.m0; externalRmsValues(rx)=gMet.extRMS; pointCounts(rx)=gNum;
    bestSelections{rx}=gSel; cpuTimes(rx)=toc;
end

gbest = max(objectiveValues);  bestParams = [];

R = packResults('results_abc','ABC','ABC (hybrid geoid)', maxFes, nRun, ...
    objectiveValues, rmseValues, mseValues, maeValues, corrValues, r2Values, ...
    m0Values, externalRmsValues, pointCounts, cpuTimes, bestSelections, traces);
saveRunResults(R);
end
