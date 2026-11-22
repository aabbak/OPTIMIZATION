function [gbest, bestParams] = pso(nRun, n, D, maxFes, GEO)
% PSO  Particle Swarm Optimisation for hybrid-geoid reference-point selection.
% The swarm searches a continuous space in [-5,5]^D; binaryTransfer maps each
% particle to a 0/1 selection vector evaluated by objectiveFunction (maximised).
% Results are written in the shared format via saveRunResults.

fun  = @(so) objectiveFunction(so, GEO);
dmin = -5;  dmax = 5;
c1 = 1.5;   c2 = 2.5;
emptyMetrics = struct('RMSE',Inf,'MSE',Inf,'MAE',Inf,'R',0,'R2',0,'m0',Inf,'extRMS',Inf);

traces = struct([]);
objectiveValues=zeros(1,nRun); pointCounts=zeros(1,nRun);
rmseValues=zeros(1,nRun); mseValues=zeros(1,nRun); maeValues=zeros(1,nRun);
corrValues=zeros(1,nRun); r2Values=zeros(1,nRun); m0Values=zeros(1,nRun);
externalRmsValues=zeros(1,nRun); bestSelections=cell(1,nRun); cpuTimes=zeros(1,nRun);

for rx = 1:nRun
    rng('shuffle');  tic;

    particle = dmin + (dmax-dmin)*rand(n,D);
    pbestPos = particle;
    pbest    = -inf(1,n);
    gbest    = -inf;  bestParams = zeros(1,D);
    gMet = emptyMetrics;  gNum = D;  gSel = [];
    vel = 0.1*particle;
    fes = 0;

    % --- Initial swarm ---
    for i = 1:n
        [met,obj,sel] = fun(binaryTransfer(particle(i,:)));
        pbest(i) = obj;  fes = fes + 1;
        if obj > gbest
            gbest=obj; bestParams=particle(i,:); gMet=met; gNum=numel(sel); gSel=sel;
        end
        traces(rx,fes).best=gbest; traces(rx,fes).rmse=gMet.RMSE; traces(rx,fes).num=gNum;
    end

    % --- Main loop ---
    iter = 1;
    while fes < maxFes
        w = 0.1 + rand*0.4;          % inertia weight
        for i = 1:n
            for j = 1:D
                vel(i,j) = w*vel(i,j) + c1*rand*(pbestPos(i,j)-particle(i,j)) ...
                                      + c2*rand*(bestParams(j)-particle(i,j));
                particle(i,j) = particle(i,j) + vel(i,j);
            end
        end
        particle = max(min(particle,dmax),dmin);

        for i = 1:n
            [met,obj,sel] = fun(binaryTransfer(particle(i,:)));
            if obj > pbest(i), pbest(i)=obj; pbestPos(i,:)=particle(i,:); end
            fes = fes + 1;
            if fes <= maxFes
                if obj > gbest
                    gbest=obj; bestParams=particle(i,:); gMet=met; gNum=numel(sel); gSel=sel;
                end
                traces(rx,fes).best=gbest; traces(rx,fes).rmse=gMet.RMSE; traces(rx,fes).num=gNum;
            end
        end
        fprintf('PSO run=%d  FEs=%d  iter=%d  f=%g  RMSE=%.4fcm  points=%d\n', ...
                 rx,fes,iter,gbest,gMet.RMSE*100,gNum);
        iter = iter + 1;
    end

    objectiveValues(rx)=gbest; rmseValues(rx)=gMet.RMSE; mseValues(rx)=gMet.MSE;
    maeValues(rx)=gMet.MAE; corrValues(rx)=gMet.R; r2Values(rx)=gMet.R2;
    m0Values(rx)=gMet.m0; externalRmsValues(rx)=gMet.extRMS; pointCounts(rx)=gNum;
    bestSelections{rx}=gSel; cpuTimes(rx)=toc;
end

gbest = max(objectiveValues);

R = packResults('results_pso','PSO','PSO (hybrid geoid)', maxFes, nRun, ...
    objectiveValues, rmseValues, mseValues, maeValues, corrValues, r2Values, ...
    m0Values, externalRmsValues, pointCounts, cpuTimes, bestSelections, traces);
saveRunResults(R);
end
