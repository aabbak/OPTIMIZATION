function R = packResults(dir, prefix, titlePrefix, maxFes, nRun, ...
    objectiveValues, rmseValues, mseValues, maeValues, corrValues, r2Values, ...
    m0Values, externalRmsValues, pointCounts, cpuTimes, bestSelections, traces)
% PACKRESULTS  Collect per-run arrays into the struct expected by saveRunResults.
R.dir=dir; R.prefix=prefix; R.titlePrefix=titlePrefix; R.maxFes=maxFes; R.nRun=nRun;
R.objectiveValues=objectiveValues; R.rmseValues=rmseValues; R.mseValues=mseValues;
R.maeValues=maeValues; R.corrValues=corrValues; R.r2Values=r2Values;
R.m0Values=m0Values; R.externalRmsValues=externalRmsValues; R.pointCounts=pointCounts;
R.cpuTimes=cpuTimes; R.bestSelections=bestSelections; R.traces=traces;
end
