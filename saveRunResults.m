function saveRunResults(R)
% SAVERUNRESULTS  Shared result writer used by every optimiser (pso/gwo/abc/tsa),
% so the output files and layout are identical across algorithms.
%
% R fields:
%   .dir .prefix .titlePrefix .maxFes .nRun
%   .objectiveValues .rmseValues .mseValues .maeValues .corrValues .r2Values
%   .m0Values .externalRmsValues .pointCounts .cpuTimes   (column vectors)
%   .bestSelections (cell) .traces (struct array, fields .best .rmse .num)

dir=R.dir; pre=R.prefix; maxFes=R.maxFes; nRun=R.nRun;
if ~exist(dir,'dir'), mkdir(dir); end

%% --- Mean convergence traces (FEs axis; comparable across algorithms) ---
objMat  = zeros(nRun,maxFes);
ptMat   = zeros(nRun,maxFes);
rmseMat = zeros(nRun,maxFes);
for i = 1:nRun
    for j = 1:maxFes
        objMat(i,j)  = R.traces(i,j).best;
        ptMat(i,j)   = R.traces(i,j).num;
        rmseMat(i,j) = R.traces(i,j).rmse;
    end
end
objConvergence   = mean(objMat,1);
pointConvergence = mean(ptMat,1);
rmseConvergence  = mean(rmseMat,1);

objectiveValues=R.objectiveValues(:); rmseValues=R.rmseValues(:); mseValues=R.mseValues(:);
maeValues=R.maeValues(:); corrValues=R.corrValues(:); r2Values=R.r2Values(:);
m0Values=R.m0Values(:); externalRmsValues=R.externalRmsValues(:); pointCounts=R.pointCounts(:);
cpuTimes=R.cpuTimes(:); bestSelections=R.bestSelections; traces=R.traces;

p = @(name) fullfile(dir,[pre '_' name]);   % filename helper

%% --- Raw arrays ---
save(p('convergenceTraces.mat'), 'traces');
save(p('objConvergence.mat'),    'objConvergence');
save(p('pointConvergence.mat'),  'pointConvergence');
save(p('rmseConvergence.mat'),   'rmseConvergence');
save(p('objectiveValues.mat'),   'objectiveValues');
save(p('rmseValues.mat'),        'rmseValues');
save(p('mseValues.mat'),         'mseValues');
save(p('maeValues.mat'),         'maeValues');
save(p('corrValues.mat'),        'corrValues');
save(p('r2Values.mat'),          'r2Values');
save(p('m0Values.mat'),          'm0Values');
save(p('externalRmsValues.mat'), 'externalRmsValues');
save(p('pointCounts.mat'),       'pointCounts');
save(p('bestSelections.mat'),    'bestSelections');

%% --- Objective statistics [mean std median best(max) worst(min) meanCPU] ---
objectiveSummary = [mean(objectiveValues), std(objectiveValues), median(objectiveValues), ...
                    max(objectiveValues),  min(objectiveValues), mean(cpuTimes)];
save(p('objectiveSummary.mat'), 'objectiveSummary');

%% --- Point-count statistics [best(min) worst(max) median mean] ---
pointStats = [min(pointCounts), max(pointCounts), median(pointCounts), mean(pointCounts)];
save(p('pointStats.mat'), 'pointStats');

%% --- Metric statistics (unit: metres) ---
rmseStats = [mean(rmseValues), std(rmseValues), median(rmseValues), min(rmseValues), max(rmseValues)];
mseStats  = [mean(mseValues),  std(mseValues),  median(mseValues),  min(mseValues),  max(mseValues)];
maeStats  = [mean(maeValues),  std(maeValues),  median(maeValues),  min(maeValues),  max(maeValues)];
corrStats = [mean(corrValues), std(corrValues), median(corrValues), max(corrValues), min(corrValues)];
r2Stats   = [mean(r2Values),   std(r2Values),   median(r2Values),   max(r2Values),   min(r2Values)];
m0Stats   = [mean(m0Values),   std(m0Values),   median(m0Values),   min(m0Values),   max(m0Values)];
extRmsStats = [mean(externalRmsValues), std(externalRmsValues), median(externalRmsValues), ...
               min(externalRmsValues),  max(externalRmsValues)];
save(p('metricStats.mat'), 'rmseStats','mseStats','maeStats','corrStats','r2Stats','m0Stats','extRmsStats');

%% --- Best run ---
[~, biIdx] = max(objectiveValues);
bestSelection   = bestSelections{biIdx};
bestPointCount  = pointCounts(biIdx);
bestObjective   = objectiveValues(biIdx);
bestRMSE        = rmseValues(biIdx);
bestMAE         = maeValues(biIdx);
bestCorr        = corrValues(biIdx);
bestR2          = r2Values(biIdx);
bestM0          = m0Values(biIdx);
bestExternalRMS = externalRmsValues(biIdx);
save(p('bestSolution.mat'), ...
     'bestSelection','bestPointCount','bestObjective', ...
     'bestRMSE','bestMAE','bestCorr','bestR2','bestM0','bestExternalRMS','biIdx');

%% --- Summary CSV (RMSE/MAE/external-RMS/m0 in cm) ---
fid = fopen(fullfile(dir,[pre '_summary.csv']),'w');
fprintf(fid,'Metric,Value\n');
fprintf(fid,'Objective_Mean,%.6f\n',   mean(objectiveValues));
fprintf(fid,'Objective_Std,%.6f\n',    std(objectiveValues));
fprintf(fid,'Objective_Best(max),%.6f\n', max(objectiveValues));
fprintf(fid,'Objective_Worst(min),%.6f\n',min(objectiveValues));
fprintf(fid,'CPU_s_mean,%.3f\n',  mean(cpuTimes));
fprintf(fid,'Points_Best(min),%d\n',  round(min(pointCounts)));
fprintf(fid,'Points_Worst(max),%d\n', round(max(pointCounts)));
fprintf(fid,'Points_Median,%.1f\n',   median(pointCounts));
fprintf(fid,'Points_Mean,%.2f\n',     mean(pointCounts));
fprintf(fid,'RMSE_Mean_cm,%.4f\n', mean(rmseValues)*100);
fprintf(fid,'RMSE_Best_cm,%.4f\n', min(rmseValues)*100);
fprintf(fid,'ExtRMS_Mean_cm,%.4f\n', mean(externalRmsValues)*100);
fprintf(fid,'ExtRMS_Best_cm,%.4f\n', min(externalRmsValues)*100);
fprintf(fid,'MAE_Mean_cm,%.4f\n',  mean(maeValues)*100);
fprintf(fid,'R_Mean,%.6f\n',   mean(corrValues));
fprintf(fid,'R2_Mean,%.6f\n',  mean(r2Values));
fprintf(fid,'m0_Mean_cm,%.4f\n', mean(m0Values)*100);
fprintf(fid,'BestRun_RMSE_cm,%.4f\n', bestRMSE*100);
fprintf(fid,'BestRun_ExtRMS_cm,%.4f\n', bestExternalRMS*100);
fprintf(fid,'BestRun_Points,%d\n', round(bestPointCount));
fclose(fid);

%% --- Per-algorithm convergence figures (English) ---
try
    figure('Visible','off');
    plot(1:maxFes, objConvergence,'LineWidth',1.5);
    xlabel('Number of Function Evaluations (FEs)'); ylabel('Objective Function Value');
    title([R.titlePrefix ' - Objective Convergence']); grid on;
    saveas(gcf, p('obj_convergence.png'));

    figure('Visible','off');
    plot(1:maxFes, rmseConvergence*100,'LineWidth',1.5);
    xlabel('Number of Function Evaluations (FEs)'); ylabel('RMSE (cm)');
    title([R.titlePrefix ' - RMSE Convergence']); grid on;
    saveas(gcf, p('rmse_convergence.png'));

    figure('Visible','off');
    plot(1:maxFes, pointConvergence,'LineWidth',1.5);
    xlabel('Number of Function Evaluations (FEs)'); ylabel('Number of Reference Points');
    title([R.titlePrefix ' - Reference Point Number Convergence']); grid on;
    saveas(gcf, p('point_convergence.png'));
    close all;
catch
end

fprintf('\n[%s] best run #%d | Obj=%.4f | RMSE=%.4f cm | MAE=%.4f cm | R2=%.4f | Points=%d\n', ...
         pre, biIdx, bestObjective, bestRMSE*100, bestMAE*100, bestR2, bestPointCount);
end
