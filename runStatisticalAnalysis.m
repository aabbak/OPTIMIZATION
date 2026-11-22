function runStatisticalAnalysis()
% RUNSTATISTICALANALYSIS  Post-processing: build the statistical tables and the
% publication figures from the per-run results of the four algorithms.
%
% USAGE: place this file (and its helpers) in the folder that contains the four
% result folders and run it:
%     results_pso  results_gwo  results_abc  results_tsa
%
% OUTPUT
%   results_statistics/  -> statistical_results.xlsx (multi-sheet Excel)
%                           Sheets: Summary | Descriptive | Wilcoxon_RMS |
%                           Wilcoxon_Objective | Friedman.
%                           (If Excel writing is unavailable, e.g. Octave, the
%                            same tables are written as .csv instead.)
%   results_graphs/      -> figures (.png + .pdf/.eps + .fig)
%
% Helper functions required alongside this file:
%   wilcoxonRankSumTest.m  friedmanTestManual.m  simpleBoxplot.m
% No toolbox required.

algos = {'pso','gwo','abc','tsa'};
PRE   = {'PSO','GWO','ABC','TSA'};
nalg  = numel(algos);
alpha = 0.05;

% --- Locate and validate the result folders ---
fols = cell(1,nalg); missing={};
for a=1:nalg
    fols{a} = ['results_' algos{a}];
    if ~exist(fols{a},'dir'), missing{end+1}=fols{a}; end %#ok<AGROW>
end
if ~isempty(missing)
    error(['The following result folders were not found in this directory:\n  %s\n' ...
           'Place this script next to the four result folders and run it.'], ...
           strjoin(missing, '\n  '));
end

statdir  = 'results_statistics';
graphdir = 'results_graphs';
if ~exist(statdir,'dir'),  mkdir(statdir);  end
if ~exist(graphdir,'dir'), mkdir(graphdir); end

% --- Load the per-run data (30 runs x 4 algorithms) ---
EXT=[]; RMSE=[]; MAE=[]; OBJ=[]; NPT=[]; M0=[]; CPU=zeros(1,nalg);
cOBJ=[]; cRMSE=[]; cNPT=[];
for a = 1:nalg
    fol=fols{a}; p=PRE{a};
    EXT (:,a) = getv(fol,p,'externalRmsValues')*100;  % cm
    RMSE(:,a) = getv(fol,p,'rmseValues')*100;         % cm
    MAE (:,a) = getv(fol,p,'maeValues')*100;          % cm
    OBJ (:,a) = getv(fol,p,'objectiveValues');
    NPT (:,a) = getv(fol,p,'pointCounts');
    M0  (:,a) = getv(fol,p,'m0Values')*100;           % cm
    summ      = getv(fol,p,'objectiveSummary');   CPU(a)=summ(6);
    cOBJ (a,:) = getrow(fol,p,'objConvergence');
    cRMSE(a,:) = getrow(fol,p,'rmseConvergence')*100; % cm
    cNPT (a,:) = getrow(fol,p,'pointConvergence');
end
maxfes = size(cOBJ,2);
fprintf('Loaded 4 algorithms (%d runs). FEs=%d\n', size(EXT,1), maxfes);

%% ================= EXCEL (multi-sheet) =================
xlsx = fullfile(statdir,'statistical_results.xlsx');
if exist(xlsx,'file'), delete(xlsx); end

Summary = summaryCell(EXT, NPT, OBJ, PRE);
Descr   = descriptiveCell(EXT, RMSE, MAE, NPT, OBJ, M0, CPU, PRE);
WilRMS  = wilcoxonCell(EXT, PRE, 1, alpha);
WilObj  = wilcoxonCell(OBJ, PRE, 0, alpha);
Fried   = friedmanCell({EXT,RMSE,NPT,OBJ}, {'External RMS','RMSE','Points','Objective'}, [1 1 1 0], PRE);

sheets = {'Summary',Summary; 'Descriptive',Descr; 'Wilcoxon_RMS',WilRMS; ...
          'Wilcoxon_Objective',WilObj; 'Friedman',Fried};

okx = true;
for s=1:size(sheets,1)
    if ~writeCellSafe(xlsx, sheets{s,1}, sheets{s,2}), okx=false; break; end
end
if okx
    fprintf('Tables -> %s\n', xlsx);
else
    if exist(xlsx,'file'), delete(xlsx); end
    for s=1:size(sheets,1)
        cellToCSV(sheets{s,2}, fullfile(statdir, [sheets{s,1} '.csv']));
    end
    fprintf('Excel writing unavailable; tables written as CSV into %s/\n', statdir);
end

%% ================= FIGURES =================
sty={'-','--','-.',':'};
col=[0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];
x=1:maxfes;

f=figure('Visible','off','Position',[100 100 720 480]); hold on;
for a=1:nalg, plot(x,cOBJ(a,:),sty{a},'Color',col(a,:),'LineWidth',1.8); end
xlabel('Number of Function Evaluations (FEs)'); ylabel('Objective Function Value');
legend(PRE,'Location','southeast'); grid on; box on;
title('Objective Function Convergence');
savefigs(f, fullfile(graphdir,'Fig1_objective_convergence'));

f=figure('Visible','off','Position',[100 100 720 480]); hold on;
for a=1:nalg, plot(x,cRMSE(a,:),sty{a},'Color',col(a,:),'LineWidth',1.8); end
xlabel('Number of Function Evaluations (FEs)'); ylabel('RMSE (cm)');
legend(PRE,'Location','northeast'); grid on; box on;
title('RMSE Convergence');
savefigs(f, fullfile(graphdir,'Fig2_rmse_convergence'));

f=figure('Visible','off','Position',[100 100 720 480]); hold on;
for a=1:nalg, plot(x,cNPT(a,:),sty{a},'Color',col(a,:),'LineWidth',1.8); end
xlabel('Number of Function Evaluations (FEs)'); ylabel('Number of Reference Points');
legend(PRE,'Location','northeast'); grid on; box on;
title('Convergence of Selected Reference Point Number');
savefigs(f, fullfile(graphdir,'Fig3_point_convergence'));

f=figure('Visible','off','Position',[100 100 640 460]);
simpleBoxplot(EXT, PRE, 'External RMS (cm)');
title('Distribution of External Accuracy (RMS) over 30 Runs'); box on;
savefigs(f, fullfile(graphdir,'Fig4_boxplot_external_rms'));

f=figure('Visible','off','Position',[100 100 640 460]);
simpleBoxplot(OBJ, PRE, 'Objective Function Value');
title('Distribution of Objective Function over 30 Runs'); box on;
savefigs(f, fullfile(graphdir,'Fig5_boxplot_objective'));

close all;
fprintf('Figures -> %s/\n', graphdir);
fprintf('=== Analysis finished. ===\n');
end

% ======================= DATA / CELL BUILDERS =======================
function v=getv(fol,pre,name)
    s=load(fullfile(fol,[pre '_' name '.mat'])); v=s.(name); v=v(:);
end
function r=getrow(fol,pre,name)
    s=load(fullfile(fol,[pre '_' name '.mat'])); r=s.(name); r=r(:)';
end

function C = descriptiveCell(EXT,RMSE,MAE,NPT,OBJ,M0,CPU,PRE)
    names={'External RMS (cm)','RMSE (cm)','MAE (cm)','Points','Objective','m0 (cm)'};
    dat ={EXT,RMSE,MAE,NPT,OBJ,M0};  low=[1 1 1 1 0 1];
    stats={'Mean','Std','Median','Best'};
    C=[{'Metric','Statistic'} PRE];
    for mI=1:numel(names)
        X=dat{mI};
        for sI=1:numel(stats)
            row=cell(1,2+numel(PRE));
            if sI==1, row{1}=names{mI}; else, row{1}=''; end
            row{2}=stats{sI};
            for a=1:numel(PRE)
                col=X(:,a);
                switch stats{sI}
                    case 'Mean',   val=mean(col);
                    case 'Std',    val=std(col);
                    case 'Median', val=median(col);
                    case 'Best',   if low(mI), val=min(col); else, val=max(col); end
                end
                row{2+a}=rnd(val,4);
            end
            C=[C; row]; %#ok<AGROW>
        end
    end
    crow=[{'CPU (s)','Mean'} num2cell(rnd(CPU,3))];
    C=[C; crow];
end

function C = wilcoxonCell(X, PRE, lowerBetter, alpha)
    nalg=numel(PRE); mn=mean(X,1);
    if lowerBetter, [~,ctrl]=min(mn); dir='lower is better'; else, [~,ctrl]=max(mn); dir='higher is better'; end
    C={ sprintf('Wilcoxon rank-sum p-value matrix (%s)',dir), '', '', '', '' };
    C=[C; [{''} PRE]];
    for i=1:nalg
        row=cell(1,1+nalg); row{1}=PRE{i};
        for j=1:nalg
            if i==j, row{1+j}='-';
            else, row{1+j}=pval(wilcoxonRankSumTest(X(:,i),X(:,j))); end
        end
        C=[C; row]; %#ok<AGROW>
    end
    C=[C; repmat({''},1,1+nalg)];
    C=[C; [{sprintf('Control algorithm (best mean): %s',PRE{ctrl})} repmat({''},1,nalg)]];
    C=[C; [{'Comparison','p-value','Significant (0.05)','Direction'} repmat({''},1,nalg-3)]];
    for j=1:nalg
        if j==ctrl, continue; end
        [p,~,~,mc,mj]=wilcoxonRankSumTest(X(:,ctrl),X(:,j));
        if lowerBetter, better=mc<mj; else, better=mc>mj; end
        if p<alpha
            if better, d='control better (+)'; else, d='control worse (-)'; end
            sig='YES';
        else
            d='not significant (=)'; sig='NO';
        end
        row=[{sprintf('%s vs %s',PRE{ctrl},PRE{j})}, {pval(p)}, {sig}, {d}, repmat({''},1,nalg-3)];
        C=[C; row]; %#ok<AGROW>
    end
end

function C = friedmanCell(datas, names, lowFlags, PRE)
    C=[{'Friedman test - mean ranks (smaller rank = better)'} repmat({''},1,numel(PRE)+4)];
    C=[C; [{'Metric'} PRE {'chi2','p','ImanDavenportF','pF'}]];
    for mI=1:numel(datas)
        [p,chi2,mr,Fid,pF]=friedmanTestManual(datas{mI}, lowFlags(mI));
        row=[{names{mI}}, num2cell(rnd(mr,3)), {rnd(chi2,4), pval(p), rnd(Fid,4), pval(pF)}];
        C=[C; row]; %#ok<AGROW>
    end
end

function C = summaryCell(EXT, NPT, OBJ, PRE)
    [~,~,mr]=friedmanTestManual(OBJ,0); [~,ord]=sort(mr);    % objective ranking (best->worst)
    [~,bestEXT]=min(mean(EXT,1)); [~,bestNPT]=min(mean(NPT,1));
    ranking=strjoin(PRE(ord),' > ');
    C={
     'SUMMARY AND HOW TO READ THESE RESULTS';
     '';
     'This workbook summarises 30 independent runs of four algorithms (PSO, GWO, ABC, TSA).';
     'Sheets: Descriptive | Wilcoxon_RMS | Wilcoxon_Objective | Friedman.';
     '';
     'HOW TO READ:';
     '- External RMS (cm): SMALLER is better (accuracy at the independent test points).';
     '- Points: SMALLER is better (solution uses fewer reference points).';
     '- Objective: LARGER is better (accuracy and sparsity combined).';
     '- Wilcoxon p<0.05 means the two algorithms differ significantly; otherwise not.';
     '- In Friedman, a smaller mean rank means a generally better algorithm.';
     '';
     'RESULT FOR THIS SETTING:';
     sprintf('- Objective-based Friedman ranking (best to worst): %s.', ranking);
     sprintf('- Best mean external accuracy: %s (%.3f cm).', PRE{bestEXT}, min(mean(EXT,1)));
     sprintf('- Fewest points on average: %s (%.2f points).', PRE{bestNPT}, min(mean(NPT,1)));
     '- GWO and ABC may be statistically tied on accuracy; see the Wilcoxon_RMS sheet.';
    };
end

% ======================= WRITE HELPERS =======================
function s=pval(p)
    if p>=1e-4, s=rnd(p,5); else, s=sprintf('%.2e',p); end
end
function y=rnd(x,n)
    f=10^n; y=round(x*f)/f;
end

function ok=writeCellSafe(file, sheet, C)
    ok=true;
    try
        writecell(C, file, 'Sheet', sheet);
    catch
        try, xlswrite(file, C, sheet); catch, ok=false; end
    end
end

function cellToCSV(C, path)
    fid=fopen(path,'w');
    for i=1:size(C,1)
        for j=1:size(C,2)
            v=C{i,j};
            if isempty(v), str='';
            elseif ischar(v), str=v;
            else, str=num2str(v);
            end
            if any(str==','), str=['"' str '"']; end
            if j>1, fprintf(fid,','); end
            fprintf(fid,'%s',str);
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
end

function savefigs(f, pathNoExt)
    try, print(f,[pathNoExt '.png'],'-dpng','-r300'); catch, saveas(f,[pathNoExt '.png']); end
    try
        if exist('exportgraphics','file')
            exportgraphics(f,[pathNoExt '.pdf'],'ContentType','vector');
        else
            print(f,[pathNoExt '.eps'],'-depsc');
        end
    catch
    end
    try, savefig(f,[pathNoExt '.fig']); catch, try, hgsave(f,[pathNoExt '.fig']); catch, end; end
end
