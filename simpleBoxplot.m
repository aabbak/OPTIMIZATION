function simpleBoxplot(M, labels, ylab)
% SIMPLEBOXPLOT  Minimal Tukey boxplot without any toolbox.
%   M : N x k matrix (columns = groups)
% Draws median, Q1-Q3 box, 1.5*IQR whiskers and outliers for each column.
[~,k]=size(M);
hold on;  w=0.3;
for g=1:k
    x=M(:,g); x=x(~isnan(x));
    q1=pct(x,25); q2=median(x); q3=pct(x,75);
    iqr=q3-q1; lo=q1-1.5*iqr; hi=q3+1.5*iqr;
    inw=x(x>=lo & x<=hi); wlo=min(inw); whi=max(inw);
    fill([g-w g+w g+w g-w],[q1 q1 q3 q3],[0.85 0.9 0.98],'EdgeColor',[0.1 0.1 0.4]);
    plot([g-w g+w],[q2 q2],'r-','LineWidth',1.6);
    plot([g g],[q3 whi],'k-'); plot([g g],[wlo q1],'k-');
    plot([g-w/2 g+w/2],[whi whi],'k-'); plot([g-w/2 g+w/2],[wlo wlo],'k-');
    out=x(x<lo | x>hi);
    if ~isempty(out), plot(g*ones(size(out)),out,'k+','MarkerSize',4); end
end
hold off;
set(gca,'XTick',1:k,'XTickLabel',labels);
xlim([0.5 k+0.5]); grid on;
if nargin>=3 && ~isempty(ylab), ylabel(ylab); end
end

function q=pct(x,p)
x=sort(x(:)); n=numel(x);
if n==1, q=x; return; end
pos=(p/100)*(n-1)+1; lo=floor(pos); hi=ceil(pos); fr=pos-lo;
q=x(lo)+(x(hi)-x(lo))*fr;
end
