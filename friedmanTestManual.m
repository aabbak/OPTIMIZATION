function [p, chi2, meanRanks, Fid, pF] = friedmanTestManual(data, lowerBetter)
% FRIEDMANTESTMANUAL  Friedman test (non-parametric, k algorithms x N blocks).
% Tie-corrected mean ranks; chi-square and Iman-Davenport (F) p-values.
% No toolbox required (gammainc/betainc are base functions).
%   data        : N x k matrix (rows = blocks/runs, columns = algorithms)
%   lowerBetter : true if smaller is better (RMSE/ext-RMS), false if larger (objective)
%   meanRanks   : mean rank per algorithm (smaller = better)
if nargin<2, lowerBetter=true; end
[N,k]=size(data);
R=zeros(N,k);
for i=1:N
    row=data(i,:);
    if ~lowerBetter, row=-row; end     % larger-is-better -> flip so rank 1 = best
    R(i,:)=tiedRank(row);
end
meanRanks=mean(R,1);
Rj=sum(R,1);
chi2=(12/(N*k*(k+1)))*sum(Rj.^2) - 3*N*(k+1);
df=k-1;
p = 1 - gammainc(chi2/2, df/2);            % chi2cdf = gammainc(x/2, df/2)
den=(N*(k-1)-chi2);
if den<=0
    Fid=Inf; pF=0;
else
    Fid=((N-1)*chi2)/den;
    df1=k-1; df2=(k-1)*(N-1);
    x=df2/(df2+df1*Fid);
    pF = betainc(x, df2/2, df1/2);
end
p=min(max(p,0),1);
end

function r=tiedRank(v)
v=v(:)'; n=numel(v); r=zeros(1,n);
[s,ord]=sort(v); rr=zeros(1,n); i=1;
while i<=n
    j=i; while j<n && s(j+1)==s(i), j=j+1; end
    rr(i:j)=mean(i:j); i=j+1;
end
r(ord)=rr;
end
