function [p, z, W, med1, med2] = wilcoxonRankSumTest(x, y)
% WILCOXONRANKSUMTEST  Wilcoxon rank-sum (Mann-Whitney U) test for two
% independent samples. Tie-corrected normal approximation with continuity
% correction; two-sided p-value. No toolbox required.
%   x,y : the 30-run metric values of two algorithms
%   p   : two-sided p-value
%   z   : standardised statistic
%   W   : rank sum of x
% Note: rank-sum is appropriate for independent runs (different seeds); for a
% paired design the signed-rank test would be used instead.
x=x(:); y=y(:); n1=numel(x); n2=numel(y); N=n1+n2;
med1=median(x); med2=median(y);
pooled=[x; y];
[s,ord]=sort(pooled); ranks=zeros(N,1); i=1;
while i<=N
    j=i; while j<N && s(j+1)==s(i), j=j+1; end
    ranks(i:j)=mean(i:j); i=j+1;
end
rnk=zeros(N,1); rnk(ord)=ranks;
W = sum(rnk(1:n1));
mu = n1*(N+1)/2;
[~,~,ic]=unique(pooled); tie=accumarray(ic,1);
sigma2 = (n1*n2/(12*N*(N-1))) * (N^3 - N - sum(tie.^3 - tie));
if sigma2<=0, p=1; z=0; return; end
sigma=sqrt(sigma2);
zc = W - mu;
z  = (zc - sign(zc)*0.5)/sigma;            % continuity correction
p  = 2*(1 - 0.5*erfc(-abs(z)/sqrt(2)));    % two-sided (normcdf = 0.5*erfc(-z/sqrt2))
p  = min(max(p,0),1);
end
