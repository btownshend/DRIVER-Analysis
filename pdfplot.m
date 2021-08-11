% pdfplot - plot PDF of data
function [h,py]=pdfplot(x,popts,dolog,nbins)
if nargin<2 || length(popts)==0
  popts='-';
end
if nargin<3
  dolog=0;
end
if nargin<4
  nbins=max(2,min([ceil(length(x)/10), 200]));
end
if (mean(x<=0) > 0.01 & dolog)
  fprintf('pdfplot Warning: %.1f%% (%d points) of data is non-positive and log plotting was requested.\n', 100*mean(x<=0),sum(x<=0));
end
if dolog
  x=x(x>0);
end
xs=sort(x);
first=1;
px=[]; py=[];
for i=1:nbins
  last=round(i/nbins*length(xs));
  if xs(last)-xs(first) < eps
    continue;
  end
  pval=(last-first+1)/length(xs)/(xs(last)-xs(first));
  px=[px,xs(first),xs(last)];
  if dolog
    py=[py,pval*mean(xs(first:last)),pval*mean(xs(first:last))];
  else
    py=[py,pval,pval];
  end
  first=last+1;
end
nlarge=sum(~isfinite(py) | py>1000);
if nlarge>0
  fprintf('pdfplot: %d samples are not finite or >1000\n', nlarge);
end
if dolog
  h=semilogx(px,py,popts);
  ylabel('PDF (scaled by X)');
else
  h=plot(px,py,popts);
  ylabel('PDF');
end
