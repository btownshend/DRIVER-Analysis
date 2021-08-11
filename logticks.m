% Add logarithmic ticks to an axis
function logticks(xlog,ylog,clog)
if nargin<1
  xlog=true;
end
if nargin<2
  ylog=xlog;
end
if nargin<3
  clog=false;
end
for ax=1:3
  if ax==1 && ~xlog
    continue;
  end
  if ax==2 && ~ylog
    continue;
  end
  if ax==3 && ~clog
    continue;
  end
  ticks=[];
  ticklabels={};
  if ax==1
    c=axis;
    cax=c([1,2]);
  elseif ax==2
    c=axis;
    cax=c([3,4]);
  elseif ax==3
    cax=10.^caxis;
  end
  if cax(1)<=0
    if cax()<=0
      fprintf('Can''t add log ticks since axes bounds are <= 0');
      continue;
    end
    cax(1)=cax(2)/1000;
  end
  ndecades=log10(cax(2))-log10(cax(1));
  for i=floor(log10(cax(1))):ceil(log10(cax(2)))
    if i<0
      fmt=sprintf('%%.%df',-i);
    elseif i>=4
      fmt='%.1g';
    else
      fmt='%.0f';
    end
    for j=1:9
      tval=j*10^i;
      if tval<cax(1) || tval>cax(2)
        continue;
      end
      ticks(end+1)=tval;
      if j==1 || (j==2 && ndecades<5) || (j==5 && ndecades <3) || ndecades<1.2
        if abs(tval)>=1 && abs(tval)<10000
          ticklabels{end+1}=sprintf('%d',tval);
        elseif abs(tval)>=0.1 && abs(tval)<1
          ticklabels{end+1}=sprintf('%.1f',tval);
        elseif abs(tval)>=0.01 && abs(tval)<0.1
          ticklabels{end+1}=sprintf('%.2f',tval);
        elseif tval/10^i == 1
          ticklabels{end+1}=sprintf('10^{%d}',i);
        else
          ticklabels{end+1}=sprintf('%dx10^{%d}',tval/10^i,i);
        end
      else
        ticklabels{end+1}='';
      end
    end
  end
  if ax==1
    set(gca,'XTick',ticks);
    set(gca,'XTickLabel',ticklabels);
  elseif ax==2
    set(gca,'YTick',ticks);
    set(gca,'YTickLabel',ticklabels);
  elseif ax==3
    % Find the colorbar
    ch=get(gcf,'Children');
    h=[];
    for i=1:length(ch)
      if isa(ch(i),'matlab.graphics.illustration.ColorBar')
        h=ch(i);
        break;
      end
    end
    if isempty(h)
      error('logticks: Unable to locate colorbar');
    end
    set(h,'Ticks',log10(ticks));
    set(h,'TickLabels',ticklabels);
  end
end
