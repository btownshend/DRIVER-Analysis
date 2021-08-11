% setfig - Create a new figure or reuse one with given name
function setfig(s)
global figlist setfig_tiling;
if strcmpi(s,'BEGIN')
  setfig_tiling=1;
  tiledlayout('flow');
  return;
elseif strcmpi(s,'END')
  setfig_tiling=0;
  return;
end
if ~exist('setfig_tiling')
  setfig_tiling=0;
end
if setfig_tiling
  fprintf('setfig: tiling on\n');
  nexttile
  return;
end
undef=0;
eval('figlist.name;','undef=1;');
if undef==1
  figlist.name={s};
  figlist.fignum=20;
  figure(figlist.fignum);
  set(figlist.fignum(end),'Name',figlist.name{1});
  return;
end
for i=1:length(figlist.name)
  if strcmp(figlist.name(i),s) == 1
    figure(figlist.fignum(i));
    set(figlist.fignum(i),'Name',figlist.name{i});
    return;
  end
end
figlist.name{end+1}=s;
figure(max(figlist.fignum)+1);
if isnumeric(gcf)
  figlist.fignum=[figlist.fignum,gcf];
else
  figlist.fignum=[figlist.fignum,get(gcf,'Number')];
end  
set(figlist.fignum(end),'Name',figlist.name{end});
return;
