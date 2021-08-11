% Save a figure for use in powerpoint or prezi
function mkinsert(fnum,varargin)
  defaults=struct('name',[],'aspect',16/9,'dpi',72,'fontsize',16,'markersize',20,'width',1280);
  args=processargs(defaults,varargin);

  if nargin<1
    fnum=gcf;
  end

  if strcmp(fnum,'all')
    global figlist;
    for i=1:length(figlist.fignum)
      try   % In case figure no longer exists
        fnum=figlist.fignum(i);
        mkinsert(fnum,'aspect',args.aspect,'dpi',args.dpi,'fontsize',args.fontsize,'markersize',args.markersize);
      catch me
      end
    end
    return;
  end

  if isempty(args.name)
    args.name=get(fnum,'Name');
    if isempty(args.name)
      args.name=sprintf('Fig%d',get(fnum,'Number'));
    end
    args.name=[args.name,'.png'];
    args.name=strrep(args.name,'/','_');
  end

  figure(fnum);
  width=args.width;   
  height=width/args.aspect;
  width=width-0.5;% Always end up with 1 extra pixel 
  height=height-0.5;

  ps=get(gcf,'position');
  ps(3)=width;
  ps(4)=height;
  set(gcf,'position',ps);
  set(gcf,'PaperUnits','inches');
  resolution=72;
  set(gcf,'PaperPosition',[0 0 width height]/resolution);
  set(gcf,'PaperSize',[width height]/resolution);
  axes=get(gcf,'children');
  for j=1:length(axes)
    props=get(axes(j));
    if isfield(props,'FontSize')
      set(axes(j),'FontSize',args.fontsize);
    end
    if isfield(props,'XLabel')
      set(get(axes(j),'XLabel'),'FontSize',args.fontsize,'FontWeight','bold');
    end
    if isfield(props,'YLabel')
      set(get(axes(j),'YLabel'),'FontSize',args.fontsize,'FontWeight','bold');
    end
    if isfield(props,'Title')
      set(get(axes(j),'Title'),'FontSize',args.fontsize,'FontWeight','bold');
    end
    c=get(axes(j),'Children');
    for i=1:length(c)
      try
        if isfield(get(c(i)),'MarkerSize')
          set(c(i),'MarkerSize',args.markersize);
        end
        if isfield(get(c(i)),'FontSize')
          set(c(i),'FontSize',args.fontsize);
        end
      catch me
        fprintf('Ignoring exception\n');
      end
    end
  end
  print(gcf,'-dpng',['-r',num2str(args.dpi)],args.name);
  fprintf('Saved figure to %s\n', args.name);
end
