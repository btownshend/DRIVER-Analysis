% Shadow of database compounds.Compounds
classdef Compounds < handle
  properties
    c;    % array of individual compounds each with:
          % compound,name,inchikey,smiles,formula,mass
    domains;  % Names of alias domains
  end

  methods(Access=private)
    function obj=Compounds()
      obj.load();
    end
  end
  
  methods(Static)
    function obj=instance()
      persistent theInstance
      if isempty(theInstance)
        theInstance=Compounds();
      end
      obj=theInstance;
    end
  end
  
  methods
    function load(obj)
    % Load from database
      NGSDatabase.open();
      fprintf('Loading compounds from database...');
      [dpk,name]=mysql('select domain,name from compounds.domains');
      for i=1:length(dpk)
        obj.domains{dpk(i)}=name{i};
      end
      [compound,name,inchikey,smiles,formula,mass]=mysql('SELECT compound,name,inchikey,smiles,formula,monoisotopicMass FROM compounds.compounds');
      obj.c=struct('compound',num2cell(compound),'name',name,'inchikey',inchikey,'smiles',smiles,'formula',formula,'mass',num2cell(mass),'aliases',{cell(1,length(dpk))});
      [compound,domain,name]=mysql('select compound,domain,name from compounds.names');
      fprintf('%d aliases...',length(compound));
      for i=1:length(obj.c)
        if mod(i,1000)==0
          fprintf('%d...',i);
        end
        sel=obj.c(i).compound==compound;
        obj.c(i).aliases(domain(sel))=name(sel);
      end
      fprintf('done\n');
    end
    
    function c=get(obj,pk)
      c=[];
      for i=1:length(pk)
        c=[c,obj.c([obj.c.compound]==pk(i))];
      end
      if length(c)~=length(pk)
        error('Only found %d/%d requested compounds', length(c), length(pk));
      end
    end

    function c=find(obj,name,domain)
      c=[];
      if ischar(name)
        name={name};
      end
      if nargin<3
        for i=1:length(name)
          c=[c,obj.c(strcmp({obj.c.name},name{i}))];
        end
      else
        % Search in domain
        if ischar(domain)
          domain=find(strcmp(domain,obj.domains));
          if isempty(domain)
            error('Bad domain: %s', domain);
          end
        else
          if domain<1 || domain>length(obj.domains)
            error('Bad domain: %d', domain);
          end
        end
        for i=1:length(name)
          sel=find(arrayfun(@(z) strcmp(z.aliases{domain},name{i}),obj.c));
          if length(sel)==0
            error('Unable to locate %s in domain %s', name{i}, obj.domains{domain});
          elseif length(sel)>1
            error('Ambiguous: %s in domain %s', name{i}, obj.domains{domain});
          end
          c=[c,obj.c(sel)];
        end
      end
      if length(c)~=length(name)
        error('Only found %d/%d requested compounds', length(c), length(name));
      end
    end

    function getinfo(obj,pk)
      if ischar(pk)
        c=find(obj,pk);
      else
        c=get(obj,pk);
      end
      res=c;
      mixtures=Mixtures.instance();
      mixes={};
      for i=1:length(mixtures.m)
        if ismember(c.compound,mixtures.m(i).contents)
          mixes{end+1}=mixtures.m(i).name;
        end
      end
      if ~isempty(mixes)
        res.mixtures=strjoin(mixes,', ');
      end
      details(res);
    end
    
    function n=getalias(obj,domain,pklist)
    % Get aliases for given pklist in domain
      query=sprintf('SELECT compound,name FROM compounds.names WHERE domain=(SELECT d.domain FROM compounds.domains d WHERE d.name=''%s'') AND compound IN (%s0)',domain,sprintf('%d,',pklist));
      [pk,name]=mysql(query);
      if length(unique(pklist))~=length(pk)
        error('Only found %d/%d unique aliases under domain %s',length(pk),length(unique(pklist)),domain);
      end
      n=cell(length(pklist),1);
      for i=1:length(pklist)
        n{i}=name{pk==pklist(i)};
      end
    end
    
    function savehtml(obj,dirname,pk,varargin)
    % Save images for compounds with given pk list (as grid of images), or just one
    % Title each group with given title(s) or index if title not provided
    % If images given, position those under each title
    % Caption each with compound name, or captions, if given, captions instead
    % Pk list can be a nested cell array, in which case captions and title should have the same dimension as the top-level cell array
      defaults=struct('captions','','title','','images',[],'maintitle','');
      args=processargs(defaults,varargin);
      if ~isdir(dirname)
        mkdir(dirname);
      end
      htmlfile=[dirname,'/index.html'];
      fd=fopen(htmlfile,'w');
      if fd<0
        error('Unable to open %s for writing.', htmlfile);
      end
      header=['<!DOCTYPE HTML>\n',...
             '<html>\n',...
              '<head>\n',...
              '<style>\n',...
              'figure{display: inline-block;text-align: center;}\n',...
              'figcaption{text-align: center;}\n',...
              '@media print { h1 {page-break-before:auto;} div{break-inside: avoid;}};\n',...
              '</style>\n',...
              '</head>\n',...
              '<body>\n'];
      fprintf(fd,header);
      if ~isempty(args.maintitle)
        fprintf(fd,'<h1>%s</h1>\n',args.maintitle);
      end
      
      if ~iscell(pk)
        pk={pk};
        args.captions={args.captions};
        args.title={args.title};
      end
      for j=1:length(pk)
        fprintf(fd,'<div class="mygroup">\n');
        if ~isempty(args.title) && ~isempty(args.title{j})
          fprintf(fd,'<h2>%s</h2>\n', args.title{j});
        else
          fprintf(fd,'<h2>Group %d</h2>\n',j);
        end
        if ~isempty(args.images)
          im=args.images{j};
          imname=sprintf('im%d.png',j);
          imwrite(im,[dirname,'/',imname]);
          fprintf(fd,'<IMG src="%s"><br>\n',imname);
        end
        c=obj.get(pk{j});
        for i=1:length(c)
          imgfile=sprintf('%d.svg',c(i).compound);
          fprintf(fd,'<figure>\n');
          fprintf(fd,' <IMG src="%s" alt="%s">\n',imgfile,c(i).name);
          fprintf(fd,' <figcaption>');
          fprintf(fd,'<center>%s</center><br>\n',c(i).name);
          if ~isempty(args.captions) && ~isempty(args.captions{j})
            fprintf(fd,'%s\n',args.captions{j}{i});
          end
          fprintf(fd,' </figcaption>\n');
          fprintf(fd,'</figure>\n');
          obj.saveimage([dirname,'/',imgfile],c(i).compound,'caption','');
        end
        fprintf(fd,'</div>\n');
      end
      fprintf(fd,'</body>');
      fprintf(fd,'</html>');
      fclose(fd);
      fprintf('Saved compound plots in %s/index.html\n', dirname);
    end
    
    function saveimage(obj,filename,pk,varargin)
    % Save image for compounds with given pk list as file 
    % Caption with compound name, or captions, if given, captions instead
    % Format can be svg, eps, or png
    % Width in pixels for png only
      defaults=struct('format','','caption','','width',1200);
      args=processargs(defaults,varargin);

      if isempty(args.format)
        sp=strsplit(filename,'.');
        args.format=sp{end};
      end
      c=obj.get(pk);
      base='/tmp/cplot';
      canfile=[base,'.can'];
      fd=fopen(canfile,'w');
      if isempty(args.caption)
        nm=c.name;
      else
        nm=args.caption;
      end
      fprintf(fd,'%s %s\n',c.smiles,nm);
      fclose(fd);
      if strcmp(args.format,'svg')
        system(sprintf('/opt/local/bin/obabel -i can %s -o svg -O %s --title "" 2>&1 | grep -v "molecule converted"',canfile,filename));
      else
        svgfile=[base,'.svg'];
        system(sprintf('/opt/local/bin/obabel -i can %s -o svg -O %s --title "" 2>/dev/null',canfile,svgfile));
        if strcmp(args.format,'eps')
          system(sprintf('/Applications/Inkscape.app/Contents/MacOS/inkscape %s -o %s --export-ignore-filters',svgfile,filename));
        elseif strcmp(args.format,'png')
          system(sprintf('/Applications/Inkscape.app/Contents/MacOS/inkscape -w %d %s -o %s',args.width,svgfile,filename));
        else
          error('Bad format: %s',args.format);
        end
      end
    end
    
  end
end