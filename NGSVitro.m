% Class for NGS in vitro data analysis
classdef NGSVitro < matlab.mixin.Copyable
  properties
    run;
    created;   % Date this record was original created/loaded
    updated;   % Date of last update from DB
    subruns; % Information on the subruns themselves
    seqs;    % NGSSeq data structure
    subrun; % subrun for each entry
    naseq;  % sequence ID for each entry
    cnt;  % count of each entry
    trpexpts;    % TRP setup
    rundata; 	  % Sequencing statistics
    basespace;	  % Basespace run data
  end
  
  methods(Access = protected)
     % Override copyElement method:
     function cpObj = copyElement(obj)
     % Make a shallow copy of all the properties
       cpObj = copyElement@matlab.mixin.Copyable(obj);
       % Make a deep copy of the NASeqs, TRPExpts
       cpObj.seqs = copy(obj.seqs);
       if ~isempty(obj.trpexpts)
         cpObj.trpexpts=copy(obj.trpexpts);
       end
     end
   end
   
  methods(Static)
    function arc(x,y,col,wid)
    % Draw an arc through the 3 points given
      xs=linspace(x(1),x(3));
      ys=(xs-x(1)).*(xs-x(3))*y(2)./((x(2)-x(1)).*(x(2)-x(3)));
      plot(xs,ys,col,'LineWidth',wid);
    end
    
    function alld=loadall(runs)
    % Load multiple datasets.  e.g. alld=NGSVitro.loadall({'NGS22','NGS24'})
      basedir=sprintf('%s/Dropbox/Synbio/NGS',getenv('HOME'));

      if nargin<1
        z=dir([basedir,'/*/analyze/data.mat']);
        runs={};
        for i=1:length(z)
          slash=find(z(i).folder=='/');
          runs{i}=z(i).folder(slash(end-1)+1:slash(end)-1);
        end
      end
      fname={};
      fprintf('Loading runs: ');
      for i=1:length(runs)
        fname{i}=sprintf('%s/%s/analyze/data.mat',basedir,runs{i});
        if ~exist(fname{i})
          error('Not found: %s',fname{i});
        end
        fprintf('%s ', runs{i});
      end
      fprintf('\n');
      alld={};
      for i=1:length(runs)
        alld{i}=NGSVitro(runs{i});
        alld{i}.load(fname{i});
      end
    end
    
    function showallstats(d,naseq)
    % Can do something like NGSVitro.showallstats({ngs22,ngs24},naseq)
      for i=1:length(d)
        fprintf('%s.showstats(%d)\n', d{i}.run, naseq);
        d{i}.showstats(naseq);
      end
    end
    
    function res=dballruns(varargin)
    % Query database for all subruns with descr like regexp
    % Make sure NGSSummary is up to date
    % To modify a previous entry, need to do a manual SQL: delete from NGSSummary.subruncnts where run='xx';
      defaults=struct('minrun',1,'descr','%','code','%','minreads',0);
      args=processargs(defaults,varargin);

      NGSDatabase.open();
      allruns=mysql('SHOW DATABASES');
      allruns=allruns(strncmp(allruns,'NGS',3));
      entered=mysql('SELECT DISTINCT run FROM NGSSummary.subruncnts');
      allruns=setdiff(allruns,entered);
      [val,ord]=sort(cellfun(@(z) sscanf(['0',z(4:end)],'%d'),allruns),'desc');
      allruns=allruns(ord(val>0));
      fprintf('Have %d/%d NGS runs to update\n', length(allruns), length(allruns)+length(entered));
      for i=1:length(allruns)
        fprintf('Updating %s...',allruns{i});
        if val(i)>=55
          hasRC='rc,';
        else
          hasRC='';
        end
        cmd=sprintf('insert into NGSSummary.subruncnts(run,subrun,descr,code,indexseq,%sleftcontext,rightcontext,seqs,rds) select ''%s'',subrun,descr,code,indexseq,%sleftcontext,rightcontext,seqs,rds FROM %s.v_subruncnts WHERE rds IS NOT NULL',hasRC, allruns{i},hasRC, allruns{i});
        tic;mysql(cmd);el=toc;
        fprintf('done in %.1f sec\n', el);
      end
      % Now query the data
      cmd=sprintf('select run,descr,code,subrun,seqs,rds from NGSSummary.subruncnts WHERE code like ''%s'' and descr like ''%s'' AND CONVERT(SUBSTRING(run, 4), SIGNED INTEGER)>=%d AND rds>=%d',args.code,args.descr,args.minrun,args.minreads)
      [run,descr,code,subrun,seqs,rds]=mysql(cmd);
      res=struct('run',{run},'descr',{descr},'code',{code},'subrun',subrun,'seqs',int32(seqs),'rds',int32(rds));
    end
  end

  methods
    function obj=NGSVitro(r)
      obj.run=r;
      obj.created=date();
      obj.seqs=NGSSeq(r);
    end

    function load(obj,file)
    % Backward compatibility with old data.mat files
      if nargin<2
        file='data.mat';
      end
      fprintf('Loading data from %s...',file);
      d=load(file);
      fprintf('done\n');
      fn=fieldnames(d);
      if length(fn)==1
        d=d.(fn{1});
        fn=fieldnames(d);
      end
      for i=1:length(fn)
        obj.(fn{i})=d.(fn{i});
      end
      obj.cleanup();
    end
    
    function save(obj,file)
    % Save data in a .mat file
      if nargin<2
        file='data.mat';
      end
      tic
      %      data=rmfield(struct(obj),'segments');
      data=struct(obj);
      seqs=struct(data.seqs); data=rmfield(data,'seqs');
      save(file,'data','seqs');
      toc
    end
    
    function loadrundata(obj,jsonfile)
    % Load run data dumped using basespace CLI (bs get run --name xx )
      if nargin<2 || isempty(jsonfile)
        jsonfile='rundata.json';
      end
      fd=fopen(jsonfile,'r');
      data=fread(fd,inf,'uint8=>char');
      fclose(fd);
      j=jsondecode(data);
      obj.basespace=j;
      % Check sequencing stats match rundata
      st=obj.basespace.SequencingStats;
      rd=obj.rundata;
      if (rd.readsTotal ~= st.ReadsTotal/3)
        fprintf('rundata.reads=%d, stats=%d\n', rd.readsTotal, st.ReadsTotal/3);
      end
      if (rd.readsPF ~= st.ReadsPfTotal/3)
        fprintf('rundata.readsPF=%d, stats=%d\n', rd.readsPF, st.ReadsPfTotal/3);
      end
    end
    
    function cleanup(obj)
    % Add computed fields where missing
      if ~isfield(obj.subruns,'codes')
        % Old version(eg NGS19) didn't have codes field
        fprintf('Adding empty subruns.codes\n');
        obj.subruns.codes=obj.subruns.desc;
        for i=1:length(obj.subruns.desc)
          obj.subruns.codes{i}='?';
        end
      end
      if ~isfield(obj.subruns,'rds')
        % Old versions didn't read total rds from objbase
        fprintf('Adding subruns.rds\n');
        NGSDatabase.open(obj.run);
        [sr,rds]=mysql(sprintf('SELECT subrun,rds FROM v_subruncnts WHERE descr is not null ORDER BY descr'));
        obj.subruns.rds=nan(size(obj.subruns.subrun));
        for i=1:length(sr)
          obj.subruns.rds(obj.subruns.subrun==sr(i))=rds(i);
        end
      end
      if ~isa(obj.trpexpts,'TRPExpt')
        % Old files had a struct for trpexpts
        te=[];
        for i=1:length(obj.trpexpts)
          te=[te,TRPExpt(obj.trpexpts(i))];
          te(i).run=obj.run;
        end
        obj.trpexpts=te;
      end
      obj.sortbysubrun();
      obj.seqs.cleanup();
    end
    
    function dbload(data,mincnt)
    % Load data from database
      if nargin<2
        mincnt=10;
        fprintf('Loading sequences that have a total of at least %d reads\n',mincnt);
      end
      if isempty(data.seqs) || isempty(data.seqs.naseq)
        data.seqs=NGSSeq(data.run);
        data.seqs.dbload(mincnt);
      else
        fprintf('Skipping reload of seqs\n');
      end
      % mysql db still open...
      NGSDatabase.open(data.run);
      
      % Run data
      data.dbloadrundata();
      
      % subruns for libraries
      fprintf('Get subruns...');
      sel='descr IS NOT NULL and descr!=''PhiX''';
      [sr,libnames,codes,rds]=mysql(sprintf('SELECT subrun,descr,code,rds FROM v_subruncnts WHERE %s ORDER BY descr',sel));
      fprintf('done (%d unique libraries)\n',length(unique(libnames)));
      
      if mincnt>1
        % Build list of naseqs to load
        fprintf('Find sequences with >= %d reads in at least one subrun ...',mincnt);
        mysql(sprintf('CREATE TEMPORARY TABLE ngs.t_seqs AS SELECT DISTINCT naseq FROM ngsentries WHERE cnt>=%d', mincnt));
        fprintf('done');
        numseq=mysql('SELECT COUNT(*) FROM ngs.t_seqs');
        fprintf('(%d)\n',numseq);
        fprintf('Get num entries...');
        numEntries = mysql(sprintf('SELECT COUNT(*) FROM ngsentries WHERE subrun IN (SELECT subrun FROM subruns WHERE %s) AND naseq IN (SELECT naseq FROM ngs.t_seqs)',sel));
        fprintf('%d\n',numEntries);
      
        cmd=sprintf(['SELECT subrun, n.naseq, cnt '...
                     'FROM ngsentries n '...
                     'WHERE subrun IN (SELECT subrun FROM subruns WHERE %s) '...
                     'AND n.naseq IN (SELECT naseq FROM ngs.t_seqs)'],sel);
      else
        fprintf('Get num entries...');
        numEntries = mysql(sprintf('SELECT COUNT(*) FROM ngsentries WHERE subrun IN (SELECT subrun FROM subruns WHERE %s)',sel));
        fprintf('%d\n',numEntries);
      
        cmd=sprintf(['SELECT subrun, n.naseq, cnt '...
                     'FROM ngsentries n '...
                     'WHERE subrun IN (SELECT subrun FROM subruns WHERE %s) '],sel);
      end
      
      fprintf('Loading %d libs with %d entries...',length(libnames),numEntries);
      [subrun,naseq,cnt]=mysql(cmd);
      fprintf('done\n');

      % TODO: merge subruns with same descr and code
      % would eliminate need to merge in database
      % but would make setup of trpexpts in database a little trickier (maybe only one entry pointing to minimum subrun? )
      % Ideally subrun would cover multiple contexts for purposes of loading to matlab and trpexpts, but would be a subsubrun in terms of contexts
      % Or maybe make (subrun,leftcontext,rightcontext) the primary key in subruns inread of just subrun (but what about ngsentries)
      
      
      data.dbloadtrpexpts();

      assert(all(ismember(unique(subrun),sr)));

      data.subruns=struct('subrun',sr,'desc',{libnames},'codes',{codes},'rds',{rds});
      data.subrun=subrun;
      data.naseq=naseq;
      data.cnt=cnt;
      data.cleanup();
      data.dbloadnormalization();
      data.created=date();
      data.updated=data.created();
    end
    
    function dbupdate(data)
    % Update information from database
      data.seqs.dbupdate();
      if strcmp(data.run,'ngs')
        % When not connected to a specific database can't run things that depend upon ngsentries
        fprintf('Skipping updates that depend on a specific NGS run\n');
      else
        data.dbloadrundata();
        data.dbloadnormalization();
      end
      data.updated=date();
    end
    
    function tbls=dbschemastats(data)
    % Find information on all schemas (runs) in db
      NGSDatabase.open();
      db=mysql('SHOW DATABASES LIKE ''NGS%''');
      tbls=table('RowNames',db);
      for i=1:length(db)
        t=mysql(sprintf('SHOW TABLES IN %s',db{i}));
        for j=1:length(t)
          tbls(i,t{j})={1};
        end
      end
    end
      
    function dbloadrundata(data)
    % Load run data from database
      NGSDatabase.open(data.run);
      fprintf('Loading run data...');
      data.rundata=[];
      [run,descr,sequenced,readsTotal,readsPF,density]=mysql('select run,descr,sequenced,readsTotal,readsPF,density from runs');
      for i=1:length(run)
        paired=mysql(sprintf('select sum(cnt) from subruns s, ngsentries e where e.subrun=s.subrun and s.run=%d and indexbarcode not in (select indexbarcode from indexbarcodes where seq in (''S0'',''S97'') OR seq LIKE ''%%U'' )',run(i)));
        [readnum,ncycles,q30frac,aligned,errorRate,intensity]=mysql(sprintf('select readnum,ncycles,q30frac,aligned,errorRate,intensity from runreads where run=%d',run(i)));
        rr=[];
        for k=1:length(readnum)
          rr=[rr,struct('readnum',readnum(k),'ncycles',ncycles(k),'q30frac',q30frac(k),'aligned',aligned(k),'errorRate',errorRate(k),'intensity',intensity(k))];
        end
        [indexbarcode,readsIdentified]=mysql(sprintf(['select ' ...
                            'b.seq,i.readsIdentified from ' ...
                            'indexstats i, indexbarcodes b where ' ...
                            'i.indexbarcode=b.indexbarcode and run=%d'],run(i)));
        indexstats=struct('indexbarcode',indexbarcode,'readsIdentified',num2cell(readsIdentified));
        data.rundata=[data.rundata,struct('run',run(i),'descr',descr{i},'readsTotal',readsTotal(i),'readsPF',readsPF(i),'readsPaired',paired,'density',density(i),'runreads',rr,'indexstats',indexstats)];
      end
      fprintf('done\n');
    end
    
    function x=dbpearstatistics(data)
      % Compare index data from database
      NGSDatabase.open(data.run);
      [indexseq,rds,rdsid]=mysql('SELECT b.seq,sum(cnt),st.readsidentified FROM indexbarcodes b, subruns s, ngsentries e, indexstats st WHERE b.indexbarcode=st.indexbarcode AND b.indexbarcode=s.indexbarcode AND e.subrun=s.subrun GROUP BY b.seq,readsidentified;');
      rds(end+1)=sum(rds);
      rdsid(end+1)=sum(rdsid);
      indexseq{end+1}='Total';
      x=table;
      x.IndexSeq=indexseq;
      x.ReadsIndexed=rdsid;
      x.ReadsPaired=rds;
      x.PctPaired=rds*100./rdsid;
    end
      
    function dbloadnormalization(data)
      % Add normalization
      NGSDatabase.open(data.run);
      fprintf('Loading normalization...');
      sel='descr IS NOT NULL';
      cmd=sprintf(['SELECT n.subrun, n.naseq, n.concentration, e.cnt '...
                   'FROM ngsentries e, normalizations n '...
                   'WHERE n.naseq=e.naseq '...
                   'AND n.subrun=e.subrun '...
                   'AND n.subrun IN (SELECT subrun FROM subruns WHERE %s)'],sel);
      [norm_subrun,norm_naseq,norm_conc,norm_cnt]=mysql(cmd);
      fprintf('%d...done\n',length(norm_subrun));
      norm=cell(length(data.subruns.subrun),3);
      data.subruns.normalization=[];
      for i=1:length(data.subruns.subrun)
        data.subruns.normalization=[data.subruns.normalization,struct('naseq',[],'conc',[],'cnt',[],'predil',1.0)];
      end
      for n=1:length(norm_subrun)
        i=find(norm_subrun(n)==data.subruns.subrun);
        if isempty(i)
          fprintf('Have normalization data for subrun %d, which is not in data.subruns\n', norm_subrun(n));
        else
          data.subruns.normalization(i).naseq(end+1)=norm_naseq(n);
          data.subruns.normalization(i).conc(end+1)=norm_conc(n);
          data.subruns.normalization(i).cnt(end+1)=norm_cnt(n);
        end
      end

      % Load predil values if available (migrated into DB at NGS59)
      try 
        [predil_sr,predil_val]=mysql('SELECT subrun, predil FROM subruns');
        for i=1:length(data.subruns.normalization)
          sel=predil_sr==data.subruns.subrun(i);
          data.subruns.normalization(i).predil=predil_val(sel);
        end
      catch me
        % no predil field, ignore
        fprintf('Warning: No predil field in %s.subruns, assuming 1.0\n',data.run);
      end
    end

    function x=findrefs(data,varargin)
    % Locate reference sequences in data
      defaults=struct('filename',[]);
      args=processargs(defaults,varargin);
      x=data.showcounts(data.refseqs(0),'filename',args.filename);
    end
    
    function x=showcounts(data,naseqs,varargin)
    % Show table of counts of given naseqs (columns) vs subruns (rows)
      defaults=struct('filename',[],'subruns',[],'sortrows',true);
      args=processargs(defaults,varargin);

      if isempty(args.subruns)
        subruns=unique(data.subrun(ismember(data.naseq,naseqs)));
      else
        subruns=args.subruns;
      end
      fprintf('Building %dx%d table\n', length(subruns), length(naseqs));
      x=table();
      x.subrun=subruns(:);
      x.totalreads=int32(arrayfun(@(z) sum(data.subruns.rds(data.subruns.subrun==z)), x.subrun));
      for j=1:length(subruns)
        ind=data.subruns.subrun==subruns(j);
        x.Properties.RowNames{j}=sprintf('%s-%s',data.subruns.desc{ind},data.subruns.codes{ind});
      end
      for i=1:length(naseqs)
        z=nan(size(x,1),1);
        for j=1:length(subruns)
          ind=data.subrun==subruns(j) & data.naseq==naseqs(i);
          z(j)=sum(data.cnt(ind));
        end
        x(:,end+1)=num2cell(z);
        name=data.seqs.getname(naseqs(i));
        if ~isempty(name)
          x.Properties.VariableNames{width(x)}=strrep(name,'-','_');
        else
          x.Properties.VariableNames{width(x)}=sprintf('N%d',naseqs(i));
        end
        x.Properties.VariableDescriptions{width(x)}=sprintf('%d %s',naseqs(i),data.seqs.getlabels(naseqs(i)));
      end
      [~,ord]=sort(x.Properties.VariableNames(2:end));
      x=x(:,[1,ord+1]);
      if args.sortrows
        x=sortrows(x,'RowNames');
      end
      % Would be nice to add a row with descriptions, but gives an error...
      %x(end+1,2:end)=x.Properties.VariableDescriptions(2:end);
      %x=x([end,1:end-1],:);
      x(end+1,:)=num2cell(sum(table2array(x)));
      x(end,1)={nan};
      x.Properties.RowNames{end}='Total';
      x.Properties.DimensionNames={'Subrun','NASeq'};
      % Remove columns with 0 count
      x=x(:,table2array(x(end,:))~=0);
      if ~isempty(args.filename)
        writetable(x,args.filename,'WriteRowNames',true);
        fprintf('Wrote table to %s\n', args.filename);
      end
    end
    
    function dbgetcodelengths(data)
    % Augment data.subruns with lengths of barcodes
      NGSDatabase.open(data.run);
      cmd='select subrun,length(leftcontext)+length(rightcontext) from v_subruns';
      [sr,len]=mysql(cmd);
      bclen=nan(size(data.subruns.subrun));
      for i=1:length(sr)
        sel=data.subruns.subrun==sr(i);
        if isnan(len(i))
          bclen(sel)=0;
        else
          bclen(sel)=len(i);
        end
      end
      data.subruns.bclen=bclen;
    end
    
    function dbloadtrpexpts(data)
      % Get trp experiment setups from database
      NGSDatabase.open(data.run);
      fprintf('Loading TRP experiments...');
      cmd=sprintf(['SELECT trpexpt, descr, cleaved, uncleaved, insr, bulkcleavage, target, targetconc, concunits, duration, idgroup, robotrun, mixture '...
                   'FROM trpexpts']);
      [trpexpt, descr, cleaved, uncleaved, insr, bulkcleavage, target, targetconc, concunits, duration, idgroup, robotrun, mixture]=mysql(cmd);
      data.trpexpts=[];
      for i=1:length(trpexpt)
        if isfinite(idgroup(i))
          idsr=mysql(sprintf('SELECT subrun FROM ids WHERE idgroup=%d',idgroup(i)));
        else
          idsr=[];
        end
        t=TRPExpt();
        t.run=data.run;
        t.trpexpt=trpexpt(i);
        t.descr=descr{i};
        t.cleavedsr=cleaved(i);
        t.uncleavedsr=uncleaved(i);
        t.insr=insr(i);
        t.idsr=idsr;
        t.bulkcleavage=bulkcleavage(i);
        t.target=target{i};
        t.targetconc=targetconc(i);
        t.concunits=concunits{i};
        t.duration=duration(i);
        if isnan(t.targetconc)
          t.targetconc=[];
        end
        t.robotrun=robotrun(i);
        t.mixture=mixture(i);
        data.trpexpts=[data.trpexpts,t];
      end
      fprintf('done\n');
    end
    
    function sortbysubrun(data)
    % Sort data by subrun description
    % Put leading 0 in Rn  descriptions
      dtmp=data.subruns.desc;
      for i=1:length(dtmp)
        if length(dtmp{i})==2
          dtmp{i}=[dtmp{i}(1),'0',dtmp{i}(2)];
        end
      end

      n=length(data.subruns.subrun);
      [~,ord]=sort(dtmp);
      fn=fieldnames(data.subruns);
      for i=1:length(fn)
        if size(data.subruns.(fn{i}),1)==n
          data.subruns.(fn{i})=data.subruns.(fn{i})(ord,:);
        end
      end
    end

    function segdistribution(data)
      % Display frequent mismatches
      for i=1:size(data.seqs.segments,2)
        segments=data.segments();
        sel=cellfun(@(z) ischar(z), segments(:,i)) & data.cnt>1;
        cnt=data.cnt(sel);
        missing=sum(data.cnt(~sel));
        [useq,ia,ic]=unique(segments(sel,i));
        fprintf('Group %d sequences - %d distinct\n', i,length(useq));
        unseq=hist(ic,1:max(ic));
        ucnt=nan(length(useq),1);
        for j=1:length(useq)
          ucnt(j)=sum(cnt(ic==j));
        end
        [ucnt,ord]=sort(ucnt,'descend');
        useq=useq(ord);
        unseq=unseq(ord);
        fprintf('    Cnt     NSeq Frac Seq\n');
        fprintf('%7d %7s %4.1f%% Missing\n', missing,'',missing/sum(data.cnt)*100);
        for j=1:min(20,length(useq))
          fprintf('%7d %7d %4.1f%% %s', ucnt(j),unseq(j), ucnt(j)/sum(data.cnt)*100,useq{j});
          if ismember(useq{j},NGSSeq.pts{i})
            fprintf(' *');
          end
          fprintf('\n');
        end
      end
    end
    
    function [conc,cnt]=getconc(data,subruns)
    % Get concentration of given subruns in M
      conc=nan(size(subruns));
      for i=1:length(subruns)
        sr=subruns(i);
        srsel=find(data.subruns.subrun==sr);
        if isempty(srsel)
          error('Subrun %d not found\n', sr);
        end
        rcnt=sum(data.subruns.normalization(srsel).cnt);
        sel=data.subrun==sr & ~ismember(data.naseq,data.subruns.normalization(srsel).naseq);
        cnt(i)=sum(data.cnt(sel));
        conc(i)=sum(data.subruns.normalization(srsel).conc)/rcnt*(data.subruns.rds(srsel)-rcnt);
      end
    end
    
    function x=overallstats(data,ind,varargin)
    % Show overall stats for the entire data structure or some subpart
      defaults=struct('increfonly',false,'oldway',false,'showall',false,'matches',[],'excludecodes',{{'BW/X','BZ/X'}});
      args=processargs(defaults,varargin);

      res=[];
      totals.desc='Total';
      if args.oldway
        fprintf('Ind Subrun Descr                     Seqs    AllRds  Ld-Ref  RefRds Conc(nM) Len\n');
      end
      [~,ord]=sort(data.subruns.desc);
      rtotal=0;
      conctotal=0;
      singleSeq=false;   % True if running against a single sequence
      if sum(data.subruns.rds) > sum(data.cnt)
        fprintf('Not all reads loaded; using all reads data in .subruns.rds for counts\n');
        useAllReads=true;
      else
        useAllReads=false;
      end
      if nargin<2
        ind=[];
      end
      if ~isempty(args.matches)
        if isempty(ind)
          ind=true(size(data.naseq));
        end
        ind=ind&ismember(data.subrun,data.subruns.subrun(cellfun(@(z) ~isempty(z),regexp(data.subruns.desc,args.matches))));
        if ~any(ind)
          error('No match');
        end
      end
      if ~isempty(ind) 
        fprintf('Concs are for listed sequence counts\n');
        useAllReads=false;
        if length(unique(data.naseq(ind)))==1
          singleSeq=true;
        end
      end
      if ~isempty(args.excludecodes)
        excsr=data.subruns.subrun(ismember(data.subruns.codes,args.excludecodes));
        nexclude=sum(data.cnt(ismember(data.subrun,excsr)));
        if nexclude>0
          fprintf('Excluding %d reads for codes [%s]\n', nexclude, strjoin(args.excludecodes,','));
        end
      end
      for ii=1:length(data.subruns.subrun)
        i=ord(ii);
        if ismember(data.subruns.codes{i},args.excludecodes)
          continue;
        end
        sr=data.subruns.subrun(i);
        if ~isempty(ind) && ~any(data.subrun(ind)==sr)
          continue;
        end
        rcnt=sum(data.subruns.normalization(i).cnt);
        sharedSubruns=strcmp(data.subruns.desc,data.subruns.desc{i});  % Same name; assume they are from same barcoding PCR, but with different prefixes
        sharedCnt=sum([data.subruns.normalization(sharedSubruns).cnt]);
        sharedConc=sum([data.subruns.normalization(sharedSubruns).conc]);
        sel=data.subrun==sr;
        if ~singleSeq
          % Don't count references (or any sequences in same cluster as refs) unless we're looking at a single sequence
          sel = sel & ~ismember(data.naseq,data.refseqs(1));
          % Note that neither 'norm' reads nor 'ref' reads include mutants (same cluster) of refs
        end
        if nargin>=2 && ~isempty(ind)
          sel=sel&ind;
        end
        cnt=sum(data.cnt(sel));
        
        if rcnt<sharedCnt/10 && rcnt<20
          fprintf('Not enough reference counts (%d) for %s-%s, using %d shared ones over %d subruns\n', rcnt, data.subruns.desc{i}, data.subruns.codes{i}, sharedCnt, sum(sharedSubruns));
          rcnt=nan;
          if ~useAllReads
            conc=sharedConc/sharedCnt*cnt*data.subruns.normalization(i).predil;
          else
            conc=sharedConc/sharedCnt*data.subruns.rds(i)*data.subruns.normalization(i).predil;
          end
        else
          if ~useAllReads
            conc=sum(data.subruns.normalization(i).conc)/rcnt*cnt*data.subruns.normalization(i).predil;
          else
            conc=sum(data.subruns.normalization(i).conc)/rcnt*(data.subruns.rds(i)-rcnt)*data.subruns.normalization(i).predil;
          end
          rtotal=rtotal+rcnt;
        end
        conctotal=conctotal+conc;
        len=sum(data.cnt(sel).*data.length(sel))/cnt;
        
        if (args.increfonly && data.subruns.rds(i)>0) || cnt>0 || ~isfinite(data.subruns.rds(i))
          % Find res entry with this name
          if ~isempty(res)
            resind=find(strcmp({res.desc},data.subruns.desc{i}));
          else
            resind=[];
          end
          if isempty(resind)
            res(end+1).desc=data.subruns.desc{i};
            resind=length(res);
          end
          res(resind).predil=sprintf('%5.0f',data.subruns.normalization(i).predil);
          if ~isfield(res(resind),'total') || isempty(res(resind).total)
            res(resind).total=0;
          end

          codename=strrep(strrep(strrep(data.subruns.codes{i},'/X1',''),'/X2',''),'/X','');
          codename=strrep(codename,'/','_');
          %          res(resind).(codename)=struct('subrun',sr,'nseqs',sum(sel), 'nreads', data.subruns.rds(i), 'nonrefrds',cnt,'nrefrds', rcnt,'conc',conc*1e9,'len',len);
          if singleSeq
            res(resind).(codename)=sprintf('%4d %8d %7.0f %8.3g',sr, cnt,rcnt,conc*1e9);
          else
            res(resind).(codename)=sprintf('%4d %6d %8d %7.0f %8.3f %3.0f',sr, sum(sel), cnt,rcnt,conc*1e9,len);
          end
          res(resind).total=res(resind).total+conc;
          if ~isfield(totals,codename)
            totals(1).(codename)=cnt;
            totals(2).(codename)=nansum([0,rcnt]);
          else
            totals(1).(codename)=totals(1).(codename)+cnt;
            totals(2).(codename)=nansum([totals(2).(codename),rcnt]);
          end
          reslen=length(res(resind).(codename));  % For filling in empty entries below
          if args.oldway
            fprintf('%3d %4d %-25.25s %6d %8d %8d %7.0f %8.3g %3.0f\n',i, sr, [data.subruns.desc{i},'-',data.subruns.codes{i}], sum(sel), data.subruns.rds(i), cnt,rcnt,conc*1e9,len);
          end
        end
      end
      if length(res)>1
        res(end+1).desc='Total';
        fn=fieldnames(totals);
        for i=1:length(fn)
          if strcmp(fn{i},'desc') || strcmp(fn{i},'predil') || strcmp(fn{i},'total')
            continue;
          end
          if singleSeq
            res(end).(fn{i})=sprintf('     %8d %7.0f         ', [totals.(fn{i})]);
          else
            res(end).(fn{i})=sprintf('            %8d %7.0f             ',[totals.(fn{i})]);
          end
        end
      end
      if args.oldway
        fprintf('Total    %-23.23s %8d %8d %8d %7.0f %8.0f\n','Total', length(unique(data.naseq)), sum(data.subruns.rds), sum(data.cnt),rtotal,conctotal*1e9);
        fprintf('References are %.02f%% of total reads\n', 100*rtotal/sum(data.subruns.rds));
      end

      % Convert total to strings
      for j=1:length(res)
        res(j).total=sprintf('%8.3g',res(j).total*1e9);
      end
      % Clean up result struct
      fn=fieldnames(res);
      for i=1:length(fn)
        for j=1:length(res)
          if isempty(res(j).(fn{i}))
            res(j).(fn{i})=blanks(reslen);
          end
        end
      end
      
      % Build header
      header=res(1); 
      for i=1:length(fn)
        if ~strcmp(fn{i},'desc') && ~strcmp(fn{i},'predil') && ~strcmp(fn{i},'total')
          if singleSeq
            header.(fn{i})='  SR      Rds     Ref     Conc';
          else
            header.(fn{i})='  SR   Seqs     Norm     Ref     Conc Len';
          end
        end
      end
      header.total='   [Str]';
      header.desc='-';
      x=struct2table([header,res]);
      x.Properties.RowNames=x.desc;
      x=removevars(x,'desc');
      if ~args.showall
        x=removevars(x,intersect(x.Properties.VariableNames,{'BW1','BZ1'}));
      end
      x=movevars(x,sort(x.Properties.VariableNames(3:end)),'After',2);
      packtable(x);

      if nargout==0   % don't assign to x if nargout==0 to supress display of 'ans'
        clear x;
      end
    end
    
    function ngslosses(data)
    % Show data losses in processing
      for i=1:length(data.rundata)
        stages={}; reads=[]; labels={};
        rr=data.rundata(i);
        stages{end+1}='Raw Reads'; reads(end+1)=rr.readsTotal;
        stages{end+1}='PF'; labels{end+1}='No PF';reads(end+1)=rr.readsPF;
        if length(rr.indexstats)>1
          isphix=strcmp({rr.indexstats.indexbarcode},'S97');
          matched=~strcmp({rr.indexstats.indexbarcode},'S0');
          identified=sum([rr.indexstats(matched).readsIdentified]);
          notphix=sum([rr.indexstats(matched&~isphix).readsIdentified]);
        else
          isphix=[];
          identified=sum([rr.indexstats.readsIdentified]);  
        end
        if ~any(isphix) && any(isfinite([rr.runreads([1,end]).aligned]))
          fprintf('Missing S97 data -- using aligned frac from runreads\n');
          stages{end+1}='Not PhiX'; labels{end+1}='PhiX';reads(end+1)=(1-mean([rr.runreads([1,3]).aligned]))*rr.readsPF;
        end
        
          
        stages{end+1}='Identified'; labels{end+1}='Unidentified';reads(end+1)=identified;
        if any(isphix)
          stages{end+1}='Not PhiX'; labels{end+1}='PhiX';reads(end+1)=notphix;
        end
        stages{end+1}='Paired'; labels{end+1}='Unpaired'; reads(end+1)=rr.readsPaired;
        stages{end+1}='Matched';labels{end+1}='Unmatched'; reads(end+1)=sum(data.cnt);
        stages{end+1}='Len>=48';labels{end+1}='Len<48';reads(end+1)=sum(data.cnt(data.length>=48));
        stages{end+1}='Cnt>=10';labels{end+1}='Cnt<10';reads(end+1)=sum(data.cnt(data.length>48&data.cnt>=10));
        labels{end+1}='Usable';
        reads=reads';
        relprior=nan(size(reads));
        relprior(2:end)=reads(2:end)./reads(1:end-1)*100;
        reltotal=reads./reads(1)*100;
        x=array2table([reads/1e6,relprior,reltotal],'VariableNames',{'Reads','RelPrior','RelTotal'},'RowNames',stages);
        x
        setfig('NGS Losses');clf;
        pie(-diff([reads;0]),labels);
        title(sprintf('%s: %.1fM Usable/%.1fM Raw Reads = %.0f%%',rr.descr,reads(end)/1e6,rr.readsTotal/1e6,reads(end)/reads(1)*100));
      end
    end
    
    function data=merge(obj,obj2,varargin)
    % Merge two different NGSVitro structures into one
      defaults=struct('offset',[],'prefix','','subruns',[]);
      args=processargs(defaults,varargin);

      % Make deep copies of the objects so that the merged one is distinct
      data=copy(obj);
      data2=copy(obj2);

      data.run=strjoin({data.run,data2.run},'+');
      % Append data2 information to data
      if isempty(args.prefix)
        args.prefix=[data2.run,'.'];
      end
      if isempty(args.offset)
        if isempty(data.subruns)
          args.offset=0;
        else
          args.offset=ceil(max(data.subruns.subrun)/1000)*1000;
        end
      end
      if isempty(args.subruns)
        ssel=true(size(data2.subruns.subrun));
        sel=true(size(data2.subrun));
      else
        sel=ismember(data2.subrun,args.subruns);
        ssel=ismember(data2.subruns.subrun,args.subruns);
      end
      fn=setdiff(fieldnames(data),'subrun');
      for i=1:length(fn)
        if size(data2.(fn{i}),1)==length(data2.subrun) && ~strcmp(fn{i},'subruns') && ~strcmp(fn{i},'trpexpts')
          data.(fn{i})=[data.(fn{i});data2.(fn{i})(sel)];
        end
      end
      data.subrun=[data.subrun;data2.subrun(sel)+args.offset];
      t1=length(data.trpexpts);
      data.trpexpts=[data.trpexpts,data2.trpexpts];
      for i=t1+1:length(data.trpexpts)
        data.trpexpts(i).descr=strcat(args.prefix,data.trpexpts(i).descr);
        data.trpexpts(i).cleavedsr=data.trpexpts(i).cleavedsr+args.offset;
        data.trpexpts(i).uncleavedsr=data.trpexpts(i).uncleavedsr+args.offset;
        if isfinite(data.trpexpts(i).insr)
          data.trpexpts(i).insr=data.trpexpts(i).insr+args.offset;
        end
        if ~isempty(data.trpexpts(i).analysis)
          data.trpexpts(i).analysis.subruns=data.trpexpts(i).analysis.subruns+args.offset;
        end
      end
      if isempty(data.subruns)
        % Handle empty starting point
        data.subruns=struct('subrun',[],'desc',{},'codes',{},'rds',[],'normalization',[]);
      end
      data.subruns(1).subrun=[data.subruns.subrun;data2.subruns.subrun(ssel)+args.offset];
      data.subruns.desc=[data.subruns.desc;strcat(args.prefix,data2.subruns.desc(ssel))];
      data.subruns.codes=[data.subruns.codes;data2.subruns.codes(ssel)];
      data.subruns.rds=[data.subruns.rds;data2.subruns.rds(ssel)];
      data.subruns.normalization=[data.subruns.normalization,data2.subruns.normalization(ssel)];
      if isfield(data.subruns,'bclen') && isfield(data2.subruns,'bclen')
        data.subruns.bclen=[data.subruns.bclen,data2.subruns.bclen(ssel)];
      end
      data.seqs.merge(data2.seqs);
    end
    
    function showstats(data,naseq,varargin)
    % Show stats for a particular naseq
      defaults=struct('details',0,'maxchanges',5);
      args=processargs(defaults,varargin);

      if ischar(naseq)
        name=naseq;
        naseq=data.seqs.findname(name);
        if isempty(naseq)
          fprintf('Unable to find naseq for %s\n', name);
          return;
        end
      end
      
      if data.seqs.getclusterindex(naseq)==0
        % Make sure its clustered
        data.clusterseqs(naseq);
      end
      data.seqs.showstats(naseq,'maxchanges',args.maxchanges);
      sel=find(data.naseq==naseq & data.cnt>0);
      if isempty(sel)
        return;
      end
      % Check for mutants
      data.checkifmutant(naseq,[],'silent',true,'minratio',6);
      
      udesc=unique(data.subruns.desc(ismember(data.subruns.subrun,data.subrun(sel))));
      [~,ord]=sort(cellfun(@(z) str2double(['0',z(ismember(z,'0123456789'))]), udesc,'UniformOutput',true));
      udesc=udesc(ord);
      prev=0;
      fprintf('Total reads: %d\n',sum(data.cnt(data.naseq==naseq)));
      for i=1:length(data.trpexpts)
        t=data.trpexpts(i);
        if isempty(t.analysis)
          continue;
        end
        tsel=find(t.analysis.naseq==naseq);
        if isempty(tsel)
          continue;
        end
        fprintf('TRP %2d: %-20.20s cleavage=%4.1f%%, ratio=%6.2f, ret=%5.1f%%, cnt: In=%4d Clv=%4d Unclv=%4d\n', i, t.descr, 100*t.analysis.cleavage(tsel), t.analysis.ratio(tsel), t.analysis.retention(tsel)*100, t.analysis.cnt(tsel,:));
      end

      if ~args.details
        return;
      end
      
      data.overallstats(data.naseq==naseq);
      return;
      
      for i=1:length(udesc)
        fprintf('%-20.20s ', udesc{i});
        ind=strcmp(data.subruns.desc,udesc{i});
        ucodes=unique(data.subruns.codes(ind));
        cnts=[];
        refcnts=[];
        for j=1:length(ucodes)
          ind2=strcmp(data.subruns.codes,ucodes{j})&ind;
          subrun=data.subruns.subrun(ind2);
          jj=intersect(sel,find(data.subrun==subrun));
          if isempty(jj)
            cnts(j)=0;
          else
            cnts(j)=data.cnt(jj);
          end
          refcnts(j)=nansum(data.subruns.normalization(ind2).cnt);
          refconc(j)=nansum(data.subruns.normalization(ind2).conc);
          fprintf('%5.5s %4d:%6d/%-6.0f %6.0f fM %5.0f ppm ', ucodes{j}, subrun, cnts(j), refcnts(j),cnts(j)/refcnts(j)*refconc(j)*1e15,cnts(j)/sum(data.cnt(data.subrun==subrun))*1e6);
        end
        fprintf('\n');
      end
    end
    
    function [cnt,info]=plotlengthpdf(data,sr,varargin)
      defaults=struct('mincnt',1,'desc',[],'newfig',true,'normalize',false,'alldata',false,'contexts',false);
      args=processargs(defaults,varargin);

      if args.alldata
        if ~isempty(args.desc)
          desc=args.desc;
        else
          desc='All data';
        end
        fprintf('Querying database...');
        NGSDatabase.open(data.run);
        if args.contexts
          [len,lcnt]=mysql('select length(e.seq)+length(s.leftcontext)+length(s.rightcontext) len,sum(cnt) total from v_entries e, v_subruns s where s.subrun=e.subrun group by len');
          desc=[desc,' (including contexts)'];
        else
          [len,lcnt]=mysql('select length(seq) len,sum(cnt) total from v_entries group by len');
        end
        cnt=zeros(max(len),1);
        cnt(len)=lcnt;
        fprintf('done\n');
        info=sprintf('%s N(rds)=%d, mean length=%.0fbp',desc,sum(cnt),dot(cnt,1:length(cnt))/sum(cnt));
      else
        sel=data.subrun==sr & data.cnt>=args.mincnt & ~ismember(data.naseq,data.refseqs(1))
        nreads=sum(data.cnt(data.subrun==sr));
        if ~isempty(args.desc)
          desc=args.desc;
        else
          desc=[data.subruns.desc{data.subruns.subrun==sr},'-',data.subruns.codes{data.subruns.subrun==sr}];
        end
        
        refcnt=[data.subruns.normalization(data.subruns.subrun==sr).cnt];
        info=sprintf('%s N(seq)=%d N(rds)=%d N(Ref)=%d',desc,sum(sel),nreads,sum(refcnt));
        for i=1:max(data.length(sel))
          cnt(i)=sum(data.cnt(sel&data.length==i));
        end
      end
      
      if args.newfig
        setfig(desc);clf;
      end
      fprintf('%s\n',info);
      if args.normalize
        cnt=cnt/sum(cnt);
      end
      bar(cnt);
      edges=[30+4-0.5,30+8+0.5,60+4-0.5,60+8+0.6]+41;
      hold on;
      c=axis;
      if ~args.contexts
        for i=1:length(edges)
          plot(edges(i)*[1,1],c(3:4),'r:');
        end
      end
      xlabel('Length (nt)');
      if ~args.normalize
        ylabel('Cnt');
      else
        ylabel('Frac');
      end

      title(info,'Interpreter','none');
    end
    
    function cnt=plotlengthcompare(data,sr1,sr2,varargin)
    % Compare length distribution over 2 different subruns
      defaults=struct('desc',[],'newfig',true,'mincnt',100);
      args=processargs(defaults,varargin);

      [c1,i1]=data.plotlengthpdf(sr1);
      [c2,i2]=data.plotlengthpdf(sr2);
      if length(c1)>length(c2)
        c2(length(c1))=0;
      elseif length(c2)>length(c1)
        c1(length(c2))=0;
      end
      f1=c1/sum(c1);
      f2=c2/sum(c2);
      
      setfig('Length Compare');clf;
      subplot(211);
      bar([f1;f2]');
      ylabel('PDF');
      xlabel('Length');
      legend(i1,i2);
      subplot(212);
      r=c2./c1;
      r(c1<=args.mincnt)=nan;
      semilogy(1:length(c1),r);
      xlabel('Length');
      ylabel('Ratio');
      yyaxis right;
      semilogy(1:length(c1),(c1+c2)/2);
      hold on;
      plot([1,length(c1)],args.mincnt*[1,1],':g');
      ylabel('Mean Cnt');
    end
    
    function newsr=mergesubruns(data,subruns)
    % Merge given subruns into a new subrun
      if length(unique(subruns))==1
        % Nothing to merge
        newsr= subruns(1);
        return;
      end
      newsr=max(data.subruns.subrun)+1;
      oldsel=[];
      desc='';
      for i=1:length(subruns)
        oldsel(i)=find(data.subruns.subrun==subruns(i));
        if i>1
          desc=[desc,'+'];
        end
        desc=[desc,data.subruns.desc{oldsel(i)}];
        if i==1
          n=data.subruns.normalization(oldsel(i));
        else
          if n.naseq~=data.subruns.normalization(oldsel(i)).naseq
            error('Cannot merge subruns: subrun %d uses reference %d whereas others use %d\n', subruns(i),data.subruns.normalization(oldsel(i)).naseq, n.naseq);
          end
          n.cnt=n.cnt+data.subruns.normalization(oldsel(i)).cnt;
        end
      end
      data.subruns.subrun(end+1)=newsr;
      data.subruns.desc{end+1}=desc;
      data.subruns.codes{end+1}=strjoin(unique(data.subruns.codes(oldsel)),',');
      data.subruns.normalization(end+1)=n;
      data.subruns.rds(end+1)=sum(data.subruns.rds(oldsel));
      x=data.commonseqs(subruns);
      n=length(x.naseq);
      data.subrun=[data.subrun;ones(n,1)*newsr];
      data.naseq=[data.naseq;x.naseq];
      data.cnt=[data.cnt;sum(x.cnt,2)];
    end
    
    function newtrp=mergetrp(data,ind,descr,varargin)
    % Merge all the underlying data of TRP experiments
      defaults=struct('target',[],'targetconc',[],'concunits','');
      args=processargs(defaults,varargin);

      if islogical(ind)
        ind=find(ind);
      end
      
      if any(strcmp(descr,{data.trpexpts.descr}))
        error('mergetrp:  %s already exists.',descr);
      end
      if all(isnan([data.trpexpts(ind).insr]))
        insr=nan;
      else
        insr=unique([data.trpexpts(ind).insr]);
      end
      if length(insr)>1
        fprintf('Warning - input subruns are different - using union\n');
        %insr=nan;   
      end
      if ~isempty(args.target)
        tgt=args.target;
      else
        tgt=unique({data.trpexpts(ind).target});
        if length(tgt)~=1
          fprintf('Warning merging runs with different targets: %s\n', strjoin(tgt,','));
          tgt=strjoin(tgt,',');
        else
          tgt=tgt{1};
        end
      end
      if isempty(tgt)
        tgtconc=[];
      elseif ~isempty(args.targetconc)
        tgtconc=args.targetconc;
      else
        tgtconc=unique([data.trpexpts(ind).targetconc]);
        tgtconc=tgtconc(isfinite(tgtconc));
        if isempty(tgtconc)
          tgtconc=nan;
        elseif length(tgtconc)~=1
          fprintf('Merging runs with different target concs: %s\n', sprintf('%f ',tgtconc));
          tgtconc=nan;
        end
      end
      if ~isempty(args.concunits)
        concunits=args.concunits;
      else
        concunits=unique({data.trpexpts(ind).concunits});
        if length(concunits)~=1
          fprintf('Merging runs with different conc units: %s\n', strjoin(concunits,','));
          concunits='';
        else
          concunits=concunits{1};
        end
      end

      newc=[data.trpexpts(ind).cleavedsr];
      newu=[data.trpexpts(ind).uncleavedsr];
      newtrp=max([data.trpexpts.trpexpt])+1;
      s=TRPExpt();
      s.run=data.trpexpts(ind(1)).run;
      s.trpexpt=newtrp;
      s.descr=descr;
      s.cleavedsr=newc;
      s.uncleavedsr=newu;
      s.insr=insr;
      s.bulkcleavage=nan;
      s.target=tgt;
      s.targetconc=tgtconc;
      s.concunits=concunits;
      s.duration=unique([data.trpexpts(ind).duration]);
      s.idsr=data.trpexpts(ind(1)).idsr;
      s.analysis=[];
      data.trpexpts(end+1)=s;

      fprintf('Completed merge of %s as trpexpt %d\n', descr, data.trpexpts(end).trpexpt);
    end
    
    function [f,c]=namedfraction(data,sel)
    % Return fraction of counts that are attributed to named sequences
      only=[data.seqs.seqnames.naseq];
      only=setdiff(only,union(data.ampliconseqs(0),data.refseqs(0)));
      c=data.cnt(sel & ismember(data.naseq,only));
      f=sum(c)/sum(data.cnt(sel));
    end
    
    function [f,nc,n,allreads]=fraction(data,subruns,ntop)
    % Return fraction of counts that are attributed to ntop top sequences (other than amplicons and refs)
      ntop=round(ntop);
      sel=ismember(data.subrun,subruns);
      fprintf('Computing fractions for ');
      for i=1:length(subruns)
        ind=data.subruns.subrun==subruns(i);
        fprintf('%d-%s-%s ',subruns(i), data.subruns.desc{ind}, data.subruns.codes{ind});
      end
      fprintf('\n');
      skip=union(data.refseqs(1),data.ampliconseqs(1));
      c=data.cnt(sel & ~ismember(data.naseq,skip));
      n=data.naseq(sel & ~ismember(data.naseq,skip));
      [n,ia,ic]=unique(n);
      nc=nan(size(n));
      for i=1:length(n)
        nc(i)=sum(c(ic==i));
      end
      [nc,ord]=sort(nc,'desc');
      n=n(ord);
      n=n(1:min(ntop,end));
      nc=nc(1:min(ntop,end));
      allreads=sum(data.subruns.rds(ismember(data.subruns.subrun,subruns)));
      f=min(nc)*length(nc)/allreads;
    end
    
    
    function [cnt,naseq]=plotcntcdf(data,sel,varargin)
    % Plot number of sequences with at least cnt reads vs cnt
      defaults=struct('newfig',true,'minreads',2,'showreads',[10,20,50]);
      args=processargs(defaults,varargin);
      if nargin<2 || isempty(sel)
        sel=data.cnt>=args.minreads;  % All data with >1 cnt
      end
      if args.newfig
        setfig('cntcdf');clf;
      end
      % Group cnt by naseq
      naseq=data.naseq(sel);
      selcnt=data.cnt(sel);
      fprintf('Total reads = %d\n', sum(selcnt));
      [naseq,ia,ic]=unique(naseq,'stable');
      cnt=zeros(size(naseq));
      for i=1:length(ic)
        cnt(ic(i))=cnt(ic(i))+selcnt(i);
      end
      naseq=naseq(cnt>=args.minreads);
      cnt=cnt(cnt>=args.minreads);
      [cnt,ord]=sort(cnt,'descend');
      naseq=naseq(ord);

      nseq=1:length(cnt);
      loglog(cnt,nseq);
      xlabel('Minimum Reads/Seq');
      ylabel('Num Seq with >= Minimum Reads');

      % Now compute how many reads were attributed to sequence with < reads/cnt
      yyaxis right
      ucnt=unique(cnt);
      avgrds=nan(size(ucnt));
      for i=1:length(ucnt)
        avgrds(i)=median(cnt(cnt>=ucnt(i)));
      end
      loglog(ucnt,avgrds);
      ylabel('Median Reads/Seq');
      
      % Label some points
      if ~isempty(args.showreads)
        hold on;
        for x=args.showreads
          y=sum(cnt>=x);
          yyaxis left
          plot(x*[1,1],[1,y],'r:');
          text(x,y,sprintf('%d',y),'VerticalAlignment','bottom');
          yyaxis right
          y2=min(avgrds(ucnt>=x));
          text(x,y2,sprintf('%.0f',y2),'VerticalAlignment','bottom','HorizontalAlignment','right');
        end
      end

      ti='Reads/Seq';
      % Check if selector is a set of subruns
      subruns=unique(data.subrun(sel));
      newsel=ismember(data.subrun,subruns);
      if all(newsel==sel) && length(subruns)<10
        srlist=arrayfun(@(z) sprintf('%d',z),subruns,'UniformOutput',false);
        ti=[ti,' (',strjoin(srlist,','),') ',strjoin(data.subruns.desc(ismember(data.subruns.subrun,subruns)),',')]
      end
      title(ti,'Interpreter','none');
    end
      
    function acc=plotlengthdist(data,sr,varargin)
    % Plot length distribution for a given subrun
      defaults=struct('mincnt',1,'desc',[],'newfig',true,'noutliers',20,'strsv',false,'maxsmallloop',9,'maxlongloop',65,'refconc',.002e-9,'includerefs',false);
      args=processargs(defaults,varargin);

      srsel = ismember(data.subrun,sr);
      sel=srsel & data.cnt>=args.mincnt;
      if ~args.includerefs
        sel=sel & ~ismember(data.naseq,data.refseqs(1));
      end
      nreads=sum(data.cnt(sel));
      if ~isempty(args.desc)
        desc=args.desc;
      else
        if length(sr)==1
          desc=[data.subruns.desc{data.subruns.subrun==sr},'-',data.subruns.codes{data.subruns.subrun==sr}];
        elseif all(ismember(data.subruns.subrun,sr))
          desc=sprintf('All Data');
        else
          desc=sprintf('%d subruns',length(sr));
        end
      end
      
      refcnt=sum([data.subruns.normalization(ismember(data.subruns.subrun,sr)).cnt]);
      info=sprintf('%s N(seq)=%d N(rds)=%d N(Ref)=%d',desc,sum(sel),nreads,refcnt);
      
      args.desc=info;
      if args.newfig
        setfig(desc);clf;
      end
      fprintf('%s\n',info);
      fsel=find(sel);
      %lib=struct('seq',{data.seq(fsel)},'cnt',data.cnt(fsel));
      l1=data.loop1len(fsel)+data.helix1len(fsel)*2;
      l1(isnan(l1))=0;
      l2=data.loop2len(fsel)+data.helix2len(fsel)*2;
      l2(isnan(l2))=0;
      args.maxlongloop=min([args.maxlongloop,max(l1),max(l2)]);
      l1(l1>args.maxlongloop)=args.maxlongloop;
      l2(l2>args.maxlongloop)=args.maxlongloop;
      cnt=data.cnt(fsel);
      d={};
      NsN60=[];
      N60Ns=[];
      NsN30=[];
      N30Ns=[];
      NsNs=[];
      N7N4=[];
      N30N30=[];
      other=[];
      ll=zeros(args.maxlongloop+1,args.maxlongloop+1);
      for i=1:length(fsel)
        v1=l1(i);
        v2=l2(i);
        ll(v2+1,v1+1)=ll(v2+1,v1+1)+cnt(i);
        d{i}=sprintf('N%dN%d',v1,v2);
        if v1==7+12 && v2==4+8
          N7N4(end+1)=cnt(i);
        elseif v1>=4+12 && v1<=args.maxsmallloop+12 && abs(v2-(60+8))<=1
          NsN60(end+1)=cnt(i);
        elseif v2>=4+8 && v2<=args.maxsmallloop+8 && abs(v1-(60+12))<=1
          N60Ns(end+1)=cnt(i);
        elseif v1>=4+12 && v1<=args.maxsmallloop+12 && abs(v2-(30+8))<=1
          NsN30(end+1)=cnt(i);
        elseif v2>=4+8 && v2<=args.maxsmallloop+8 && abs(v1-(30+12))<=1
          N30Ns(end+1)=cnt(i);
        elseif abs(v2-30-8)<=1 && abs(v1-30-12)<=1
          N30N30(end+1)=cnt(i);
        elseif v1<12+12 && v2<12+8
          NsNs(end+1)=cnt(i);
        else
          other(end+1)=cnt(i);
        end
      end

      ll(end+1,:)=nan;
      ll(:,end+1)=nan;
      pcolor(log10(ll));colorbar;
      axis equal
      axis([1,size(ll,1),1,size(ll,2)]);
      shading faceted;
      colormap(hot);
      xlabel('Stem 1 Length');
      ylabel('Stem 2 Length');
      title(info,'Interpreter','none');
      xt=unique([0,(4:2:8)+12,30+12,60+12],'sorted')
      set(gca,'XTick',xt+1);
      tl={};
      for i=1:length(xt)
        if xt(i)<12+3
          tl{i}=sprintf('%d',xt(i));
        else
          tl{i}=sprintf('%d+12',xt(i)-12);
        end
      end
      set(gca,'XTickLabel',tl);
      set(gca,'XTickLabelRotation',45);
      yt=unique([0,(4:2:8)+8,30+8,60+8],'sorted');
      set(gca,'YTick',yt+1.5);
      tl={};
      for i=1:length(yt)
        if yt(i)<8+3
          tl{i}=sprintf('%d',yt(i));
        else
          tl{i}=sprintf('%d+8',yt(i)-8);
        end
      end
      set(gca,'YTickLabel',tl);
      
      [ud,ia,ic]=unique(d);
      fprintf(' Seq      Cnt /Seq  Std Desc\n');
      for i=1:length(ud)
        nseq=sum(ic==i);
        nread=sum(cnt(ic==i));
        if nread>0.01*sum(cnt)
          fprintf('%5d %7d %4.0f %4.0f %s\n', nseq, nread, nread/nseq, std(cnt(ic==i)), ud{i});
        end
      end
      fprintf('%5d %7d %4.0f %4.0f N*N60\n',length(NsN60),sum(NsN60),mean(NsN60),std(NsN60));
      fprintf('%5d %7d %4.0f %4.0f N60N*\n',length(N60Ns),sum(N60Ns),mean(N60Ns),std(N60Ns));
      fprintf('%5d %7d %4.0f %4.0f N*N30\n',length(NsN30),sum(NsN30),mean(NsN30),std(NsN30));
      fprintf('%5d %7d %4.0f %4.0f N30N*\n',length(N30Ns),sum(N30Ns),mean(N30Ns),std(N30Ns));
      fprintf('%5d %7d %4.0f %4.0f N30N30\n',length(N30N30),sum(N30N30),mean(N30N30),std(N30N30));
      fprintf('%5d %7d %4.0f %4.0f N<12N<12\n',length(NsNs),sum(NsNs),mean(NsNs),std(NsNs));
      fprintf('%5d %7d %4.0f %4.0f Other\n',length(other),sum(other),mean(other),std(other));
      asdesigned=[sum(NsN60),sum(N60Ns),sum(NsN30),sum(N30Ns),max(N7N4) ];
      if args.strsv
        acc=sum(asdesigned)/sum(cnt)*100;
      else
        acc=sum(asdesigned(1:4))/sum(cnt)*100;
      end
      fprintf('%.1f%% as designed\n',acc);
      fprintf('\n');
      ax=axis;
      text(ax(2)-10,ax(4)-5,sprintf('%.0f%% Correct',acc),'Color','w');
      if (args.newfig)
        mkinsert(gcf,'aspect',1.05);
      end
      % Check large count outliers
      l1=data.loop1len(fsel);
      l2=data.loop2len(fsel);
      shortloop=min(l1,l2);
      longloop=max(l1,l2);
      outlier=shortloop<4 | shortloop>args.maxsmallloop | (abs(longloop-30)>1 & abs(longloop-60)>1);
      fprintf('Unexpected sequences with high counts:\n');
      data.listsequences(fsel(outlier),'maxlist',args.noutliers,'mincnt',10);
    end

    function labelled(data,sel,varargin)
    % List labelled sequences
      defaults=struct('addtags',{{}},'skiptags',{{}}, 'sort','cnt');
      args=processargs(defaults,varargin);

      naseqs=[data.seqs.seqnames.naseq];
      if ~isempty(args.addtags)
        fprintf('Including tags: %s',strjoin(args.addtags,','));
        naseqs=union(naseqs,[data.seqs.tags(ismember({data.seqs.tags.name},args.addtags)).naseq]);
        args.skiptags=setdiff(args.skiptags,args.addtags);  % Don't skip added ones
      end
      if ~isempty(args.skiptags)
        fprintf(' Skipping tags: %s\n', strjoin(args.skiptags,','));
        naseqs=setdiff(naseqs,[data.seqs.tags(ismember({data.seqs.tags.name},args.skiptags)).naseq]);
      end
      fprintf('\n');
      if nargin<2 || isempty(sel)
        sel=ismember(data.naseq,naseqs);
      else
        sel=sel&ismember(data.naseq,naseqs);
      end
      data.listsequences(sel,'maxlist',1000,'sort',args.sort);
    end
    
    function x=listcleavage(data,trpexpt,varargin)
      defaults=struct('mincnt',50,'bycluster',false);
      args=processargs(defaults,varargin);

      data.trpanalyze('mincnt',args.mincnt);
      a=data.trpexpts(trpexpt).analysis;
      cnt=sum(a.cnt(:,2:3),2);
      sel=find(cnt>=args.mincnt);
      extra=containers.Map('KeyType','double','ValueType','char');
      for ii=1:length(sel)
        i=sel(ii);
        extra(a.naseq(i))=sprintf('%5.1f',a.cleavage(i)*100);
      end
      fprintf('Cleavage for %s for sequences with >= %d reads\n',data.trpexpts(trpexpt).descr,args.mincnt);
      data.listsequences(ismember(data.naseq,a.naseq(sel))&ismember(data.subrun,a.subruns(2:3)),'sort','name','maxlist',length(sel),'extra',extra,'extraheader',' Clv ','bycluster',args.bycluster,'includerefs',false);
      x=table();
      x.naseq=int32(a.naseq(sel));
      x.cnt=int32(cnt(sel));
      x.cleavage=a.cleavage(sel);
      x.name=arrayfun(@(z) data.seqs.getname(z), x.naseq,'UniformOutput',false);
    end
    
    function listsequences(data,sel,varargin)
    % List sequences selected
    % Extra is other data as a map indexed by naseq
      defaults=struct('maxlist',30,'extra',containers.Map(),'extraheader','','sort','cnt','showcnt',true,'maxpergroup',10,'bycluster',false,'includerefs',true,'ignoresegs',[],'mincnt',10,'hitdetails',false,'maxclusterbuild',20);
      args=processargs(defaults,varargin);

      if nargin<2 || isempty(sel)
        sel=data.cnt>=args.mincnt;
        args.showcnt=false;
        fprintf('Listing all sequences with cnt>=%d in at least 1 subrun\n',args.mincnt);
      end
      if islogical(sel)
        sel=find(sel);
      end
      if isempty(sel)
        fprintf('No sequences selected\n');
        return;
      end
      if ~args.includerefs
        sel=sel(~ismember(data.naseq(sel),data.refseqs(1)));
      end
      sr=unique(data.subrun(sel));
      if length(sr)==1
        fprintf('Subrun %d, %s - %s\n', sr, data.subruns.desc{data.subruns.subrun==sr}, data.subruns.codes{data.subruns.subrun==sr});
        cnt=data.cnt(sel);
        naseq=data.naseq(sel);
      else
        % Only keep one of each naseq
        if length(sel) >20000
          thresh=prctile(data.cnt(sel),100*(1-20000/length(data.cnt(sel))));
          nsel=sel(data.cnt(sel)>thresh);
          if length(nsel)<length(sel) && length(nsel)>args.maxlist
            fprintf('Would take a while to compute total counts for %d entries; only considering the %d with >%.0f counts\n', length(sel),length(nsel),thresh);
            sel=nsel;
          end
        end
        c1=data.cnt(sel);
        [naseq,subsel,pick]=unique(data.naseq(sel),'stable');
        sel=sel(subsel);
        cnt=nan(size(naseq));
        for i=1:length(naseq)
          cnt(i)=sum(c1(pick==i));
        end
      end
      if ~isempty(args.ignoresegs)
        % Ignore sequence in given segments and regroup
        seqs={};
        for i=1:length(naseq)
          segs=data.seqs.getsegments(naseq(i));
          seqs{i}=strjoin(segs(setdiff(1:length(segs),args.ignoresegs)),'-');
        end
        [useqs,ia,ic]=unique(seqs);
        gnaseq=naseq(ia);
        for i=1:length(gnaseq)
          gcnt(i)=sum(cnt(ic==i));
        end
        assert(sum(cnt)==sum(gcnt));
        sel=sel(ia);
        naseq=gnaseq;
        cnt=gcnt;
      end
      
      if ~args.bycluster        
        if strcmp(class(args.sort),'containers.Map')
          keys=arrayfun(@(z) args.sort(z), naseq);
          [~,ord]=sort(keys);
        elseif strcmp(args.sort,'cnt')
          [~,ord]=sort(cnt,'descend');
        elseif strcmp(args.sort,'label') || strcmp(args.sort,'name')
          flags={};
          for i=1:length(naseq)
            flags{i}=data.seqs.getlabels(naseq(i),'hitdetails',args.hitdetails);
            if isempty(flags{i})
              flags{i}='~';  % Default to end of sort order
            end
          end
          [~,ord]=sort(flags);
        elseif strcmp(args.sort,'fixed')
          ord=1:length(cnt);
        elseif strcmp(args.sort,'seq')
          [~,ord]=sort(data.seqs.getseq(naseq));
        else
          error('Unknown sort key: %s (choices are cnt,label,name,fixed,seq)\n', args.sort);
        end

        if args.maxlist<length(ord)
          fprintf('Showing only first %d/%d sequences\n', args.maxlist, length(ord));
          ord=ord(1:args.maxlist);
        end
        
        sel=sel(ord);
        cnt=cnt(ord);
        naseq=naseq(ord);
        nprint=min(args.maxlist,length(sel));
      else
        cluster=data.getclusterindex(naseq);
        if any(cluster==0)
          csel=cluster==0 & cnt>=args.mincnt;
          rebuild=naseq(csel);
          [~,ord]=sort(cnt(csel),'desc');
          rebuild=rebuild(ord);
          if ~isempty(rebuild)
            fprintf('Building clusters for %d/%d sequences\n',min(args.maxclusterbuild,length(rebuild)),length(rebuild));
            data.clusterseqs(rebuild(1:min(args.maxclusterbuild,end)),'minaddcnt',1,'minreads',1,'mincnt',1);   % Add to existing ones, small ones will be pruned at save time
            cluster=data.getclusterindex(naseq);
          end
        end
          
        % Build sort order:  first by total cluster reads, then by sequence reads
        ucluster=unique(cluster);
        ccnt=[];
        for i=1:length(ucluster)
          ccnt(i)=sum(cnt(cluster==ucluster(i)));
        end
        [ccnt,ord]=sort(ccnt,'desc');
        ucluster=ucluster(ord);
        % Now keep only the first maxpergroup of each cluster, ordered by cnt desc
        reorder=[];
        lastcluster=-1;   % Last cluster to print to keep number of printed lines < maxlist
        nkept=0;   % Number of entries kept so far
        for i=1:length(ucluster)
          s=find(cluster==ucluster(i));
          c=cnt(s);
          [c,ord]=sort(c,'desc');
          s=s(ord);
          reorder=[reorder;s];
          nkept=nkept+min(length(s),args.maxpergroup);
          if nkept>=args.maxlist && lastcluster==-1
            lastcluster=ucluster(i);
          end
        end
        cluster=cluster(reorder);
        sel=sel(reorder);
        cnt=cnt(reorder);
        naseq=naseq(reorder);
        fprintf('Have %d unique clusters\n', length(unique(cluster)));
        if lastcluster>0
          nprint=find(cluster==lastcluster,1,'last');
        else
          nprint=length(sel);
        end
      end

      fseq=data.seqs.naseqformat(data.naseq(sel(1:nprint)),'ignoresegs',args.ignoresegs);
      if nprint==1
        fseq={fseq};
      end
      cat=data.categorize(sel(1:nprint));
      if args.showcnt
        fprintf('    Cnt');
      end
      fprintf('  AllCnt     NASeq     Loops   Type   Seq%s%s\n',blanks(length(fseq{1})-3),args.extraheader);
      hasflags=false;
      fcnt=0;
      for ii=1:nprint
        i=sel(ii);
            
        if args.bycluster
          if ii==1 || cluster(ii)~=cluster(ii-1)
            fcnt=0;
            if cluster(ii)==0
              if ii==1 || cluster(ii-1) ~= 0
                fprintf('No Cluster: %d reads over %d sequences\n',sum(cnt(cluster==0)),sum(cluster==0));
              end
            else
              fprintf('Cluster: %d %s: %d reads over %d sequences.\n',data.seqs.clusters(cluster(ii)).cluster,data.seqs.clusters(cluster(ii)).name,sum(cnt(cluster==cluster(ii))),sum(cluster==cluster(ii)));
            end
          end
          fcnt=fcnt+1;
          if fcnt>args.maxpergroup 
            if fcnt==args.maxpergroup+1
              fprintf('...\n');
            end
            if isempty(data.seqs.getname(data.naseq(i)))
              continue;  % Skip the rest of this cluster (except for named seqs)
            end
          end
        end

        flags=data.seqs.getlabels(data.naseq(i),'hitdetails',args.hitdetails);
        if args.showcnt
          fprintf('%7d ', cnt(ii));
        end
        fprintf('%7d %9d ', sum(data.cnt(data.naseq==data.naseq(i))), data.naseq(i));
        fprintf('%9.9s %6.6s %s ',cat{ii},data.seqs.seqclasses{data.getclass(i)}{1}, fseq{ii});
        if args.extra.isKey(data.naseq(i))
          fprintf('%s ',args.extra(data.naseq(i)));
        end
        fprintf('%s\n',flags);
        if ~isempty(flags)
          hasflags=true;
        end
      end
      if length(sel)>nprint
        fprintf('... and %d more sequences with a total of %d reads',length(sel)-nprint, sum(cnt(nprint+1:end)));
        if args.bycluster
          fprintf(' in %d clusters\n',length(unique(cluster(nprint+1:end))));
        else
          fprintf('\n');
        end
      end
      % if hasflags
      %   fprintf('Flags:  [M]utant, [C]lose to a common seq\n');
      % end
    end

    function [total,name]=mixfrac(data,mixsteps,varargin)
    % Show percentages of reads (or mass) following a nested cell array, mixsteps
    % For cell elements that have more than one column, 2nd column is relative number of reads planned (otherwise 1)
      defaults=struct('bymass',false,'planned',[],'ind',[],'names',{{}},'prodname',{{}});
      args=processargs(defaults,varargin);
      if nargin<2 || isempty(mixsteps)
        mixsteps=unique(data.subruns.desc);
      end
      t=[];
      names={};
      subruns=[];
      if ~iscell(mixsteps)
          mixsteps=num2cell(mixsteps);
      end
      if isempty(args.ind)
        args.ind=true(size(data.cnt));
      end
      if size(mixsteps,2) >= 2
        args.planned=cell2mat(mixsteps(:,2));
        args.planned=args.planned/sum(args.planned);
      end
      if size(mixsteps,2) >= 3
        args.names=mixsteps(:,3);
      end
      mixsteps=mixsteps(:,1);
      for i=1:length(mixsteps)
        if ischar(mixsteps{i})
          t(i)=sum(data.cnt(args.ind & ismember(data.subrun,data.subruns.subrun(strcmp(data.subruns.desc,mixsteps{i})))));
          subruns(i)=nan;
          names{i}=mixsteps{i};
        elseif length(mixsteps{i}) > 1  || iscell(mixsteps{i})
          if isempty(args.names)
            if length(mixsteps{i})==2
              prodname=sprintf('%s+%s',mixsteps{i}{1},mixsteps{i}{2});
            else
              prodname=sprintf('Prod%d',i);
            end
          else
            prodname=args.names{i};
          end
          [t(i),names{i}]=mixfrac(data,mixsteps{i},'bymass',args.bymass,'ind',args.ind,'prodname',prodname);
          subruns(i)=nan;
        else
          t(i)=sum(data.cnt(args.ind & data.subrun==mixsteps{i}));
          subruns(i)=mixsteps{i};
          names{i}=sprintf('%s-%s',data.subruns.desc{data.subruns.subrun==mixsteps{i}},data.subruns.codes{data.subruns.subrun==mixsteps{i}});
        end
      end
      if ~isempty(args.names)
        assert(length(names)==length(args.names));
        names=args.names;
      end
      if isempty(args.planned)
        args.planned=repmat(1/length(mixsteps),length(mixsteps),1);
      else
        args.planned=args.planned/sum(args.planned);
      end
      for i=1:length(mixsteps)
        if sum(t)>0 && args.planned(i)>0
          stars=[repmat('*',1,round(t(i)/sum(t)/args.planned(i)*10)),blanks(10)];
        else
          stars=blanks(10);
        end
        stars(10)='|';
        fprintf('%4d %-35.35s %7d  %5.1f%%  %5.1f%% %4.2f %s\n', subruns(i), names{i}, t(i), t(i)*100/sum(t),args.planned(i)*100,t(i)/sum(t)/args.planned(i),stars);
      end
      total=sum(t);
      if isempty(args.prodname)
        name=sprintf('Total(%d,%s...)',length(names),names{i});
      else
        name=args.prodname;
      end
      fprintf('     %-35.35s %7d  %5.1f%% %7.0f/part\n\n', name, total, 100,total/length(names))
    end
    
    function [c,sel]=bccnt(data,bclist)
    % Get total number of reads over given barcode list
      allbc=arrayfun(@(z) str2double(z.indexbarcode(2:end)),data.rundata.indexstats);
      sel=ismember(allbc,bclist);
      c=sum([data.rundata.indexstats(sel).readsIdentified]);
    end
    
    %
    % Delegates to NGSSeq
    %  Ones which start with get* expect an array of naseq as argument, others expects a logical selector or indices into naseq
    %
    function seq=getseq(data,naseq)
      seq=data.seqs.getseq(naseq);
    end
    
    function s=seq(data,sel)
      if nargin>=2
        s=data.getseq(data.naseq(sel));
      else
        s=data.getseq(data.naseq);
      end
    end
    
    function len=getlength(data,naseq)
      len=data.seqs.getlength(naseq);
    end

    % Get cluster indices
    function c=getclusterindex(data,naseq)
      c=data.seqs.getclusterindex(naseq);
    end

    % Get cluste
    function [descr,c]=getcluster(data,naseq)
      [descr,c]=data.seqs.getcluster(naseq);
    end

    function len=length(data,sel)
      if nargin>=2
        len=data.getlength(data.naseq(sel));
      else
        len=data.getlength(data.naseq);
      end
    end
    
    function c=clusterindex(data,sel)
      if nargin>=2
        c=data.getclusterindex(data.naseq(sel));
      else
        c=data.getclusterindex(data.naseq);
      end
    end
    
    function n=getloop1len(data,naseq)
      n=data.seqs.getloop1len(naseq);
    end
      
    function n=gethelix1len(data,naseq)
      n=data.seqs.gethelix1len(naseq);
    end
      
    function n=loop1len(data,sel)
      if nargin>=2
        n=data.getloop1len(data.naseq(sel));
      else
        n=data.getloop1len(data.naseq);
      end
    end

    function n=helix1len(data,sel)
      if nargin>=2
        n=data.gethelix1len(data.naseq(sel));
      else
        n=data.gethelix1len(data.naseq);
      end
    end

    function n=getloop2len(data,naseq)
      n=data.seqs.getloop2len(naseq);
    end
      
    function n=gethelix2len(data,naseq)
      n=data.seqs.gethelix2len(naseq);
    end
      
    function n=loop2len(data,sel)
      if nargin>=2
        n=data.getloop2len(data.naseq(sel));
      else
        n=data.getloop2len(data.naseq);
      end
    end

    function n=helix2len(data,sel)
      if nargin>=2
        n=data.gethelix2len(data.naseq(sel));
      else
        n=data.gethelix2len(data.naseq);
      end
    end

    function l1=loop1(data,sel)
      if nargin>=2
        l1=data.seqs.getloop1(data.naseq(sel));
      else
        l1=data.seqs.getloop1(data.naseq);
      end
    end

    function l2=loop2(data,sel)
      if nargin>=2
        l2=data.seqs.getloop2(data.naseq(sel));
      else
        l2=data.seqs.getloop2(data.naseq);
      end
    end
    
    function n=seglen(data,segnum,sel)
      n=data.seqs.getseglen(segnum,data.naseq(sel));
    end

    function res=seqformat(data,inds)
      res=data.seqs.naseqformat(data.naseq(inds));
    end

    function id=segmentid(data,segnum,sel)
      if nargin>=3
        id=data.seqs.getsegmentid(segnum,data.naseq(sel));
      else
        id=data.seqs.getsegmentid(segnum,data.naseq);
      end
    end    

    function c=getclass(data,sel)
      if nargin>=2
        c=data.seqs.getclass(data.naseq(sel));
      else
        c=data.seqs.getclass(data.naseq);
      end
    end

    function s=segments(data,sel)
      if nargin>=2
        s=data.seqs.getsegments(data.naseq(sel));
      else
        s=data.seqs.getsegments(data.naseq);
      end
    end
    
    %
    % End of delegates to NGSSeq
    %
    

    function s=categorize(data,sel)
    % Return short string categorization of selected sequences
      s1=data.helix1len(sel);
      l1=data.loop1len(sel);
      s2=data.helix2len(sel);
      l2=data.loop2len(sel);
      s=cell(length(l1),1);
      for i=1:length(l1)
        s{i}=sprintf('%d.%d/%d.%d',s1(i),l1(i),s2(i),l2(i));
      end
      s=strrep(s,'NaN','?');
    end
    
    function getmutclassdist(data,ind)
      if nargin<2
        ind=true(size(data.cnt));
      end
      
      totalc=sum(data.cnt(ind));
      mainclasses=find(cellfun(@(z) z{5}, data.seqs.seqclasses));
      totalp=sum(data.cnt(ind & ismember(data.getclass,mainclasses)));
      fprintf('Fraction of library as planned: %.1f%%\n', totalp/totalc*100);
      for i=1:length(data.seqs.seqclasses)
        c=sum(data.cnt(data.getclass==i & ind));
        if c==0
          continue;
        end
        nseq=length(unique(data.naseq(data.getclass==i & ind)));
        fprintf('%2d %-12.12s: %6d %7d  %5.1f%%', i, data.seqs.seqclasses{i}{1}, nseq, c, c*100/totalc);
        if ismember(i,mainclasses)
          fprintf(' %5.1f%%', c*100/totalp);
        end
        fprintf('\n');
      end
    end
    
    function barstack(data,vu,desc,classnames,varargin)
      defaults=struct('cm',[],'sort',[]);
      args=processargs(defaults,varargin);

      if isempty(args.cm)
        args.cm=colorcube(64);   % Colormap to use
        args.cm=args.cm(mod((0:7:7*size(args.cm,1)-1),size(args.cm,1))+1,:);
      end

      if ~isempty(args.sort)
        [~,ord]=sort(vu(end,args.sort),'descend');
        vu(:,args.sort)=vu(:,args.sort(ord));
        classnames(args.sort)=classnames(args.sort(ord));
      end

      b=bar(vu,'stacked');
      for k = 1:length(b)
        b(k).FaceColor = args.cm(k,:);
      end
      % Label sections
      hold on;
      yt=[]; y1=[];
      for i=1:size(vu,2)
        yt(i)=(i+0.5)*100/(size(vu,2)+1);
        y1(i)=sum(vu(end,1:i-1))+vu(end,i)/2;
      end
      % Adjust to optimal spacing
      yt=[-1,yt,102];
      y1=[-1,y1,102];
      minsep=min(2,100/size(vu,2)/3);
      for iter=1:50
        for i=1:length(yt)-1
          if yt(i)<y1(i) && yt(i+1)-yt(i)>minsep
            yt(i)=min(y1(i),yt(i+1)-minsep);
          end
        end
        for i=2:length(yt)
          if yt(i)>y1(i) && yt(i)-yt(i-1)>minsep
            yt(i)=max(y1(i),yt(i-1)+minsep);
          end
        end
      end
      y1=y1(2:end-1);
      yt=yt(2:end-1);
      textx=(size(vu,1)+0.4)*1.05;
      for i=1:size(vu,2)
        text(textx,yt(i),sprintf('%s (%.2f%%)',classnames{i},vu(end,i)),'Interpreter','none');
        plot([size(vu,1)+0.4,textx-.02],[y1(i),yt(i)],'-k');
      end
      
      set(gca,'XTick',1:length(desc));
      set(gca,'XTickLabel',desc);
      set(gca,'XTickLabelRotation',45);
      set(gca,'TickLabelInterpreter','None');
      ylabel('Percent of reads');
      c=axis;
      c(4)=101;
      c(1)=0.5;
      c(2)=c(2)*1.15;
      axis(c);
  end
  
    function [v1,v,desc,classes]=graphstruct(data,subruns,varargin)
    % Make a graph of the structural types vs subrun
      defaults=struct('cm',[],'thresh',0.5,'allparts',false,'sel',[]);  
      args=processargs(defaults,varargin);

      setfig('Structure of sequences'); clf;
      minlength=0;
      
      if nargin<2 || isempty(subruns)
        subruns=data.subruns.subrun;
      end
      
      v=[];
      l1len=data.loop1len();
      l2len=data.loop2len();
      desc={};
      class=data.getclass;
      for i=1:length(subruns)
        sel=data.subrun==subruns(i);
        if ~isempty(args.sel)
          sel=sel&args.sel;
        end
        t=sum(data.cnt(sel));
        for k=1:length(data.seqs.seqclasses)
          v(i,k)=sum(data.cnt(sel&class==k));
        end
        desc{i}=data.subruns.desc{data.subruns.subrun==subruns(i)};
      end
      % Collapse with same desc
      udesc=unique(desc);
      vc=[];
      for i=1:length(udesc)
        vc(i,:)=sum(v(strcmp(desc,udesc{i}),:),1);
      end
      v=vc;
      desc=udesc;
      
      % Total
      sel=ismember(data.subrun,subruns);
      if ~isempty(args.sel)
        sel=sel&args.sel;
      end
      totalindex=size(v,1)+1;
      for k=1:length(data.seqs.seqclasses)
        v(totalindex,k)=sum(data.cnt(sel&data.getclass==k));
      end
      desc{end+1}='Total';

      for i=1:size(v,1)
        v(i,:)=v(i,:)/sum(v(i,:));
      end
      v=v*100;
      used=any(v>args.thresh);
      classes=data.seqs.seqclasses(used);
      data.barstack(v(:,used),desc,cellfun(@(z) z{1},classes,'UniformOutput',false),'sort',3:sum(used)-1);
      
      v1=v(:,used);
      title('Structure of sequences');
      if ~args.allparts
        return;
      end
      
      pause(0.1);
      
      
      % Show breakdown of prefix, core, suffix
      v={};
      segs=[3,5,6,7,9];
      titles={'S1a','S1b','Core','S2a','S2b'};
      for k=1:length(segs)
        fprintf('Processing %s...',titles{k})
        id=data.segmentid(segs(k));
        sel=ismember(data.subrun,subruns);
        if ~isempty(args.sel)
          sel=sel&args.sel;
        end
        ind=data.seqs.naseq2ind(unique(data.naseq(sel)));
        naseq=data.seqs.naseq(ind);
        seq=data.seqs.segments(ind,segs(k));
        [useq,ia,ic]=unique(seq);
        fprintf('%d unique seqs...',length(useq));
        for i=1:length(subruns)
          sel=data.subrun==subruns(i);
          if ~isempty(args.sel)
            sel=sel&args.sel;
          end
          for j=1:length(useq)
            v{k}(i,j)=sum(data.cnt(sel&ismember(data.naseq,naseq(ic==j))));
          end
        end
        % Total
        v{k}(end+1,:)=sum(v{k},1);
        % Normalize and convert to percent
        for i=1:size(v{k},1)
          v{k}(i,:)=v{k}(i,:)/sum(v{k}(i,:))*100;
        end
        setfig(titles{k});clf;
        sel=any(v{k}>args.thresh,1);
        data.barstack(v{k}(:,sel),desc,useq(sel),'sort',1:sum(sel));
        title(sprintf('Structure of %s',titles{k}));
        fprintf('done\n');
        pause(0.1);
      end

    end

    function deletesubrun(data,sr)
      sel=~ismember(data.subruns.subrun,sr);
      data.subruns.subrun=data.subruns.subrun(sel);
      data.subruns.desc=data.subruns.desc(sel);
      data.subruns.codes=data.subruns.codes(sel);
      data.subruns.normalization=data.subruns.normalization(sel);
      sel=~ismember(data.subrun,sr);
      data.subrun=data.subrun(sel);
      data.naseq=data.naseq(sel);
      data.cnt=data.cnt(sel);
    end
    
    function s=commonseqs(data,subruns,mincnt)
    % Get counts of sequences in common for a given set of subruns
      if nargin<3
        mincnt=1;
      end
      subruns=subruns(:);
      naseqs=unique(data.naseq(ismember(data.subrun,subruns)));
      desc=cell(1,length(subruns));
      codes=desc;
      
      for j=1:length(subruns)
        if subruns(j)>0
          sel=data.subruns.subrun==subruns(j);
          if sum(sel)==0
            error('Subrun %d not found\n', subruns(j));
          end
          desc{j}=data.subruns.desc{sel};
          codes{j}=data.subruns.codes{sel};
        end
      end

      all=nan(length(naseqs),length(subruns));
      for j=1:length(subruns)
        nn=data.naseq(data.subrun==subruns(j));
        cc=data.cnt(data.subrun==subruns(j));
        [~,~,ib]=intersect(nn,naseqs,'stable');
        all(ib,j)=cc;
      end
      all(isnan(all))=0;
      all(:,isnan(subruns))=nan;

      [~,~,ib]=intersect(naseqs,data.naseq,'stable');
      minsel=nansum(all')>=mincnt;
      ib=ib(minsel);
      all=all(minsel,:);
      s=struct('subruns',subruns,'desc',{desc},'codes',{codes});
      fn=setdiff(fieldnames(data),{'run','subruns','seqs','trpexpts','subrun'});
      for i=1:length(fn)
        if ~isfield(s,fn{i}) && size(data.(fn{i}),1)==length(data.subrun)
          s.(fn{i})=data.(fn{i})(ib,:,:);
        end
      end
      s.loop1len=data.loop1len(ib);
      s.loop2len=data.loop2len(ib);
      s.cnt=all;
      s.selector=ib;
    end
    
    function addtrp(obj,descr,insr,cleavedsr,uncleavedsr,target,targetconc,concunits)
      if nargin<6
        target='';
      end
      t=struct('trpexpt',length(obj.trpexpts)+1,'descr',descr,'cleavedsr',cleavedsr,'uncleavedsr',uncleavedsr,'insr',insr,'bulkcleavage',nan,'target',target,'targetconc',targetconc,'concunits',concunits,'duration',nan,'idsr',[],'analysis',[]);
      obj.trpexpts=[obj.trpexpts,t];
      obj.trpanalyze();
    end
    
    function reflist=refseqs(obj,andclusters)
    % Return naseqs for references (include any in same cluster if second arg is present and true)
      if nargin<2
        andclusters=false;
      end
      reflist=union([obj.subruns.normalization.naseq],[obj.seqs.tags(strcmp({obj.seqs.tags.name},'ref')).naseq]);
      if andclusters
        reflist=union(reflist,abs([obj.seqs.clusters(ismember([obj.seqs.clusters.root],reflist)).members]));
      end
    end
    
    function ampliconlist=ampliconseqs(obj,andclusters)
    % Return naseqs for amplicons (include any in same cluster if second arg is present and true)
      if nargin<2
        andclusters=false;
      end
      ampliconlist=[obj.seqs.tags(strcmp({obj.seqs.tags.name},'amplicon')).naseq];
      if andclusters
        ampliconlist=union(ampliconlist,abs([obj.seqs.clusters(ismember([obj.seqs.clusters.root],ampliconlist)).members]));
      end
    end
    
    function hitlist=hitseqs(obj, varargin)
    % Return naseqs for hits (include any in same cluster if second arg is present and true)
      defaults=struct('minfold',1,'match','','clusters',false);
      args=processargs(defaults,varargin);

      sel=[obj.seqs.hits.fold]>=args.minfold;
      if ~isempty(args.match)
        sel=sel&cellfun(@(z) ~isempty(z), regexp({obj.seqs.hits.name},args.match));
      end
      hitlist=[obj.seqs.hits(sel).naseq];
      if args.clusters
        hitclusters=arrayfun(@(z) any(ismember([z.root,z.members],hitlist)),obj.seqs.clusters);
        hitlist=union(hitlist,[obj.seqs.clusters(hitclusters).members]);
      end
    end

    
    function trpanalyze(obj,varargin)
    % Analyze all trp experiments 
      defaults=struct('mincnt',10,'minlength',48,'nboot',100,'removeendmutants',true,'force',false);
      args=processargs(defaults,varargin);

      for j=1:length(obj.trpexpts)
        t=obj.trpexpts(j);
        if ~isprop(t,'analysis')
          fprintf('Adding missing analysis field to trpexpts(%d)',j);
          t.analysis=[];
        end
        if ~isempty(t.analysis) && t.analysis.args.mincnt<=args.mincnt && t.analysis.args.minlength==args.minlength && ~args.force
          continue;
        end
        ntotal=sum(obj.cnt(ismember(obj.subrun,[t.cleavedsr,t.uncleavedsr]))>=args.mincnt);
        if ntotal<5
          continue;
        end
        fprintf('Analyzing %s\n',t.descr);
        x=obj.commonseqs([t.insr,t.cleavedsr,t.uncleavedsr],args.mincnt);
        x.groups={1:length(t.insr),length(t.insr)+(1:length(t.cleavedsr)),length(t.insr)+length(t.cleavedsr)+(1:length(t.uncleavedsr))};
        x.srcnt=x.cnt;   % Count by subrun
        x.cnt=nan(size(x.cnt,1),3);
        x.refcnts=nan(1,3);
        for i=1:length(x.groups)
          x.cnt(:,i)=sum(x.srcnt(:,x.groups{i}),2);
          x.refcnts(i)=nansum([obj.subruns.normalization(ismember(obj.subruns.subrun,x.subruns(x.groups{i}))).cnt]);
        end
        if any(isnan(x.refcnts(2:end)))||any(x.refcnts(2:end)==0)
          fprintf('Refcnts zero: run without normalization\n');
          x.refcnts(:)=1;
        elseif isnan(x.refcnts(1))||x.refcnts(1)==0
          fprintf('Input refcnts zero: run without input normalization\n');
          x.refcnts(1)=1;
        end
        c=x.cnt(:,2);
        u=x.cnt(:,3);
        refc=x.refcnts(2);
        refu=x.refcnts(3);
        seqs=obj.getseq(x.naseq);
        x.endmutant=~strncmp(seqs,'GCTGTC',6) | ~strncmp(cellfun(@(z) z(end:-1:1),seqs,'UniformOutput',false),'CGACAAAGCAG',11);
        isref=ismember(x.naseq,[obj.subruns.normalization(ismember(obj.subruns.subrun,x.subruns)).naseq]);
        keep=(c+u)>=args.mincnt&~isref;  % commonseqs call included inputs as well
        if args.removeendmutants && ~isempty(keep&x.endmutant)
          keep=keep&~x.endmutant;
        end
        if ~any(keep)
          fprintf('No data to keep\n');
          continue;
        end
        fn=fieldnames(x);
        %fprintf('Keeping %d/%d\n', sum(keep), length(keep));
        for i=1:length(fn)
          if size(x.(fn{i}),1)==length(keep) && ~ismember(fn{i},{'subruns'})
            x.(fn{i})=x.(fn{i})(keep,:,:);
          end
        end
        c=c(keep);
        u=u(keep);
        
        lratio=nan(length(c),1);
        x.stdlratio=nan(length(c),1);
        x.ratioboot=nan(length(c),args.nboot);
        for i=1:length(c)
          [lratio(i),x.stdlratio(i),x.ratioboot(i,:)]=obj.bootcleave([c(i),u(i)],x.refcnts(2:3),'nboot',args.nboot);
        end
        x.ratio=exp(lratio);
        
        x.length=obj.getlength(x.naseq);
        x.ratio(x.length<args.minlength)=nan;   % Nor short products
        refratio=sum(refc)./sum(refu);
        if refratio < 0.09 || refratio > 11
          fprintf('TRP Experiment %s has RefC=%d, RefU=%d -- seems like missing data\n', t.descr, sum(x.refcnts(:,2:3),1));
          %x.ratio(:)=nan;
        end
        x.cleavage=1./(1+1./x.ratio);
        x.args=args;
        % Compute ratio of input to output
        x.retscale=x.refcnts(1)./[refc,refu];
        x.retention=((c+1)/(refc+1)+(u+1)/(refu+1))./((x.cnt(:,1)+1)/(x.refcnts(1)+1));   % Use an estimator that avoid zero divides
        x.retention(isnan(x.ratio))=nan;
        % Would be nice to have the absolute retention, but this would require knowledge of all the intervening steps
        x.retscale=x.retscale/median(x.retention(isfinite(x.retention)));
        x.retention=x.retention/median(x.retention(isfinite(x.retention)));
        % Compute fraction of pool
        % Group together frac by insr,cleavedsr,uncleavedsr
        x.frac=x.cnt;
        for i=1:length(x.groups)
          x.frac(:,i)=x.cnt(:,i)/sum(obj.cnt(ismember(obj.subrun,x.subruns(x.groups{i}))));   
        end
        x.fitness=x.frac(:,2).*x.frac(:,3)./(x.frac(:,1).^2);
        x.loop1len=obj.getloop1len(x.naseq);
        x.loop2len=obj.getloop2len(x.naseq);
        % Determine sources
        if ~isempty(t.idsr)
          src=obj.commonseqs(t.idsr);
          [c,ia,ib]=intersect(x.naseq,src.naseq);
          x.srccnt=nan(length(x.naseq),length(t.idsr));
          for k=1:length(ia)
            x.srccnt(ia(k),:)=src.cnt(ib(k),:);
          end
          [m,o]=max(x.srccnt');
          x.src=t.idsr(o);
          x.src(isnan(m))=nan;
        end
        obj.trpexpts(j).analysis=x;
      end
    end

    function [y1,y2]=zstattocnt(obj, z, x, sx, sy)
    % Find the values of y that correspond to zstat statistical difference from x using two-proportion z-test
      a=-sx^2;
      b=sx*sy*(z^2+2*x);
      c=sy*(sx*x*z^2-x.^2*sy);
      y1=(-b + sqrt(b.^2 - 4*a*c))/(2*a);
      y2=(-b - sqrt(b.^2 - 4*a*c))/(2*a);
      %      y=(-(sx*sy)*(z^2+2*x)+sqrt(sx^2*sy^2*(z^2+2*x).^2-4*sx^2*sy*(sx*x*z^2-x.^2*sy)))/(2*sx^2);
      %q1=(a*y1.^2 + b.*y1 + c)/a
      %q2=(a*y2.^2 + b.*y2 + c)/a
    end
    
    function retentionplot(obj, trpexpt,varargin)
      defaults=struct('maxlist',30,'mincnt',10,'zthresh',3,'errthresh',1.5,'plotcleavage',true,'minlength',70,'bycluster',false,'bistable',[],'trunc',[],'maxlength',110,'optimize',true);
      args=processargs(defaults,varargin);

    % Plot retention
      if nargin<2 || isempty(trpexpt)
        trpexpt=1:length(obj.trpexpts);
      end
      obj.trpanalyze('mincnt',args.mincnt);
      setfig(sprintf('Retention (%s)',sprintf('%d ',trpexpt)));clf;
      pnum=1;
      nvalid=0;
      for i=1:length(trpexpt)  
        t=obj.trpexpts(trpexpt(i)).analysis;
        if sum(isfinite(t.retention))>0
          nvalid=nvalid+1;
        end
      end
        
      for i=1:length(trpexpt)  
        t=obj.trpexpts(trpexpt(i)).analysis;
        if sum(isfinite(t.retention))==0
          continue;
        end
        if length(trpexpt)>1
          subplot(2,ceil(nvalid/2),pnum);pnum=pnum+1;
        elseif args.plotcleavage
          subplot(1,2,1);
        end
        keep=find(~ismember(t.naseq,obj.refseqs(1)) & ~ismember(t.naseq,obj.ampliconseqs(1)) & sum(t.cnt,2)>=args.mincnt & t.length>=args.minlength & t.length <= args.maxlength);
        if ~isempty(args.bistable)
          isbis = arrayfun(@(z) obj.seqs.isbistable(z), t.naseq(keep));
          if args.bistable
            keep = keep(isbis);
          else
            keep = keep(~isbis);
          end
        end
        if ~isempty(args.trunc)
          istrunc = obj.seqs.istruncation(t.naseq(keep));
          if args.trunc
            keep = keep(istrunc);
          else
            keep = keep(~istrunc);
          end
        end
        fprintf('%s - N=%d, ',obj.trpexpts(trpexpt(i)).descr,length(keep));
        fprintf('Retention = (%.4f*C + %.4f*U)/IN\n', t.retscale);
        naseq=t.naseq(keep);

        % Check for loss of long sequences on input
        maxoutlength=max(obj.seqs.getlength(naseq));
        maxinlength=max(obj.seqs.getlength(naseq(t.cnt(keep,1)>0)));
        fprintf('Max input length = %d, max output length = %d\n', maxinlength, maxoutlength);
        if maxinlength<maxoutlength
          error('Sequences greater than %d long seem to be missing from input subrun, whereas outputs have up to %d long - should set maxlength', maxinlength,maxoutlength);
        end
        x=max(0.5,t.cnt(keep,1));
        y=t.retention(keep).*x;
        if args.optimize
          y0=y;
          for refit=1:2
            fitsel=keep(y./x >0.3 & y./x <2);
            fit=robustfit(t.cnt(fitsel,2:3),max(0.5,t.cnt(fitsel,1)),[],[],'off');
            fprintf('Fitted retention model using %d points: retention = (%.4f*C + %.4f*U)/IN\n',  length(fitsel), fit);
            y=t.cnt(keep,2:3)*fit;
          end
        else
          % Scale such that output count is approx. equal to cleaved+uncleaved reads on average
          scaling=sum(sum(t.cnt(keep,2:3)))/nansum(y);
          y=y*scaling;
        end
        retention=y./x;
        % Add noise to fuzz points for better scatter plotting
        xn=x+0.2*(rand(size(x))-0.5);
        yn=y+0.2*(rand(size(y))-0.5);
        isnamed=ismember(naseq,[obj.seqs.seqnames.naseq]);
        loglog(xn(~isnamed),yn(~isnamed),'.');
        hold on;
        loglog(xn(isnamed),yn(isnamed),'g.');
        c=axis;
        c(1)=max(c(1),0.3);
        c(3)=max(c(3),0.9);
        axis(c);
        plot(c(1:2),c(1:2)*nansum(y)/nansum(x),'-g');
        % Compute z-statistic for testing whether these are consistent
        py=y/nansum(y); px=x/nansum(x);p=(x+y)/(nansum(x)+nansum(y));
        zstat=(px-py)./sqrt(p.*(1-p)*(1/nansum(x)+1/nansum(y)));
        zstat(isnan(zstat))=0;   % So Nan's will be zero and at end of sort
        if isnan(args.zthresh)
          % Use a zthresh for a p-value of 0.01 based on Bonferroni correction
          pval=0.01;
          zthresh=-norminv(.01/2/length(keep));
          fprintf('Using p-value=%.2g -> zthresh = %.2f (Bonferroni corrected)\n', pval, zthresh);
        else
          zthresh=args.zthresh;
        end
        % Plot contours
        xz=logspace(0,log10(c(2)));
        [yz1,yz2]=obj.zstattocnt(zthresh,xz,nansum(x),nansum(y));
        plot(xz,yz1,'r:');
        plot(xz,yz2,'r:');

        % Highlight ones with |zstat|>zthresh and (ratio>errthresh or ratio<1/errthresh)
        highlight=abs(zstat)>=zthresh & (retention>args.errthresh | retention<1/args.errthresh);
        loglog(xn(highlight),yn(highlight),'r.');
        extra=containers.Map('KeyType','double','ValueType','char');
        sortkey=containers.Map('KeyType','double','ValueType','double');
        for kk=1:length(x)
          k=keep(kk);
          extra(naseq(kk))=sprintf('%3.0f -> %3.0f+%3.0f ret=%4.2f, zstat=%5.2f',t.cnt(k,:),retention(kk),zstat(kk));
          sortkey(naseq(kk))=-abs(log(retention(kk)));
        end
        fprintf('Low retention:\n');
        obj.listsequences(ismember(obj.naseq,naseq(highlight & zstat>0))&ismember(obj.subrun,t.subruns),'extra',extra,'extraheader','retention','maxlist',args.maxlist,'sort',sortkey,'bycluster',args.bycluster);
        fprintf('High retention:\n');
        obj.listsequences(ismember(obj.naseq,naseq(highlight & zstat<0))&ismember(obj.subrun,t.subruns),'extra',extra,'extraheader','retention','maxlist',args.maxlist,'sort',sortkey,'bycluster',args.bycluster);
        fprintf('\n');
        xlabel('Input Reads');
        ylabel('Output Reads');
        hleg=legend(sprintf('N=%d',length(keep)),'Location','NorthWest');
        set(hleg,'FontSize',10);
        title(sprintf('Retention %s',obj.trpexpts(trpexpt(i)).descr),'Interpreter','none');
        if length(trpexpt)==1 && args.plotcleavage
          % Plot retention vs cleavage
          subplot(122);
          loglog(t.ratio(keep),retention,'.');
          cleaveticks(1,0,'pct',true);
          logticks(0,1);
          xlabel('Fraction Cleaved');
          ylabel('Retention');
        end
      end
    end

    function cleaveplot(obj,trpexpt, mincnt)
      if nargin<3
        mincnt=10;
      end
      if nargin<2 || isempty(trpexpt)
        setfig('Cleavage');clf;
        for i=1:length(obj.trpexpts)
          subplot(2,ceil(length(obj.trpexpts)/2),i);
          obj.cleaveplot(i,mincnt);
        end
        return;
      end
      
      obj.trpanalyze('mincnt',mincnt);
      x=obj.trpexpts(trpexpt).analysis;
      sel=sum(x.cnt(:,2:3),2)>=mincnt & isfinite(x.ratio);
      if sum(isfinite(x.ratio(sel)))>100
        pdfplot(x.ratio(sel),[],1);
      end
      cleaveticks(1,0);
      title(sprintf('Cleavage %s N=%d mincnt=%d',obj.trpexpts(trpexpt).descr,sum(isfinite(x.ratio(sel))),mincnt),'Interpreter','none');
    end

    function cleavecdf(obj,trpexpt, mincnt,idsr)
      if nargin<3
        mincnt=30;
      end
      obj.trpanalyze('mincnt',mincnt);
      x=obj.trpexpts(trpexpt).analysis;
      %sel=x.cnt(:,1)>=mincnt/2 & x.cnt(:,2)+x.cnt(:,3)>=mincnt;
      sel=x.cnt(:,2)+x.cnt(:,3)>=mincnt & ~ismember(x.naseq,obj.refseqs(1));
      ti=sprintf('CleaveCDF(%d)',trpexpt);
      if nargin>=4
        sel=sel&ismember(x.naseq,obj.naseq(obj.subrun==idsr));
        ti=sprintf('%s (only %s)',ti,obj.subruns.desc{obj.subruns.subrun==idsr});
      end
      setfig(ti);clf;
      cdfplot(x.ratio(sel));
      leg={sprintf('All N=%d',sum(sel))};
      hold on;
      sel2=sel&abs(x.loop1len-30)<=1&abs(x.loop2len-6)<=3;
      if sum(sel2)>=10
        leg{end+1}=sprintf('N30N* N=%d',sum(sel2));
        cdfplot(x.ratio(sel2));
      end
      sel2=sel&abs(x.loop2len-30)<=1&abs(x.loop1len-6)<=3;
      if sum(sel2)>=10
        leg{end+1}=sprintf('N*N30 N=%d',sum(sel2));
        cdfplot(x.ratio(sel2));
      end
      sel2=sel&abs(x.loop1len-60)<=1&abs(x.loop2len-6)<=3;
      if sum(sel2)>=10
        leg{end+1}=sprintf('N60N* N=%d',sum(sel2));
        cdfplot(x.ratio(sel2));
      end
      sel2=sel&abs(x.loop2len-60)<=1&abs(x.loop1len-6)<=3;
      if sum(sel2)>=10
        leg{end+1}=sprintf('N*N60 N=%d',sum(sel2));
        cdfplot(x.ratio(sel2));
      end
      seq=obj.seqs.getseq(x.naseq);
      bistable=cellfun(@(z) ~isempty(z), regexp(seq,'GAA[AG][CT][AG]G[CT].*GAA[AG][CT][AG]G[CT]'));
      ishit=ismember(x.naseq,[obj.seqs.hits.naseq]);
      leg{end+1}=sprintf('Bistable non-hit N=%d',sum(sel&~ishit&bistable));
      cdfplot(x.ratio(sel&bistable&~ishit));
      leg{end+1}=sprintf('Not Bistable non-hit N=%d',sum(sel&~ishit&~bistable));
      cdfplot(x.ratio(sel&~bistable&~ishit));
      % leg{end+1}=sprintf('Bistable hit N=%d',sum(sel&ishit&bistable));
      % cdfplot(x.ratio(sel&bistable&ishit));
      % leg{end+1}=sprintf('Not Bistable hit N=%d',sum(sel&ishit&~bistable));
      % cdfplot(x.ratio(sel&~bistable&ishit));
      
      legend(leg);
      set(gca,'XScale','log');
      cleaveticks(1,0);
      title(sprintf('%s %s mincnt=%d',ti,obj.trpexpts(trpexpt).descr,mincnt),'Interpreter','none');
    end

    function [meanlogratio,stdlogratio,ratio]=bootcleave(obj,cnts,refcnts,varargin)
    % Calculate mean log(ratio) and standard deviation(log(ratio)) 
      defaults=struct('nboot',1000,'nquant',1000);
      args=processargs(defaults,varargin);
      if nargin<3 || isempty(refcnts)
        refcnts=[1,1];
      end
      
      if any(refcnts==0)
        error('Refcnts cannot be zero');
      end
      c=cnts(1); u=cnts(2);
      
      rng=[1e-8,1-1e-8];  % Avoid extrema
      for i=1:3
        %fprintf('rng=[%g,%g]\n',rng);
        p=linspace(rng(1),rng(2),args.nquant);
        loglike=log(p).*c + log(1-p).*u;
        loglike=loglike-max(loglike(isfinite(loglike)));   % Shift to avoid bounds problems
        like=exp(loglike);
        like=like/sum(like);   % Convert to prob
        cumlike=cumsum(like);
        sel=[find(cumlike>0.1/args.nboot,1),find(1.0-cumlike<0.1/args.nboot,1)];
        % Redo over new interval
        rng=[p(max(1,sel(1)-1)),p(min(sel(end)+1,end))];
      end
      sel=find(like>0);
      p=p(sel);
      like=like(sel);

      cumlike=cumsum(like);
      r=p./(1-p)*refcnts(2)/refcnts(1);
      ratio=interp1(cumlike,r,linspace(cumlike(2),cumlike(end-1),args.nboot),'pchip');
      lratio=log(ratio);
      meanlogratio=mean(lratio);
      stdlogratio=std(lratio);
    end
    
    function [meanlogfold,stdlogfold,fold]=bootfold(obj,r1,r2,varargin)
    % Calculate mean log(fold change) and standard deviation(log(fold)) going from (c1,u1) to (c2,u2)  (i.e if cleavage1 > cleavage1, fold is >1)
      defaults=struct('nboot',[],'metric','ratio_of_ratio');
      args=processargs(defaults,varargin);
      % Sample each individual 
      if isempty(args.nboot)
        args.nboot=max([100,2*length(r1),2*length(r2)]);
      end
      br1=r1(randi(length(r1),args.nboot,1));
      br2=r2(randi(length(r2),args.nboot,1));
      bc1=br1./(1+br1);
      bc2=br2./(1+br2);
      bu1=1-bc1;
      bu2=1-bc2;
      if strcmp(args.metric,'ratio_of_ratio')
        fold=sqrt(br2./br1);
      elseif strcmp(args.metric,'ratio_of_cleavage')
        fold=(bc2)./(bc2+bu2)./((bc1)./(bc1+bu1));
      elseif strcmp(args.metric,'ratio_of_uncleavage')
        fold=(bu1./(bc1+bu1))./(bu2./(bc2+bu2));
      elseif strcmp(args.metric,'fitness')
        fold=2*sqrt(bu1./(bu1+bc1).*(bc2./(bc2+bu2))); 
      else
        error('Unknown fold-change metric: %s (valid choices: ratio_of_ratio, ratio_of_cleavage, ratio_of_uncleavage\n', args.metric);
      end
      lfold=log(fold);
      meanlogfold=mean(lfold);
      stdlogfold=std(lfold);
    end
    
    function expectedCI(obj,c1,c2,varargin)
    % Get statistics for fold change going from c1 to c2 as a function of number of reads
      defaults=struct('nboot',1000,'readratio',0.5,'metric','ratio_of_uncleavage','ci',[5,95],'tgtfold',3);
      args=processargs(defaults,varargin);
      alln=unique(round(logspace(log10(10),log10(500),20)));
      ci=[];allmean=[];phigh=[];phighexact=[];
      for i=1:length(alln)
        n=alln(i);
        nc1=binornd(n,c1,args.nboot,1);
        nc2=binornd(n,c2,args.nboot,1);
        r1=[];r2=[];
        for j=1:length(nc1)
          logr=obj.bootcleave([nc1(j),n-nc1(j)],[1,1]);
          r1(j)=exp(logr);
          logr=obj.bootcleave([nc2(j),n-nc2(j)],[1,1]);
          r2(j)=exp(logr);
        end
        [m,s,f]=obj.bootfold(r1,r2,'metric',args.metric);
        ci(i,:)=prctile(f,args.ci);
        phigh(i)=mean(f>args.tgtfold);
        allmean(i)=exp(m);

        % Exact method
        nc=0:n;
        p1=binopdf(nc,n,c1);
        p2=binopdf(nc,n,c2);
        r=[];
        for j=1:length(nc)
          r(j)=obj.bootcleave([nc(j),n-nc(j)],[1,1]);
        end
        r=exp(r);
        u=1./(1+r);
        p=0;
        for j=1:length(nc)
          p=p+p1(j)*sum(p2(u(j)>args.tgtfold*u));
        end
        phighexact(i)=p;
      end
      setfig('expectedCI');clf;
      subplot(211);
      loglog(alln,ci,'k:');
      hold on;
      plot(alln,allmean);
      logticks(0,1);
      xlabel('N');
      ylabel(args.metric,'Interpreter','none');
      subplot(212)
      loglog(alln,phigh);
      hold on;
      loglog(alln,phighexact);
      logticks(1,1);
      sel=find(phighexact<1e-6,1);
      axis([alln([1,sel]),phighexact(sel),1]);
      legend('Bootstrap','Analytic ratio of uncleavage');
      xlabel('N');
      ylabel(sprintf('P(fold>%.1f)',args.tgtfold));
      suptitle(sprintf('c1=%.2f, c2=%.02f, nboot=%d',c1,c2,args.nboot));
    end
    
    function x=switchscatter(obj,expts, varargin)
      defaults=struct('mincnt',30,'nscale',1,'minlength',48,'title',[],'minfold',3.0,'naseqs',[],'rescale',1,'metric','ratio_of_ratio','showoff',false,'highlight',[],'showerr',1,'showcleave',1,'alphalabel',true,'maxstd',2,'publication',false,'maxlist',100,'bycluster',true,'dbdescr','','hitdetails',false);
      args=processargs(defaults,varargin);
      if strcmp(args.metric,'fitness') && args.rescale==1
        error('Rescaling of fitness metric by fitness not supported');
      end
      if ~ismember(args.rescale,[0,1,2,3])
        error('Unsupported rescale option: %d;  expected 0-no rescale, 1-rescale by metric, 2-rescale by ratio, 3-rescale by ratio of NSRef seqs',args.rescale);
      end
    % Compare first experiment in list with each of the others
      if nargin<2 || isempty(expts)
        expts=1:length(obj.trpexpts)
      end
      if isempty(args.title)
        args.title=sprintf('switchscatter(%s)',sprintf('%d ',expts));
      end
      obj.trpanalyze('mincnt',args.mincnt,'minlength',args.minlength);
      pcntr=1;
      if ~args.showcleave && ~args.showerr && nargout==0
        error('Nothing to show -- both showcleave and showerr are false and no output\n');
      elseif args.showcleave && args.showerr
        nc=2;
        nr=length(expts)-1;
      else
        nc=max(1,floor(sqrt(3/4*(length(expts)-1))));
        nr=ceil((length(expts)-1)/nc);
      end
      if args.showcleave || args.showerr
        if ~args.showcleave
          setfig([args.title,'-err']);cla;
        else
          setfig(args.title);cla;
        end
      end

      for i=1:1 % length(expts)
        t1=obj.trpexpts(expts(i)).analysis;
        for j=i+1:length(expts)
          t2=obj.trpexpts(expts(j)).analysis;
          sel1=t1.cnt(:,2)+t1.cnt(:,3)>=args.mincnt & isfinite(t1.ratio);
          sel2=t2.cnt(:,2)+t2.cnt(:,3)>=args.mincnt & isfinite(t2.ratio);
          naseqs=intersect(t1.naseq(sel1),t2.naseq(sel2));
          naseqs=setdiff(naseqs,obj.refseqs(1));
          if ~isempty(args.naseqs)
            onlyseqs=ismember(naseqs,args.naseqs);
            naseqs=naseqs(onlyseqs);
          end
          [~,~,ia]=intersect(naseqs,t1.naseq,'stable');
          [~,~,ib]=intersect(naseqs,t2.naseq,'stable');
          n1=sum(t1.cnt(ia,2:3),2);
          n2=sum(t2.cnt(ib,2:3),2);
          c1=t1.cleavage(ia);
          c2=t2.cleavage(ib);
          if strcmp(args.metric,'ratio_of_ratio') || strcmp(args.metric,'fitness')
            r1=t1.ratio(ia);
            r2=t2.ratio(ib);
          elseif strcmp(args.metric,'ratio_of_cleavage')
            r1=t1.cleavage(ia);
            r2=t2.cleavage(ib);
          elseif strcmp(args.metric,'ratio_of_uncleavage')
            r1=-(1-t1.cleavage(ia));
            r2=-(1-t2.cleavage(ib));
          else
            error('Unsupported metric: %s',args.metric);
          end

          meanlogfold=[]; stdlogfold=[];
          if args.rescale==2 || args.rescale==3
            rsel=ones(size(ia));
            if args.rescale==3
              rescaleseqs=[obj.seqs.seqnames(strncmp({obj.seqs.seqnames.name},'NSRef',5)).naseq];
              rsel=ismember(naseqs,rescaleseqs);
              if length(rsel)>=4
                fprintf('Rescaling using %d NSRef sequences\n',sum(rsel));
              else
                fprintf('Only %d NSRef sequences; rescaling using overall median instead\n', sum(rsel));
              end
            end
            medlogfold=nanmedian(log(t1.ratio(ia(rsel))./t2.ratio(ib(rsel))));
          else
            medlogfold=0;
          end
          for k=1:length(ia)
            [meanlogfold(k),stdlogfold(k)]=obj.bootfold(t2.ratioboot(ib(k),:)*exp(medlogfold),t1.ratioboot(ia(k),:),'metric',args.metric);
          end
          if args.rescale==1
            medlogfold=nanmedian(meanlogfold);
            meanlogfold=meanlogfold-medlogfold;
          end
          % Find extreme points
          npts=length(meanlogfold);
          alpha=1/(npts*args.nscale);
          z=-norminv(alpha);
          % Find ones where logfold-medlogfold is > z*stdlogfold
          extremeon=find(meanlogfold > z * stdlogfold & exp(meanlogfold)>args.minfold);
          if strcmp(args.metric,'fitness')
            extremeoff=[];	  % Would show too many low-fitness sequences
          else
            extremeoff=find(meanlogfold < -z * stdlogfold & exp(meanlogfold)<1/args.minfold);
          end
          extreme=union(extremeon,extremeoff);
          islabelled=ismember(naseqs(extreme),[obj.seqs.seqnames.naseq]);
          % How many would've been detectable at minfold?
          detectable=sum(log(args.minfold)>z*stdlogfold);
          %fprintf('Of %d total sequences, %d had enough data to detect a %.2f fold change\n', length(stdlogfold), sum(detectable), args.minfold);

          if ~isempty(args.dbdescr)
            % Create hit entry in database
            assert(strcmp(obj.trpexpts(expts(i)).run,obj.trpexpts(expts(j)).run));   % Both from same run (in case of merge)
            NGSDatabase.open();
            [hitrun,descr]=mysql(sprintf('select hitrun,descr from ngs.hitruns where run=''%s'' and minus=''%s'' and plus=''%s''',obj.trpexpts(expts(i)).run,obj.trpexpts(expts(i)).descr,obj.trpexpts(expts(j)).descr));
            if ~isempty(hitrun)
              descr=descr{1};
              fprintf('Inserting into existing hitrun %d (%s)\n',hitrun,descr);
              if ~strcmp(descr,args.dbdescr)
                fprintf('Changing name of hitrun %d from %s to %s\n',hitrun, descr, args.dbdescr);
                mysql(sprintf('update ngs.hitruns set descr=''%s'' where hitrun=%d',args.dbdescr,hitrun));
              end
              nold=mysql(sprintf('select count(*) from ngs.hitseqs where hitrun=%d',hitrun));
              if nold>0
                fprintf('Deleting %d old hitseqs\n', nold);
                mysql(sprintf('DELETE from ngs.hitseqs WHERE hitrun=%d',hitrun));
              end
            else
              cmd=sprintf('insert into ngs.hitruns(run,minus,plus,descr) values (''%s'',''%s'',''%s'',''%s'')',obj.trpexpts(expts(i)).run,obj.trpexpts(expts(i)).descr,obj.trpexpts(expts(j)).descr,args.dbdescr)
              nins=mysql(cmd);
              assert(nins==1);
              hitrun=mysql('select LAST_INSERT_ID()');
              fprintf('Inserting into new hitrun %d\n',hitrun);
            end
            if ~isempty(extremeon)
              cmd='INSERT INTO ngs.hitseqs(hitrun,naseq,fold,foldstd,clvminus,clvplus) VALUES ';
              for eii=1:length(extremeon)
                ei=extremeon(eii);
                cmd=sprintf('%s\n(%d,%d,%.3f,%.3f,%.3f,%.3f),',cmd,hitrun,naseqs(ei),exp(meanlogfold(ei)),exp(stdlogfold(ei)),c1(ei),c2(ei));
              end
              cmd(end)=';';
              mysql(cmd);
            end
            obj.seqs.dbloadhits();
          end
          
          % Highlights
          hsel=ismember(naseqs,args.highlight);
          nval=sprintf('N=%d',length(naseqs));
          ti=sprintf('%s vs %s N=%d, N(det)=%d',obj.trpexpts(expts(i)).descr,obj.trpexpts(expts(j)).descr,length(naseqs),detectable);
          
          %%%%%%%%%%%%%%%%%
          % Plot stderr vs fold
          %%%%%%%%%%%%%%%%%
          if args.showerr
          subplot(nr,nc,nc*(j-2)+1);
          loglog(exp(meanlogfold),exp(stdlogfold),'.');
          hold on;
          % Add highlights
          if ~isempty(hsel)
            loglog(exp(meanlogfold(hsel)),exp(stdlogfold(hsel)),'sm','MarkerSize',10,'LineWidth',2);
          end
          % Add scaling references (NSRef)
          if exist('rescaleseqs','var') && ~isempty(rescaleseqs)
            rssel=ismember(naseqs,rescaleseqs);
            loglog(exp(meanlogfold(rssel)),exp(stdlogfold(rssel)),'+m');
          end
          
          if args.publication
            xlabel('Fold Change');
            ylabel('Error Estimate');
          else
            xlabel(args.metric,'Interpreter','none');
            ylabel('std (multiplicative)');
          end
          % Add error markers
          ax=axis;
          ax(4)=min(ax(4),args.maxstd);
          axis(ax);
          plot(0*[1,1],ax(3:4),':k');
          plot(args.minfold*[1,1],ax(3:4),'r--');
          if ~strcmp(args.metric,'fitness') && args.showoff
            plot(1/args.minfold*[1,1],ax(3:4),'r--');
          end
          % Plot detectability threshold
          % plot(ax(1:2),exp(log(args.minfold)/z)*[1,1],'r--');
          
          ax=axis;
          ax(2)=ax(2)*1.05;
          fvals=linspace(log(ax(1)/2),log(ax(2)*2));
          if strcmp(args.metric,'fitness')
            fvals=fvals(fvals>=0);
          end
          for alpha=1./(npts.*[args.nscale]) % ,0.1,1,10])
            z=-norminv(alpha);
            if alpha==1./(npts*args.nscale)
              col=':r';
            else
              col=':k';
            end
            plot(exp(fvals),exp(abs(fvals)./z),col);
            if alpha*npts>=1
              astr=sprintf('alpha=%.0f/N',alpha*npts);
            else
              astr=sprintf('alpha=1/%.0fN',1/(alpha*npts));
            end
            if args.alphalabel
              text(ax(2)*1.01,exp(abs(log(ax(2)))./z),astr);
            end
          end
          axis(ax);

          hold on;
          for k=1:length(extreme)
            if islabelled(k)
              %col='g';
              col='r';   % Do this for simplifying the figures in the paper
            else
              col='r';
            end
            plot(exp(meanlogfold(extreme(k))),exp(stdlogfold(extreme(k))),['o',col]);
          end
          logticks(1,0);
          if ~args.publication
            title(ti,'Interpreter','none');
          end
          end
          
          %%%%%%%%%%%%%%%%%
          % Plot cleavages
          %%%%%%%%%%%%%%%%%
          if args.showcleave
          if args.showerr
            subplot(nr,nc,nc*(j-2)+nc);
          elseif nr*nc>1
            subplot(nr,nc,j-1);
          end
          loglog(r1,r2,'.');
          hold on;
          % Plot highlights
          if ~isempty(hsel)
            loglog(r1(hsel),r2(hsel),'ms','MarkerSize',10,'LineWidth',2);
          end
          % Plot rescale refs (NSRefs)
          if exist('rescaleseqs','var') && ~isempty(rescaleseqs)
            rssel=ismember(naseqs,rescaleseqs);
            loglog(r1(rssel),r2(rssel),'+m');
          end
          for k=1:length(extreme)
            if islabelled(k)
              %col='g';
              col='r';   % Do this for simplifying the figures in the paper
            else
              col='r';
            end
            loglog(r1(extreme(k)),r2(extreme(k)),['o',col]); %,'MarkerSize',20*exp(stdlogfold(extreme(k))));
          end
          ax=axis;
          % Plot minratio curve
          neg=0.001:0.001:0.999;
          if strcmp(args.metric,'fitness')
            for fit=0.5:0.1:1.0
              pos=1-(fit.^2./neg);
              pos(pos>1 | pos<0)=nan;
              % Fitness is plotted on ratio_of_ratio axes
              plot(neg./(1-neg),pos./(1-pos)/exp(medlogfold),'r--');
            end
          else
            for fold=[1/args.minfold,1,args.minfold]
              if fold<1 && ~args.showoff
                continue;
              end
              if strcmp(args.metric,'ratio_of_ratio')
                negratio=neg./(1-neg);
                posratio=negratio/(fold^2);
                pos=posratio./(1+posratio);
                %pos=sqrt(fold./neg.*(1-(1-neg)./fold));
              elseif strcmp(args.metric,'ratio_of_cleavage')
                pos=neg./fold;
              elseif strcmp(args.metric,'ratio_of_uncleavage')
                pos=1-fold*(1-neg);
              else
                error('Unknown fold-change metric: %s (valid choices: ratio_of_ratio, ratio_of_cleavage, ratio_of_uncleavage\n', args.metric);
              end
              if args.rescale
                if strcmp(args.metric,'ratio_of_ratio') || args.rescale==2
                  posratio=pos./(1-pos);
                  posratio=posratio/exp(medlogfold);
                  pos=posratio./(1+posratio);
                elseif strcmp(args.metric,'ratio_of_cleavage')
                  pos=neg./exp(medlogfold);
                elseif strcmp(args.metric,'ratio_of_uncleavage')
                  pos=1-exp(medlogfold)*(1-neg);
                end
              end
              pos(pos<=0|pos>=1)=nan;   % Kill impossible situations
              if fold==1
                sym=':k';
              else
                sym='r--';
              end
              if strcmp(args.metric,'ratio_of_ratio')
                plot(neg./(1-neg),pos./(1-pos),sym);
              elseif strcmp(args.metric,'ratio_of_cleavage')
                plot(neg,pos,sym);
              elseif strcmp(args.metric,'ratio_of_uncleavage')
                plot(neg-1,pos-1,sym);
              end
            end
          end

          cleaveticks(1,1,'metric',args.metric);
          xlabel(sprintf('%s',obj.trpexpts(expts(i)).descr),'Interpreter','none');
          ylabel(sprintf('%s',obj.trpexpts(expts(j)).descr),'Interpreter','none');
          if args.rescale==1
            ti=sprintf('%s median fold=%.2fx',ti,exp(medlogfold));
          elseif args.rescale==2
            ti=sprintf('%s median ratio-fold=%.2fx',ti,exp(medlogfold));
          end
          if ~args.publication
            title(ti,'Interpreter','none');
          end
          end
          
          %%%%%%%%%%%%%%%%%
          % Summarize data
          %%%%%%%%%%%%%%%%%
          fprintf('\n%s vs %s: %s',obj.trpexpts(expts([i,j])).descr,args.metric);
          if args.rescale==2
            fprintf(', rescaling around median ratio-fold change of %.2f', exp(medlogfold));
          elseif args.rescale
            fprintf(', rescaling around median fold change of %.2f', exp(medlogfold));
          end
          fprintf(' N=%d',npts);
          fprintf('\n');

          % Results to return
          x=table();
          x.naseq=naseqs;
          x.logfold=meanlogfold';
          x.stdlogfold=stdlogfold';
          x.cleavage=[c1,c2];
          x.extreme=false(size(naseqs));
          x.extreme(extreme)=true;
          
          extra=containers.Map('KeyType','double','ValueType','char');
          ratios=containers.Map('KeyType','double','ValueType','double');
          natmp=union(extreme,find(ismember(naseqs,args.highlight)));

          for k=1:length(natmp)
            kk=natmp(k);
            %seq=t1.seq{t1.naseq==naseqs(kk)};
            extra(naseqs(kk))=sprintf('%5d %5d %3.0f%%%3.0f%% %5.2f %6.2f',n1(kk),n2(kk),c1(kk)*100, 100*c2(kk),exp(meanlogfold(kk)),exp(stdlogfold(kk)));
            ratios(naseqs(kk))=exp(meanlogfold(kk));
          end
          selon=find(ismember(obj.naseq,naseqs(extremeon)) & ismember(obj.subrun,t1.subruns));
          seloff=find(ismember(obj.naseq,naseqs(extremeoff)) & ismember(obj.subrun,t1.subruns));
          if strcmp(args.metric,'fitness')
            extraheader='  n1    n2   c1  c2  fit    std* Flags';
          else
            extraheader='  n1    n2   c1  c2  fold   std* Flags';
          end
          if sum(selon)>0
            [~,ord]=sort(arrayfun(@(z) ratios(z),obj.naseq(selon)),'descend');
            obj.listsequences(selon(ord),'extra',extra,'showcnt',false,'extraheader',extraheader,'sort','fixed','maxlist',args.maxlist,'bycluster',args.bycluster,'hitdetails',args.hitdetails);
          end
          if sum(seloff)>0  && args.showoff
            [~,ord]=sort(arrayfun(@(z) ratios(z),obj.naseq(seloff)),'ascend');
            obj.listsequences(seloff(ord),'extra',extra,'showcnt',false,'extraheader',extraheader,'sort','fixed','maxlist',args.maxlist,'bycluster',args.bycluster,'hitdetails',args.hitdetails);
          end
          if ~isempty(args.highlight)
            fprintf('Highlighted Sequences:\n');
            hsel=find(ismember(obj.naseq,args.highlight) & ismember(obj.subrun,t1.subruns([t1.groups{2:3}])));
            obj.listsequences(hsel,'extra',extra,'showcnt',true,'extraheader',extraheader,'sort','name','maxlist',length(args.highlight),'bycluster',args.bycluster,'hitdetails',args.hitdetails);
          end
        end
      end
      if args.nscale>1
        astr=sprintf('1/%.0fN',args.nscale);
      else
        astr=sprintf('%.0f/N',1/args.nscale);
      end
      stitle=sprintf('%s mincnt=%d minfold=%.1f alpha=%s',args.title,args.mincnt,args.minfold, astr);
      if ~args.publication && (args.showerr || args.showcleave)
        suptitle(stitle);
      end
      if nargout==0
        % Suppress ans display
        clear x;
      end
    end

    function trpstats(obj,varargin)
      defaults=struct('plotcnts',false,'matches','','sort',false);
      args=processargs(defaults,varargin);

      if isempty(obj.trpexpts)
        fprintf('No TRP data\n');
        return;
      end
      
      fprintf('Counts ignore amplicons and references\n');
      fprintf('%2s %2s %-25.25s %3s %10s %8s %8s/%-6s =   nM %8s/%-6s =   nM %8s/%-6s =   nM %6s\n','#','T#','Description','Run','Target','Conc','In','Ref','Clvd','Ref','Unclvd','Ref','Cleavage RdBalance Amplicon ICU');
      all=nan(length(obj.trpexpts),4);
      descr={obj.trpexpts.descr};
      if args.sort
        [descr,ord]=sort(descr);
      else
        ord=1:length(descr);
      end
      for ii=1:length(obj.trpexpts)
        i=ord(ii);
        t=obj.trpexpts(i);
        if ~isempty(args.matches)
          if isempty(regexp(t.descr,args.matches))
            continue;
          end
        end
        selc=ismember(obj.subrun,t.cleavedsr);
        selu=ismember(obj.subrun,t.uncleavedsr);
        selin=ismember(obj.subrun,t.insr);
        %total=[sum(obj.cnt(selin)), sum(obj.cnt(selc)), sum(obj.cnt(selu))];
        %total=[sum(obj.subruns.rds(ismember(obj.subruns.subrun,t.insr))),sum(obj.subruns.rds(ismember(obj.subruns.subrun,t.cleavedsr))),sum(obj.subruns.rds(ismember(obj.subruns.subrun,t.uncleavedsr)))];
        isamplicon= ismember(obj.naseq,obj.ampliconseqs(1));
        isref = ismember(obj.naseq,obj.refseqs(1));
        isshort = obj.getclass()==1;
        total=[sum(obj.cnt(selin & ~isshort & ~isamplicon & ~isref)), sum(obj.cnt(selc & ~isshort & ~isamplicon & ~isref)), sum(obj.cnt(selu & ~isshort & ~isamplicon & ~isref))];
        sr={t.insr,t.cleavedsr,t.uncleavedsr};
        refcnt=nan(size(sr));
        refconc=nan(size(sr));
        predil=nan(size(sr));
        for k=1:length(sr)
          refcnt(k)=nansum([obj.subruns.normalization(ismember(obj.subruns.subrun,sr{k})).cnt]);
          refconc(k)=nansum([obj.subruns.normalization(ismember(obj.subruns.subrun,sr{k})).conc]);
          if isfield(obj.subruns.normalization,'predil')
            predil(k)=nansum([obj.subruns.normalization(ismember(obj.subruns.subrun,sr{k})).predil]);
          else
            predil(k)=1;
          end
        end
        scaled=(total)./refcnt.*refconc.*predil;
        clv=scaled(2)/(scaled(2)+scaled(3));
        bal=total./scaled;
        if t.concunits=='M'
          cfmt=concfmt(t.targetconc);
        else
          cfmt=sprintf('%g%s',t.targetconc,t.concunits);
        end
        ampliconFrac=sum(obj.cnt((selc|selu)&(isamplicon|isshort)&~isref))/sum(obj.cnt((selc|selu)&~isref));
        ccode=obj.subruns.codes{obj.subruns.subrun==t.cleavedsr(1)};
        ucode=obj.subruns.codes{obj.subruns.subrun==t.uncleavedsr(1)};
        fprintf('%2d %2d %-25.25s %3d %10.10s %8.8s %8d/%-6d %6.1f %8d/%-6d %6.1f %8d/%-6d %6.1f %6.2f%%  %5.1f     %4.1f%% %s,%s%c,%s%c\n',i, t.trpexpt, t.descr, t.robotrun, t.target, cfmt,[total;refcnt;scaled*1e9],clv*100,bal(2)/bal(3),ampliconFrac*100,intlist(t.insr),intlist(t.cleavedsr),ccode(1),intlist(t.uncleavedsr),ucode(1));
        all(i,:)=[total(2:3),refcnt(2:3)];
        all(i,:)=all(i,:)/sum(all(i,:));
      end
      if args.plotcnts
        setfig('trpstats');cla;
        bar(all,'stacked');
        set(gca,'XTick',1:length(obj.trpexpts));
        set(gca,'XTickLabel',descr);
        set(gca,'XTickLabelRotation',45);
        set(gca,'TickLabelInterpreter','None');
        ylabel('Fraction');
        legend('Clvd','Unclvd','RefC','RefU');
      end
    end
    
    function plotbulkbyid(obj,expts)
    % Plot bulk cleavage (based on NGS) by ID
      if nargin<2
        expts=1:length(obj.trpexpts);
      end
      setfig(sprintf('bulkbyid(%s)',sprintf('%d ',expts)));clf;
      for i=1:length(expts)
        if length(expts)>1
          subplot(2,ceil(length(expts)/2),i);
        end
        t=obj.trpexpts(expts(i));
        a=t.analysis;
        bulk=[];
        descr={};
        for i=1:length(t.idsr)
          bulk(i)=nansum(a.cleavage.*a.srccnt(:,i))/nansum(a.srccnt(:,i));
          descr{i}=obj.subruns.desc{obj.subruns.subrun==t.idsr(i)};
          if length(descr{i})==2 && descr{i}(1)=='R'
            descr{i}=[descr{i}(1),'0',descr{i}(2)];
          end
        end
        [descr,o]=sort(descr);
        bulk=bulk(o);
        bulk(end+1)=nanmean(a.cleavage);
        descr{end+1}=['Total'];
        bar(bulk);
        set(gca,'XTick',1:length(bulk));
        set(gca,'XTickLabel',descr);
        set(gca,'XTickLabelRotation',45);
        set(gca,'TickLabelInterpreter','None');
        title(sprintf('Bulk Cleavage %s',t.descr),'Interpreter','none');
        c=axis;
        c(3)=0; c(4)=1.01;
        axis(c);
      end
    end

    function cleavebycond(obj,naseqs,varargin)
    % For each of the naseqs given, produce a bar graph across the conditions for which we have data
      defaults=struct('mincnt',10,'trpexpts',[],'subplots',[1,1],'rescale',1,'ci',[5,95]);
      args=processargs(defaults,varargin);
      if isempty(args.trpexpts)
        args.trpexpts=1:length(obj.trpexpts);
      end
      
      obj.trpanalyze('mincnt',args.mincnt);

      frame=1; sp=1;
      t1=obj.trpexpts(args.trpexpts(1));
      for i=1:length(naseqs)
        naseq=naseqs(i);
        conds={};
        ratio=[];
        ci=[];
        for ii=1:length(args.trpexpts)
          i=args.trpexpts(ii);
          t=obj.trpexpts(i);
          if args.rescale
            [c,ia,ib]=intersect(t1.analysis.naseq,t.analysis.naseq);
            scale=exp(nanmedian(log(t1.analysis.ratio(ia))-log(t.analysis.ratio(ib))));
            if isnan(scale)
              fprintf('Unable to scale %s to %s\n',t.descr, t1.descr);
              continue;
            end
          else
            scale=1;
          end
          ind=find(naseq==t.analysis.naseq);
          if ~isempty(ind) && sum(t.analysis.cnt(ind,2:3))>=args.mincnt && isfinite(t.analysis.cleavage(ind))
            if scale==1
              conds{end+1}=t.descr;
            else
              conds{end+1}=sprintf('%sx%.1f',t.descr,scale);
            end
            ratio(end+1)=t.analysis.ratio(ind)*scale;
            ci(end+1,:)=prctile(t.analysis.ratioboot(ind,:)*scale,args.ci);
          end
        end
        if length(ratio)>1
          if sp==1
            setfig(sprintf('ByCond%d',frame));clf;
            suptitle(sprintf('ByCond%d',frame));
          end
          if prod(args.subplots)>1
            subplot(args.subplots(1),args.subplots(2),sp);
          end
          sp=sp+1;
          if sp>prod(args.subplots)
            sp=1;
            frame=frame+1;
          end
          barh(ratio);
          hold on;
          for i=1:length(ratio)
            plot(ci(i,:),i*[1,1],'r-');
          end
          ax=axis;
          plot(median(ratio)*[1,1],ax(3:4),'g:')
          set(gca,'XScale','log');
          cleaveticks(1,0);
          set(gca,'YTick',1:length(conds));
          set(gca,'YTickLabel',conds);
          %set(gca,'YTickLabelRotation',45);
          set(gca,'TickLabelInterpreter','None');
          xlabel('Cleavage');
          title(sprintf('%d %s CI=[%.0f,%.0f]',naseq,obj.seqs.getlabels(naseq),args.ci));
        end
      end
    end
       
    function mutstats(obj, subrun)
    % Compute mutations statistics for given subrun
      sel=find(obj.subrun==subrun & ~ismember(obj.naseq,union(obj.ampliconseqs(1),obj.refseqs(1))));
      naseq=obj.naseq(sel);
      cnt=obj.cnt(sel);
      len=obj.length(sel);
      seq=obj.seq(sel);
      [basecnt,ind]=max(cnt);
      basenaseq=naseq(ind);
      baseseq=seq{ind};
      baselen=len(ind);
      fprintf('Mutations assuming base = %d\n', basenaseq);
      onedel=sum(cnt(len==baselen-1));
      oneins=sum(cnt(len==baselen+1));
      subsel=find(len==baselen);
      nsubs=zeros(baselen+1,1);
      for ii=1:length(subsel)
        i=subsel(ii);
        nsub=sum(seq{i}~=baseseq);
        nsubs(nsub+1)=nsubs(nsub+1)+cnt(i);
      end
      short=sum(cnt(len<48));
      other=sum(cnt)-onedel-oneins-sum(nsubs)-short;
      data={'Exact',basecnt,'S1',nsubs(2),'S2',nsubs(3),'S>2',sum(nsubs(4:end)),'D1',onedel,'I1',oneins,'Short',short,'Other',other};
      for i=1:2:length(data)
        fprintf('%-20.20s\t%.1f%%\n',data{i},data{i+1}/sum(cnt)*100);
      end
    end
    
    
    function ratio=mutspectrum(obj,naseq,trpexpts,varargin)
    % Show mutation spectrum for given naseq over trpexpts
      defaults=struct('showstruct',true,'showbar',false,'showpcolor',true,'metric','ratio_of_ratio','nlist',5,'varna',true);
      args=processargs(defaults,varargin);

      ratio=[];
      obj.trpanalyze();
      t=obj.trpexpts(trpexpts(1));
      ta=t.analysis;
      if length(trpexpts)==1
        sel=ta.naseq==naseq;
        if sum(sel)==0
          fprintf('No data for %d in trpexpt %d\n', naseq,trpexpts(1));
          return;
        end
        baseline=ta.ratio(ta.naseq==naseq);
        ti=sprintf('%s cleavage - %s',obj.seqs.getname(naseq),t.descr);
      elseif length(trpexpts)==2
        t2=obj.trpexpts(trpexpts(2));
        ta2=t2.analysis;
        if ~ismember(naseq,ta.naseq) || ~ismember(naseq,ta2.naseq)
          fprintf('No data for %d in trpexpt %d and/or %d\n', naseq,trpexpts);
          return;
        end
          
        baseline=exp(obj.bootfold(ta2.ratioboot(ta2.naseq==naseq,:),ta.ratioboot(ta.naseq==naseq,:),'metric',args.metric));
        ti=sprintf('%s switching - %s vs. %s',obj.seqs.getname(naseq),t.descr,t2.descr);
      else
        error('Expected 1 or 2 trpexpts, got %d\n', length(trpexpts));
      end

      seq=obj.getseq(naseq);
      seq=seq{1};

      actg='ACTG';
      acug='ACUG';


      ratio=ones(length(seq),length(actg))*baseline;
      allseq=nan(length(seq),length(actg));
      for i=1:length(seq)
        mutseq=seq;
        for j=1:4
          if actg(j)==seq(i)
            continue;
          end
          mutseq(i)=actg(j);
          mutnaseq=obj.seqs.findseq(mutseq);
          if mutnaseq>0
            ind=find(ta.naseq==mutnaseq);
            if isempty(ind)
              %fprintf('No data for naseq %d (%c%d%c) in %s\n', mutnaseq, seq(i), i, actg(j), t.descr);
              ratio(i,j)=nan;
              continue;
            end
            allseq(i,j)=mutnaseq;
            if length(trpexpts)==1
              ratio(i,j)=ta.ratio(ind);
            else
              ind2=find(ta2.naseq==mutnaseq);
              if isempty(ind2)
                fprintf('No data for naseq %d (%c%d%c) in %s\n', mutnaseq, seq(i), i, actg(j), t2.descr);
                ratio(i,j)=nan;
                continue;
              end
              ratio(i,j)=exp(obj.bootfold(ta2.ratioboot(ind2,:),ta.ratioboot(ind,:),'metric',args.metric));
            end
          end
        end
      end
      
      if args.showbar
        setfig([ti,'-bar']);clf;
        h=[];
        colors='rbgm';
        for i=1:size(ratio,2)
          h(i)=bar(ratio(:,i),colors(i),'BaseValue',baseline);
          alpha(0.3);
          hold on;
        end 
        for i=1:size(ratio,1)
          for j=1:size(ratio,2)
            if isfinite(ratio(i,j))
              if ratio(i,j)==max(ratio(i,:))
                align='bottom';
              elseif ratio(i,j)==min(ratio(i,:))
                align='top';
              else
                align='middle';
              end
              text(i,ratio(i,j),acug(j),'HorizontalAlignment','center','VerticalAlignment',align,'FontSize',12);
            end
          end
        end
        set(gca,'YScale','log');
        alpha(0.2);
        ax=axis;
        ax(3)=min(ratio(:));
        ax(4)=max(ratio(:));
        r0=(ax(4)/ax(3))^(1/10);
        ax(4)=ax(4)*r0;
        ax(3)=ax(3)/r0;
        axis(ax);
        if length(trpexpts)==1
          cleaveticks(0,1);
          ax2=axis;
          ax2(3)=min(ax(3),ax2(3));
          ax2(4)=max(ax(4),ax2(4));
          axis(ax2);
          baseval=baseline/(1+baseline);  % Convert to cleavage
        else
          logticks(0,1);
          baseval=baseline;
        end
        text(1,baseline,sprintf('%.2f',baseval),'HorizontalAlignment','left','VerticalAlignment','bottom');

        % Add vertical marks at places that we usually break seq strings
        seqfmt=obj.seqs.naseqformat(naseq);
        k=1;
        ax=axis;
        for i=1:length(seqfmt)
          if seqfmt(i)==' '
            plot(k*[1,1]-0.5,ax(3:4),'k');
          else
            k=k+1;
          end
        end
        set(gca,'TickLength',[0.002,0.002]);
        for i=1:length(seq)
          text(i,ax(3),strrep(seq(i),'T','U'),'HorizontalAlignment','center','VerticalAlignment','bottom');
        end
        legend(h,{'A','C','U','G'});
        if length(trpexpts)==1
          ylabel('Fraction Cleaved');
        else
          ylabel(args.metric,'Interpreter','none');
        end
        title(ti);
        set(gca,'XGrid','on');
        set(gca,'XMinorGrid','on');
        %axis tight;
      end

      if length(trpexpts)==1
        cleavage=ratio./(1+ratio);
      end
      
      if args.showpcolor
        setfig([ti,'-pcolor']);clf;
        pc=ratio';
        if length(trpexpts)==1
          pc=log10(ratio');
        else
          pc=log10(ratio');
        end
        pc(end+1,:)=nan;
        pc(:,end+1)=nan;
        pcolor(pc);
        hold on;
        h=colorbar;
        if length(trpexpts)==1
          set(get(h,'Label'),'string','Fraction Cleaved')
          cleaveticks(0,0,'colorbar',h);
        elseif strcmp(args.metric,'ratio_of_uncleavage')
          set(get(h,'Label'),'string','Uncleaved Ratio')
        elseif strcmp(args.metric,'ratio_of_ratio')
          cax=caxis;
          cax(1)=0;
          caxis(cax);
          logticks(0,0,1);
          set(get(h,'Label'),'string','Fold Change of Ratio')
        else
          set(get(h,'Label'),'string',args.metric,'Interpreter','none');
        end
        cax=caxis;
        set(gca,'YTick',1.5:1:4.5);
        set(gca,'YTickLabel',{'A','C','U','G'});
        shading flat
        
        % Add vertical marks at loop boundaries
        segext=obj.seqs.getsegextents(naseq);
        k=1;
        ax=axis;
        markpos=segext([4,8],:);
        for k=1:size(markpos,1)
          plot(markpos(k,1)*[1,1],ax(3:4),'k','Color','red','LineWidth',3);
          plot(markpos(k,2)*[1,1],ax(3:4),'k','Color','red','LineWidth',3);
          plot(markpos(k,:),ax(3)*[1,1],'k','Color','red','LineWidth',3);
          plot(markpos(k,:),ax(4)*[1,1],'k','Color','red','LineWidth',3);
        end

        % Shorten ticks
        set(gca,'TickLength',[0.002,0.002]);

        % Mark unmutated seq with +
        for i=1:length(seq)
          ypos=find(seq(i)==actg);
          plot(i+0.5,ypos+0.5,'+w');
          %text(i+0.5,ax(3),strrep(seq(i),'T','U'),'HorizontalAlignment','center','VerticalAlignment','bottom');
        end

        title(ti);
        set(gca,'XGrid','on');
        set(gca,'XMinorGrid','on');
        set(gca,'TickDir','out');
        lbls=get(gca,'XTickLabel');
        set(gca,'XTick',get(gca,'XTick')+0.5);  % Shift ticks
        set(gca,'XTickLabel',lbls);  % So we use old values
        xlabel('Position');
        ylabel('Mutation');
        %axis tight;
      end
      
      sortratio=sort(ratio(isfinite(ratio(:))),'descend');
      if length(trpexpts)==1
        fprintf('Listing sequences with cleavage>=%.2f\n', 1./(1+1./sortratio(min(end,args.nlist))));
      else
        fprintf('Listing sequences with ratio>=%.2f\n', sortratio(min(end,args.nlist)));
      end
      for i=1:size(ratio,1)
        for j=1:size(ratio,2)
          if isfinite(ratio(i,j))
            if ratio(i,j)>=sortratio(min(end,args.nlist))
              if length(trpexpts)==1
                fprintf('%5.2f %d %c%d%c ',1./(1+1./ratio(i,j)),allseq(i,j),seq(i),i,acug(j));
              else
                fprintf('%5.2f %d %c%d%c',ratio(i,j),allseq(i,j),seq(i),i,acug(j));
              end
              fprintf('%s\n', obj.seqs.getname(seq(i)));
            end
          end
        end
      end
      
      if ~args.showstruct
        return;
      end
      
      setfig([ti,'-structure']);clf;

      if length(trpexpts)==1
        rtmp=log10(ratio);
      else
        rtmp=log10(ratio);
        %rtmp(rtmp==baseline)=nan;   % Blank out the original so we plot the effect of the best mutation
      end

      if ~args.varna
        val=max(rtmp,[],2);
        rng=[min(val),max(val)];
        cm=colormap();
        ind=round(interp1(rng,[1,size(cm,1)],val));
        cm=[1,1,1;cm];   % White for missing data
        ind=ind+1; % (2:n+1)
        ind(isnan(ind))=1;  % Maps to white
        col=cm(ind,:);
        ribodraw(seq,obj.seqs.getloop1(naseq),obj.seqs.getloop2(naseq),col);
        caxis(rng);
      else
        % Map effects onto 2d structure
        numnucs=4;
        imgc={};
        if ~exist('cax','var')
          % Warning: cax not set
          cax=[0,4];
        end
        for nuc=1:numnucs
          imgc{nuc}=obj.seqs.pcolorstruct(naseq,rtmp(:,nuc),'caxis',cax);
        end
        img=[imgc{1},imgc{2};imgc{3},imgc{4}];
        imshow(img);
        if numnucs>1
          locs=[0.25,0
               .755,0
                0.25,0.5
                .75,0.5];
          for nuc=1:4
            text(locs(nuc,1)*size(img,2),locs(nuc,2)*size(img,1),sprintf('Mutations to %c',acug(nuc)),'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','top');
          end
        end
        h=colorbar;
        caxis(cax);   % Force the colorbar to match
        if length(trpexpts)==1
          set(get(h,'Label'),'string','Fraction Cleaved')
          cleaveticks(0,0,'colorbar',h);
        elseif strcmp(args.metric,'ratio_of_uncleavage')
          set(get(h,'Label'),'string','Uncleaved Ratio')
        elseif strcmp(args.metric,'ratio_of_ratio')
          logticks(0,0,1);
          set(get(h,'Label'),'string','Fold Change of Ratio')
        else
          set(get(h,'Label'),'string',args.metric,'Interpreter','none');
        end
        title(ti);
        pause(0.1);
      end
    end
    
    function x=getsubrunname(obj,sr)
      sel=obj.subruns.subrun==sr;
      assert (sum(sel)==1);
      x=sprintf('%s-%s',obj.subruns.desc{sel},obj.subruns.codes{sel});
      if isfield(obj.subruns,'bclen')
        x=sprintf('%s(%d)',x,obj.subruns.bclen(sel));
      end
    end

    function sr=getsubrunsbydesc(obj,re,code)
    % Get all subruns whose name match desc RE
      s=regexp(obj.subruns.desc,re);
      match=cellfun(@(z) ~isempty(z), s);
      if nargin==3
        s2=regexp(obj.subruns.codes,code);
        match=match & cellfun(@(z) ~isempty(z), s2);
      end
      sr=obj.subruns.subrun(match);
    end
    
    function checkifmutant(obj,naseq,subrun,varargin)
    % Check if given naseq might be a mutant in given subrun
      defaults=struct('minratio',50,'silent',false);
      args=processargs(defaults,varargin);
      if nargin<3 || isempty(subrun)
        % Use the subrun with the highest count for naseq
        sel=find(obj.naseq==naseq);
        cnts=obj.cnt(sel);
        [~,ind]=max(cnts);
        subrun=obj.subrun(sel(ind));
        if ~args.silent
          fprintf('Checking against subrun %d (%s) cnt=%d\n', subrun, obj.getsubrunname(subrun),cnts(ind));
        end
      end
      cnt=obj.cnt(obj.naseq==naseq & obj.subrun==subrun);
      if isempty(cnt)
        cnt=0;
      end
      thresh=cnt*args.minratio;
      others=obj.naseq(obj.cnt>=thresh & obj.subrun==subrun);
      if isempty(others)
        return;
      end
      other=obj.seqs.compare(naseq, others,'silent',args.silent);
      if isempty(other)
        if ~args.silent
          fprintf('No sequences with >=%d counts similar to %d\n', thresh, naseq);
        end
      else
        ocnt=obj.cnt(obj.naseq==other & obj.subrun==subrun);
        fprintf('%d (cnt=%d) may be a mutant of %d (cnt=%d) - ratio = %.0f in subrun %d\n', naseq, cnt, other, ocnt, ocnt/cnt,subrun);
      end
    end
    
    function x=comparefreq(obj,s1,s2,ind,varargin)
    % Compare frequency of sequences in subruns s1 and s2
    % Optionally, only look at entries given in ind
      defaults=struct('mincnt',10,'zthresh',3.0,'lengthplot',true,'maxlist',30,'minlength',0,'minratio',0,'maxlength',0);
      args=processargs(defaults,varargin);
      x=obj.commonseqs([s1,s2],args.mincnt);
      x.length=obj.getlength(x.naseq);

      % naseqs to ignore
      ignores=union(obj.refseqs(1), obj.ampliconseqs(1));
      sel=x.length>30  & ~ismember(x.naseq,ignores);
      if nargin>3 && ~isempty(ind)
        sel=sel & ismember(x.naseq,obj.naseq(ind));
      end
      if args.minlength>0
        sel=sel & x.length>=args.minlength;
      end
      if args.maxlength>0
        sel=sel & x.length<=args.maxlength;
      end
      x.cnt=x.cnt(sel,:);
      x.cntfuzzed=x.cnt+(rand(size(x.cnt))-0.5).*2.*min(x.cnt/50,0.5);  % Fuzz the counts for plotting separate points
      x.cntfuzzed=max(x.cntfuzzed,0.5);   % Plot 0 as 0.5 so it shows up on log plots
      x.naseq=x.naseq(sel);
      x.selector=x.selector(sel);
      x.length=x.length(sel);
      x.loop1len=x.loop1len(sel);
      x.loop2len=x.loop2len(sel);
      %rdratio=obj.subruns.rds(obj.subruns.subrun==s2)/obj.subruns.rds(obj.subruns.subrun==s1);
      rdratio=nanmedian((x.cnt(:,2)+1)./(x.cnt(:,1)+1));
      ti=sprintf('comparefreq(%s,%s) N=%d',obj.getsubrunname(s1),obj.getsubrunname(s2),length(x.naseq));
      if args.minlength>0
        ti=sprintf('%s, minlen=%d',ti,args.minlength);
      end
      setfig(ti);cla;
      if args.lengthplot
        subplot(211);
      end
      islong = x.length>91;
      loglog(x.cntfuzzed(~islong,1),x.cntfuzzed(~islong,2),'.b');
      % Plot long ones with a o
      hold on;
      loglog(x.cnt(islong,1),x.cnt(islong,2),'ob');
      
      xlabel(obj.getsubrunname(s1),'Interpreter','none');
      ylabel(obj.getsubrunname(s2),'Interpreter','none');
      
      % Calculate z-statistics for testing whether fractions are different for each seq
      zstat=[];
      totalcnt=sum(x.cnt);
      for i=1:length(x.naseq)
        p1=x.cnt(i,1)/totalcnt(1);
        p2=x.cnt(i,2)/totalcnt(2);
        p=sum(x.cnt(i,:))/sum(totalcnt);
        zstat(i)=(p2-p1)/sqrt(p*(1-p)*(1/totalcnt(1)+1/totalcnt(2)));
      end
      x.zstat=zstat;
      
      [~,ord]=sort(abs(zstat),'desc');

      c=axis;
      hold on;
      plot(c(1:2),c(1:2)*rdratio,':g');

      for k=1:2
        if k==1
          fprintf('%s Read ratio=%.3f, z>=%.1f\n',ti, rdratio,args.zthresh);
          kord=ord(zstat(ord)>=args.zthresh);
        else
          fprintf('%s Read ratio=%.3f, z<= -%.1f\n',ti, rdratio,args.zthresh);
          kord=ord(zstat(ord)<=-args.zthresh);
        end
        ratio=(x.cnt(:,2)+1)./(x.cnt(:,1)+1)/rdratio;
        kord=kord(ratio(kord) > args.minratio | ratio(kord)<1/args.minratio);
        if length(kord)>0
          extra=containers.Map('KeyType','double','ValueType','char');
          sortkey=containers.Map('KeyType','double','ValueType','double');
          extraheader='ratio   n1    n2 zstat';
          for ii=1:length(kord)
            i=kord(ii);
            extra(x.naseq(i))=sprintf('%5.2f %5d %5d %4.1f', ratio(i), x.cnt(i,:),zstat(i));
            sortkey(x.naseq(i))=-abs(zstat(i));
          end
          obj.listsequences(ismember(obj.naseq,x.naseq(kord)),'extra',extra,'extraheader',extraheader,'sort',sortkey,'maxlist',args.maxlist);
          plot(x.cntfuzzed(kord,1),x.cntfuzzed(kord,2),'or');
        end
      end
      logticks(1,1);
      lbls=get(gca,'XTickLabels');lbls=strrep(lbls,'0.5','0');set(gca,'XTickLabels',lbls);
      lbls=get(gca,'YTickLabels');lbls=strrep(lbls,'0.5','0');set(gca,'YTickLabels',lbls);
      
      if args.lengthplot
        % Plot dependence on length
        subplot(212);
        if length(x.length)>100
          sym='.';
        else
          sym='o';
        end
        semilogy(x.length+(rand(size(x.length))-0.5)/2,x.cnt(:,2)./x.cnt(:,1),sym);
        hold on;
        sel=abs(x.zstat)>args.zthresh;
        if sum(sel)>0
          semilogy(x.length(sel)+(rand(size(x.length(sel)))-0.5)/2,x.cnt(sel,2)./x.cnt(sel,1),'ro');
        end
        logticks(0,1);
        xlabel('Length (nt)');
        ylabel('Count Ratio');
        suptitle(sprintf('Count Compare (N=%d, mincnt=%d)',size(x.cnt,1),args.mincnt));
      else
        title(sprintf('Count Compare (N=%d, mincnt=%d)',size(x.cnt,1),args.mincnt));
      end
      
    end
    
    
    function [frac,fracs,naseq]=comparefraction(obj,naseq,subruns,varargin)
    % Compare fraction of naseqs over given subruns (for which it is expected that the fraction is constant, such as with replicates)
      defaults=struct('mincnt',100);
      args=processargs(defaults,varargin);
      % naseqs to ignore
      ignores=union(obj.refseqs(1),obj.ampliconseqs(1));
      ignores=union(ignores,obj.naseq(obj.segmentid(10)~=1 | obj.segmentid(2)~=1));  % Only ones that start with GCTGTC and end with GAAACAGC
      if isempty(naseq)
        % Use any sequences which have adequate counts
        naseq=obj.naseq(obj.subrun==subruns(1) & obj.cnt>=args.mincnt);
        naseq=setdiff(naseq,ignores);
        fprintf('Analyzing %d sequences...\n', length(naseq));
      end

      for i=1:length(naseq)
        names{i}=obj.seqs.getname(naseq(i));
        if isempty(names{i})
          names{i}=obj.seqs.getlabels(naseq(i));
        end
        if isempty(names{i})
          names{i}=sprintf('%d',naseq(i));
        end
        names{i}=sprintf('L%03d-%d',obj.seqs.getlength(naseq(i)),naseq(i));
      end

      [names,ord]=sort(names);
      naseq=naseq(ord);
      
      ti=sprintf('Fraction Comparison [%s]',strjoin(arrayfun(@(z) sprintf('%d',z),subruns,'UniformOutput',false),','));
      fprintf('Comparing: ');
      for i=1:length(subruns)
        fprintf(' %s',obj.getsubrunname(subruns(i)));
      end
      fprintf('\n');
      
      enrich=zeros(length(naseq),1);
      total=[];lbls={};ignore=[];
      for i=1:length(subruns)
        total(i)=obj.subruns.rds(obj.subruns.subrun==subruns(i));
        ignore(i)=sum(obj.cnt(obj.subrun==subruns(i) & ismember(obj.naseq,ignores)));
        lbls{i}=obj.getsubrunname(subruns(i));
      end
      fprintf('Ignored %d/%d (%.1f%%)\n', sum(ignore),sum(total),sum(ignore)/sum(total)*100);
      total=total-ignore;

      frac=[]; fracs=[];
      for j=1:length(naseq)
        match(j,:)=arrayfun(@(z) sum(obj.cnt(obj.subrun==z & obj.naseq==naseq(j))), subruns);
        frac(j,:)=match(j,:)./total;
        fracs(j,:)=frac(j,:)/median(frac(j,:));
        fracerr(j,:)=sqrt(match(j,:))./total;
      end

      setfig(ti);clf;
      boxplot(fracs');
      set(gca,'XTick',1:length(naseq));
      set(gca,'TickLabelInterpreter','None');
      set(gca,'XTickLabel',names);
      set(gca,'XTickLabelRotation',45);
      set(gca,'YScale','log');
      logticks(0,1);
      ylabel('Fraction');
      title(ti,'Interpreter','none');
      
      setfig([ti,' by subrun']);clf;
      boxplot(fracs);
      set(gca,'XTick',1:length(subruns));
      set(gca,'TickLabelInterpreter','None');
      set(gca,'XTickLabel',lbls);
      set(gca,'XTickLabelRotation',45);
      set(gca,'YScale','log');
      logticks(0,1);
      ylabel('Fraction');
      title([ti,' by subrun'],'Interpreter','none');
      
    end
      
    function [frac,naseq]=plotenrichment(obj,naseq,subruns,varargin)
    % Plot enrichment (fraction of total) of given naseq over subruns (or all if naseq not given)
    % 'named' is 0 (just go by enrichment order), 1 (also include all named seqs even if poor enrichment), 2 (only named ones)
      defaults=struct('print',true,'mincnt',100,'named',1,'maxplot',20,'minenrich',0,'plotstep',true,'plotdist',true,'lastround',nan);
      args=processargs(defaults,varargin);
      % naseqs to ignore
      ignores=union(obj.refseqs(1),obj.ampliconseqs(1));
      if isempty(naseq)
        % Use any sequences which have adequate counts
        naseq=union(obj.naseq(obj.subrun==subruns(1) & obj.cnt>=args.mincnt), obj.naseq(obj.subrun==subruns(end) & obj.cnt>=args.mincnt));
        naseq=setdiff(naseq,ignores);
        if args.named==2
          naseq=intersect(naseq,[obj.seqs.seqnames.naseq]);
        end
        fprintf('Analyzing %d sequences...\n', length(naseq));
      end
      if length(naseq)==1
        ti=sprintf('Enrichment of %.20s',obj.seqs.getlabels(naseq));
      elseif args.named==2
        ti='Enrichment of Named Sequences';
      elseif args.named==1
        ti=sprintf('Named and %d Highest-enriched Sequences',args.maxplot);
      else
        ti=sprintf('%d Highest-enriched Sequences',args.maxplot);
      end
      ti=sprintf('%s %s-%s',ti,obj.subruns.desc{obj.subruns.subrun==subruns(1)},obj.subruns.desc{obj.subruns.subrun==subruns(end)});
      enrich=zeros(length(naseq),1);
      total=[];lbls={};ignore=[];
      for i=1:length(subruns)
        total(i)=obj.subruns.rds(obj.subruns.subrun==subruns(i));
        ignore(i)=sum(obj.cnt(obj.subrun==subruns(i) & ismember(obj.naseq,ignores)));
        lbls{i}=obj.subruns.desc{obj.subruns.subrun==subruns(i)};
      end
      fprintf('Ignored %d/%d (%.1f%%)\n', sum(ignore),sum(total),sum(ignore)/sum(total)*100);
      total=total-ignore;

      for j=1:length(naseq)
        match(j,:)=arrayfun(@(z) sum(obj.cnt(obj.subrun==z & obj.naseq==naseq(j))), subruns);
        frac(j,:)=match(j,:)./total;
        fracerr(j,:)=sqrt(match(j,:))./total;
      end

      enrich=(frac(:,end)./frac(:,1)).^(1/(length(subruns)-1));
      f0=log10(frac(:,end))-log10(enrich).*args.lastround;
      setfig([ti,'-ebar']);clf;
      [~,ord]=sort(enrich,'descend');
      h=[];legs={};
      for ii=1:length(ord)
        i=ord(ii);
        if enrich(i)<args.minenrich
          break;
        end
        nm=obj.seqs.getname(naseq(i));
        if ii>args.maxplot && (isempty(nm) || args.named==0)
          continue;
        end
        h(end+1)=errorbar(frac(i,:),fracerr(i,:),'o-');
        hold on;
        if isempty(nm)
          nm=sprintf('%d',naseq(i));
        else
          set(h(end),'LineWidth',get(h(end),'LineWidth')*4);
        end
        fprintf('%9d %-12.12s %6.0f  %6.4f\n', naseq(i), nm,nanmean(match(i)),enrich(i));
        legs{end+1}=sprintf('%s \\eta=%4.2f f_0=10^{%.0f}',nm,enrich(i),f0(i));
      end
      set(gca,'XTick',1:length(subruns));
      set(gca,'TickLabelInterpreter','None');
      set(gca,'XTickLabel',lbls);
      set(gca,'XTickLabelRotation',45);
      set(gca,'YScale','log');
      ylabel('Fraction');
      legend(h,legs,'Location','EastOutside');
      title(ti,'Interpreter','none');
      if length(enrich)>20 && args.plotdist
        setfig([ti,'-dist']);clf;
        histogram(enrich);
        xlabel('Enrichment');
        title(sprintf('%s (mean=%.2f)', ti, mean(enrich)),'Interpreter','none');
      end
      if length(naseq)~=1 && args.plotstep
        % Plot step enrichment
        ti=sprintf('Enrich([%s],mincnt=%d)',sprintf('%d,',subruns),args.mincnt);
        setfig([ti,'-track']);clf;
        pnum=1;
        for i=1:length(subruns)-2
          for j=i+1:length(subruns)-1
            subplot(length(subruns)-2,length(subruns)-2,(i-1)*(length(subruns)-2)+j-1)
            loglog(frac(:,j+1)./frac(:,j),frac(:,i+1)./frac(:,i),'.');
            hold on;
            sel1=enrich>1.1;loglog(frac(sel1,j+1)./frac(sel1,j),frac(sel1,i+1)./frac(sel1,i),'.g');
            sel2=enrich<0.9;loglog(frac(sel2,j+1)./frac(sel2,j),frac(sel2,i+1)./frac(sel2,i),'.r');
            xlabel(sprintf('%s->%s',obj.subruns.desc{obj.subruns.subrun==subruns(j)},obj.subruns.desc{obj.subruns.subrun==subruns(j+1)}),'Interpreter','none');
            ylabel(sprintf('%s->%s',obj.subruns.desc{obj.subruns.subrun==subruns(i)},obj.subruns.desc{obj.subruns.subrun==subruns(i+1)}),'Interpreter','none');
            axis([0.5,2.0,0.5,2.0]);
          end
        end
        suptitle(ti);
      end
    end
      
    function x=fractionbyround(obj,naseq,subruns,rounds,varargin)
    % Plot enrichment (fraction of total) of given naseqs over subruns, which refer to round numbers given in rounds
    % If subruns is a cell array of strings, they should be of the format 'RUN.descr.code' (e.g. 'NGS108.R87A.T7W/X')
    %    and in which case use the database directly
    % If args.trpexpts is not empty, it is a map from round number to trpexpt for that round
      defaults=struct('other',[],'mincnt',100,'maxplot',20,'maxlegend',20,'trpexpts',[],'minclvcnt',30,'groupbycluster',true,'seqonly',false);
      args=processargs(defaults,varargin);

      if iscell(naseq)
        % Using names instead of ids
        names=naseq;
        naseq=nan(size(names));
        for i=1:length(names)
          x=obj.seqs.findname(names{i});
          if isempty(x)
            error('Seq name "%s" not found\n', names{i});
          end
          naseq(i)=x;
        end
      end
      
      % naseqs to ignore
      amplicons=obj.ampliconseqs(1);
      refs=obj.refseqs(1);

      if isempty(naseq)
        if iscell(subruns)
          error('naseqs must be specified is subruns is a cell array'); 
        end
        % Use any sequences which have adequate counts
        naseq=unique(obj.naseq(ismember(obj.subrun,subruns) & obj.cnt>=args.mincnt));
        naseq=setdiff(naseq,union(amplicons,refs));
        fprintf('Analyzing %d sequences...\n', length(naseq));
      end

      ti='Enrichment';
      ti=sprintf('%s R%d-R%d',ti,min(rounds),max(rounds));

      urnds=unique(rounds,'sorted');
      total=[];lbls={};ignore=[];
      if iscell(subruns)
        % Direct database access
        % Get total rds, amplicons, refs per round
        naseqlist=strjoin(arrayfun(@(z) sprintf('%.0f',z),naseq,'Unif',false),',');
        ampliconlist=strjoin(arrayfun(@(z) sprintf('%.0f',z),amplicons,'Unif',false),',');
        reflist=strjoin(arrayfun(@(z) sprintf('%.0f',z),refs,'Unif',false),',');

        cmd1='';cmd2='';cmd3='';cmd4='';
        for i=1:length(subruns)
          s=strsplit(subruns{i},'.');
          q=sprintf('SELECT %d rnd,sum(cnt) FROM %s.ngsentries e WHERE e.subrun=(SELECT s.subrun FROM %s.v_subruns s WHERE s.subrun AND s.descr=''%s'' AND s.code=''%s'')',rounds(i), s{1}, s{1}, s{2}, s{3});
          q2=sprintf('%s AND naseq IN (%s)',q, ampliconlist);
          q3=sprintf('%s AND naseq IN (%s)',q, reflist);
          q4=sprintf('SELECT %d rnd,naseq,cnt FROM %s.ngsentries e WHERE e.subrun=(SELECT s.subrun FROM %s.v_subruns s WHERE s.subrun AND s.descr=''%s'' AND s.code=''%s'')',rounds(i), s{1}, s{1}, s{2}, s{3});
          q4=sprintf('%s AND naseq IN (%s)',q4, naseqlist);
          if i>1
            cmd1=sprintf('%s\nUNION ',cmd1);
            cmd2=sprintf('%s\nUNION ',cmd2);
            cmd3=sprintf('%s\nUNION ',cmd3);
            cmd4=sprintf('%s\nUNION ',cmd4);
          end
          cmd1=[cmd1,q];
          cmd2=[cmd2,q2];
          cmd3=[cmd3,q3];
          cmd4=[cmd4,q4];
        end
        fprintf('Executing total query...');
        [r1,c1]=mysql(cmd1);
        fprintf('done\nExecuting amplicon query...');
        [r2,c2]=mysql(cmd2);
        fprintf('done\nExecuting ref query...');
        [r3,c3]=mysql(cmd3);
        fprintf('done\nExecuting naseq-specifc query...');
        [r4,n4,c4]=mysql(cmd4);
        fprintf('done\n');
        
        for i=1:length(urnds)
          total(i)=nansum(c1(r1==urnds(i)));
          ampliconcnt(i)=nansum(c2(r2==urnds(i)));
          refcnt(i)=nansum(c3(r3==urnds(i)));
        end
      else
        for i=1:length(urnds)
          sr=subruns(rounds==urnds(i));
          total(i)=sum(obj.subruns.rds(ismember(obj.subruns.subrun,sr)));
          ampliconcnt(i)=sum(obj.cnt(ismember(obj.subrun,sr) & ismember(obj.naseq,amplicons)));
          refcnt(i)=sum(obj.cnt(ismember(obj.subrun,sr) & ismember(obj.naseq,refs)));
        end
      end
      fprintf('Ignored %d refs, %d amplicons out of %d total reads (%.1f%%)\n', nansum(refcnt), nansum(ampliconcnt),nansum(total),nansum(ampliconcnt+refcnt)/nansum(total)*100);
      total=total-ampliconcnt-refcnt;

      cnt=nan(length(naseq),length(urnds));
      frac=nan(length(naseq),length(urnds));
      fracerr=nan(length(naseq),length(urnds));
      for j=1:length(naseq)
        for k=1:length(urnds)
          if iscell(subruns)
            cnt(j,k)=nansum(c4(r4==urnds(k) & n4==naseq(j)));
          else
            cnt(j,k)=sum(obj.cnt(ismember(obj.subrun,subruns(rounds==urnds(k))) & obj.naseq==naseq(j)));
          end
        end
        frac(j,:)=(cnt(j,:)+1)./(total+2);
        fracerr(j,:)=sqrt(cnt(j,:)+1-1e-10)./(total+2);
      end
      fracerr(frac(:)==0)=nan;
      frac(frac(:)==0)=nan;
      rounds=urnds;   % In sorted order now
      
      if ~isempty(args.other)
        for i=1:length(args.other)
          x=args.other(i);
          for j=1:length(x.naseq)
            jj=find(x.naseq(j)==naseq);
            if isempty(jj)
              jj=length(naseq)+1;
              naseq(jj)=x.naseq(j);
            end
            for k=1:length(x.rounds)
              kk=find(x.rounds(k)==rounds);
              if isempty(kk)
                kk=length(rounds)+1;
                rounds(kk)=x.rounds(k);
                ampliconcnt(1,kk)=x.ampliconcnt(k);
                total(1,kk)=x.total(k);
              else
                fprintf('Warning not overwriting ampliconcnt(1,%d) or total(1,%d) with other\n', kk, kk);
              end
              frac(jj,kk)=x.frac(j,k);
              fracerr(jj,kk)=x.fracerr(j,k);
              cnt(jj,kk)=x.cnt(j,k);
            end
          end
        end
        % Sort by round
        [rounds,ord]=sort(rounds);
        cnt=cnt(:,ord);
        frac=frac(:,ord);
        fracerr=fracerr(:,ord);
      end
      
      % Only keep maxplot
      if size(frac,1)>args.maxplot
        totalcnt=sum(cnt,2);
        [~,ord]=sort(totalcnt,'desc');
        sel=ord(1:args.maxplot);
        fprintf('Plotting %d/%d (maxplot=%d)\n', length(sel), length(ord),args.maxplot);
        frac=frac(sel,:);
        fracerr=fracerr(sel,:);
        naseq=naseq(sel);
        cnt=cnt(sel,:);
      elseif size(frac,1)<args.maxplot
        fprintf('Only have %d/%d to plot -- might want to decrease mincnt from %d\n', size(frac,1), args.maxplot, args.mincnt);
      end
      
      % Labels
      for i=1:length(rounds)
        if mod(rounds(i),2)==0  || length(rounds)<20
          lbls{i}=sprintf('%d',rounds(i));
        else
          lbls{i}='';
        end
      end
      
      % Plot
      setfig(ti);clf;
      h=[];legs={};
      for i=1:length(naseq)
        nm=obj.seqs.getname(naseq(i));
        h(end+1)=errorbar(frac(i,:),fracerr(i,:),'o');
        hold on;
        for k=1:size(frac,2)-1
          line(k:k+1,frac(i,k:k+1),'Color',get(h(end),'Color'));
        end
        if isempty(nm)
          nm=sprintf('%d',naseq(i));
        else
          set(h(end),'LineWidth',get(h(end),'LineWidth')*4);
        end
        legs{end+1}=nm;
      end
      jumps=find(rounds(1:end-1)+1 ~= rounds(2:end));   % Jump in round number
      
      set(gca,'XTick',1:length(subruns));
      set(gca,'TickLabelInterpreter','None');
      set(gca,'XTickLabel',lbls);
      set(gca,'XTickLabelRotation',45);
      set(gca,'YScale','log');
      % Separate jumps in round number
      ax=axis;
      for i=1:length(jumps);
        r=jumps(i);
        plot((r+0.5)*[1,1],ax(3:4),':');
      end

      xlabel('Selection Round');
      ylabel('Relative Abundance');
      if args.maxlegend>0
        legend(h(1:min(end,args.maxlegend)),legs(1:min(end,args.maxlegend)),'Location','EastOutside');
      end
      title(ti,'Interpreter','none');
      
      % Sort by peak position (which round has maximum frac for each naseq)
      peakpos=[];
      for i=1:size(frac,1)
        sel=cnt(i,:)>0;
        f=frac(i,:).*sel;
        peakpos(i)=nansum(f.*(1:size(frac,2)))/nansum(f);
      end
      [~,ord]=sort(peakpos);
              
      % Group together ones from same cluster
      cluster=obj.getclusterindex(naseq);
      if args.groupbycluster
        for ii=1:length(naseq)-1
          if cluster(ord(ii))==cluster(ord(ii+1))
            continue;
          end
          i=ord(ii);
          for jj=ii+2:length(naseq)
            j=ord(jj);
            if cluster(i)==cluster(j)
              ord=ord([1:ii,jj,ii+1:jj-1,jj+1:end]);
              break;
            end
          end
        end
      end
      naseq=naseq(ord);
      cluster=cluster(ord);
      frac=frac(ord,:);
      fracerr=fracerr(ord,:);
      cnt=cnt(ord,:);
      
      % Pcolor plot
      setfig([ti,'-pc']);clf;
      data=frac;
      data(end+1,:)=ampliconcnt./total;
      data(end+1,:)=1./total;
      data(end+1,:)=nan; data(:,end+1)=nan;
      pcolor(log10(data));
      h1=gca;
      caxis(log10([nanmin(frac(:)),nanmax(frac(:))]));
      h=colorbar;
      %set(h,'TickLabels',arrayfun(@(z) sprintf('%.1e',10^z),get(h,'Ticks'),'Unif',false));
      set(gca,'XTick',(1:length(rounds))+0.5);
      set(gca,'TickLabelInterpreter','None');
      set(gca,'XTickLabel',arrayfun(@(z) sprintf('%d',z),rounds,'Unif',false));
      set(gca,'YTick',(1:length(naseq)+2)+0.5)
      ylbl={};
      hitseqs=obj.hitseqs('clusters',true);
      for i=1:length(naseq)
        if args.seqonly
          ylbl{i}=sprintf('%d',naseq(i));
        else
        nm=obj.seqs.getname(naseq(i));
        if isempty(nm)
          ylbl{i}=sprintf('   %9d',naseq(i));
        else
          ylbl{i}=sprintf('   %-9.9s',nm);
        end
        if ismember(naseq(i),[obj.seqs.hits.naseq])
          ylbl{i}(1)='H';
        elseif ismember(naseq(i),hitseqs)
          ylbl{i}(1)='h';
        end
        if obj.seqs.isbistable(naseq(i))
          ylbl{i}(2)='B';
        end
        if i==1 || cluster(i)~=cluster(i-1)
          ylbl{i}=sprintf('%s %4d',ylbl{i},cluster(i));
        else
          ylbl{i}=sprintf('%s     ',ylbl{i});
        end
        end
      end
      ylbl{end+1}='Amplicons';
      ylbl{end+1}='MinDetect';
      set(get(gca,'YAxis'),'FontSize',8,'FontName','Andale Mono');
      set(gca,'YTickLabel',ylbl);
      xlabel('Round');
      ylabel('Seq');
      title('Fraction by Round');
      
      if ~isempty(args.trpexpts)
        % Todo pcolor plot of cleavage

        setfig([ti,'-clv']);clf;
        clv=nan(size(frac));
        for i=1:length(rounds)
          if isKey(args.trpexpts,rounds(i))
            a=obj.trpexpts(args.trpexpts(rounds(i))).analysis;
            for j=1:length(naseq)
              ind=a.naseq==naseq(j);
              if sum(ind)==1 && sum(a.cnt(ind,2:3))>=args.minclvcnt
                clv(j,i)=a.cleavage(ind);
              end
            end
          end
        end
        data=clv;
        data(end+1,:)=nan; data(:,end+1)=nan;
        pcolor(data);
        h2=gca;
        linkaxes([h1,h2]);
        caxis([0.1,0.9]);
        h=colorbar;
        set(gca,'XTick',(1:length(rounds))+0.5);
        set(gca,'TickLabelInterpreter','None');
        set(gca,'XTickLabel',arrayfun(@(z) sprintf('%d',z),rounds,'Unif',false));

        set(get(gca,'YAxis'),'FontSize',8,'FontName','Andale Mono');
        set(gca,'YTick',(1:length(naseq))+0.5);
        set(gca,'YTickLabel',ylbl(1:length(naseq)));
        xlabel('Round');
        ylabel('Seq');
        title('Cleavage by Round');
      else
        clv=[];
      end
      
      x=struct('frac',frac,'fracerr',fracerr,'cnt',cnt,'ampliconcnt',ampliconcnt,'total',total,'rounds',rounds,'naseq',naseq,'cleavage',clv,'cluster',cluster);
    end
      
    function [r_naseq,r_seqstr]=similarseqs(obj,naseq,varargin) 
      defaults=struct('maxerrors',1,'print',true,'mincnt',3);
      args=processargs(defaults,varargin);
      sel=ismember(obj.seqs.naseq,obj.naseq(obj.cnt>=args.mincnt));
      [r_naseq,r_seqstr]=obj.seqs.similarseqs(naseq,'maxerrors',args.maxerrors,'print',args.print,'sel',sel);
    end

    function x=adjustfixed(obj,subrun)
    % Show steps to adjust fixed mix so that all labelled sequences, that are not invalid, and are not 'notmade' are included and equalized
    % Use 'subrun' as the list of subrun with the current counts
      conc=100;
      NGSDatabase.open(obj.run);
      [naseqs,names]=mysql('select naseq,name from ngs.naseqnames where naseq not in (select naseq from ngs.naseqtags where tagtype in (4,5,8,9))');
      cnts=zeros(size(naseqs));
      id=repmat('',length(cnts),1);
      for i=1:length(naseqs)
        c=sum(obj.cnt(ismember(obj.subrun,subrun) & obj.naseq==naseqs(i)));
        if ~isempty(c)
          cnts(i)=c;
        end
        tagsel=[obj.seqs.tags.naseq]==naseqs(i) & strcmp({obj.seqs.tags.name},'id');
        if sum(tagsel)==1
          id{i}=obj.seqs.tags(tagsel).value;
        end
      end
      totalcnt=sum(obj.cnt(ismember(obj.subrun,subrun)));
      add=zeros(size(naseqs));
      target=median(cnts(cnts>10));
      target=median(cnts(cnts>=target/4));
      target=median(cnts(cnts>=target/4));
      for i=1:length(naseqs)
        if cnts(i)<target/2
          add(i)=target-cnts(i);
        end
      end
      [~,ord]=sort(add-cnts/1000);
      fprintf('%9s %20s %20s %5s %5s\n','naseq','name','id','cnt','add');
      for ii=1:length(naseqs)
        i=ord(ii);
        fprintf('%9d %20.20s %20.20s %5d %5d\n', naseqs(i), names{i}, id{i}, cnts(i),add(i));
      end
      fprintf('%9s %20.20s %20s%6d\n', '','Other','',totalcnt-sum(cnts));
      matchcnt=sum(cnts);
      newtotal=matchcnt+sum(add);
      fprintf('%9s %20.20s %20s %5.1f%%%5.1f%%\n', '','Total','',matchcnt/newtotal*100,sum(add)/newtotal*100);
      fprintf('After additions, total reads/median over the %d desired switches = %.1f\n', length(cnts), newtotal/median(add+cnts));
      x=table();
      x.naseq=int32([naseqs(ord);nan]);
      x.name=[names(ord);'other'];
      x.id={id{ord},''}';
      x.left=arrayfun(@(z) '',x.naseq,'UniformOutput',false);
      x.right=x.left;
      for i=1:length(x.id)
        if isempty(x.id{i})
          continue;
        end
        idsplit1=strsplit(x.id{i},'->');
        idsplit2=strsplit(idsplit1{1},'+');
        if length(idsplit1)==2
          x.id{i}=idsplit1{2};
        elseif length(idsplit2)==2
          x.id{i}='';
        end
        if length(idsplit2)==2
          x.left{i}=idsplit2{1};
          x.right{i}=idsplit2{2};
        end
      end
      x.cnt=[cnts(ord);totalcnt-sum(cnts)];
      x.relative=round(x.cnt/target,2);
      x.add=[round(add(ord)/target,2);0];
      ss=obj.seqs.naseqformat(x.naseq,'squeeze',true);
      x.seq=ss;
      writetable(x,'fixedadj.csv');
    end
    
    function tree=switchphylo(obj,target)
    % Create phylogenetic tree of switches for each target
      x=obj.seqs.switchstats();
      if nargin<2
        targets=unique(cellfun(@(z) z(1:min(end,4)),x.name,'UniformOutput',false));
        for i=1:length(targets)
          obj.switchphylo(targets{i});
        end
        return;
      end
      
      ignoreseqs=[obj.seqs.tags(ismember({obj.seqs.tags.name},{'invalid'})).naseq];
      sel=strncmp(x.name,target,length(target)) & ~ismember(x.naseq,ignoreseqs);
      if sum(sel)<2
        fprintf('Not enough switch names match %s-* (found %d)\n',target,sum(sel));
        return;
      end
      seqs=x.seq(sel);
      meanlen=mean(cellfun(@(z) sum(z~= ' '), seqs));
      dist=seqpdist(seqs,'alphabet','nt','method','p-distance')*meanlen;
      tree=seqlinkage(dist,'average',x.name(sel));
      % setfig(target);clf;
      plot(tree);
      title(target);
    end
    
    function x=compareseqs(obj,naseqs,varargin)
    % Compare given sequences across all trpexpts
      defaults=struct('mincnt',10);
      args=processargs(defaults,varargin);
      x=struct();
      x.trpexpts=1:length(obj.trpexpts);
      x.name={obj.trpexpts(x.trpexpts).descr};
      x.cnt=nan(length(x.trpexpts),length(naseqs),3);
      x.cleavage=nan(length(x.trpexpts),length(naseqs));
      x.keep=false(size(x.trpexpts));
      for i=x.trpexpts
        for k=1:length(naseqs)
          t=obj.trpexpts(i).analysis;
          sel=t.naseq==naseqs(k);
          if sum(sel)==1
            x.cnt(i,k,:)=t.cnt(sel,:);
            x.cleavage(i,k)=t.cleavage(sel);
          end
        end
        if sum(sum(x.cnt(i,:,2:3),3)>=args.mincnt)>=2
          keep(i)=true;
        end
      end
      setfig('compareseqs');clf;
      bar(x.cleavage(keep,:));
      set(gca,'XTick',1:sum(keep));
      set(gca,'XTickLabel',x.name(keep));
      set(gca,'XTickLabelRotation',45);
      legend(arrayfun(@(z) obj.seqs.getname(z),naseqs,'UniformOutput',false));
      ylabel('Cleavage');
      if length(naseqs)==2
        setfig('compareseqs-scatter');clf;
        loglog(-(1-x.cleavage(keep,1)),-(1-x.cleavage(keep,2)),'o');
        hold on;
        loglog(-(1-[.01,.99]),-(1-[.01,.99]),':');
        fk=find(keep);
        for i=1:length(fk)
          text(-(1-x.cleavage(fk(i),1)),-(1-x.cleavage(fk(i),2)),['  ',x.name(fk(i))],'HorizontalAlignment','left','VerticalAlignment','middle');
        end
        cleaveticks(1,1,'metric','ratio_of_uncleavage');
        xlabel(sprintf('Cleavage of %s',obj.seqs.getname(naseqs(1))));
        ylabel(sprintf('Cleavage of %s',obj.seqs.getname(naseqs(2))));
      end
    end

    function x=gelcompare(obj,expt,varargin)
    % Compare gel cleavage with values based on the given expt
      defaults=struct('usegain',false);
      args=processargs(defaults,varargin);

      obj.trpanalyze();
      geltags=strcmp({obj.seqs.tags.name},'gel') & ~strcmp({obj.seqs.tags.value},'TBD');
      fprintf('Have %d sequences with gel cleavage data\n', sum(geltags));
      naseqs=[obj.seqs.tags(geltags).naseq]';
      gelcleavage=str2double({obj.seqs.tags(geltags).value})';
      ngscleavage=nan(size(gelcleavage));
      retention=nan(size(gelcleavage));
      ugain=nan(size(gelcleavage));
      cnt=nan(length(gelcleavage),3);
      t=obj.trpexpts(expt).analysis;
      name=cell(length(naseqs),1);
      for i=1:length(naseqs)
        ind=find(naseqs(i)==t.naseq);
        if isempty(ind)
          fprintf('No NGS data for naseq %d\n', naseqs(i));
          continue;
        end
        ngscleavage(i)=t.cleavage(ind);
        ugain(i)=t.cnt(ind,3)./t.cnt(ind,1);
        cnt(i,:)=t.cnt(ind,:);
        name{i}=obj.seqs.getname(naseqs(i));
        retention(i)=t.retention(ind);
      end
      r=corrcoef(gelcleavage,ngscleavage*100);
      rsqd=r(1,2)^2;
      setfig(sprintf('gelcompare %s',obj.trpexpts(expt).descr));clf;
      if args.usegain
        semilogy(ugain,gelcleavage,'o');
        xlabel('CleaveSeq Uncleaved/Input Ratio');
      else
        plot(ngscleavage*100,gelcleavage,'o');
        xlabel('CleaveSeq Cleavage');
      end
      %axis([0,100,0,100]);
      %axis equal
      text(5,90,sprintf('R^2=%.2f',rsqd));
      ylabel('Gel-based Cleavage');
      x=table();
      x.naseq=int32(naseqs);
      x.name=name;
      x.gel=gelcleavage;
      x.ngs=ngscleavage*100;
      x.ugain=ugain;
      x.incnt=int32(cnt(:,1));
      x.ccnt=int32(cnt(:,2));
      x.ucnt=int32(cnt(:,3));
      x.retention=retention;
      x=sortrows(x,'gel');
    end
  
    function x=arcompare(obj,expts,varargin)
    % Compare ar cleavage with values based on the given 2 expts
      defaults=struct('usegain',false,'metric','ratio_of_uncleavage','target',[]);
      args=processargs(defaults,varargin);

      if length(expts)~=2
        error('Usage: arcompare([trp1,trp2],options,...)');
      end
      
      obj.trpanalyze();
      artags=find(strcmp({obj.seqs.tags.name},'yeast') & ~strcmp({obj.seqs.tags.value},'TBD'));
      fprintf('Have %d sequences with yeast AR data\n', length(artags));
      naseqs=[]; ar=[];
      for i=1:length(artags)
        tag=obj.seqs.tags(artags(i));
        seqname=obj.seqs.getname(tag.naseq);
        if isempty(args.target) || strncmp(seqname,args.target,length(args.target))
          naseqs(end+1)=tag.naseq;
          val=tag.value;
          % Parse string of form 'ar=2.2;1.89'
          start=find(val==';' | val=='=',1,'last');
          ar(end+1)=str2double(val(start+1:end));
        end
      end
      ngsswitching=nan(size(ar));
      ta=obj.trpexpts(expts(1)).analysis;
      ta2=obj.trpexpts(expts(2)).analysis;
      name=cell(length(naseqs),1);
      for i=1:length(naseqs)
        name{i}=obj.seqs.getname(naseqs(i));
        ind1=ta.naseq==naseqs(i);
        ind2=ta2.naseq==naseqs(i);
        if sum(ind1)==1 && sum(ind2)==1
          ngsswitching(i)=exp(obj.bootfold(ta2.ratioboot(ind2,:),ta.ratioboot(ind1,:),'metric',args.metric));
        end
      end
      setfig(sprintf('arcompare %s/%s',obj.trpexpts(expts(1)).descr,obj.trpexpts(expts(2)).descr));clf;
      loglog(ngsswitching,ar,'o');
      xlabel(args.metric,'Interpreter','none');
      ylabel('Activity Ratio');
      logticks(0,1);
      
      x=table();
      x.naseq=int32(naseqs)';
      x.name=name;
      x.ar=ar';
      x.ngs=ngsswitching';
      x=sortrows(x,'ar');
    end

    function x=refcheck(obj,subruns,varargin)
    % Compare reference counts in the given subruns
      defaults=struct('mincnt',100,'dopc',false,'alpha',0.1,'subplots',false);
      args=processargs(defaults,varargin);

      if nargin<2 || isempty(subruns)
        subruns=obj.subruns.subrun;
      end
      sel=strcmp({obj.seqs.tags.name},'ref');
      prefix={'A','W','Z','T7W','T7Z'};
      naseqs=[];names={};
      for i=1:length(prefix)
        ns=[obj.seqs.tags(sel & strcmp({obj.seqs.tags.value},prefix{i})).naseq];
        nm=arrayfun(@(z) [prefix{i},'_',obj.seqs.getname(z)],ns,'UniformOutput',false);
        naseqs=[naseqs,ns];
        names={names{:},nm{:}};
      end
      fprintf('refs={%s}\n',strjoin(names,','));
      cnt=zeros(length(naseqs),length(subruns));
      for i=1:length(subruns)
        for j=1:length(naseqs)
          ind=obj.subrun==subruns(i) & obj.naseq==naseqs(j);
          if any(ind)
            cnt(j,i)=sum(obj.cnt(ind));
          end
        end
      end
      data=cnt;
      sel=sum(cnt)>args.mincnt;
      fprintf('Keeping %d/%d subruns with mincnt>=%d\n', sum(sel),length(sel),args.mincnt);
      cnt=cnt(:,sel);
      subruns=subruns(sel);
      codes=arrayfun(@(z) obj.subruns.codes{z==obj.subruns.subrun}, subruns,'UniformOutput',false);
      ci=nan(size(cnt,1),size(cnt,2),2);
      data=nan(size(cnt));
      sumcnt=sum(cnt,1);
      for i=1:size(cnt,2)
        for j=1:size(cnt,1)
          [data(j,i),ci(j,i,:)]=binofit(cnt(j,i),sumcnt(i),args.alpha);
        end
      end
      % Only keep active prefixes
      setfig('refcheck');clf;
      keep=false(size(prefix));
      % Only keep prefixes in use
      for i=1:length(prefix)
        sel=strncmp(codes,prefix{i},length(prefix{i})) & (sum(data>0)>0)';
        keep(i)=any(sel);
      end
      fprintf('Keeping %d/%d prefixes\n', sum(keep), length(keep));
      prefix=prefix(keep);
      if args.subplots
        firsti=1;
      else
        firsti=length(prefix)+1;
      end
      
      for i=firsti:length(prefix)+1
        if args.subplots
          subplot(length(prefix)+1,1,i);
        end
        if i>length(prefix)
          sel=(sum(data>0)>0)';
        else
          sel=strncmp(codes,prefix{i},1) & (sum(data>0)>0)';
        end
       
        if ~any(sel)
          continue;
        end
        dd=data(:,sel);
        de=ci(:,sel,:);
        sel2=any(dd>0.01,2)';
        dd=dd(sel2,:);
        de=de(sel2,:,:);
        srnames=arrayfun(@(z,c) sprintf('%s N=%d',obj.getsubrunname(z),c), subruns(sel),sumcnt(sel)','UniformOutput',false);
        if args.dopc
          dd(end+1,:)=nan;dd(:,end+1)=nan;
          pcolor(dd);
          set(gca,'YTick',(1:size(dd,1))+0.5);
          set(gca,'YTickLabel',names(sel2));
          set(gca,'XTick',(1:length(srnames))+0.5);
          set(gca,'XTickLabel',srnames);
          set(gca,'XTickLabelRotation',30);
          colorbar();
        else
          for k=1:size(dd,2)
            errorbar((1:size(dd,1))+(k/size(dd,2)-0.5)/4,dd(:,k),dd(:,k)-de(:,k,1),de(:,k,2)-dd(:,k));
            hold on;
          end
          %set(gca,'YScale','log');
          set(gca,'XTick',(1:size(dd,1)));
          set(gca,'XTickLabel',names(sel2));
          set(gca,'XTickLabelRotation',45);
          if length(srnames)<=20
            legend(srnames,'location','best','Interpreter','none');
          end
        end
        set(gca,'TickLabelInterpreter','None');
        if i>length(prefix)
          title('All');
        else
          title(sprintf('Prefix=%s',prefix{i}));
        end
      end
      if nargout>0
        sel=sum(cnt,2)>0;
        x=struct('names',{names(sel)},'cnt',cnt(sel,:),'srnames',{srnames},'subruns',subruns);
      end
    end
    
    function refleakage(obj,varargin)
      defaults=struct('subruns',[],'mincnt',100,'thresh',0.001);
      args=processargs(defaults,varargin);

      if isempty(args.subruns)
        args.subruns=obj.subruns.subrun;
      end

      fprintf('Mincnt: %d, Thresh: %.1f%%\n', args.mincnt, args.thresh*100);
      
      % Find all refs
      allrefs=[obj.seqs.tags(strcmp({obj.seqs.tags.name},'ref')).naseq];
      [tgts,ord]=sort({obj.seqs.tags(strcmp({obj.seqs.tags.name},'ref')).value});
      allrefs=allrefs(ord);
      
      codes=unique(obj.subruns.codes);
      cnts=nan(length(allrefs),length(codes));
      for i=1:length(codes)
        srs=intersect(obj.subruns.subrun(strcmp(obj.subruns.codes,codes{i})),args.subruns);
        for j=1:length(allrefs)
          sel=ismember(obj.subrun,srs) & obj.naseq==allrefs(j);
          cnts(j,i)=sum(obj.cnt(sel));
        end
      end
      sel=sum(cnts,2)>=args.mincnt;
      sel2=sum(cnts,1)>=args.mincnt;
      allrefs=allrefs(sel);
      tgts=tgts(sel);
      cnts=cnts(sel,sel2);
      codes=codes(sel2);
      
      % Merge all the X variants
      for i=1:length(codes)
        codes{i}=regexprep(codes{i},'/X.*','/X');
        codes{i}=regexprep(codes{i},'NNNN_W/','W/');
        codes{i}=regexprep(codes{i},'GCGT_Z/','Z/');
      end
      ucodes=unique(codes);
      for i=1:length(ucodes)
        m=strcmp(codes,ucodes{i});
        ucnts(:,i)=sum(cnts(:,m),2);
      end
      cnts=ucnts;
      codes=ucodes;
      
      dcnt=[];
      for i=1:length(allrefs)
        desired=strncmp(codes,tgts{i},length(tgts{i}));
        dcnt(i)=nansum(cnts(i,desired));
      end
      
      for i=1:length(allrefs)
        desired=strncmp(codes,tgts{i},length(tgts{i}));
        fprintf('%9d %-25.25s tgt=%3s %4.2fx ',allrefs(i), obj.seqs.getname(allrefs(i)), tgts{i}, dcnt(i)/nanmedian(dcnt(dcnt>args.mincnt)));
        if any(desired)
          fprintf('%8.8s:%-7d ', codes{desired}, dcnt(i));
        else
          fprintf('%8.8s:%-7d ', '',dcnt(i));
        end
        sel=find(cnts(i,:)'>args.thresh*dcnt(i) & ~desired);
        for kk=1:length(sel)
          k=sel(kk);
          fprintf('%8.8s:%d=%4.1f%% ',codes{k},cnts(i,k),cnts(i,k)/dcnt(i)*100);
        end
        fprintf('\n');
      end
      
    end
    
    function x=refcompare2(obj,varargin)
    % Compare refs across subrun descs (ie. summed over all codes)
      defaults=struct('mincnt',100);
      args=processargs(defaults,varargin);
    
      naseqs=[obj.seqs.tags(strcmp({obj.seqs.tags.name},'ref')).naseq];
      udesc=unique(obj.subruns.desc);
      for j=1:length(naseqs)
        total(j)=sum(obj.cnt(obj.naseq==naseqs(j)));
      end
      naseqs=naseqs(total>=args.mincnt);
      
      cnt=zeros(length(udesc),length(naseqs));
      for i=1:length(udesc)
        subruns=obj.subruns.subrun(strcmp(obj.subruns.desc,udesc{i}));
        for j=1:length(naseqs)
          cnt(i,j)=sum(obj.cnt(obj.naseq==naseqs(j) & ismember(obj.subrun,subruns)));
        end
      end
      keepdesc=max(cnt,[],2)>=args.mincnt;
      udesc=udesc(keepdesc);
      keepnaseq=max(cnt,[],1)>=args.mincnt;
      naseqs=naseqs(keepnaseq);
      cnt=cnt(keepdesc,keepnaseq);
      % Convert to relative
      rel=zeros(size(cnt));
      for i=1:size(cnt,1)
        rel(i,:)=cnt(i,:)/sum(cnt(i,:));
      end
      setfig('refcompare2');clf;
      plot(100*rel','o-');
      ylabel('Fraction (%)');
      legend(udesc,'Location','best');
      names=arrayfun(@(z) obj.seqs.getname(z), naseqs,'Unif',false);
      set(gca,'XTick',1:length(names));
      set(gca,'XTickLabel',names);
      set(gca,'XTickLabelRotation',45);
      set(gca,'TickLabelInterpreter','None');
      
      x=array2table(round(rel*100),'VariableNames',names,'RowNames',udesc);
    end
    
    
    function [c,ratio,subruns,refs,prefix]=refcompare(obj, varargin)
    % Compare given refs across subruns
      defaults=struct('mincnt',100,'doplots',true,'subruns',[]);
      args=processargs(defaults,varargin);
      if isempty(args.subruns)
        subruns=obj.subruns.subrun;
      else
        subruns=args.subruns;
      end
      
      if ~isfield(obj.subruns,'bclen')
        obj.dbgetcodelengths();
      end
        
      refs={'R6R30','R7R4','R25R6','R35R6','R60R5'};
      prefix={'A','W','Z','T7W','T7Z'};
      for i=1:length(refs)
        hasref=cellfun(@(z) ~isempty(strfind(z,refs{i})),{obj.seqs.seqnames.name});
        n1=[obj.seqs.seqnames(hasref).naseq];
        for j=1:length(prefix)
          n2=[obj.seqs.tags(strcmp({obj.seqs.tags.name},'ref') & strcmp({obj.seqs.tags.value},prefix{j})).naseq];
          n=intersect(n1,n2);
          rds=[];
          for k=1:length(n)
            rds(k)=sum(obj.cnt(obj.naseq==n(k)));
          end
          fprintf('Have %d references for %s-%s with [%s] rds\n',length(n),prefix{j},refs{i},sprintf('%d ',rds));
          if length(n)~=1
            if length(n)>1
              n=n(1);
            else
              n=nan;%assert(false);
            end
          end
          naseqs(i,j)=n;
        end
      end
      
      ratio=nan(size(naseqs,1),size(naseqs,2),length(subruns));
      for i=1:size(naseqs,1)
        for j=1:size(naseqs,2)
          for k=1:length(subruns)
            c(i,j,k)=sum(obj.cnt(obj.naseq==naseqs(i,j) & obj.subrun==subruns(k)));
          end
          ratio(i,j,:)=c(i,j,:)./c(1,j,:);
          lowcnt=c(1,j,:)<args.mincnt;
          ratio(i,j,lowcnt)=nan;
          ratio(i,j,:)=ratio(i,j,:)./nanmedian(ratio(i,j,:));  % Normalize around 1.0
        end
      end
      if ~args.doplots
        return;
      end

      for i=1:length(prefix)
        t=ratio(2:end,i,:);
        keep(i)=any(isfinite(t(:)));
      end
      prefix=prefix(keep);
      ratio=ratio(:,keep,:);
      naseqs=naseqs(:,keep);
      
      setfig('refcompare');clf;
      desc={};
      for i=1:length(subruns)
        sel=obj.subruns.subrun==subruns(i);
        desc{i}=sprintf('%s-%s L=%d',obj.subruns.desc{sel},obj.subruns.codes{sel},obj.subruns.bclen(sel));
      end
        
      for j=1:length(prefix)
        subplot(length(prefix),1,j);
        r=squeeze(ratio(2:end,j,:));
        sel=find(any(isfinite(r)));
        if ~any(sel)
          continue;
        end
        [~,ord]=sort(obj.subruns.bclen(sel));   % Sort by barcode length
        sel=sel(ord);
        r=r(:,sel);
        r(end+1,:)=nan;
        r(:,end+1)=nan;
        pcolor(log10(r));
        set(gca,'XTick',(1:length(sel))+0.5);
        set(gca,'XTickLabel',desc(sel));
        set(gca,'XTickLabelRotation',45);
        set(gca,'YTick',(1:(length(refs)-1))+0.5);
        set(gca,'YTickLabel',refs(2:end));
        set(gca,'TickLabelInterpreter','None');
        ylabel(sprintf('Relative to %s',refs{1}));
        title(sprintf('Prefix=%c mincnt=%d',prefix{j},args.mincnt));
        colorbar();
        caxis(log10([nanmin(ratio(:)),nanmax(ratio(:))]));
      end
      
      setfig('refdotplot');clf;
      for i=1:length(prefix)
        subplot(2,ceil(length(prefix)/2),i);
        cc=squeeze(c(:,i,:));
        h=[];
        for j=2:length(refs)
          sel=cc(1,:)>=args.mincnt;
          if ~any(sel)
            continue;
          end
          h(end+1)=cdfplot(cc(j,sel)./cc(1,sel));
          set(gca,'XScale','log');
          hold on;
          hupper=cdfplot((cc(j,sel)+sqrt(cc(j,sel)))./(cc(1,sel)-sqrt(cc(1,sel))));
          set(hupper,'Color',get(h(end),'Color'));
          set(hupper,'LineStyle',':');
          hlower=cdfplot((cc(j,sel)-sqrt(cc(j,sel)))./(cc(1,sel)+sqrt(cc(1,sel))));
          set(hlower,'Color',get(h(end),'Color'));
          set(hlower,'LineStyle',':');
          grid off;
        end
        ylabel('CDF');
        xlabel('Ref Ratio');
        legend(h,refs(2:end));
        title(sprintf('%s-Prefix',prefix{i}));
      end
    end
    
    function data=vitrovivo(obj,expt,vivo,varargin)
    % Plot in vitro vs in vivo data
      defaults=struct('mincnt',50,'ivexpt',5,'minvivocnt',30,'sel',[],'medlen',20);
      args=processargs(defaults,varargin);

      % Locate all the finite in vivo data points
      ivsel=isfinite(vivo.invivo(args.ivexpt).ratio) & sum(vivo.invivo(args.ivexpt).hist,2)>=args.minvivocnt;
      ivseq={vivo.tmpl.elib(ivsel).seq};
      ratio=vivo.invivo(args.ivexpt).ratio(ivsel);
      fprintf('Have %d in vivo measurements\n', sum(ivsel));
      setfig('vitrovivo');clf;
      t=obj.trpexpts(expt).analysis;
      sel=sum(t.cnt(:,2:3),2)>=args.mincnt & isfinite(t.cleavage);
      if ~isempty(args.sel)
        sel=sel&ismember(t.naseq,obj.naseq(args.sel));
      end
      fprintf('Have %d in vitro measurements for %s\n', sum(sel),obj.trpexpts(expt).descr);
      cleavage=t.cleavage(sel);
      naseq=t.naseq(sel);
      vitroseqs=obj.seqs.getseq(naseq);
      [c,ia,ib]=intersect(ivseq,vitroseqs);
      fprintf('Have %d measurements with both\n', length(c));
      cleavage=cleavage(ib);
      naseq=naseq(ib);
      ratio=ratio(ia);
      loglog(-(1-cleavage),ratio,'.');
      hold on;
      [sc,ord]=sort(cleavage);
      sr=ratio(ord);
      plot(-(1-sc),medfilt1(sr,args.medlen,'truncate'),':g','LineWidth',3);
      %logticks(0,1);
      cleaveticks(1,0,'metric','ratio_of_uncleavage','pct',true);
      markseq=[54029,202280];
      for i=1:length(markseq)
        sel=find(naseq==markseq(i));
        if sum(sel)==0
          continue;
        end
        fprintf('%d %s cleavage=%f, ratio=%f\n', naseq(sel), obj.seqs.getname(naseq(sel)), cleavage(sel), ratio(sel));
        plot(cleavage(sel)*100,ratio(sel),'or');
      end
      xlabel('CleaveSeq cleavage');
      ylabel('GFP/mCherry in vivo');
      logticks(0,1);
      data=table();
      data.naseq=naseq;
      data.cleavage=cleavage;
      data.ratio=ratio;
    end

    function revertcompare(obj,revset,varargin)
      defaults=struct('mincnt',100);
      args=processargs(defaults,varargin);
      setfig('reversion');clf;
      lbls={};
      for i=1:length(revset)
        r=revset{i};
        concs=r{4};
        for e=1:length(r{3})
          t=obj.trpexpts(r{3}(e)).analysis;
          sel1=t.naseq==r{1}; 
          sel2=t.naseq==r{2};
          if ~any(sel1) || ~any(sel2)
            fprintf('No data for revset %d, expt %d\n', i, e);
            continue;
          end
          cnt=sum(t.cnt(sel1|sel2,2:3),2);
          if min(cnt)<args.mincnt
            fprintf('Skipping revset %d, expt %d -- only have %d data points\n', i, e, min(cnt));
            continue;
          end
          r1=t.ratio(sel1);
          r2=t.ratio(sel2);
          if e==1
            col='r';
          elseif e==length(r{3})
            col='g';
          else
            col='k';
          end
          semilogy(i*2+[-1,0],[r1,r2],[col,'o-']);
          hold on;
          % Display a legend of concs
          text(i*2,r2,sprintf('  %.0f\\muM',concs(e)),'HorizontalAlignment','Left','VerticalAlignment','middle');
        end
        lbls{end+1}=obj.seqs.getname(r{1});
        lbls{end+1}=obj.seqs.getname(r{2});
      end
      cleaveticks(0,1);
      % Add dividing lines
      c=axis; c=c+[-0.5 0.9 0 0]; axis(c);
      for i=2:2:c(2)-1
        plot((i+0.82)*[1,1],c(3:4),':k');
      end
      set(gca,'XTick',1:length(lbls));
      set(gca,'XTickLabel',lbls);
      set(gca,'XTickLabelRotation',45);
      ylabel('Cleavage');
    end
  
    % Compute the Kullbakc-Liebler divergence between the distributions of sequences in subruns s1 and
    % Can compare across NGS runs by specifying 'other'
    function d=kldivergence(obj,s1,s2,varargin)
      defaults=struct('minlength',70,'other',[]);
      args=processargs(defaults,varargin);

      if isempty(args.other)
        args.other=obj;
      end
      
      sel1=obj.subrun==s1 & obj.length>=args.minlength & ~ismember(obj.naseq,[obj.seqs.tags.naseq]);
      naseqs1=obj.naseq(sel1);
      cnt1=obj.cnt(sel1);
      sel2=args.other.subrun==s2 & args.other.length>=args.minlength;
      naseqs2=args.other.naseq(sel2);
      cnt2=args.other.cnt(sel2);
      allseqs=union(naseqs1,naseqs2);
      p=zeros(size(allseqs));
      q=p;
      [~,ia,ib]=intersect(naseqs1,allseqs);
      p(ib)=cnt1(ia);
      [~,ia,ib]=intersect(naseqs2,allseqs);
      q(ib)=cnt2(ia);
      k=.01;
      p=(p+k)/(sum(cnt1)+length(allseqs)*k);
      q=(q+k)/(sum(cnt2)+length(allseqs)*k);
      d=sum(q.*log(q./p));
    end
  
    function [d,xp,lbl]=klplot(obj,objs,subruns,xpos)
      assert(iscell(objs));
      assert(iscell(subruns));
      assert(iscell(xpos));
      
      d=[];
      lbl={};
      xp=[];
      ipos=0;
      for i1=1:length(subruns)
        for i=1:length(subruns{i1})
          ipos=ipos+1;
          lbl{ipos}=objs{i1}.subruns.desc{objs{i1}.subruns.subrun==subruns{i1}(i)};
          xp(ipos)=xpos{i1}(i);
          d(ipos,:)=nan;
          jpos=0;
          for j1=1:length(subruns)
            for j=1:length(subruns{j1})
              jpos=jpos+1;
              if xpos{j1}(j)>xp(ipos) && (xpos{j1}(j)<xp(ipos)+30 || ismember(xp(ipos),[2,57,63,70]))
                d(ipos,jpos)=objs{i1}.kldivergence(subruns{i1}(i),subruns{j1}(j),'other',objs{j1});
              else
                d(ipos,jpos)=nan;
              end
            end
          end
        end
      end
      % Reorder by xpos
      [xp,ord]=sort(xp);
      d=d(ord,ord);
      lbl=lbl(ord);

      d(end+1,:)=nan;
      d(:,end+1)=nan;
      setfig('klplot');clf;
      pcolor(d);
      set(gca,'XTick',(1:length(lbl))+0.5);
      set(gca,'XTickLabel',lbl);
      set(gca,'XTickLabelRotation',45);
      set(gca,'YTick',(1:length(lbl))+0.5);
      set(gca,'YTickLabel',lbl);
      title('KL Divergence');
      colorbar;
    end
    
    function [frac,desc]=trackseqs(obj,naseqs,subruns,varargin)
    % Track abundance of given naseqs over given subruns
      defaults=struct('doplot',true,'otherfrac',[],'otherdesc',[],'xpos',[],'highlight',[]);
      args=processargs(defaults,varargin);

      ti=obj.run;
      frac=nan(length(naseqs),length(subruns));
      desc={};
      totalcnt=[];
      for i=1:length(subruns)
        srsel=find(obj.subruns.subrun==subruns(i));
        if isempty(srsel)
          error('Subrun %d does not exist',subruns(i));
        end
        rds=nan(length(naseqs),1);
        for j=1:length(naseqs)
          rds(j)=sum(obj.cnt(obj.naseq==naseqs(j) & obj.subrun==subruns(i)));
        end
        totalcnt(i)=sum(obj.cnt(obj.subrun==subruns(i)));
        frac(:,i)=(rds)/(totalcnt(i));
        if ~isempty(obj.trpexpts) && ismember(subruns(i),union([obj.trpexpts.cleavedsr],[obj.trpexpts.uncleavedsr]))
          desc{i}=['TRP(',obj.subruns.desc{srsel},')->',obj.subruns.codes{srsel}(1)];
        else
          desc{i}=obj.subruns.desc{srsel};
        end
      end
      frac(end+1,:)=1./(totalcnt);

      if ~isempty(args.otherfrac)
        frac=[args.otherfrac,frac];
        desc=[args.otherdesc,desc];
      end
      
      if args.doplot
        setfig([ti,'-trackseqs']);clf;
        if isempty(args.xpos)
          xp=1:length(desc);
        else
          xp=args.xpos+(1:length(args.xpos))/length(args.xpos)/10;   % Make monotonic and jittered
        end
        semilogy(xp,frac(1:end-1,:)','o-');
        hold on;
        semilogy(xp,frac(end,:)',':k');
        ylabel('Fraction');
        set(gca,'XTick',xp);
        set(gca,'XTickLabel',desc);
        set(gca,'XTickLabelRotation',45);
        set(gca,'TickLabelInterpreter','None');
        lbls=arrayfun(@(z) num2str(z), naseqs,'UniformOutput',false);
        for i=1:length(naseqs)
          nm=obj.seqs.getlabels(naseqs(i),'includetags',false,'includemuts',false);
          if ~isempty(nm)
            lbls{i}=[lbls{i},' ',nm];
          end
        end
        lbls{end+1}='floor';
        legend(lbls);
        if length(naseqs)>20
          % Hide legend (but displaynames still set)
          legend('off');
        end
        % Highlights
        for i=1:length(args.highlight)
          sel=find(naseqs==args.highlight(i));
          h=get(gca,'Children');
          set(h(length(h)+1-sel),'LineWidth',3);
          %plot(args.xpos,frac(sel,:),'o-','LineWidth',5);
          last=find(frac(sel,:)>0,1,'last');
          id=obj.seqs.getname(args.highlight(i));
          if isempty(id)
            id=sprintf('%d',args.highlight(i));
          end
          text(args.xpos(last)+0.5,frac(sel,last),id);
        end
      end
    end
    
    function fixtargets(obj)
    % Was previously saving the target name (abbrev.) and conc in desc and/or targets
    % Now, set 'target' to full name of target, 'targetconc' to concentration in M (or x), 'concunits' to 'M' or 'x'
      if isempty(obj.trpexpts)
        return;
      end
      if isfield(obj.trpexpts,'targetconc')   % FIXME: obj.trpexpts is a class member, so isfield() is always false
        % Already fixed, but make sure we don't have nan's in targetconc
        for i=1:length(obj.trpexpts)
          if isempty([obj.trpexpts(i).targetconc]) || isnan(obj.trpexpts(i).targetconc)
            obj.trpexpts(i).targetconc=[];
            obj.trpexpts(i).concunits='';
          end
        end
        return;
      end
      
        
      obj.trpexpts(1).targetconc=[];
      obj.trpexpts(1).concunits='';
      
      for i=1:length(obj.trpexpts)
        t=obj.trpexpts(i);
        if isempty(t.target) || ~isempty(t.targetconc)
          if isnan(t.targetconc)
            t.targetconc=[];
            obj.trpexpts(i)=t;
          end
          continue;
        end
        at=find(t.target=='@');
        if ~isempty(at)
          t.targetconc=str2double(t.target(at+1:end))*1e-6;
          t.concunits='M';
          t.target=t.target(1:at-1);
        else
          at=find(t.descr=='@');
          if isempty(at)
            t.targetconc=1;
            t.concunits='x';
          else
            t.targetconc=str2double(t.descr(at+1:end))*1e-6;
            t.concunits='M';
          end
        end
        abbrevs=containers.Map();
        abbrevs('Acic')='Aciclovir';
        abbrevs('Loxo')='Loxoribine';
        abbrevs('S-Ret')='S-Reticuline';
        abbrevs('SRet')='S-Reticuline';
        abbrevs('TZea')='Trans-Zeatin';
        abbrevs('Nosc')='Noscapine';
        abbrevs('Gard')='Gardiquimod';
        abbrevs('Imiq')='Imiquimod';
        abbrevs('Theo')='Theophylline';
        if abbrevs.isKey(t.target)
          t.target=abbrevs(t.target);
        end
        if strcmp(t.target,'Theophylline') || strcmp(t.target,'Loxoribine')
          t.targetconc=t.targetconc*1000;   % In mM
        end
        obj.trpexpts(i)=t;
      end
    end
    

    function findamplicons(obj,varargin)
    % Find sequences which a TRP shows output (>=mincnt), but no input
      defaults=struct('mincnt',10);
      args=processargs(defaults,varargin);
      for i=1:length(obj.trpexpts)
        t=obj.trpexpts(i);
        a=t.analysis;
        if sum(a.cnt(:,1))==0
          continue;   % No input subrun
        end
        sel=find(a.cnt(:,1)==0 & (sum(a.cnt(:,2:3),2)>20));
        if length(sel)==0
          continue;
        end
        naseq=a.naseq(sel);
        fprintf('%s:\n',t.descr);
        extra=containers.Map('KeyType','double','ValueType','char');
        for i=1:length(sel)
          extra(a.naseq(sel(i)))=sprintf('%d->%d+%d',a.cnt(sel(i),:));
        end
        obj.listsequences(ismember(obj.naseq,naseq) & ismember(obj.subrun,a.subruns),'extra',extra);
        fprintf('\n');
      end
    end
    
    function processcheck(obj,stages,varargin)
    % Various diagnostics of a TRP run that has data for Input, RT, Ext, and PCR products
    % Use obj.trpexpts for the Ext stage to figure out input and output prefixes
    % stages is cell array of names (==subruns.desc) of stages In,RT,Ext,PCR
    %   - empty strings if stage not acquired
      defaults=struct('ind',[],'class',[],'prefixin','W','prefixout','Z');
      args=processargs(defaults,varargin);
      % ind - only use data based on this selector (logical(length(obj.cnt)))
      if ~isempty(args.class)
        args.ind=ismember(obj.getclass,args.class);
        fprintf('Including only sequences from class(es): ');
        for i=1:length(args.class)
          fprintf('%s ',obj.seqs.seqclasses{args.class(i)}{1});
        end
        fprintf('\n');
      end
      
      if isempty(args.ind)
        args.ind=true(size(obj.naseq));
      end
      sr=nan(length(stages),4);   % Second index is for each of: inprefix,outprefix,t7inprefix,t7outprefix
      for i=1:length(stages)
        trpexpt=find(strcmp({obj.trpexpts.descr},stages{i}),1);
        if ~isempty(trpexpt)
          break;
        end
      end
      if isempty(trpexpt)
        error('Unable to find any TRP expt with name same as any stage (should always include Ext stage)\n');
      end
      
      incode=obj.subruns.codes{obj.subruns.subrun==obj.trpexpts(trpexpt).uncleavedsr};
      outcode=obj.subruns.codes{obj.subruns.subrun==obj.trpexpts(trpexpt).cleavedsr};
      fprintf('Basing process check on trpexpt %d: %s->%s\n', trpexpt, incode, outcode);
      codes={incode,outcode,['T7',incode],['T7',outcode]};
      for i=1:length(stages)
        for j=1:length(codes)
          sel=strcmp(obj.subruns.desc,stages{i}) & strcmp(obj.subruns.codes,codes{j});
          if sum(sel)==1
            sr(i,j)=obj.subruns.subrun(sel);
          elseif sum(sel)==0
            fprintf('Unable to locate %s-%s\n', stages{i}, codes{j});
          elseif sum(sel)>1
            fprintf('Multiple possible subruns for %s-%s\n', stages{i}, codes{j});
          end
        end
      end
      cnt=nan(size(sr));
      refcnt=nan(size(sr));
      refconc=nan(size(sr));
      for i=1:length(sr(:))
        if isnan(sr(i))
          continue;
        end
        norm=obj.subruns.normalization(obj.subruns.subrun==sr(i));
        refcnt(i)=nansum(norm.cnt);
        refconc(i)=nansum(norm.conc)*norm.predil;
        cnt(i)=sum(obj.cnt(obj.subrun==sr(i) & args.ind & ~ismember(obj.naseq,norm.naseq)));
      end
      if nansum(nansum(refconc(:,3:4)))==0
        % No references for T7 prefixes, use those without
        refcnt(:,3:4)=refcnt(:,1:2);
        refconc(:,3:4)=refconc(:,1:2);
      end
      conc=cnt./refcnt.*refconc;
      fprintf('%20.20s %7s  %7s  %7s  %7s\n', '',codes{:});
      for i=1:length(stages)
        fprintf('%-20.20s %7.2f  %7.2f  %7.2f  %7.2f\n', stages{i}, conc(i,:)*1e9);
        if i>1
          fprintf('%-20.20s ',' Gain');
          for k=1:size(conc,2)
            g=conc(i,k)./conc(i-1,k);
            if g>=1000
              fprintf('%7.0fx ', g);
            elseif g>1
              fprintf('%7.2fx ', g);
            else
              fprintf('%7.2f/ ', 1/g);
            end
          end
          fprintf('\n\n');
        end
      end
    end
    
    function downsample(obj,subruns,factor)
    % Downsample the given subruns using the given factor
      assert(factor>1);
      sel=ismember(obj.subrun,subruns);
      cnt=obj.cnt(sel);
      total=sum(cnt);
      fprintf('Before downsampling, these subruns have %d reads over %d entries\n', total, length(cnt));
      newcnt=zeros(size(cnt));
      for i=1:length(cnt)
        newcnt(i)=binornd(cnt(i),1/factor);
      end
      fprintf('After downsampling, these subruns have %d reads\n', sum(newcnt));
      obj.cnt(sel)=newcnt;
      zcnt=obj.cnt==0;
      obj.subrun=obj.subrun(~zcnt);
      obj.naseq=obj.naseq(~zcnt);
      obj.cnt=obj.cnt(~zcnt);
    end
    
    function clusterseqs(obj,naseq,varargin)
    % Cluster given sequences 
    % If naseq is empty, then cluster all with at least minreads in some subrun, and at least mincnt total reads
      defaults=struct('maxdist',5,'mincnt',10,'minreads',3,'minaddcnt',[],'onlycluster',[]);
      args=processargs(defaults,varargin);
      if nargin<2 || isempty(naseq)
        % Find the unique sequences
        useqs=unique(obj.naseq(obj.cnt>=args.minreads),'stable');   % Only naseqs that have minreads in some subrun
      else
        useqs=unique(naseq);
      end
      if ~isempty(obj.seqs.clusters)
        useqs=setdiff(useqs,abs([obj.seqs.clusters.members]));
      end
      if nargin<2 || isempty(naseq)
        fprintf('Have %d unique unclustered seqs that have at last %d reads in a single subrun\n', length(useqs),args.minreads);
      else
        fprintf('Have %d unique unclustered seqs\n', length(useqs));
      end
      ucnt=zeros(length(useqs),1);
      % Compute total reads (over all subruns) for each unique one
      fprintf('Computing total reads...');
      cntsel=ismember(obj.naseq,useqs);
      cntnaseq=obj.naseq(cntsel);
      cntcnt=obj.cnt(cntsel);
      for i=1:length(useqs)
        if mod(i,1000)==0
          fprintf('%d...',i);
        end
        ucnt(i)=sum(cntcnt(cntnaseq==useqs(i)));
      end
      fprintf('done\n');
      obj.seqs.clusterseqs(useqs,ucnt,'maxdist',args.maxdist,'mincnt',args.mincnt,'minaddcnt',args.minaddcnt,'onlycluster',args.onlycluster);
      % And do short ones
      [useqs,ia]=setdiff(useqs,abs([obj.seqs.clusters.members]));
      ucnt=ucnt(ia);
      if length(useqs)>1
        fprintf('Have %d unclustered (short?) sequences to cluster with maxdist=%d\n', length(useqs),min(args.maxdist,2));
        obj.seqs.clusterseqs(useqs,ucnt,'maxdist',min(args.maxdist,2),'mincnt',args.mincnt,'minlength',5,'minaddcnt',args.minaddcnt,'onlycluster',args.onlycluster);
      end
    end
    
    function [naseq,divfrac,total]=dbuniqueseqs(obj,varargin)
    % Determine how many sequences in run were unique and not in any cluster 
    % Use database in case we haven't loaded everything
      defaults=struct('descr','','nsample',100);
      args=processargs(defaults,varargin);
      
      if isempty(args.descr)
        subrunsel='descr is not null and descr != ''PhiX''';
      else
        subrunsel=sprintf('descr=''%s''',args.descr);
      end
      NGSDatabase.open(obj.run);

      ignoreseqs=sprintf('%d,',union(obj.refseqs(1),obj.ampliconseqs(1))); ignoreseqs=ignoreseqs(1:end-1);
      lensel='((LENGTH(s.seq)>=75 AND LENGTH(s.seq<=79)) OR (LENGTH(s.seq)>=105 AND LENGTH(s.seq)<=109))';
      cmd=sprintf('SELECT e.naseq FROM ngsentries e, ngs.naseqs s WHERE s.naseq=e.naseq AND %s AND subrun IN (SELECT subrun FROM subruns where %s) and e.naseq NOT IN (SELECT naseq FROM ngs.clustermembers) and e.naseq NOT IN (%s) group by naseq having sum(cnt)=1',lensel, subrunsel,ignoreseqs);
      NGSDatabase.open(obj.run);
      fprintf('Query: %s\n', cmd);
      naseq=mysql(cmd);
      if isempty(naseq)
        error('No sequences loaded;  query: %s', cmd);
      end
      total=mysql(sprintf('SELECT SUM(cnt) FROM ngsentries e, ngs.naseqs s WHERE s.naseq=e.naseq AND %s AND subrun IN (SELECT subrun FROM subruns where %s) AND e.naseq NOT IN (%s)',lensel,subrunsel,ignoreseqs));
      loaded=ismember(naseq,obj.naseq);
      fracUncl1=length(naseq)/total;
      fprintf('Have %d/%d (%.1f%%) unclustered, singleton reads with correct length of which %d are loaded\n', length(naseq),total, fracUncl1*100, sum(loaded));
      
      % Do the rest on a subsample for speed
      if args.nsample*20<length(naseq)
        sel=randperm(length(naseq),args.nsample*20);
      else
        sel=true(size(naseq));
      end
      snaseq=naseq(sel);
      
      if length(snaseq)>args.nsample
        snaseq=snaseq(1:args.nsample);
      end
      % Check to see how many would cluster on average
      fprintf('Checking clustering of %d random sequences sampled from singletons...',length(snaseq));
      obj.seqs.clusterseqs(snaseq,ones(length(snaseq),1),'mincnt',1,'minaddcnt',2);
      fprintf('done\n');
      ncl=sum(ismember(snaseq,abs([obj.seqs.clusters.members])));
      fracUncl2=1-ncl/length(snaseq);
      fprintf('%d/%d (%.1f%%) of the singletons clustered with existing sequences\n', ncl, length(snaseq), (1-fracUncl2)*100);
      divfrac=fracUncl1*fracUncl2;
      fprintf('Estimated diversity is <= %.1f%%(uncl)  = %.3f%%: %.0f of %.0f total reads\n',fracUncl1*fracUncl2*100,divfrac*100, [divfrac,1]*total);
    end
    
    function clusterstats(obj,ind)
    % Display stats about clustering
      if nargin<2 || isempty(ind)
        ind=true(size(obj.cnt));
      end
      totalcnt=sum(obj.cnt(ind));
      ccnt=zeros(length(obj.seqs.clusters),1);
      fprintf('Summing %d cluster cnts...',length(obj.seqs.clusters));
      for i=1:length(obj.seqs.clusters)
        ccnt(i)=sum(obj.cnt(ind & ismember(obj.naseq,obj.seqs.clusters(i).members)));
      end
      fprintf('done\n');
      nmembers=arrayfun(@(z) length(z.members), obj.seqs.clusters);
      fprintf('Have %d/%d clusters covering %d/%d (%.1f%%) of selected reads\n', sum(ccnt>0),length(ccnt), sum(ccnt), totalcnt, sum(ccnt)/totalcnt*100);
      singles=obj.naseq(ind & obj.cnt==1 & ~ismember(obj.naseq,[obj.seqs.clusters.members]));
      h=hist(singles,unique(singles));
      nsingle=sum(h==1);   % 1 read in only 1 subrun
      fprintf('Have %d (%.1f%%) singleton selected reads that are not in any cluster\n', nsingle, nsingle/totalcnt*100);
      fprintf('Overall, have %d/%d clusters with only 1 member\n', sum(nmembers==1), length(nmembers));
    end
    
    function listclusters(obj,ind,varargin)
      defaults=struct('maxlist',50);
      args=processargs(defaults,varargin);

      if nargin<2
        ind=true(size(obj.naseq));
      end
      extra=containers.Map('KeyType','double','ValueType','char');
      sortkey=containers.Map('KeyType','double','ValueType','double');
      c=obj.seqs.clusters;
      for i=1:length(c);
        cnt=sum(obj.cnt(ind&ismember(obj.naseq,c(i).members)));
        extra(c(i).root)=sprintf('%6d',cnt);
        sortkey(c(i).root)=-cnt;
      end
      obj.listsequences(ind&ismember(obj.naseq,[c.root]),'extra',extra,'sort',sortkey,'extraheader','ClCnt','maxlist',args.maxlist);
    end
    
    function dumpcluster(obj,cluster,varargin)
      defaults=struct('listfrac',0.99,'maxlist',20);
      args=processargs(defaults,varargin);

      if ischar(cluster)
        c=obj.seqs.clusters(strcmp({obj.seqs.clusters.name},cluster));
        if isempty(c)
          error('Cluster %s not found\n', cluster);
        end
      else
        c=obj.seqs.clusters([obj.seqs.clusters.cluster]==cluster);
        if isempty(c)
          error('Cluster %d not found\n', cluster);
        end
      end
      fprintf('Cluster %d',c.cluster);
      if ~isempty(c.name)
        fprintf(' (%s)',c.name);
      end
      cnt=zeros(size(c.members));
      for i=1:length(c.members)
        cnt(i)=sum(obj.cnt(obj.naseq==c.members(i)));
      end
      total=sum(cnt);
      fprintf(' %d members, %d total cnt\n', length(c.members), total);
      ord=sort(cnt,'descend');
      nlist=min(args.maxlist,find(cumsum(cnt)>=total*args.listfrac,1,'first')+1);
      obj.listsequences(ismember(obj.naseq,c.members),'maxlist',nlist);
    end
    
    function findparent(obj,ind,seq,varargin)
    % Find parent of given sequence by searching for substrings in all other seqs
      defaults=struct('minlength',6,'maxlength',10,'rc',false);
      args=processargs(defaults,varargin);
      if iscell(seq)
        seq=seq{1};
      end
      if ischar(seq)
        naseq=obj.seqs.findseq(seq);
      else
        naseq=seq;
        seq=obj.getseq(naseq);
        seq=seq{1};
      end

      otherseqs=obj.seq(ind);
      othercnt=obj.cnt(ind);
      same=strcmp(otherseqs,seq);
      otherseqs=otherseqs(~same);
      othercnt=othercnt(~same);

      if args.rc
        seq=rc(seq);
      end
      maxlen=min(length(seq),args.maxlength);
      cnt=nan(maxlen,length(seq));
      fprintf('Searching %d sequences for substrings of length ', length(otherseqs));
      for n=args.minlength:maxlen
        fprintf('%d...',n);
        for pos=1:length(seq)-n+1
          s=seq(pos:pos+n-1);
          sel=cellfun(@(z) ~isempty(strfind(z,s)),otherseqs);
          cnt(n,pos)=sum(othercnt(sel));
        end
      end
      fprintf('done\n');
      % Remove all the longer substring cnts from the next smaller ones
      for n=1:size(cnt,1)-1
        cnt(n,:)=cnt(n,:)-cnt(n+1,:);
      end
      % Compute entropy of observations
      entropy=nan(size(cnt));
      totalcnt=sum(othercnt);
      for n=1:size(entropy,1)
        entropy(n,:)=-cnt(n,:)/totalcnt*log2(4^-n);
      end
      
      if isfinite(naseq) && naseq>0
        if args.rc
          ti=sprintf('FindParent(RC(%d))',naseq);
        else
          ti=sprintf('FindParent(%d)',naseq);
        end
      else
        ti=sprintf('FindParent(%s)',seq);
      end
      setfig(ti);
      valid=any(isfinite(cnt'));
      plot(entropy(valid,:)');
      ax=axis;
      ax(2)=length(seq);
      axis(ax);
      set(gca,'XTick',1:length(seq));
      set(gca,'XTickLabel',num2cell(seq));
      ylabel('Entropy (bits)');
      legend(arrayfun(@(z) sprintf('N=%d',z),find(valid),'UniformOutput',false));
      title(sprintf('%s %d reads over %d seqs',ti,totalcnt,length(otherseqs)));
    end
    
    function res=prefixcnt(obj,prefix,varargin)
    % Breakdown of fraction of reads with given prefix by subrun
      defaults=struct('mincnt',100,'code','','subruns',[],'rounds',[],'minlength',48);
      args=processargs(defaults,varargin);

      if isempty(args.subruns)
        args.subruns=obj.subruns.subrun;
      end
      if ~isempty(args.rounds) 
        assert(length(args.rounds)==length(args.subruns));
      end
      res=[];
      for ii=1:length(args.subruns)
        sr=args.subruns(ii);
        if length(args.rounds)==length(args.subruns)
          rnd=args.rounds(ii);
        else
          rnd=nan;
        end
        i=find(obj.subruns.subrun==sr);
        if ~isempty(args.code) && ~ismember(obj.subruns.codes{i},args.code)
          continue;
        end
        sel=obj.subrun==sr & ~ismember(obj.naseq,obj.refseqs(1)) & obj.length>=args.minlength;
        cnt=obj.cnt(sel);
        total=sum(cnt);
        if total<args.mincnt
          continue;
        end
        seq=obj.seq(sel);
        match=cellfun(@(z) strncmp(z,prefix,length(prefix)),seq);
        mcnt=sum(cnt(match));
        %fprintf('Round %3d Subrun %3d:  %-20.20s %d/%d = %.1f%%\n', rnd, sr, sprintf('%s-%s',obj.subruns.desc{i},obj.subruns.codes{i}), mcnt, total, mcnt/total*100);
        res=[res,struct('round',rnd, 'run',obj.run, 'subrun',sr,'desc',obj.subruns.desc{i},'codes',obj.subruns.codes{i},'mcnt',mcnt,'total',total,'pct',mcnt/total*100)];
      end
    end

    function x=liststructs(obj,sr,varargin)
    % Count sequences with loops replaced with N's
      defaults=struct('mincnt',10,'maxlist',50,'foldN',false);
      args=processargs(defaults,varargin);
      allcnt=[];
      allseqs={};
      refs=obj.refseqs(true);
      for k=1:length(sr)
        sel=find(obj.subrun==sr(k) & ~ismember(obj.naseq,refs));
        seqs=obj.seq(sel);
        cnts=obj.cnt(sel);
        segextents=nan(length(cnts),11,2);
        [~,ia,ib]=intersect(obj.seqs.naseq,obj.naseq(sel));
        assert(length(ia)==length(sel));
        segextents(ib,:,:)=obj.seqs.segextents(ia,:,:);

        % Convert loop seqs to NaN
        fmtseqs={};
        for i=1:length(seqs)
          f='';
          sei=squeeze(segextents(i,:,:));
          for s=1:size(segextents,2)
            se=squeeze(sei(s,:));
            if all(isfinite(se))
              if s==4 || s==8
                % Loop1 or Loop2
                if args.foldN
                  nstring='N';
                else
                  nstring=sprintf('N%d',se(2)-se(1)+1);
                end
                f=[f,nstring];
              else
                f=[f,seqs{i}(se(1):se(2))];
              end
            end
            if s>1 && length(f)<(s-1)*8
              f=[f,blanks((s-1)*8-length(f))];
            end
          end
          fmtseqs{i}=f;
        end
        [useqs,ia,ic]=unique(fmtseqs);
        ucnt=zeros(size(useqs));
        for i=1:length(ic)
          ucnt(ic(i))=ucnt(ic(i))+cnts(i);
        end
        allseqs=union(allseqs,useqs,'stable');
        for i=1:length(ucnt)
          allcnt(strcmp(allseqs,useqs{i}),k)=ucnt(i);
        end
      end

      totalcnt=sum(allcnt,1);
      zstat=[];
      for i=1:size(allcnt,1)
        p1=allcnt(i,1)/totalcnt(1);
        for j=2:size(allcnt,2)
          pj=allcnt(i,j)/totalcnt(j);
          p=sum(allcnt(i,[1,j]))/sum(totalcnt([1,j]));
          zstat(i,j)=(pj-p1)/sqrt(p*(1-p)*(1/totalcnt(1)+1/totalcnt(j)));
        end
      end

      zmax=max(abs(zstat),[],2);
      [zmax,ord]=sort(zmax,'desc');
      allseqs=allseqs(ord);
      allcnt=allcnt(ord,:);
      zstat=zstat(ord,:);
      srnames=arrayfun(@(z) obj.getsubrunname(z), sr,'UniformOutput',false);

      for i=1:length(srnames)
        fprintf('%10.10s ',srnames{i});
      end
      for i=2:length(sr)
        fprintf(' zstat ');
      end
      for i=2:length(sr)
        fprintf('enrich ');
      end
      fprintf('\n');
      fprintf('%s%s Subruns\n', sprintf('%10d ',sr), sprintf('%6s ',''));
      fprintf('%s%s Total\n', sprintf('%10d ',sum(allcnt,1)), sprintf('%6s ',''));
      for i=1:args.maxlist
        if sum(allcnt(i,:))<args.mincnt
          continue;
        end
        frac=allcnt(i,:)./totalcnt;
        delta=(frac-frac(1))/frac(1);
        fprintf('%s%s %s %s\n', sprintf('%10d ',allcnt(i,:)), sprintf('%6.2f ',zstat(i,2:end)),sprintf('%6.1f ',100*delta(2:end)),allseqs{i});
      end
      x=struct('subruns',sr,'desc',{srnames},'seq',{allseqs},'cnt',allcnt,'zstat',zstat);
      if length(sr)==2
        setfig('liststructs');clf;
        loglog(x.cnt(:,1),x.cnt(:,2),'o');
        hold on;
        sel=zmax>2;
        loglog(x.cnt(sel,1),x.cnt(sel,2),'or');
        
        xlabel(x.desc{1},'Interpreter','none');
        ylabel(x.desc{2},'Interpreter','none');
      end
    end

    function res=runsummary(obj)
      udesc=unique(obj.subruns.desc);
      isref=ismember(obj.naseq,obj.refseqs(1));
      isamplicon=ismember(obj.naseq,obj.ampliconseqs(1)) & ~isref;
      pooledseqs=[obj.seqs.tags(strcmp({obj.seqs.tags.name},'pool')).naseq];
      % Sequences in same clusters as pooled seqs
      pooledclusters=unique(obj.seqs.getclusterindex(pooledseqs));
      pooledclusters=pooledclusters(pooledclusters>0);
      pooledclusterseqs=abs(union([obj.seqs.clusters(pooledclusters).members],[obj.seqs.clusters(pooledclusters).root]));
      inpool=ismember(obj.naseq,pooledseqs)&~isref&~isamplicon;
      inpoolc=ismember(obj.naseq,pooledclusterseqs)&~inpool&~isref&~isamplicon;
      taggedseqs=[obj.seqs.tags.naseq];
      taggedseqs=union(taggedseqs,abs([obj.seqs.clusters(ismember([obj.seqs.clusters.root],taggedseqs)).members]));
      hastags=ismember(obj.naseq,taggedseqs)&~isref&~isamplicon&~inpool;
      isshort=obj.getclass()==1 & ~hastags&~isref&~isamplicon&~inpool;
      issingle=obj.cnt==1 & ~isshort &~hastags&~isref&~isamplicon&~inpool;

      fprintf('%-25.25s %7s %5s %7s %5s %5s %5s %5s %5s %5s %5s\n','Name','total','ref','nrtotal','ampl','pool','poolc','tag','short','1-rd','other');
      res=[];
      for i=1:length(udesc)
        sr=obj.subruns.subrun(strcmp(obj.subruns.desc,udesc{i}));
        sel=ismember(obj.subrun,sr);
        total=sum(obj.cnt(sel));
        nrtotal=sum(obj.cnt(sel&~isref));
        other=sel&~isref&~isamplicon&~inpool&~inpoolc&~hastags&~isshort&~issingle;
        fprintf('%-25.25s %7d %5.2f %7d %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f ',udesc{i},total,...
                sum(obj.cnt(sel&isref))*100/total,nrtotal,...
                sum(obj.cnt(sel&isamplicon))*100/nrtotal,...
                sum(obj.cnt(sel&inpool))*100/nrtotal,...
                sum(obj.cnt(sel&inpoolc))*100/nrtotal,...
                sum(obj.cnt(sel&hastags))*100/nrtotal,...
                sum(obj.cnt(sel&isshort))*100/nrtotal,...
                sum(obj.cnt(sel&issingle))*100/nrtotal,...
                sum(obj.cnt(other)*100/nrtotal));
        % Non-ref fractions by prefix
        for j=1:length(sr)
          frac=sum(obj.cnt(sel&~isref&obj.subrun==sr(j)))/nrtotal;
          if frac>0.01
            fprintf('%s(%d)=%.2f ',obj.subruns.codes{obj.subruns.subrun==sr(j)},sr(j),frac*100);
          end
        end
        fprintf('\n');
        res=[res,struct('desc',udesc{i},'subruns',sr,'total',total,'all',sel,'other',other)];
      end
      if nargout==0
        clear res;
      end
    end
    
    function x=dumptrp(obj,file,trpexpt,varargin)
    % Dump CleaveSeq results to a CSV file for use outside this environment
      defaults=struct('mincnt',10,'naseq',[],'expr','');
      args=processargs(defaults,varargin);
      t=obj.trpexpts(trpexpt);
      a=t.analysis;
      sel=sum(a.cnt(:,2:3),2)>=args.mincnt;
      if ~isempty(args.naseq)
        sel=sel & ismember(a.naseq,args.naseq);
      end
      if ~isempty(args.expr)
        seqs=obj.seqs.getorloadseq(a.naseq(sel));
        match=cellfun(@(z) ~isempty(z),regexp(seqs,args.expr,'forcecelloutput'));
        fsel=find(sel);
        sel(fsel(~match))=false;
      end
      if sum(sel)==0
        error('No data to dump\n');
      end
      fprintf('Dumping %d/%d entries from %s\n', sum(sel), length(sel), t.descr);
      x=table();
      x.naseq=int32(a.naseq(sel));
      x.seq=obj.seqs.getorloadseq(x.naseq);
      x.incnt=a.cnt(sel,1);
      x.clvcnt=a.cnt(sel,2);
      x.unclvcnt=a.cnt(sel,3);
      x.cleavage=a.cleavage(sel);
      x.loop1=obj.seqs.getloop1(x.naseq);
      x.loop2=obj.seqs.getloop2(x.naseq);
      x=sortrows(x,{'loop1','loop2'});
      writetable(x,file);
    end
    
    function res=cleaveseqwells(obj,varargin)
    % Show read count by well
      defaults=struct('link',true,'relref',false,'naseq',[]);
      % link - true to link colorbar axes of all subplots
      % relref - compute nonref,nonamplicon reads divided by ref reads
      % naseq - count only reads for the given naseqs
      args=processargs(defaults,varargin);

      rrmissing=find(isnan([obj.trpexpts.robotrun]));
      if ~isempty(rrmissing)
        fprintf('No robotrun set for %d expts: %s\n', length(rrmissing),sprintf('%d ',rrmissing));
      end
      
      NGSDatabase.open(obj.run);
      [trpexpt,well,pname]=mysql('select t.trpexpt,s.well,p.name from trpexpts t, robot.samples s, robot.programs p where s.program=(select program from robot.runs where run=t.robotrun) and p.program=s.program and t.robotrun is not NULL and s.plate=''Samples'' and t.descr LIKE CONCAT(s.name,''%'')');
      missing=setdiff([obj.trpexpts.trpexpt],trpexpt);
      if ~isempty(missing)
        fprintf('Unable to locate well position for %d expts: %s\n', length(missing),sprintf('%d ',missing));
      end
      row=[];col=[];
      for i=1:length(well)
        row(i)=well{i}(1)-'A'+1;
        col(i)=str2num(well{i}(2:end));
      end
      setfig('CleaveSeq Wells');clf;
      tiledlayout('flow');
      uruns=unique([obj.trpexpts.robotrun]);
      uruns=uruns(isfinite(uruns));
      h=[];
      ref=ismember(obj.naseq,obj.refseqs(0));
      refc=ismember(obj.naseq,obj.refseqs(1));
      amplicon=ismember(obj.naseq,obj.ampliconseqs(1));
      if ~isempty(args.naseq)
        sel=ismember(obj.naseq,args.naseq);
      else
        sel=~refc & ~amplicon;
      end
      maxcnt=0;
      res=[];
      samplename=cell(8,12);
      for i=1:length(uruns)
        nexttile;
        data=nan(9,13);
        for j=1:length(trpexpt)
          t=obj.trpexpts([obj.trpexpts.trpexpt]==trpexpt(j));
          if t.robotrun==uruns(i)
            sr=[t.cleavedsr,t.uncleavedsr];
            if args.relref
              data(row(j),col(j))=sum(obj.cnt(ismember(obj.subrun,sr)&sel))/sum(obj.cnt(ismember(obj.subrun,sr)&ref));
            else
              data(row(j),col(j))=sum(obj.cnt(ismember(obj.subrun,sr)&sel));
            end
            samplename{row(j),col(j)}=obj.subruns.desc(obj.subruns.subrun==t.uncleavedsr);
            pnamei=pname{j};
          end
        end
        pcolor(data);
        res=[res,struct('run',uruns(i),'samplename',{samplename},'reads',data(1:8,1:12))];
        h(i)=gca;
        colorbar;
        axis ij;
        set(gca,'XTick',(1:12)+0.5);
        set(gca,'XTickLabel',arrayfun(@(z) sprintf('%d',z),1:12,'Unif',false));
        set(gca,'YTick',(1:8)+0.5);
        set(gca,'YTickLabel',arrayfun(@(z) sprintf('%c',z+'A'-1),1:8,'Unif',false));;
        title(sprintf('Run %d - %s',uruns(i),pnamei));
        maxcnt=max(maxcnt,nanmax(data(:)));
      end
      if args.link
        linkprop(h,'CLim');
        set(h(1),'CLim',[0,maxcnt]);
      end
      if nargout==0
        clear res;
      end
    end
    
    
  end  % methods
end     

