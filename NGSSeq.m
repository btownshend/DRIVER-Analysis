% Class encapsulating information sequences in database
classdef NGSSeq  < matlab.mixin.Copyable
  properties
    database;   % Which database this represents (e.g. NGS34)
    naseq;
    seq;
    length;
    segextents; % (:,i,j) - segment i (1:11), j=1 is start, j=2 is end,  end==start-1 for missing segment
    class; 	% Class of each sequence (index into seqclasses)
    seqnames;	% Mapping between naseq and names of sequences
    tags;	% Tags 
    hits;	% Hits
    clusters;   % Master data for each cluster
  end

  properties (Transient)
    segments;   % (:,i) Identity of each of i=1:11 segments (pre,s3a,s1a,loop1,s1b,core,s2a,loop2,s2b,s3b,post), automatically rebuilt as needed
    otherseqs;  % Lookup for temporary seqs read from database
    clusterdists;   % Distance between cluster roots (struct:  clusters(N),dists(N,N) )
  end
  
  properties (Constant)
    % Seq classes
    % Each entry consists of name, totallengthrange, l1lenrange or regexp, l2lenrange or regexp, planned?
    % Evaluated in order
    seqclasses={{'Short',[0,40],[],[],false},
                {'NsN30',[],[1,9],[28,32],true},
                {'N0N30',[],[0,0],[28,32],false},
                {'N30N30',[],[28,32],[28,32],false},
                {'N30Ns',[],[28,32],[1,9],true},
                {'N30N0',[],[28,32],[0,0],false},
                {'NsNs',[],[0,9],[0,9],false},
    % Minimum length of any structured library is 84nt
                % {'3wj',[84,999],'NNN RYRYRYRYRY NNNNN RYRYRYRYRY NNN RYRYRYRYRY NNNNN RYRYRYRYRY NNN',[],true},
                % {'4wj',[84,999],'NN RYRYRYRY NNNN RYRYRYRY NN RYRYRYRYRY NNNN RYRYRYRYRY NN RYRYRYRY NNNN RYRYRYRY NN',[],true},
                % {'5wj',[84,999],'NN RYRYRYRY NNNN RYRYRYRY NN RYRYRYRY NNNN RYRYRYRY NN RYRYRYRY NNNN RYRYRYRY NN RYRYRYRY NNNN RYRYRYRY NN',[],true},
                % {'GR',[84,999],'NNNNNNN CGCGTGGATATGGCACGCA NNNNNNNNNN GGGCACCGTAAATGTCC NNNNNN',[],true},
                % {'GRP1d3',[84,999],'NNNN CGCGTGGATATGGCACGCA NNNNNNNNNN GGGCACCGTAAATGTCC NNN',[],true},
                % {'CG',[84,999],'NNNNNNNNN CAAACCATTCGAAAGAGTGGGACG NNNNN CCTCCGGCCTAAACCGAAAGGTAGGTAGCGGGG NNNNNNN',[],true},
                % {'CG_P1d2',[84,999],'NNNNNNN CAAACCATTCGAAAGAGTGGGACG NNNNN CCTCCGGCCTAAACCGAAAGGTAGGTAGCGGGG NNNN',[],true},
                % {'HH',[84,999],'NNNN TGGTATCCAATGAAAATGTACTACCA NNNNNNNNNN CCCAAATAGG NNNNNNN',[],true},
                % {'HH_P1d2',[84,999],'NN TGGTATCCAATGAAAATGTACTACCA NNNNNNNNNN CCCAAATAGG NNNNN',[],true},
                % {'GLY',[84,999],'NNNN GAACCGTTTAATCGGTCGC NNNNNNNN CAAGCTCTGCGCATATGCAGAGTGAAA NNNNNNN GCAAAA NNN',[],true},
                % {'TET',[84,999],'AAAA NNNNN CAGATTTCGATCTG NNNN GGTG NNNNNNNNNNN CACCT',[],true},
                % {'GLN',[84,999],'N?N?NNN TTGGCCCAGTTTATCTGGG NNNNNNNN AGGTCTTTGGCCT NNNN CAA NNN',[],true},
                % {'FMN_p1',[84,999],'NNNNCAGGGTGAAATTCCCGANNNNTGGTATAGTCCANNNAAGTATTTGCTTNNNTTTGGTGAAATTCCAAAACNNNNDGTAGAGTCHNNNNNNNNN',[],true},
                % {'FMN_p4',[84,999],'NNNTTTGGTGAAATTCCAAAACNNNNDGTAGAGTCHNNNNNNNNNGAAGAGAAATCTTCNNNNCAGGGTGAAATTCCCGANNNNTGGTATAGTCCANNN',[],true},
                % {'SAM_p1',[84,999],'NNNNNNGGTGGAGGGACTGGCCCGATGAAACCCNNNNNNNAGAAATNNNNNNAATTCCTGCAGCGGTTTCGCTGAAA',[],true},
                % {'SAM_p2',[84,999],'NNNNAATTCCTGCAGCGGTTTCGCTGAAACTNNNNNNNNNN          GGTGGAGGGACTGGCCCGATGAAACCCNNNNN',[],true},
                {'NsN60',[],[1,9],[58,62],true},
                {'N0N60',[],[0,0],[58,62],false},
                {'N60Ns',[],[58,62],[1,9],true},
                {'N60N0',[],[58,62],[0,0],false},
                {'N60N60',[],[58,62],[58,62],false},
                {'NsNm',[],[0,9],[10,29],false},
                {'NsNl',[],[0,9],[30,999],false},
                {'NmNs',[],[10,29],[0,9],false},
                {'NmNm',[],[10,29],[10,29],false},
                {'NmNl',[],[10,29],[30,999],false},
                {'NlNs',[],[30,999],[0,9],false},
                {'NlNm',[],[30,999],[10,29],false},
                {'NlNl',[],[30,999],[30,999],false},
                {'Unclassifiable',[],[0,100],[0,100],false},
                {'Unparseable',[],[],[],false}};

    % Ribozyme parts (empty sets are wildcards)
    pts={{},
         {'GCTGTC','ACTGTC','TCTGTC','GCTTTC','GCTTTTC','GCTGTA','GCTATC','GTTGTC','GATGTC','GCAGTC','GCCGTC','GCTGTT','GCTGCC','GCTGGC','GCTGAC','GCTCTC','GCGGTC','GCTGTG','GGTGTC','CTGTC','TGTC'},
         {'ACCGGA','ACAGGA','ACTGGA','AACGAA','ACGGGA','ACCGAA','ATCGGA','ACCAGA','ACAGAA'},
         {},
         {'TCCGGT','TCCAGT','TCTGGT','TCCTGT','TCCGAT','TTCTGT','TCTTGT','TCTAGT','TCCCGT','TACTGT','TCATGT','TTCGGT','TCTGT'},
         {'CTGATGA','CTGACGA','CCGATGA','CTGAAGA','CTGTCTGA','CTCATGA','TTGAAGA','CTGAGGA','CTGATTA','CTAATGA','ATGATGA','CTTATGA','CTGATGG','CGGATGA','CTGATGT','CTGATAA','TTGATGA','CTGGTGA','TTGATAA','CTGTTGA','CTAACGA','CAGATGA','CCGACGA','GTGATGA','CTGCTGA','TGATGA','TGACGA'},
         {'GTCC','GTCT','ATCC','TTCC','GTTC','GCCC'},
         {},
         {'GGAC','CGAC'},
         {'GAAACAGC'},
         {}};
    scoringmatrix=NGSSeq.editdistance();
  end
  
  methods(Static)
    function obj=loadobj(obj)
    % Handler for loading from file
      if isstruct(obj)
        fprintf('Warning: loading a struct -- run NGSSeq::cleanup()\n');
      else
        obj.cleanup();
      end
    end
    
    function sm=editdistance()
    % Scoring matrix that counts edit distance
    % Use with gapopen=1
      sm=double(nuc44()>0);
      % N's match anything
      sm(end,:)=1;
      sm(:,end)=1;
      sm=sm-1;
    end
    
    function [dist,al]=seqcompare(seq1,seq2)
    % Align 2 sequences using scoring based on edit differences
    % Substituions count as 1, insertion/deletions as 1
    %      [dist,al]=nwalign(seq1(seq1~='N'),seq2(seq2~='N'),'alphabet','nt','scoringmatrix',eye(4)-1,'gapopen',1);
      [dist,al]=nwalign(seq1,seq2,'alphabet','nt','scoringmatrix',NGSSeq.scoringmatrix,'gapopen',1);
      dist=-dist;
    end
  end
  
  methods
    function obj=NGSSeq(database)
      obj.database=database;
      obj.otherseqs=containers.Map('KeyType','double','ValueType','char');
      obj.seq={};
      obj.segments={};
    end
    
    function cleanup(obj)
    % Add computed fields where missing
      if isempty(obj.otherseqs)
        obj.otherseqs=containers.Map('KeyType','double','ValueType','char');
      end
      if isempty(obj.length) || length(obj.length)~=length(obj.naseq)
        fprintf('Computing lengths...');
        obj.length=cellfun(@(z) length(z), obj.seq);
        fprintf('done\n');
      end
      if isempty(obj.segextents) || size(obj.segextents,1)~=length(obj.naseq)
        if isempty(obj.segextents)
          fprintf('Segmenting...');
        else
          fprintf('Re-segmenting because number of segments is %d, but have %d seqs...',size(obj.segextents,1),length(obj.naseq));
        end
        obj.segmentall();
        fprintf('done\n');
      end
      if strcmp(class(obj.segextents),'double')
        obj.segextents=single(obj.segextents);
      end
      if ~strcmp(class(obj.length),'uint16')
        obj.length=uint16(obj.length);
      end
      if ~strcmp(class(obj.class),'uint8')
        obj.class=uint8(obj.class);
      end
      if isempty(obj.segments) || size(obj.segments,1)~=length(obj.naseq)
        fprintf('Setting segments...\n');
        obj.setsegments();
      end
      if isempty(obj.class) || length(obj.class)~=length(obj.naseq) || max(obj.class)>length(obj.seqclasses)
        fprintf('Classifying...');
        obj.classifyseqs();
        fprintf('done\n');
      end
    end
    
    
    function dbopen(data)
    % Open a database connection using given user, close all other connections
      NGSDatabase.open(data.database);
    end
    
    function dbenablewrites(data,password)
    % Enable writes by changing which user access the database
      NGSDatabase.enablewrites(password);
    end
    
    function dbdisablewrites(data)
    % Disable writes by changing which user access the database
      NGSDatabase.disablewrites();
    end
    
    function dbload(data,mincnt)
    % Load data from database
      if nargin<2
        mincnt=1;
      end
        
      data.dbopen();
      cmd=sprintf(['SELECT naseq,seq '...
                   'FROM ngs.naseqs ' ...
                   'WHERE naseq IN (SELECT naseq FROM ngsentries WHERE cnt>=%d) '...
                   'AND naseq IN (SELECT naseq FROM ngsentries WHERE subrun IN (SELECT subrun FROM subruns WHERE descr IS NOT NULL AND descr!=''PhiX''))'],mincnt);

      fprintf('Loading sequences...');
      [naseq1,seq1]=mysql(cmd);
      fprintf('1...');
      cmd2=sprintf(['SELECT s.naseq,seq '...
                    'FROM ngs.naseqs s, ngs.naseqnames n ' ...
                    'WHERE s.naseq=n.naseq']);
      % For some reason combining cmd and cmd2 in a UNION results in very very slow execution (>10minutes), even though explain plan looks the same as separate (which takes 5sec)
      [naseq2,seq2]=mysql(cmd2);
      [naseq,ia,ib]=union(naseq1,naseq2,'stable');
      seq=[seq1(ia);seq2(ib)];

      fprintf('2...');
      cmd3=sprintf(['SELECT s.naseq,seq '...
                    'FROM ngs.naseqs s, ngs.clusters n ' ...
                    'WHERE s.naseq=n.root']);
      [naseq3,seq3]=mysql(cmd3);
      [naseq,ia,ib]=union(naseq,naseq3,'stable');
      seq=[seq(ia);seq3(ib)];

      fprintf('done\n');

      fprintf('Computing seq lengths...');
      len=cellfun(@(z) length(z), seq);
      fprintf('done\n');

      data.naseq=naseq;
      data.seq=seq;
      data.length=uint16(len);
      data.dbupdate();
    end
    
    function dbupdate(data)
      data.dbloadnames();
      data.dbloadtags();
      data.dbloadhits();
      data.dbloadclusters();
      data.cleanup();
    end
    
    function dbloadnames(data)
    % Load naseqnames
      data.dbopen();

      cmd=['SELECT n.naseq,n.name,s.seq ',...
           'FROM ngs.naseqnames n, ngs.naseqs s '...
           'WHERE n.naseq=s.naseq '...
           'ORDER BY n.naseq;'...
          ];
      fprintf('Loading names...');
      [naseq,name,seq]=mysql(cmd);
      fprintf('loaded %d names\n',length(naseq));
      data.seqnames=struct('naseq',num2cell(naseq),'name',name);
      % Don't load these naseqs as that causes the segments to be out of sync
      % for i=1:length(naseq)
      %   if ~ismember(naseq(i),data.naseq)
      %     % Need to add
      %     data.naseq(end+1)=naseq(i);
      %     data.seq{end+1}=seq{i};
      %     data.length(end+1)=length(seq{i});
      %   end
      % end
    end

    function dbloadtags(data)
    % Load naseqtags
      data.dbopen();

      cmd=['SELECT t.naseq,t.name,t.value,s.seq ',...
           'FROM ngs.v_tags t, ngs.naseqs s '...
           'WHERE t.naseq=s.naseq '...
           'ORDER BY t.naseq;'...
          ];
      fprintf('Loading tags...');
      [naseq,name,value,seq]=mysql(cmd);
      fprintf('loaded %d tags\n',length(naseq));
      data.tags=struct('naseq',num2cell(naseq),'name',name,'value',value);
      % Don't load these naseqs as that causes the segments to be out of sync
      % for i=1:length(naseq)
      %   if ~ismember(naseq(i),data.naseq)
      %     % Need to add
      %     data.naseq(end+1)=naseq(i);
      %     data.seq{end+1}=seq{i};
      %     data.length(end+1)=length(seq{i});
      %   end
      % end
    end

    function dbloadhits(data)
    % Load ngs.hitseqs
      data.dbopen();

      cmd=sprintf('SELECT s.naseq,r.descr,s.fold FROM ngs.hitruns r, ngs.hitseqs s WHERE s.hitrun=r.hitrun ORDER by r.descr');
      fprintf('Loading hits...');
      [naseq,name,fold]=mysql(cmd);
      fprintf('loaded %d hits\n',length(naseq));
      data.hits=struct('naseq',num2cell(naseq),'name',name,'fold',num2cell(fold));
    end


    function dbaddtag(data,naseq,name,value,replace)
      if nargin<5
        replace=false;   % Don't replace existing tags
      end
      cmd=sprintf(['SELECT tagtype '...
                   'FROM ngs.tagtypes '...
                   'WHERE name=''%s'' '],name);
      ttype=mysql(cmd);
      if isempty(ttype)
        fprintf('Tag type ''%s'' not found in database\n', name);
        return;
      end
      if nargin<4
        value=[];
      elseif ~ischar(value)
        if round(value)==value
          value=sprintf('%d',value);  
        else
          value=sprintf('%f',value);
        end
      end
      
      % Check how many there are for this name, naseq
      existing=mysql(sprintf('SELECT value FROM ngs.naseqtags WHERE tagtype=%d and naseq=%d',ttype,naseq));
      if replace && length(existing)>0
        if (~isempty(value) && ismember(value,existing)) || (isempty(value) && ismember('',existing))
          % Already exists
          return;
        end
        if length(existing)==1
          fprintf('Replacing existing entry of %s=%s for naseq %d\n', name, strjoin(existing,','), naseq);
          if ~isempty(value)
            % Has a value
            cmd=sprintf('UPDATE ngs.naseqtags SET value=''%s'' WHERE naseq=%d AND tagtype=%d',value,naseq,ttype);
          else
            % No value, this is a flag
            cmd=sprintf('UPDATE ngs.naseqtags SET value=NULL WHERE naseq=%d AND tagtype=%d',naseq,ttype);
          end
        else
          error('Have multiple existing values for naseq=%d, tag %s\n', naseq, name);
        end
      else
        if length(existing)>0
          fprintf('Note: Already have %d entries for naseq %.0f, tagtype %s;  adding another\n', length(existing), naseq, name);
        end
        if ~isempty(value)
          % Has a value
          cmd=sprintf(['INSERT INTO ngs.naseqtags(naseq,tagtype,value) '...
                       'VALUES(%d,%d,''%s'')'],naseq,ttype,value);
        else
          % No value, this is a flag
          cmd=sprintf(['INSERT INTO ngs.naseqtags(naseq,tagtype,value) '...
                       'VALUES(%d,%d,NULL)'],naseq,ttype);
        end
      end
      fprintf('mysql: %s\n', cmd);
      nadd=mysql(cmd);
      fprintf('Added %d tags\n',nadd);
      % If adding an ID tag, make sure to remove any 'notmade' tag
      if strcmp(name,'id')
        cmd=sprintf('DELETE FROM ngs.naseqtags WHERE naseq=%d AND tagtype=(SELECT tagtype FROM ngs.tagtypes WHERE name=''notmade'')',naseq)
        fprintf('Should run: %s\n',cmd);
      end
    end
    
    function dbaddname(data,naseq,name)
      if any([data.seqnames.naseq]==naseq & strcmp({data.seqnames.name},name))
        % Already there
        return;
      end
      data.dbopen();
      % Check how many there are for this name, naseq
      num_all=mysql(sprintf('SELECT COUNT(*) FROM ngs.naseqnames WHERE name=''%s'' and naseq=%d',name,naseq));
      if num_all==1
        % Already exists
        return;
      end
      oldnaseq=mysql(sprintf('SELECT naseq FROM ngs.naseqnames WHERE name=''%s''',name));
      if ~isempty(oldnaseq)
        error('dbaddname: Attempt to add name %s for naseq %d, but it is already bound to %d\n',name,naseq,oldnaseq);
      end
      oldname=mysql(sprintf('SELECT name FROM ngs.naseqnames WHERE naseq=%d',naseq));
      if ~isempty(oldname)
        fprintf('Updating name for %d from %s to %s\n', naseq, oldname{1}, name);
        cmd=sprintf(['UPDATE ngs.naseqnames '...
                     'SET name=''%s'' '...
                     'WHERE naseq=%d;'],name,naseq);
      else
        cmd=sprintf(['INSERT INTO ngs.naseqnames(naseq,name) '...
                     'VALUES '...
                     '(%d,''%s'');'],naseq,name);
      end
      fprintf('mysql: %s\n', cmd);
      mysql(cmd);
      %data.dbloadnames();
    end
    
    function prune(data,naseq)
      if any(~ismember(naseq,data.naseq))
        error('Some naseq passed in are not in seqs');
      end
      naseq=union(naseq,[data.seqnames.naseq]);
      naseq=union(naseq,[data.tags.naseq]);
      naseq=union(naseq,[data.clusters.root]);
      if isempty(setdiff(data.naseq,naseq))
        return;
      end
      fprintf('Keeping only %d/%d naseqs in seqs structure\n', length(naseq), length(data.naseq));
      pause(5);
      keep=ismember(data.naseq,naseq);
      fns={'naseq','seq','length','segextents','class','segments'};
      for i=1:length(fns)
        data.(fns{i})=data.(fns{i})(keep,:,:);
      end
    end
    
    function [s,ext]=segment_re(data,seq)
    % Decompose sequences into parts using regular expression (all or none)
    % Return segments and extents
      expr='^';
      for i=1:length(NGSSeq.pts)
        expr=[expr,'('];
        if isempty(NGSSeq.pts{i})
          if i>1 && i<length(NGSSeq.pts) % No padding on end s
            expr=[expr,'.*'];
          end
        else
          for j=1:length(NGSSeq.pts{i})
            if j>1
              expr=[expr,'|'];
            end
            expr=[expr,NGSSeq.pts{i}{j}];
          end
        end
        expr=[expr,')'];
      end
      expr=[expr,'$'];
      fprintf('RE=%s\n',expr);  
      [sc,extc]=regexp(seq,expr,'tokens','tokenExtents');
      s=cell(length(seq),length(NGSSeq.pts));
      ext=nan(length(seq),length(NGSSeq.pts),2,'single');
      for i=1:length(sc)
        if isempty(sc{i})
          %s(i,:)=nan;
          ext(i,:,:)=nan;
        else
          s(i,:)=sc{i}{1};
          ext(i,:,:)=extc{i}{1};
        end
      end
    end
    
    function segmentall(data)
      if isempty(data.segextents)
        [~,data.segextents]=data.segment_re(data.seq);
        fprintf('Segmented %d (%.1f%%) sequences completely\n',sum(isfinite(data.segextents(:,1,1))),100*mean(isfinite(data.segextents(:,1,1))));
      elseif size(data.segextents,1)~=length(data.naseq)
        data.segextents(end+1:length(data.naseq),:,:)=nan;
      end
      missing=any(any(isnan(data.segextents(:,:,:)),3),2);
      fprintf('Decomposing %d/%d unsegmented sequences\n', sum(missing), length(missing));
      ext=data.decompose(data.seq(missing));
      data.segextents(missing,:,:)=ext;
      fprintf('Coverage: %s\n', sprintf('%.1f ',mean(all(isfinite(ext),3))*100));
      data.fixstemhelices();
      data.retryUnparseables();
      data.setsegments();
      data.classifyseqs();
    end

    function helices=findhelices(data,seq,varargin)
    % Find all RNA helices possible in seq
    % Not currently used
      defaults=struct('minstem',3,'minloop',3);
      args=processargs(defaults,varargin);
      helices=zeros(length(seq),length(seq));
      for i=1:length(seq)
        % Start at position i
        for stemlen=args.minstem:floor((length(seq)-(i-1)-args.minloop)/2)
          for looplen=args.minloop:length(seq)-(i-1)-2*stemlen
            s1=seq(i:i+stemlen-1);
            s2=seq(i+stemlen+looplen:i+stemlen*2+looplen-1);
            if data.checkhelix(s1,s2)
              %fprintf('Helix     at %d,%d,%d: %s/%s\n',i,stemlen,looplen,s1,s2);
              helices(i,i+stemlen*2+looplen-1)=stemlen;
            else
              ;%fprintf('Not helix at %d,%d,%d: %s/%s\n',i,stemlen,looplen,s1,s2);
            end
          end
        end
      end
    end
    
    function retryUnparseables(data)
    % Take another pass through unparseables to see if they can be segmented
      unp=find(any(any(isnan(data.segextents(:,3:9,:)),3),2) & isfinite(data.segextents(:,2,2)) & isfinite(data.segextents(:,10,1))& data.length>=41);   % Missing something in the middle
      fprintf('Have %d/%d unparseable sequences\n', length(unp), size(data.segextents,1));
      found=0;
      nadebug=[472219451];   % naseqs to monitor
      for ii=1:length(unp)
        ind=unp(ii);
        debug=ismember(data.naseq(ind),nadebug);
        if debug
          fprintf('---- naseq %d -----\n', data.naseq(ind));
        end
        seq=data.seq{ind}(data.segextents(ind,2,2)+1:data.segextents(ind,10,1)-1);
        % find all possible helices
        helices={};
        minloop=3; minstem=2;
        for s1a=1:min(3,length(seq))
          for l1=s1a+minstem:length(seq)
            s1len=l1-s1a;
            for s1b=l1+minloop:length(seq)-s1len-6-2*minstem-minloop+1
              l1len=s1b-l1;
              corea=s1b+s1len;
              if ~data.checkhelix(seq(s1a:l1-1),seq(s1b:s1b+s1len-1)) || (l1len>minloop+2 && data.checkhelix(seq(l1),seq(s1b-1)))
                continue;
              end
              if s1a>1 && data.checkhelix(seq(s1a-1),seq(corea-1))
                % Would form helix with an outer nt
                continue;
              end
              for s2a=corea+6:corea+7
                for l2=s2a+minstem:length(seq)
                  s2len=l2-s2a;
                  for s2b=l2+minloop:length(seq)-s2len+1
                    l2len=s2b-l2;
                    if s2b+s2len<length(seq)-1
                      continue;   % Pushes too much onto next element
                    end
                    if ~data.checkhelix(seq(s2a:l2-1),seq(s2b:s2b+s2len-1)) || (l2len>minloop+2 && data.checkhelix(seq(l2),seq(s2b-1)))
                      continue;
                    end
                    if s2b+s2len < length(seq) && data.checkhelix(seq(s2a-1),seq(s2b+s2len))
                      % Would form a helix with outer nt
                      continue;
                    end
                    helices{end+1}=[s1a,l1,s1b,s1b+s1len,s2a,l2,s2b,s2b+s2len];
                  end
                end
              end
            end
          end
        end
      
        if length(helices)>1
          % Find best one
          if debug
            fprintf('Have %d possible segmentations:\n',length(helices));
          end
          cores={};score=[];excess=[];
          for j=1:length(helices)
            h=helices{j};
            cores{j}=seq(h(4):h(5)-1);
            score(j)=nwalign(cores{j},'CTGATGA');
            excess(j)=h(1)-1+(length(seq)+1-h(end));
            if debug
              fprintf('%s %f %d\n', cores{j}, score(j),excess(j));
            end
          end
          [~,best]=max(score-100*excess);
          bestscore=score(best);
          bestexcess=excess(best);
          % if data.naseq(ind)==429603190
          %   keyboard;
          % end
          if bestscore>10 || bestexcess==0
            helices={helices{best}};
          else
            if debug
              fprintf('Poor scores: %f, excess=%d\n', bestscore, bestexcess);
            end
          end
          %fprintf('Now have %d possible segmentations:\n',length(helices));
        end
        if length(helices)==1
          segextents=squeeze(data.segextents(ind,:,:));
          segextents(3:10,1)=helices{1}+segextents(2,2);
          segextents(2:9,2)=segextents(3:10,1)-1;
          if debug
            fprintf('%s\n',data.naseqformat(data.naseq(ind)));
          end
          data.segextents(ind,:,:)=segextents;
          if debug
            fprintf('%s\n',data.naseqformat(data.naseq(ind)));
          end
          found=found+1;
        end
      end
      fprintf('Successful on %d/%d unparseable sequences\n', found, length(unp));
    end
    
    function valid=checkhelix(data,s1,s2)
      valid=false;
      if length(s1)==length(s2)
        for i=1:length(s1)
          if (s1(i)=='A' && s2(end+1-i)~='T') || (s1(i)=='C' && s2(end+1-i)~='G') || (s1(i)=='G' && s2(end+1-i)~='C' && s2(end+1-i)~='T') || (s1(i)=='T' && s2(end+1-i)~='A' && s2(end+1-i)~='G')
            %fprintf('checkhelix(%s,%s)->%d\n', s1, s2, valid);
            return
          end
        end
        valid=true;
      end
      %fprintf('checkhelix(%s,%s)->%d\n', s1, s2, valid);
    end
    
    function extents=decompose(data, seq)
    % Decompose sequences into parts sequentially
    % Return segmentation as matrix of start positions and segmented strings
    % Return NaN in extents for unmatched parts
      p={};
      origseq=seq;
      fprintf('Decomposing %d sequences\n',length(seq));
      searchpos=ones(length(seq),1);   % Starting position of searches
      nonempty=cellfun(@(z) ~isempty(z), NGSSeq.pts);
      pts=NGSSeq.pts();
      extents=nan(length(seq),length(pts),2);
      for i=1:length(pts)
        fprintf('Group %d\n',i);
        hits=nan(length(seq),length(pts{i}));
        for j=1:length(pts{i})
          fprintf(' %s...', pts{i}{j});
          loc=strfind(seq,pts{i}{j});
          sel=cellfun(@(z) ~isempty(z), loc);   % Hits
          fprintf('%d hits (%.1f%%)\n', sum(sel),mean(sel)*100);
          hits(sel,j)=cellfun(@(z) z(1), loc(sel));
        end

        m=nan(size(hits,1),1);
        o=m;
        if i==2
          % Use hit that occurs earliest in seq
          for k=1:size(hits,2)
            sel=isfinite(hits(:,k)) & (isnan(o) | hits(:,k)<m) & hits(:,k)<5;  % Don't allow more than 4nt prefix
            o(sel)=k;
            m(sel)=hits(sel,k);
          end
        else
          % Use "best" hit (first in list)
          for k=1:size(hits,2)
            o(isnan(o)&isfinite(hits(:,k)))=k;
            m(isnan(m))=hits(isnan(m),k);
          end
        end
        fnd=isfinite(m);
        fprintf(' Covered %.1f%% of sequences\n', mean(fnd)*100);
        if i>1
          % Contiguous only
          noncontig=m>1 & isfinite(extents(:,i-1,2));
          fprintf('  Of possible matches %.1f%% were rejected for not being contiguous with prior\n', mean(noncontig&fnd)/mean(fnd)*100);
          fnd(noncontig)=false;
        end
        fnd=find(fnd);
        extents(fnd,i,1)=m(fnd)+searchpos(fnd)-1;
        len=cellfun(@(z) length(z),pts{i}(o(fnd)));
        extents(fnd,i,2)=extents(fnd,i,1)+len(:)-1;
        if i==5 || i==9
          % Check stem helix
          valid=false(size(fnd));
          for k=1:length(fnd)
            if all(isfinite(extents(fnd(k),[i,i-2],[1,2])))
              valid(k)=data.checkhelix(origseq{fnd(k)}(extents(fnd(k),i-2,1):extents(fnd(k),i-2,2)),origseq{fnd(k)}(extents(fnd(k),i,1):extents(fnd(k),i,2)));
            else
              valid(k)=false;
            end
          end
          fprintf('%d/%d of the helices are valid\n', sum(valid), length(valid));
          extents(fnd(~valid),i,:)=nan;
          fnd=fnd(valid);
        end
        searchpos(fnd)=extents(fnd,i,2)+1;
        for k=1:length(fnd)
          seq{fnd(k)}=origseq{fnd(k)}(searchpos(fnd(k)):end);
        end
      end
      for rpt=1:3  % Inter-dependent, so redo a few times
        % Fill in adjacent segments
        extents(:,1,1)=1;
        extents(:,end,2)=cellfun(@(z) length(z),origseq);
        for i=1:size(extents,2)-1
          sel=isnan(extents(:,i,2))&isfinite(extents(:,i+1,1));
          extents(sel,i,2)=extents(sel,i+1,1)-1;
        end
        for i=2:size(extents,2)
          sel=isnan(extents(:,i,1))&isfinite(extents(:,i-1,2));
          extents(sel,i,1)=extents(sel,i-1,2)+1;
        end
        % % Infer specific ranges
        % for k=[3,7]
        %   ptlen=median(cellfun(@(z) length(z), NGSSeq.pts{k}));
        %   sel=isfinite(extents(:,k,1))&~isfinite(extents(:,k,2))&isfinite(extents(:,k+2,1))&(extents(:,k+2,1)-extents(:,k,1)>=ptlen);  
        %   extents(sel,k,2)=extents(sel,k,1)+ptlen-1;
        %   fprintf('Inferred group %d length %d for %d sequences\n', k, ptlen, sum(sel));
        % end
        % for k=[5]
        %   ptlen=median(cellfun(@(z) length(z), NGSSeq.pts{k}));
        %   sel=isfinite(extents(:,k,2))&~isfinite(extents(:,k,1))&isfinite(extents(:,k-2,2))&(extents(:,k,2)-extents(:,k-2,2)>=ptlen);  
        %   extents(sel,k,1)=extents(sel,k,2)-ptlen+1;
        %   fprintf('Inferred group %d length %d for %d sequences\n', k, ptlen, sum(sel));
        % end
        % % Missing loop2 with (short) unmatched stem2a:
        % sel=isnan(extents(:,7,2)) & isnan(extents(:,8,1)) & (extents(:,8,2)-extents(:,7,1))+1<=4;
        % fprintf('Inferred short loop2 for %d sequences\n',sum(sel));
        % extents(sel,7,2)=extents(sel,8,2);
        % extents(sel,8,1)=extents(sel,8,2)+1;
        % 5/2019 - above is handled better by fixstemhelices()
      end
      
    end
    
    function hlen=splithelices(data,seq)
    % Split seq into 3 parts -- left helix, loop, right helix
    % Return extents (n,3,2) showing where each extent starts/ends
      minhelix=2;   % Mininum length of a helix
      minloop=3;    % Minimum loop length
      seqlen=cellfun(@(z) length(z), seq);
      hlen=zeros(length(seq),1);   % Length of helices
      sel=find(seqlen>=2);   % Seqs with possibly longer helix
      pos=1;
      while ~isempty(sel)
        fprintf('Checking %d sequences for length %d helices\n', length(sel), pos);
        nt1=cellfun(@(z) z(pos), seq(sel));
        nt2=cellfun(@(z) z(end-pos+1), seq(sel));
        match=((nt1=='A' & nt2=='T') | (nt1=='C' & nt2=='G') | (nt1=='G' & (nt2=='C'|nt2=='T')) | (nt1=='T' & (nt2=='A'|nt2=='G')));
        sel=sel(match);
        hlen(sel)=pos;
        pos=pos+1;
        sel=sel(seqlen(sel)>=pos*2);
      end
      for i=1:3
        remove=hlen>0 & seqlen-2*hlen<minloop;
        fprintf('Have %d helices with loop length < %d -> opening by 1 base pair (Pass %d)\n', sum(remove),minloop,i);
        hlen(remove)=hlen(remove)-1;
      end
      remove=hlen>0 & hlen<minhelix;
      fprintf('Have %d helices with length < %d -> removing\n', sum(remove),minhelix);
      hlen(remove)=0;
    end
    
    function fixstemhelices(data)
    % Adjust segmentation so the stems are true helices
      for pos=[3,7]   % Both stems
        fprintf('Fixing helices at segment %d-%d\n', pos, pos+2);
        p1=data.segextents(:,pos,1);
        p2=data.segextents(:,pos+2,2);
        sel=isfinite(p1)&isfinite(p2) & p2>p1;
        seqs=arrayfun(@(z,a,b) z{1}(a:b), data.seq(sel), p1(sel), p2(sel),'UniformOutput',false);
        hlen=data.splithelices(seqs);
        oldhlen=data.segextents(sel,pos,2)-data.segextents(sel,pos,1)+1;
        fprintf('Have different helix length for %d/%d sequences\n', sum(hlen~=oldhlen), length(hlen));
        %oldextents=data.segextents;
        data.segextents(sel,pos,2)=data.segextents(sel,pos,1)+hlen-1;
        data.segextents(sel,pos+2,1)=data.segextents(sel,pos+2,2)-hlen+1;
        data.segextents(sel,pos+1,1)=data.segextents(sel,pos,2)+1;
        data.segextents(sel,pos+1,2)=data.segextents(sel,pos+2,1)-1;
      end
      % Redo segment strings and classifications
      data.setsegments();
      data.classifyseqs();
    end
    
    function setsegments(data)
    % Fill in segments for inferred ranges (i.e. prior/next are set)
      data.segments=cell(size(data.segextents,1),size(data.segextents,2));
      [data.segments{:}]=deal('');
      
      for i=1:size(data.segextents,2)
        sel=find(isfinite(data.segextents(:,i,1))&isfinite(data.segextents(:,i,2)));
        if ~isempty(sel)
          empty=data.segextents(sel,i,1)>data.segextents(sel,i,2);
          fprintf('Setting segments for %d empty and %d non-empty sequences for part %d...', sum(empty),length(sel)-sum(empty), i);
          [data.segments{sel(empty),i}]=deal('');
          sel=sel(~empty);
          s=cell(length(sel),1);
          seqs=data.seq(sel);
          start=data.segextents(sel,i,1);
          finish=data.segextents(sel,i,2);
          for j=1:length(seqs)
            s{j}=seqs{j}(start(j):finish(j));
          end
          data.segments(sel,i)=s;
          fprintf('done\n');
        end
      end
    end

    function naseq=findseq(data,seq)
    % Find given seq, return -1 if not found
      cleanseq=upper(seq(seq~=' '));
      cleanseq(cleanseq=='U')='T';
      if any(~ismember(cleanseq,'ACGTN'))
        error('Sequence invalid: %s',seq);
      end
      ind=strcmp(data.seq,cleanseq);
      if ~any(ind)
        % Check database
        %fprintf('Checking database...');
        cmd=sprintf('select naseq from ngs.naseqs where seq=''%s''',cleanseq);
        res=mysql(cmd);
        %fprintf('done\n');
        if ~isempty(res)
          naseq=res;
        else
          naseq=-1;
        end
      else
        naseq=data.naseq(ind);
      end
    end
    
    function naseq=dbnewnaseq(data,seq) 
    % Add new seq to db and return its naseq
      addseqcmd='INSERT INTO ngs.naseqs(seq) VALUES(''%s'')';
      mysql(sprintf(addseqcmd,seq));
      naseq=mysql('select LAST_INSERT_ID()');
    end
    
    function ind=naseq2ind(data,naseq)
      [~,ind]=ismember(naseq,data.naseq);
      if any(ind==0)
        fprintf('naseq2ind: naseqs [%s] not found\n',sprintf('%d ',naseq(ind==0)));
      end
    end
    
    function res=naseqformat(data,naseqs,varargin)
      defaults=struct('squeeze',false,'ignoresegs',[]);
      args=processargs(defaults,varargin);

      indsz=data.naseq2ind(naseqs);
      inds=indsz(indsz~=0);
      % Show decomposition of seq{inds}
      segments=data.segments(inds,:);
      segextents=data.segextents(inds,:,:);
      seqs=data.seq(inds);
      minpos=ones(1,size(segments,2));
      for i=2:length(minpos)
        lens=seglen(data,i-1,inds);
        lens=lens(isfinite(lens));
        if isempty(lens)
          minpos(i)=minpos(i-1)+1;
        else
          minpos(i)=ceil(prctile(lens,95))+minpos(i-1)+1;
        end
      end

      res=cell(length(inds),1);
      for i=1:size(segments,1)
        pos=1;
        seq=seqs{i};
        s=segments(i,:);
        e=squeeze(segextents(i,:,:));
        startpos=nan; 
        r='';
        for j=1:size(e,1)
          if ismember(j,args.ignoresegs)
            continue;
          end
          if isfinite(e(j,1))
            assert(~isempty(args.ignoresegs) || isnan(startpos));
            startpos=e(j,1);
            if (j==5 || j==9) && isfinite(e(j,2)) && length(r)+(e(j,2)-startpos+2) < minpos(j+1)
              % Right justify 2nd half of stem loops
              r=[r,blanks(minpos(j+1)-length(r)-(e(j,2)-startpos+2))];
            elseif length(r)<minpos(j)
              r=[r,blanks(minpos(j)-length(r))];
            end
          end
          if isfinite(startpos) && isfinite(e(j,2))
            r=[r,seq(startpos:e(j,2)),' '];
            startpos=nan;
          end
        end
        if isempty(args.ignoresegs) && ~strcmp(strrep(r,' ',''),seq)
          fprintf('Mismatch of "%s" and "%s"\n', strrep(r,' ',''),seq);
          keyboard
        end
        res{i}=r;
      end
      if size(res,1)<length(naseqs)
        % Some missing ones
        fprintf('Appending...\n');
        resz=cell(length(naseqs),1);
        resz(indsz~=0)=res;  % Make the mapping match
        [resz{indsz==0}]=deal('?');
        res=resz;
      end

      if length(res)==1
        res=res{1};
      else
        lens=cellfun(@(z) length(z), res);
        tgtlen=round(prctile(lens,95));
        for i=1:length(res)
          if lens(i)<tgtlen
            res{i}=[res{i},blanks(tgtlen-lens(i))];
          end
        end
      end
      
      if args.squeeze
        res=regexprep(res,'  *',' ');
        res=regexprep(res,'^ ','');
        res=regexprep(res,' $','');
      end
    end
    
    function showstats(data,naseq,varargin)
      defaults=struct('maxchanges',5);
      args=processargs(defaults,varargin);

      if length(naseq)>1
        error('showstats: expected exactly 1 naseq, got %d\n', length(naseq));
      end
      if ~ismember(naseq,data.naseq)
        error('showstats: naseq %d not loaded\n', naseq);
      end
      lbls=data.getlabels(naseq);
      fprintf('%d %s (%d)', naseq, data.naseqformat(naseq), data.getlength(naseq));
      if ~isempty(lbls)
        fprintf(' %s',lbls);
      end
      fprintf('\n');
      fprintf('Class: %s, %d/%d\n', data.seqclasses{data.getclass(naseq)}{1},data.getloop1len(naseq),data.getloop2len(naseq));
      %      fprintf('prefix=%d, core=%d, suffix=%d, ref=%d, loops=(%d,%d)\n', data.hasprefix(sel(1)), data.hascore(sel(1)), data.hassuffix(sel(1)),data.isref(sel(1)),data.loop1len(sel(1)),data.loop2len(sel(2)));
      data.comparelabelled(naseq,'maxchanges',args.maxchanges);
    end
    
    function comparelabelled(data,naseq,varargin)
    % Compare given naseq with all labelled sequences
      defaults=struct('maxchanges',5);
      args=processargs(defaults,varargin);
      data.compare(naseq,unique([data.seqnames.naseq,data.tags.naseq]),'maxchanges',args.maxchanges);
    end
    
    function [bestnaseq,descr,bestchanges]=compare(data,naseq,others,varargin)
    % Perform local alignment (counting edit difference) to other naseqs and show closest match
      defaults=struct('maxchanges',5,'silent',false,'exact',false);
      args=processargs(defaults,varargin);

      seq=data.getorloadseq(naseq);
      seq=seq{1};
      bestnaseq=[]; bestchanges=args.maxchanges+1; bestalign=[];
      if isempty(seq)
        descr='';
        return;
      end
      oseqs=data.getorloadseq(others);
      for i=1:length(others)
        other=others(i);
        if ~args.exact && other==naseq
          continue;
        end
        oseq=oseqs{i};
        if isempty(oseq)
          fprintf('compare: other naseq %d not loaded\n',other);
          continue;
        end
        if abs(length(seq)-length(oseq)) > bestchanges
          % Not possibly to match best
          continue;
        end
        [changes,align]=data.seqcompare(seq,oseq);
        if changes<bestchanges
          bestnaseq=other;
          bestchanges=changes;
          bestalign=align;
        end
      end
      if isempty(bestalign)
        if ~args.silent
          fprintf('No alignment returned.\n');
        end
      else
      end
      if bestchanges>args.maxchanges
        if ~args.silent
          fprintf('Closest sequence has >%d changes\n', args.maxchanges);
        end
        descr='';
        bestnaseq=[];
        return;
      end
      % Convert to a descriptive string of mutations going from bestalign(3,:) to bestalign(1,:)
      descr='';
      for i=1:size(bestalign,2)
        if bestalign(2,i)~='|'
          descr=sprintf('%s%c%d%c ',descr,bestalign(3,i),i,bestalign(1,i));
        end
      end
      if length(descr)>0
        descr=descr(1:end-1);  % Remove trailing space
      end
      
      if args.silent || isempty(bestnaseq)
        return;
      end
      
      fprintf('Best matching seq has %d changes:\n',bestchanges);
      if ismember(naseq,data.naseq)
        s=data.naseqformat(naseq);
      else
        % Not loaded so can't format; just use raw seq
        s=seq;
      end
      % Restore and '-' in s that were in bestalign(1,:)
      dashes=find(bestalign(1,:)=='-');
      for i=1:length(dashes)
        d=dashes(i);
        lpos=find(s~=' ');
        if d<=length(lpos)   % Dash at end would be a problem
          s=[s(1:lpos(d)-1),'-',s(lpos(d):end)];
        else
          s=[s,'-'];
        end
      end
      al=blanks(length(s));
      al(2,:)=al(1,:);
      al(3,:)=al(1,:);
      al(:,s~=' ')=bestalign;
      fprintf('%9d %s %s %s\n',naseq,al(1,:),descr,data.getlabels(naseq));
      fprintf('          %s\n',al(2,:));
      fprintf('%9d %s %s\n',bestnaseq,al(3,:),data.getlabels(bestnaseq));
    end
    
    % Getters - all of form getxx()
    %  Take an array of naseqs as an arg (not a selector)
    %  Versions without the get prefix act on a selector
    function res=getseq(data,naseqs)
      inds=data.naseq2ind(naseqs);
      res=data.seq(inds);
    end
    
    % Force db load if not found
    function res=getorloadseq(data,naseqs)
      if isempty(data.otherseqs)
        data.otherseqs=containers.Map('KeyType','double','ValueType','char');
      end
      res=cell(length(naseqs),1);
      exists=ismember(naseqs,data.naseq);
      res(exists)=data.getseq(naseqs(exists));
      if any(~exists)
        for i=1:length(naseqs)
          if ~exists(i) & isKey(data.otherseqs,naseqs(i))
            res{i}=data.otherseqs(naseqs(i));
            exists(i)=true;
          end
        end
        if any(~exists)
          load=naseqs(~exists);
          for k=1:1000:length(load)
            cmd='SELECT naseq,seq FROM ngs.naseqs WHERE naseq IN (';
            nl=0;
            for i=k:min(k+999,length(load))
              cmd=sprintf('%s%d,',cmd,load(i));
              nl=nl+1;
            end
            cmd(end)=')';
            fprintf('Loading %d seqs from db (k=%d)...',nl,k);
            [nn,seq]=mysql(cmd);
            fprintf('done\n');
            for i=1:length(nn)
              ind=find(naseqs==nn(i));
              for k=1:length(ind)
                res{ind(k)}=seq{i};
              end
              data.otherseqs(nn(i))=seq{i};
            end
          end
        end
      end
    end
    
    function l1=getloop1(data,naseqs)
      l1=data.loop1(data.naseq2ind(naseqs));
    end
    
    function l1=loop1(data,sel)
      if nargin<2
        l1=data.segments(:,4);
      else
        l1=data.segments(sel,4);
      end
    end
    
    function n=getseglen(data,segnum,naseqs)
      sel=data.naseq2ind(naseqs);
      n=data.seglen(segnum,sel);
    end
      
    function n=seglen(data,segnum,sel)
      if nargin<3
        n=data.segextents(:,segnum,2)-data.segextents(:,segnum,1)+1;
      else
        n=data.segextents(sel,segnum,2)-data.segextents(sel,segnum,1)+1;
      end
    end
      
    function n=loop1len(data,sel)
      n=data.seglen(4,sel);
    end

    function n=helix1len(data,sel)
      n=data.seglen(3,sel);
    end

    function n=getloop1len(data,naseqs)
      n=data.loop1len(data.naseq2ind(naseqs));
    end
    
    function n=gethelix1len(data,naseqs)
      n=data.helix1len(data.naseq2ind(naseqs));
    end
    
    function l2=loop2(data,sel)
      if nargin<2
        l2=data.segments(:,8);
      else
        l2=data.segments(sel,8);
      end
    end
    
    function l2=getloop2(data,naseqs)
      l2=data.loop2(data.naseq2ind(naseqs));
    end
    
    function n=loop2len(data,sel)
      n=data.seglen(8,sel);
    end
    
    function n=helix2len(data,sel)
      n=data.seglen(7,sel);
    end
    
    function n=getloop2len(data,naseqs)
      n=data.loop2len(data.naseq2ind(naseqs));
    end
    
    function n=gethelix2len(data,naseqs)
      n=data.helix2len(data.naseq2ind(naseqs));
    end
    
    function n=getlength(data,naseqs)
      n=double(data.length(data.naseq2ind(naseqs)));
    end
      
    function c=getclass(data,naseqs)
      c=data.class(data.naseq2ind(naseqs));
    end
      
    function isref=getisref(data,naseqs)
      isref=false(size(naseqs));
      for j=1:length(naseqs)
        isref(j)=any([data.tags.naseq]==naseqs(j) & strcmp({data.tags.name},'ref'));
      end
    end

    function i=getsegments(data,naseqs)
      i=data.segments(data.naseq2ind(naseqs),:);
    end

    function i=getsegextents(data,naseqs)
      i=squeeze(data.segextents(data.naseq2ind(naseqs),:,:));
    end

    function id=getsegmentid(data,segnum,naseqs)
    % Get index into pts of segment
    % Return 0 if not matched
    % Only valid for segments (2,3,5,6,7,9) (which are pts(1:6))
      sel=data.naseq2ind(naseqs);
      segs=data.segments(sel,segnum);
      pts=NGSSeq.pts{segnum};
      fnd=~cellfun('isempty',segs);
      id=zeros(length(fnd),1);
      [~,id(fnd)]=ismember(segs(fnd),pts);
    end
    
      
    function id=hascore(data,naseqs)
      id=getsegmentid(data,6,data.naseq2ind(naseqs));
    end

    function merge(data,data2)
    % Append data2 information to this one
      numnaseq=size(data.naseq,1);
      [data.naseq,ia,ib]=union(data.naseq,data2.naseq,'stable');
      fn=setdiff(fieldnames(data),{'seqclasses','pts','naseq'});
      for i=1:length(fn)
        if size(data.(fn{i}),1)==numnaseq && size(data2.(fn{i}),1)==size(data2.naseq,1)
          %fprintf('Merging %d %s\n',length(data.naseq)-numnaseq,fn{i});
          data.(fn{i})(numnaseq+1:length(data.naseq),:,:)=data2.(fn{i})(ib,:,:);
        else
          ;%fprintf('Skipping %d %s\n',size(data2.(fn{i}),1),fn{i});
        end
      end
      % If starting with an empty struct, copy other fields in
      fn={'seqnames','tags','clusters','clusterdists'};
      for i=1:length(fn)
        if isempty(data.(fn{i}))
          data.(fn{i})=data2.(fn{i});
        end
      end
      fprintf('NGSSeq:merge: Note that the merge product will still access the database as %s\n',data.database);
    end
      
    function order(data,naseqs,varargin)
    % Order oligo in 2 parts
      defaults=struct('overlap',24,'prefix','W','parts',2);
      args=processargs(defaults,varargin);

      for naseq=naseqs
        ind=find(data.naseq==naseq,1);
        if length(ind) ~= 1
          error('Naseq %d not found', naseq);
        end
        if args.prefix=='W'
          prefix='GGGAAACAAACAAA';
        elseif args.prefix=='Z'
          prefix='GGGACAAAACAAAA';
        else
          error('Unknown prefix: %s\n',args.prefix);
        end
        X='TTTTATTTTTCTTTTT';
        seq=[prefix,data.seq{ind},X];
        fmtseq=[prefix,data.naseqformat(data.naseq(ind)),X];
        totallen=length(seq);
        if args.parts==1
          fprintf('%s_%d_X\t%3d\t%s\t%8.8s\n',args.prefix,naseq,length(seq),fmtseq,string2hash(seq));
        elseif args.parts==2
          p2=floor((totallen-args.overlap)/2);
          if (totallen-p2+1)>60 && (totallen-p2+1)<70
            % Adjust so len2<=60
            p2=totallen-60+1;
          end
          while seq(p2)~='C' && seq(p2)~='G'
            p2=p2+1;
          end
          p1=p2+args.overlap-1;
          while seq(p1)~='C' && seq(p1)~='G'
            p1=p1+1;
          end
          seq1=seq(1:p1);
          seq2=rc(seq(p2:end));
          overlap=seq(p2:p1);
          fprintf('%s_%d_Left\t\t\t\t\t\t\t\t\t\t\t\t%3d\t%s\t%8.8s\n',args.prefix,naseq,length(seq1),seq1,string2hash(seq1));
          fprintf('X_%d_Right-RC\t\t\t\t\t\t\t\t\t\t\t\t%3d\t%s\t%8.8s\n',naseq,length(seq2),seq2,string2hash(seq2));
          %fprintf('Overlap\t         \t%d\t%s\n', length(overlap),overlap);
          %fprintf('Product\t%d\t%s\n', length(seq), fmtseq);
        end
      end
    end
    
    function classifyseqs(data)
    % Classify the sequences 
      [naseqs,ua,ub]=unique(data.naseq);
      seqs=data.seq(ua);
      len=data.length(ua);
      l1len=data.loop1len(ua);
      l2len=data.loop2len(ua);
      s1len=data.helix1len(ua);
      s2len=data.helix2len(ua);
      % Compensate for shortened/extended stem helices
      adjl1len=max(l1len+(s1len-6)*2,0);
      adjl1len(isnan(l1len))=nan;
      adjl2len=max(l2len+(s2len-4)*2,0);
      adjl2len(isnan(l2len))=nan;

      loop1=data.loop1(ua);
      class=nan(length(ua),1);
      for i=1:length(data.seqclasses)
        sc=data.seqclasses{i};
        nm=sc{1};
        sel=isnan(class);
        if ~isempty(sc{2})
          sel=sel & len>=sc{2}(1) & len <= sc{2}(2);
        end
        if ~isempty(sc{3})
          if ischar(sc{3})
            % Regexp for loop1
            l1=sc{3};
            re=strrep(l1,'R','(G|A)');
            re=strrep(re,'Y','(T|C)');
            re=strrep(re,'D','(G|A|T)');
            re=strrep(re,'H','(A|T|C)');
            re=strrep(re,'N','.');
            re=strrep(re,' ','');
            re=['^',re,'$'];
            fprintf('Searching in %d sequences for %s: l1=(%d) %s\n', length(sel),nm,length(l1(l1~=' ')), l1);
            m=regexp(loop1(sel),re);
            sel2=cellfun(@(z) ~isempty(z), m);
            sel(sel)=sel2;
          else
            sel=sel & ((l1len>=sc{3}(1) & l1len<= sc{3}(2)) | (adjl1len>=sc{3}(1) & adjl1len<= sc{3}(2)));
          end
        end
        if ~isempty(sc{4})
          if ischar(sc{4})
            % Regexp for loop2
            assert(false);  % Not implemented
          else
            sel=sel & ( (l2len>=sc{4}(1) & l2len <= sc{4}(2)) |  (adjl2len>=sc{4}(1) & adjl2len <= sc{4}(2)) );
          end
        end
        fprintf('Identified %d sequences as %s\n', sum(sel), nm);
        class(sel)=i;
      end
      class(isnan(class))=length(data.seqclasses);
      assert(length(data.seqclasses) < 256);   % So we can store it in a uint8
      data.class=uint8(class(ub));
    end
    
    function f=getribofolding(obj,naseq)
    % Get the folding of naseq that follows the expected canonic structure
    % Assume the following folding by segment:
    %    unfolded:  pre (1), core(6), s3b(10)(first part), post(11)
    %    helix:  s3a(2):s3b(10)(last part), s1a(3):s1b(5), s2a(7):s2b(9)
    %    fold independently:  loop1(4), loop2(8)
      segs=obj.getsegments(naseq);
      assert(length(segs{2})-1 == length(segs{10})-3);  % segments 2 and 9 contain a bit more than the helix
      loop1fold=fold(struct('name','loop1','seq',strrep([segs{3},segs{4},segs{5}],'T','U')),1);
      l1f=loop1fold.ufoldings(1).folding(length(segs{3})+1:end-length(segs{5}));
      % Could have some single-stranded in planned helix
      while sum(l1f==')') > sum(l1f=='(')
        l1f(find(l1f==')',1,'last'))='.';
      end
      while sum(l1f==')') < sum(l1f=='(')
        l1f(find(l1f=='(',1))='.';
      end
      loop2fold=fold(struct('name','loop1','seq',strrep([segs{7},segs{8},segs{9}],'T','U')),1);
      l2f=loop2fold.ufoldings(1).folding(length(segs{7})+1:end-length(segs{9}));
      % Could have some single-stranded in planned helix
      while sum(l2f==')') > sum(l2f=='(')
        l2f(find(l2f==')',1,'last'))='.';
      end
      while sum(l2f==')') < sum(l2f=='(')
        l2f(find(l2f=='(',1))='.';
      end
      f=[repmat('.',1,length(segs{1})), ...
         repmat('(',1,length(segs{2})-1), ...
         '.', ...
         repmat('(',1,length(segs{3})), ...
         l1f, ...
         repmat(')',1,length(segs{5})), ...
         repmat('.',1,length(segs{6})), ...
         repmat('(',1,length(segs{7})), ...
         l2f, ...
         repmat(')',1,length(segs{9})), ...
         repmat('.',1,3), ...
         repmat(')',1,length(segs{10})-3), ...
         repmat('.',1,length(segs{11})), ...
         ];
      seq=obj.getseq(naseq); seq=seq{1};
      fprintf('%s\n%s\n',seq,f);
      assert(length(f)==obj.getlength(naseq));
    end
    
    function img=pcolorstruct(obj,naseq,x,varargin)
     % Draw structure with colormapping of each nucleotide controlled by x
      defaults=struct('colormap',colormap,'caxis',[]);
      args=processargs(defaults,varargin);
      
      if isempty(args.caxis)
        args.caxis=[nanmin(x),nanmax(x)];
      end
      
      % Make sure x is 1xn
      x=x(:)';
      
      seq=obj.getseq(naseq);seq=seq{1};
      assert(length(seq)==length(x));
      
      useq=strrep(seq,'T','U');
      folding=obj.getribofolding(naseq);

      % Scale X to 0..1 using args.caxis
      sx=(x-args.caxis(1))/(args.caxis(2)-args.caxis(1));
      
      % Convert x to a 1xNx3 RGB image using the colormap
      rgb=zeros(1,length(x),3);
      for i=1:length(x)
        rgb(1,i,:)=interp1(0:1/(size(args.colormap,1)-1):1,args.colormap,sx(i));
      end
      nstyles=49;
      % Convert rgb to indexed color using nstyles quantization levels
      [ix,smap]=rgb2ind(rgb,nstyles,'nodither');
      ix=ix+1;
      
      cargs={};
      % Make styles and apply themn
      for i=1:size(smap,1)
        if any(ix==i & isfinite(x))
          cargs{end+1}=sprintf('basesStyle%d',i);
          cargs{end+1}=['fill=#',sprintf('%02x',round(smap(i,:)*255))];
          cargs{end+1}=sprintf('applyBasesStyle%don',i);
          cargs{end+1}=strjoin(arrayfun(@(z) sprintf('%d',z), find(ix==i & isfinite(x)),'Unif',false),',');
        end
      end
      res=3;
      img=varna(useq,folding,'resolution',res,cargs{:});
    end
    
    function dist=getlocalradius(data, fold, p1, p2)
    % Find stretched-out distance between positions p1 and p2 on fold
      debug=false;
      baselength=3.4e-10;   % meters along backbone
      helixlength=15.3e-10;
      f=rnaconvert(fold);
      f=f+f';
      dist=0;
      pos=blanks(length(fold));
      pos(p1)='S';
      while (p1~=p2)
        hjump=find(f(p1,:));
        if ~isempty(hjump) && abs(hjump-p2)<abs(p1-p2)
          p1=hjump;
          pos(p1)='H';
          dist=dist+helixlength;
        else
          p1=p1+(p2-p1)/abs(p2-p1);
          pos(p1)='.';
          dist=dist+baselength;
        end
      end
      pos(p2)='E';
      if debug
        fprintf('%s %d %d\n', fold, p1, p2);
        fprintf('%s %.1f\n', pos, dist);
      end
    end
    
    function f=showstructs(data, naseq, varargin)
      defaults=struct('maxfoldings',6,'prefix',[],'stop',[],'mincoverage',0.6,'label',[],'color',true,'stochastic',true,'nstoch',100,'maxhelixdiffs',2,'ligate',false,'cut',false,'C9T',true,'plot',true,'verbose',false);
      args=processargs(defaults,varargin);


      B='CGGAAATTTCAAAGGTGCTTC';
      A='CTTTTCCGTATATCTCGCCAG';
      W='AAACAAACAAA';
      Z='ACAAAACAAAAC';
      X='AAAAAGAAAAATAAAAA';
      if args.C9T
        % New C9T stops (2018-2018)
        AStop='GACAGCCTGGCGAGATATACGGAAAAGAGGCUGTCACTGGAXUTTTCTTTTT';  % BT2018 with trailing overlap to ribo removed
        WStop='GACAGCTTTGTTTGTTTCCCAAGCUGTCACTGGAXUTTTCTTTTT'; % BT2018 with trailing overlap to ribo removed
        ZStop='GACAGCGTTTTGTTTTGTCCCACGCUGTCACTGGAXUTTTCTTTTT';  % BT2019 with trailing overlap to ribo removed
      else
        AStop='GACAGCCTGGCGAGATATACGGAAAAGAGGCUGTCACCGGATCCGGTCTGATGAGUCCTTTCTTTTT';  % BT1305 with trailing overlap to ribo removed
        ZStop='GACAGCGTTTTGTTTTGTCCCACGCUGTCACCGGATCCGGTCTGATGAGUCCTTTCTTTTT'; % BT1508 with trailing overlap to ribo removed
        WStop='GACAGCTTTGTTTGTTTCCCAAGCUGTCACCGGATCCGGTCTGATGAGUCCTTTCTTTTT'; % BT1316 with trailing overlap to ribo removed
      end
      seq=data.getorloadseq(naseq); %unique(data.seq(data.naseq==naseq));
      seq=seq{1};
      % Map to regions ids
      % From segmentation have 1:11 (front extra,s1a+c12,s2a,l2,s2b,core,s3a,l3,s3b,c31+s1b,back extra)
      % -> region 1,2,3,4,3,5,6,7,6,2,1
      % After that we have:
      % 8-GGG prefix
      % 9-prefix
      % 10-X suffix
      % 11-stop (other than prefix/suffix parts)
      % 12-Uracil cut point
      regions=ones(size(seq));  % Regions for color-mapping
      mapsegs=[1,2,3,4,3,5,6,7,6,2,1];
      if ismember(naseq,data.naseq)
        ext=data.getsegextents(naseq);
        for i=1:size(ext,1)
          if all(isfinite(ext(i,1:2)))
            regions(ext(i,1):ext(i,2))=mapsegs(i);
          end
        end
      end
          
      rnaparams='RNA37';

      if isempty(args.prefix)
        name=sprintf('%d',naseq);
      else
        name=sprintf('%s_%d_X',args.prefix,naseq);
        if strcmpi(args.prefix,'clvd')
          if ~strncmp(seq,'GCTGTC',6)
            error('Unable to find cleaved sequence for %s, does not start with GCTGTC',seq);
          end
          seq=seq(7:end);
          regions=regions(7:end);
          seq=[seq,X];
          regions=[regions,repmat(10,1,length(X))];
        else
          if strcmp(args.prefix,'A')
            prefix=A;
          elseif strcmp(args.prefix,'B')
            prefix=B;
          elseif strcmp(args.prefix,'W')
            prefix=W;
          elseif strcmp(args.prefix,'Z')
            prefix=Z;
          else
            error('Invalid prefix: ',args.prefix);
          end
          seq=['GGG',prefix,seq,X];
          regions=[repmat(8,1,3),repmat(9,1,length(prefix)),regions,repmat(10,1,length(X))];
        end
      end

      if ~isempty(args.stop)
        if isempty(args.prefix)
          error('Not valid to apply stop without a prefix');
        end
        if args.stop=='A'
          stop=AStop;
          stopprefix=A;
        elseif args.stop=='W'
          stop=WStop;
          stopprefix=W;
        elseif args.stop=='Z'
          stop=ZStop;
          stopprefix=Z;
        else
          error('Unsupported stop: %s\n', args.stop);
        end
        name=[args.stop,'Stop_',name];
        % Find overlap between stop and RNA
        stoppos=strfind(rc(seq),stop(end-8:end));
        seq=seq(1:end-stoppos-8);
        regions=regions(1:end-stoppos-8);
        
        if ~strcmp(seq(end-11:end),'GGACGAAACAGC')
          fprintf('*** Warning sequence ends with %s instead of GGAC GAA ACAGC',seq(end-11:endstr));
        end
        seq=[stop,rc(seq)];
        regions=[repmat(11,1,length(stop)),regions(end:-1:1)];
        regions(seq=='U')=12;
        % Change prefix region
        pr=strfind(seq,rc(stopprefix));
        assert(~isempty(pr));
        regions(pr:pr+length(stopprefix)-1)=9;
        if strcmp(seq(pr+length(stopprefix):pr+length(stopprefix)+2),'CCC')
          % GGG before prefix
          regions(pr+length(stopprefix):pr+length(stopprefix)+2)=8;
        end
        % And suffix region
        assert(strcmp(rc(seq(length(stop)-8:length(stop))),X(1:9)));
        regions(length(stop)-8:length(stop))=10;
        % And s3a region
        regions(1:6)=2;
        rnaparams='dna';   % Since we're adding the stop, this must be DNA, not RNA
      end
      if args.ligate
        if ~args.cut
          error('Cannot ligate if not cutting after');
        end
        if isempty(args.stop)
          error('Cannot ligate and cut if no stop specified');
        end
        keep=[find(seq=='U',1,'last')+1:length(seq),1:find(seq=='U',1)-1];
        seq=seq(keep);
        regions=regions(keep);
        name=[name,'_lig_cut'];
      elseif args.cut
        keep=[find(seq=='U',1,'last')+1:length(seq)];
        seq=seq(keep);
        regions=regions(keep);
        name=[name,'_cut'];
      end
      
      assert(length(regions)==length(seq));
      if args.verbose
        fprintf('%s\n%s\n%s\n',seq,sprintf('%1d',floor(regions/10)),sprintf('%1d',mod(regions,10)));
      end
      
      if ~args.stochastic && any(seq=='X')
        fprintf('Removing spacer from seq to allow folding\n');
        regions=regions(seq~='X');
        seq=seq(seq~='X');
      end
      
      f=struct('name',name,'seq',seq,'regions',regions,'rnaparams',rnaparams,'stochastic',args.stochastic);

      if ~args.stochastic
        f.maxFoldings = args.maxfoldings;
        f.minCoverage = args.mincoverage;
        f=fold(f,1);
        if ~isempty(args.label)
          f.ufoldings(1).label=args.label;
          for i=2:length(f.ufoldings)
            f.ufoldings(i).label=args.label;
          end
        end
        for i=1:length(f.ufoldings)
          f.ufoldings(i).isactive=true;   % Lock down so the label colors don't change
        end
      else
        % Use stochastic folding
        f.maxhelixdiffs=args.maxhelixdiffs;
        if strcmp(f.rnaparams,'dna')
          f.viennaparams='dna_mathews2004';
        else
          f.viennaparams='rna_langdon2018';
        end
        % Get stochastic folds
        stochfolds=stochastic(seq,args.nstoch,'vienna','viennaparams',f.viennaparams,'dangles',2,'material',f.rnaparams,'verbose',args.verbose);
        f.stochfolds={stochfolds.folding};
        [u,c]=groupstochastic(f.stochfolds,'maxhelixdiffs',args.maxhelixdiffs,'verbose',args.verbose);   % Group ones within 2 helices difference
        f.ufoldings=struct('folding',u,'prob',num2cell(c/args.nstoch),'dG',nan);
      end

      % Compute some stats
      if ~isempty(args.stop) && ~args.ligate && ~args.cut
        % Compute ligation distances
        for i=1:length(f.ufoldings)
          f.ufoldings(i).localradius=data.getlocalradius(f.ufoldings(i).folding,1,length(f.seq));
          vol=(4/3*pi*(f.ufoldings(i).localradius)^3)*1000;   % Volume in liters
          f.ufoldings(i).ligconc=1/6.022e23 / vol;  % in M
        end
        f.ligconc=mean([f.ufoldings.ligconc]);
      end

      if args.plot
        cargs={};
        % Make styles and apply themn
        ureg=unique(regions);
        assert(max(ureg)<=12);
        smap=colorcube(12)*0.5+0.5;   % Lightened colorcube
        for i=1:length(ureg)
          r=ureg(i);
          cargs{end+1}=sprintf('basesStyle%d',r);
          cargs{end+1}=['fill=#',sprintf('%02x',round(smap(r,:)*255))];
          cargs{end+1}=sprintf('applyBasesStyle%don',r);
          cargs{end+1}=strjoin(arrayfun(@(z) sprintf('%d',z), find(regions==r),'Unif',false),',');
        end
        f=plotstructs(f,args.maxfoldings,false,cargs);
      end
    end

    function name=getname(obj,naseq)
    % Return name of given sequence, or empty string if not found
      sel=find(naseq==[obj.seqnames.naseq]);
      if isempty(sel)
        name='';
      else
        name=obj.seqnames(sel).name;
      end
    end
    
    function naseq=findname(obj,name)
    % Get the naseq with given name or [] if not found
      sel=strcmp(name,{obj.seqnames.name});
      if ~any(sel)
        naseq=[];
      else
        naseq=obj.seqnames(sel).naseq;
      end
    end

    function tags=gettags(obj,naseq)
    % Return tags of given sequence, or empty string if not found
      sel=find(naseq==[obj.tags.naseq]);
      t=obj.tags(sel);
      [~,ord]=sort(arrayfun(@(z) [z.name,'.',z.value],t,'Unif',false));
      t=t(ord);
      tags='';
      for i=1:length(t)
        if isempty(t(i).value)
          tags=[tags,sprintf(' +%s',t(i).name)];
        else
          if i>1 && strcmp(t(i).name,t(i-1).name)
            tags=[tags,sprintf(',%s',t(i).value)];
          else
            tags=[tags,sprintf(' %s=%s',t(i).name,t(i).value)];
          end
        end
      end
      if length(tags)>0
        tags=tags(2:end);  % Get rid of leading space
      end
    end
    
    function [naseq,values]=findbytag(obj,name,value)
      sel=strcmp({obj.tags.name},name);
      if nargin>2
        sel=sel&strcmp({obj.tags.value},value);
      end
      naseq=[obj.tags(sel).naseq];
      if nargout>=2
        values={obj.tags(sel).value};
      end
    end
    
    function is=isbistable(obj,naseq)
      is=false;
      seq=obj.getorloadseq(naseq);
      if ~isempty(seq)
        %seq=seq{1};
        is=cellfun(@(z) ~isempty(z),regexp(seq,'GAA[AG][CT][AG]G[GCT].*GAA[AG][CT][AG]G[CT]'));
      end
    end
    
    function [is,matches]=istruncation(obj,naseq)
    % Is this sequence a truncation of the right end of another sequence that is present
      seq=obj.getorloadseq(naseq);
      is=false(length(naseq),1);
      matches=cell(length(naseq),1);
      for i=1:length(naseq)
        ind=obj.naseq==naseq(i);
        if ~any(ind)
          continue;   % Skip unloaded sequences
        end
        sel=find(obj.length>=obj.length(ind)+12);
        match=strncmp(seq{i},obj.seq(sel),length(seq{i})-12);  % Ignore s2b to end
        matches{i}=obj.naseq(sel(match));
        is(i)=any(match);
      end
    end
    
    function flags=getlabels(obj,naseq,varargin)
    % Return labels for given naseq or empty string
      defaults=struct('includename',true,'includetags',true,'includemuts',true,'includeclusters',true,'includebistable',true,'includehits',true,'includetruncations',false,'hitdetails',false);
      args=processargs(defaults,varargin);
      flags='';

      if args.includename
        name=obj.getname(naseq);
        if ~isempty(name)
          flags=[flags,' ',name];
        end
      end
     
      if args.includeclusters
        c=obj.getcluster(naseq);
        if ~isempty(c)
          flags=[flags,' ',c];
        end
      end
      
      if args.includehits
        sel=find([obj.hits.naseq]==naseq);
        if isempty(sel)
          ;
        elseif args.hitdetails || length(sel)==1
          flags=[flags,' Hits['];
          for i=1:length(sel)
            flags=sprintf('%s%s=%.2f ',flags,obj.hits(sel(i)).name,obj.hits(sel(i)).fold);
          end
          flags(end)=']';
        else
          [maxfold,ind]=max([obj.hits(sel).fold]);
          flags=sprintf('%s Hits[%s=%.2f+%d]',flags,obj.hits(sel(ind)).name,maxfold,length(sel)-1);
        end
      end
      
      if args.includetags
        tags=obj.gettags(naseq);
        if ~isempty(tags)
          flags=[flags,' ',tags];
        end
      end

      if args.includebistable && obj.isbistable(naseq)
        flags=[flags,' bistable=S3'];
      end

      if args.includetruncations && obj.istruncation(naseq)
        flags=[flags,' trunc'];
      end

      if length(flags)>0 && flags(1)==' '
        flags=flags(2:end);  % Get rid of leading space
      end
    end

    function c=loopinteractions(obj,naseq,varargin)
    % Find any possible loop interactions between the loops
      defaults=struct('conc',1e-6,'cutoff',1e-6,'verbose',false,'sodium',0.05,'mg',0.0005,'temp',37,'prefix','none','maxsize',2,'dangles','some','material','rna','pairthresh',0.01);
      args=processargs(defaults,varargin);

      ind=find(obj.naseq==naseq,1);
      if length(ind)~=1
        error('Naseq %d not found in data\n', naseq);
      end
      seq=obj.seq{ind};
      suffix='AAAAAGAAAAATAAAA';
      if strcmp(args.prefix,'W')
        prefix='GGGAAACAAACAAA';
      elseif strcmp(args.prefix,'Z')
        prefix='GGGACAAAACAAAAC';
      elseif strcmp(args.prefix,'none')
        prefix='';
        suffix='';
      else
        error('Unrecognized prefix code: %s\n', args.prefix);
      end
      seq=[prefix,seq,suffix];
      seq(seq=='T')='U';
      c=complexes({seq},'temp',args.temp,'concentrations',args.conc*[1],'material',args.material,'pairs',true,'sodium',args.sodium,'mg',args.mg,'cutoff',args.cutoff,'verbose',args.verbose,'maxsize',args.maxsize,'dangles',args.dangles);
      fprintf('Interactions between strands at %s\n',concfmt(args.conc));
      for i=1:length(c.ocomplex)
        oc=c.ocomplex(i);
        % Compute base pair probabilities for this complex
        p=nu_pairs({seq},oc.perm,'temp',args.temp,'cutoff',args.cutoff,'verbose',args.verbose,'sodium',args.sodium,'mg',args.mg,'material',args.material,'dangles',args.dangles,'maxsize',args.maxsize);
        c.ocomplex(i).pairs=p;	% Keep for reference
        oc=c.ocomplex(i);

        frac=oc.eqconc*sum(oc.strands)/args.conc;
        fprintf('%s%5.1f %.3g%%\n', sprintf('%d ',oc.strands), oc.dG, frac*100);
        ti1=sprintf('%d f(%d strands)',naseq, sum(oc.strands));
        ti=sprintf('%s=%.3g%% at %s',ti1, frac*100, concfmt(args.conc));
        setfig(ti1);clf;
        p=oc.pairs.pairfrac;
        n=length(seq);
        self=sum(p(1:n,1:n),2);
        unpaired=p(1:n,end);
        if oc.strands==1
          h=bar([self,unpaired],'stacked');
          leg={'Self','Unpaired','Loop1','Loop2'};
        elseif oc.strands==2
          cross=sum(p(1:n,n+1:end-1),2);
          h=bar([self,unpaired,cross],'stacked');
          leg={'Self','Unpaired','Cross','Loop1','Loop2'};
        elseif oc.strands==3
          cross=sum(p(1:n,n+1:2*n),2);
          cross2=sum(p(1:n,2*n+1:end-1),2);
          h=bar([self,unpaired,cross,cross2],'stacked');
          leg={'Self','Unpaired','Cross1','Cross2','Loop1','Loop2'};
        else
          error('Unexpected number of strands: %d',oc.strands);
        end
        h(2).FaceColor=[1,1,1];
        hold on;

        % Find most likely matched nucleotide
        matches=zeros(1,n);
        for i=1:n
          % Draw lines connecting same-strand base pairs
          for j=i+1:n
            psum=p(i,j);
            if psum>args.pairthresh
              obj.arc([i,(i+j)/2,j],[0,(j-i)/n*0.4,0],'y',psum*2);
            end
          end

          if oc.strands>=2
            % Draw lines connecting cross-strand base pairs
            for j=n+1:2*n
              psum=sum(p(i,j:n:end-1));   % Sum to any other strand
              if psum>args.pairthresh
                plot([i,j-n],[1.0,0.0],'-c','LineWidth',psum*2);
              end
            end

            [maxpair,pos]=max(p(i,n+1:end-1));
            pos=mod(pos-1,n)+1;   % Could be 3rd or higher strand
            if maxpair>args.pairthresh
              matches(i)=pos;
              if i>1 && matches(i-1)-1 ~= matches(i)  % New helix
                ht=text(i-0.5,1.02,'|','HorizontalAlignment','center');
              end
              text(i,1.02,seq(pos),'HorizontalAlignment','center');
            end
          end
        end
        
        h(end+1)=plot(length(prefix)+squeeze(obj.segextents(ind,4,:)),-0.05*[1,1],'-g','LineWidth',2,'Clipping','off');
        h(end+1)=plot(length(prefix)+squeeze(obj.segextents(ind,8,:)),-0.05*[1,1],'-r','LineWidth',2,'Clipping','off');
        legend(h,leg);
        a=axis;
        a(4)=1.04;
        a(3)=0;
        axis(a);
        set(gca,'XTick',1:n);
        set(gca,'XTickLabel',arrayfun(@(z) z, seq,'UniformOutput',false));
        set(gca,'TickLength',[0,0]);
        ylabel('Base-pairing probability');
        title(ti);
      end
    end
    
    function loopinteractions1(obj,naseq)
    % Find any possible loop interactions between the loops
      ind=find(obj.naseq==naseq,1);
      if length(ind)~=1
        error('Naseq %d not found in data\n', naseq);
      end
      l1=obj.loop1(ind);
      l2=obj.loop2(ind);
      s=[l1{1},'*',l2{1}];  % Concatenated sequence
      s=obj.naseqformat(naseq);
      s=s(s~=' ');
      s=strrep(s,'T','U');
      m=nan(length(s));
      for i=1:length(s)
        for j=1:length(s)
          m(i,j)=(s(i)=='A'&&s(j)=='U')||(s(i)=='C'&&s(j)=='G')||(s(i)=='G'&&(s(j)=='U'||s(j)=='C'))||(s(i)=='U'&&(s(j)=='G'||s(j)=='A'));
        end
      end
      for i=2:size(m,1)
        sel=find(m(i,1:end-1)==1);
        m(i,sel)=1+m(i-1,sel+1);
      end
      maxhelix=max(m(:));
      fprintf('Maximum helix = %d\n', maxhelix);
      [x,y]=find(m==maxhelix);
      src=s(end:-1:1);
      for k=1:length(x)
        i=x(k);
        j=y(k);
        fprintf('s[%d:%d]=%s, s[%d:%d]=%s\n', i-maxhelix+1,i,s(i+(-maxhelix+1:0)),j,j+maxhelix-1,s(j+(0:maxhelix-1)));
        nb=(length(s)-j+1)-i;
        nb1=max(nb,0);
        nb2=max(-nb,0);
        fprintf('%s%s\n',blanks(nb1),s);
        fprintf('%s%s\n',blanks(nb1+i-maxhelix),repmat('|',1,maxhelix));
        fprintf('%s%s\n',blanks(nb2),src);
      end
    end

    
    function [r_naseq,r_seqstr]=similarseqs(obj,naseq,varargin) 
    % Find sequences similar to given one
    % Return list of similar sequences and their sequence (with differences shown)
      defaults=struct('maxerrors',1,'usecache',false,'sel',[],'print',false,'ignore',[]);
      % If sel is set use only the selected entries
      % If usecache is set, reuses a cell2mat conversion of sequences (which is slow)
      args=processargs(defaults,varargin);

      global lastlength lastc2mat lastnaseq    % Cache the cell2mat conversion
      
      if ischar(naseq)
        %fprintf('similarseqs(%s)\n',naseq);
        fmt=naseq;
        seq=upper(naseq(naseq~=' ' & naseq~='-'));
      else
        ind=find(obj.naseq==naseq,1);
        if isempty(ind)
          error('NASeq %d not found',naseq);
        end
        seq=obj.seq{ind};
        fmt=obj.naseqformat(obj.naseq(ind));
      end
      % Find substitutions
      if ~isempty(args.sel)
        sel=args.sel & obj.length==length(seq);
      else
        sel=obj.length==length(seq);
      end
      if ~isempty(args.ignore)
        %fprintf('Ignoring %s\n',sprintf('%d ',args.ignore));
        sel=sel&~ismember(obj.naseq,args.ignore);
      end
      % Cache the cell2mat conversion
      if args.usecache && length(seq)==lastlength
        ss=lastc2mat;
        naseq=lastnaseq;
      else
        ss=cell2mat(obj.seq(sel));
        naseq=obj.naseq(sel);
        lastlength=length(seq);
        lastc2mat=ss;
        lastnaseq=naseq;
      end
      errs=zeros(length(naseq),1);
      ferrs=true(size(errs));
      for i=1:size(ss,2)
        if seq(i)~='N'
          errs=errs+(ss(:,i)~=seq(i));
        end
        if i==20
          % Prune it out
          ferrs=errs<=args.maxerrors;
          %fprintf('Reduced from %d to %d\n', length(ferrs), sum(ferrs));
          ss=ss(ferrs,:);
          naseq=naseq(ferrs);
          errs=errs(ferrs);
        end
      end
      ferrs=errs<=args.maxerrors;
      ss=ss(ferrs,:);
      naseq=naseq(ferrs);
      errs=errs(ferrs);
      r_naseq=[];
      r_seqstr={};
      for j=1:size(ss,1)
        d=lower(ss(j,:));
        d(d~=lower(seq))=upper(d(d~=lower(seq)));
        d2=fmt;
        d2(d2~=' ' & d2~='-')=d;
        r_naseq(end+1)=naseq(j);
        r_seqstr{end+1}=d2;
        if args.print
          fprintf('%9d %s %d %s\n',naseq(j), d2,errs(j),obj.getlabels(naseq(j)));
        end
      end
      % Order them 
      [r_seqstr,ord]=sort(r_seqstr);
      r_naseq=r_naseq(ord);
      if args.maxerrors>0
        % Deletions
        % if args.print
        %   fprintf('Check for deletion and %d other error:\n',args.maxerrors-1);
        % end
        for i=1:length(fmt)
          if fmt(i)==' ' || fmt(i)=='-'
            continue;
          end
          [r1,s1]=obj.similarseqs([fmt(1:i-1),'-',fmt(i+1:end)],'maxerrors',args.maxerrors-1,'usecache',true,'sel',args.sel,'print',args.print,'ignore',[args.ignore,r_naseq]);
          r_naseq=[r_naseq,r1];
          r_seqstr={r_seqstr{:},s1{:}};
        end
        % Insertions
        % if args.print
        %   fprintf('Check for insertion and %d other error:\n',args.maxerrors-1);
        % end
        for i=1:length(fmt)+1
          if i<=length(fmt) && fmt(i)==' '
            % No need to insert before a blank
            continue;
          end
          [r1,s1]=obj.similarseqs([fmt(1:i-1),'N',fmt(i:end)],'maxerrors',args.maxerrors-1,'usecache',true,'sel',args.sel,'print',args.print,'ignore',[args.ignore,r_naseq]);
          r_naseq=[r_naseq,r1];
          r_seqstr={r_seqstr{:},s1{:}};
        end
      end
    end

    function naseqs=getgold(obj)
    % Get all 'gold' naseqs
      sel=strcmp({obj.tags.name},'gold');
      naseqs=[obj.tags(sel).naseq];
    end
  
    function g=groupseqs(obj,maxsubs,allseqs)
    % Form groups of sequences that have at most 'maxsubs' substitutions
      if nargin<3
        allseqs=obj.seq;
        alllen=obj.length;
      else
        alllen=cellfun(@(z) length(z),allseqs);
      end
      lens=unique(alllen);
      lens=lens(lens>maxsubs);
      %lens=lens(lens<39);   % For testing
      g=nan(length(allseqs),1);
      for i=1:length(lens)
        sel=find(alllen==lens(i));
        if length(sel)<2*maxsubs
          continue;
        end
        fprintf('Comparing %d sequences with length %d...', length(sel), lens(i));
        seqs=allseqs(sel);
        d=seqpdist(seqs,'method','p-distance','squareform',true)*lens(i);
        % tree() reorders leaf nodes, save indices as labels so we can find the leaves later
        tree=seqlinkage(d,'complete',arrayfun(@(z) int2str(z), 1:length(sel),'UniformOutput',false));
        clust=cluster(tree,maxsubs+.01);
        % Extract the mapping of clust index to sel 
        tsel=str2double(get(tree,'LeafNames'));
        g(sel(tsel))=clust+max([0;g]);
        maxdist=0;
        for i=1:length(clust)
          for j=1:length(clust)
            if clust(i)==clust(j)
              maxdist=max(maxdist,d(tsel(i),tsel(j)));
            end
          end
        end
        
        fprintf('merged %d sequences into %d clusters (max dist=%d)\n', length(sel), length(unique(clust)),maxdist);
        if length(unique(clust))<length(sel)
          ;%keyboard;
        end
        
        % merged=0;
        % for j=1:length(sel)-1
        %   if isfinite(g(sel(j)))
        %     continue;
        %   end
        %   diffs=seqs(j+1:end,:)~=seqs(j,:);
        %   total=sum(diffs');
        %   matches=find(total<=maxsubs)+j;
        %   if lens(i)==14 && length(matches)>0
        %     fprintf('len=%d, j=%d, have %d matches\n', lens(i), j, length(matches));
        %   end
        %   matches=matches(isnan(g(sel(matches))));   % Skip already matched ones
        %   g(sel(matches))=sel(j);
        %   merged=merged+length(matches);
        % end
        % fprintf('merged %d/%d sequences, now have %d groups\n', merged, length(sel), sum(isfinite(g)));
        % if merged~=sum(isfinite(g(sel)))
        %   fprintf('merged=%d, sum=%d\n', merged, sum(isfinite(g(sel))));
        %   keyboard;
        % end
      end
    end
    
    function x=switchstats(obj,varargin)
    % Produce table of switch statistics for all labelled switches
      defaults=struct('savefile',[],'goldonly',false,'summary',[],'naseq',[]);
      args=processargs(defaults,varargin);
      x=table();
      if args.goldonly
        naseq=int32(obj.getgold())';
      else
        naseq=int32(setdiff([obj.seqnames.naseq],[obj.tags(ismember({obj.tags.name},{'invalid','amplicon','ref','aptamer'})).naseq])');
      end
      if ~isempty(args.naseq)
        naseq=intersect(naseq,args.naseq);
      end
      x.naseq=naseq;
      x.Properties.VariableDescriptions{'naseq'}='ID#';

      x.name=arrayfun(@(z) obj.getname(z), x.naseq,'UniformOutput',false);
      x.Properties.VariableDescriptions{'name'}='Name';
      x=sortrows(x,'name');
      emptycolumn=arrayfun(@(z) '', x.naseq,'UniformOutput',false);  % All blank to start

      % Grab tags as individual columns
      alltags=unique({obj.tags(ismember([obj.tags.naseq],x.naseq)).name});
      for i=1:length(alltags)
        x.(alltags{i})=emptycolumn;
        x.Properties.VariableDescriptions{alltags{i}}=alltags{i};
      end
      x.Properties.VariableDescriptions{'src'}='Source';
      % And clusters
      x.cluster=emptycolumn;
      for i=1:length(x.naseq)
        tsel=find([obj.tags.naseq]==x.naseq(i));
        for k=1:length(tsel)
          t=obj.tags(tsel(k));
          if isempty(t.value)
            x.(t.name){i}='Y';
          else
            x.(t.name){i}=t.value;
          end
        end
        cl=obj.getcluster(x.naseq(i));
        x.cluster{i}=cl;
      end
      x.Properties.VariableDescriptions{'cluster'}='Clusters';
      
      % Cleavage results
      if isempty(args.summary)
        fprintf('No summary table -- skipping\n');
      else
        y=args.summary.summarytable();
        [~,ia,ib]=intersect(x.naseq,y.naseq);
        for k=2:length(y.Properties.VariableNames)
          nm=y.Properties.VariableNames{k};
          v=y.(nm);
          vord=nan(length(x.naseq),size(v,2));
          vord(ia,:)=v(ib,:);
          x.(nm)=vord;
          x.Properties.VariableDescriptions{nm}=y.Properties.VariableDescriptions{k};
        end
      end

      % Sequence
      x.seq=obj.naseqformat(x.naseq);
      % WT Core
      x.wtcore=cellfun(@(z) ~isempty(strfind(z,'CTGATGA')),x.seq);
      if ~isempty(args.savefile)
        writetable(x,args.savefile);
      end
    end

    % ------ CLUSTERING SUPPORT ----------
    function dbclusterremovemember(obj,naseq)
    % Remove naseqs from any clusters
      cmd=sprintf('delete from ngs.clustermembers where naseq in (%s0)',sprintf('%d,',naseq))
      mysql(cmd);
      obj.dbloadclusters();
    end
    
    function clusterseqs(obj,naseq,cnt,varargin)
    % Determine clusters for given naseqs using cnt for each as guide to ordering new additions
    % Only attempt to cluster sequences with >= mincnt reads (as given in input)
    % Only create new clusters based on seqs with >=minaddcnt reads (which defaults to 'mincnt')
      defaults=struct('maxdist',5,'mincnt',10,'minlength',20,'minaddcnt',[],'debug',false,'clusterunloaded',false,'onlycluster',[],'statusinterval',200);
      args=processargs(defaults,varargin);
      
      if isempty(args.minaddcnt)
        args.minaddcnt=args.mincnt;
      end

      unloaded=~ismember(naseq,obj.naseq);
      if sum(unloaded)>0
        if args.clusterunloaded
          fprintf('WARNING: clustering %d sequences that are not loaded - could create duplicates\n', sum(unloaded));
        else
          fprintf('WARNING: not clustering %d sequences that are not loaded\n', sum(unloaded));
          naseq=naseq(~unloaded);
        end
      end
      
      % Remove any that we already have clusters for from list
      orignum=length(naseq);
      if ~isempty(obj.clusters)
        [useqs,ia]=setdiff(naseq,abs([obj.clusters.members]));
        cnt=cnt(ia);
      else
        useqs=naseq;
      end
      % Remove any with <mincnt
      useqs=useqs(cnt>=args.mincnt);
      cnt=cnt(cnt>=args.mincnt);
      
      % Order by cnt descending and form new clusters whenever dist is too high
      [cnt,ord]=sort(cnt,'desc');
      useqs=useqs(ord);
      
      fprintf('Have %d/%d new seqs to cluster with at least %d cnts', length(useqs),orignum,args.mincnt);
      if isempty(useqs)
        fprintf('\n');
        return;
      end
      fprintf('...');
      seqs=obj.getorloadseq(useqs);
      nwcnt=0;
      origclusters=length(obj.clusters);
      origmembers=length([obj.clusters.members]);
      
      if isempty(args.onlycluster)
        clist=1:length(obj.clusters);
      else
        clist=find(ismember([obj.clusters.cluster],args.onlycluster));
        if isempty(clist)
          error('None of the clusters in "onlyclusters" exist');
        end
      end
      clistnamed=clist(arrayfun(@(z) ~isempty(z.name), obj.clusters(clist)));
      clistunnamed=clist(arrayfun(@(z) isempty(z.name), obj.clusters(clist)));
      fprintf('Have %d named and %d unnamed clusters out of %d\n', length(clistnamed), length(clistunnamed), length(clist));
      clear clist   % Not needed any more
      for i=1:length(useqs)
        if mod(i,args.statusinterval)==1
          fprintf('%d(%d,%d,%d)...',i,cnt(i),length(obj.clusters),nwcnt);
          nwcnt=0;
        end
        
        seq=seqs{i};
        assert(~iscell(seq));
        slen=length(seq);
        if slen<args.minlength
          %fprintf('S');
          continue;
        end
        mindist=args.maxdist+1;
        best=-1;
        % First, look for a named one
        for cmplendiff=0:args.maxdist   % Start with same length and work away
          for j=clistnamed
            assert(~isempty(obj.clusters(j).name));
            lendiff=length(obj.clusters(j).rootseq)-slen;
            if abs(lendiff)==cmplendiff && abs(lendiff)<mindist
              % Only align after we've tried all the same length ones
              dist=obj.seqcompare(obj.clusters(j).rootseq,seq);
              nwcnt=nwcnt+1;
              if dist<mindist
                mindist=dist;
                best=j;
                if args.debug
                  fprintf('\nDist to named cluster %s(%d) is %d\nroot:%s\nseq :%s\n', obj.clusters(j).name, obj.clusters(j).cluster, dist,obj.clusters(j).rootseq,seq);
                end
              end
            end
          end
        end
        if best==-1
          % No named one found, check unnamed
          for cmplendiff=0:args.maxdist   % Start with same length and work away
            if cmplendiff>=mindist
              break;
            end
            for j=clistunnamed
              assert(isempty(obj.clusters(j).name));
              lendiff=length(obj.clusters(j).rootseq)-slen;
              if abs(lendiff)==cmplendiff && abs(lendiff)<mindist
                % Only align after we've tried all the same length ones
                dist=obj.seqcompare(obj.clusters(j).rootseq,seq);
                nwcnt=nwcnt+1;
                if dist<mindist
                  mindist=dist;
                  best=j;
                  if args.debug
                    fprintf('\nDist to unnamed cluster (%d) is %d\nroot:%s\nseq :%s\n', obj.clusters(j).cluster, dist,obj.clusters(j).rootseq,seq);
                  end
                end
              end
            end
          end
        end
        if mindist<=args.maxdist
          %obj.clusters(best).cnt=obj.clusters(best).cnt+cnt(i);
          obj.clusters(best).members(end+1)=-useqs(i);
          obj.clusters(best).mcnt(end+1)=cnt(i);
          obj.clusters(best).dist(end+1)=mindist; 
        elseif isempty(args.onlycluster) && cnt(i)>=args.minaddcnt
          % Add a new cluster
          if isempty(obj.clusters)
            cid=-1;
          else
            cid=-(max(abs([obj.clusters.cluster]))+1);
          end
          name=obj.getname(useqs(i));
          c=struct('cluster',cid,'root',useqs(i),'rootseq',seq,'members',[-useqs(i)],'mcnt',cnt(i),'dist',0,'name',name);
          assert(~iscell(c.rootseq));
          obj.clusters=concatstructs(obj.clusters,c);
          if isempty(name)
            clistunnamed(end+1)=length(obj.clusters);
          else
            clistnamed(end+1)=length(obj.clusters);
          end
        end
      end
      fprintf('%d.\nAdded %d clusters, %d members.\n', length(useqs),length(obj.clusters)-origclusters,length([obj.clusters.members])-origmembers);
    end
    
    function ind=getclusterindex(data,naseq)
    % Get the array indicees for clusters for naseq, 0 if not found
      ind=zeros(size(naseq));
      for i=1:length(data.clusters)
        ind(ismember(naseq,abs([data.clusters(i).root,data.clusters(i).members])))=i;
      end
    end
    
    
    function [descr,cluster]=getcluster(data,naseq)
    % Get the cluster that the given seq belongs to, if any, and the descriptive diffs
      assert(length(naseq)==1);
      ind=data.getclusterindex(naseq);
      if ind==0
        descr='';
        cluster=[];
      else
        cluster=data.clusters(ind);
        [~,diffs,ch]=data.compare(naseq,cluster.root,'maxchanges',99,'silent',true,'exact',true);
        if ch==99
          error('Cluster member %d does not align with cluster root %d for cluster %d',naseq,cluster.root,ind);
        end
        if naseq~=cluster.root
          assert(~isempty(diffs));
        end
        if isempty(cluster.name)
          descr=sprintf('C%d(%s)',cluster.cluster,diffs);
        else
          descr=sprintf('%s(%s)',cluster.name,diffs);
        end
      end
    end
    
    function [neighbors,dists]=getcloseseqs(data,naseq,varargin)
    % Find all seqs (that are clustered) that are within maxdist of given naseq
      defaults=struct('maxdist',5);
      args=processargs(defaults,varargin);

      seq=data.getorloadseq(naseq);
      seq=seq{1};
      [~,cl]=data.getcluster(naseq);
      mydist=data.seqcompare(cl.rootseq,seq);  % Distance from naseq to assigned cluster
      fprintf('Naseq %d is distance %d from root of cluster %d\n', naseq, mydist, cl.cluster);
      neighbors=[]; dists=[];
      for i=1:length(data.clusters)
        c=data.clusters(i);
        dist=data.seqcompare(c.rootseq,seq);
        cldist=abs(data.seqcompare(c.rootseq,cl.rootseq));
        closestpossible=max(1,cldist-5);   % Closest possible distance between naseq and elements of this cluster
        if closestpossible<=args.maxdist && length(c.members)>0
          % Check the members of this cluster
          fprintf('Cluster %5d: %3d members could be within %d...', c.cluster, length(c.members), closestpossible);
          cnt1=length(dists);
          mseq=data.getorloadseq(abs(c.members));
          for j=1:length(c.members)
            pos=cldist-c.dist(j)-mydist;
            if pos<=args.maxdist
              d=data.seqcompare(mseq{j},seq);
              if d<=args.maxdist
                neighbors(end+1)=abs(c.members(j));
                dists(end+1)=d;
              end
            end
          end
          fprintf('%d\n',length(dists)-cnt1);
        end
      end
      fprintf('Found %d sequences within %d of %d\n', length(neighbors), args.maxdist, naseq);
    end
    
    function dbloadclusters(data)
      data.dbopen();  % Make sure to open here -- have to be USEing the right db for ngsentries to correspond
      
      %cmd='select cluster,root,name from ngs.clusters where cluster in (select cluster from ngs.clustermembers where naseq in (select naseq from ngsentries))';
      cmd='select cluster,root,n.name from ngs.clusters c left join ngs.naseqnames n on n.naseq=c.root order by root';
      fprintf('Loading clusters...');
      [cluster,root,name]=mysql(cmd);
      fprintf('%d\n', length(cluster));
      if strcmp(data.database,'ngs')
        % Load all
        cmd='select c.naseq, cluster, 1 total, dist from ngs.clustermembers c group by c.naseq,cluster,dist';
      else
        % Just sequences we're using
        cmd='select c.naseq, cluster, sum(e.cnt) total, dist from ngs.clustermembers c, ngsentries e where c.naseq=e.naseq group by c.naseq,cluster,dist order by total desc';
      end
      [naseq,mcluster,cnt,dist]=mysql(cmd);
      fprintf('done\n');
      data.clusters=[];
      for i=1:length(cluster)
        c=struct('cluster',cluster(i), 'root',root(i),'name',name{i});
        c.members=naseq(mcluster==cluster(i))';
        c.dist=dist(mcluster==cluster(i))';
        c.mcnt=cnt(mcluster==cluster(i))';
        data.clusters=[data.clusters,c];
      end
      if ~isempty(data.clusters)
        rootseqs=data.getorloadseq([data.clusters.root]);
        for i=1:length(data.clusters)
          data.clusters(i).rootseq=rootseqs{i};
          % 12/9/20: shouldn't need this
          % data.clusters(i).members=union(data.clusters(i).members,data.clusters(i).root);
        end
      end
      
      fprintf('Loaded %d clusters with %d members\n', length(cluster), length(naseq));
    end
    
    function dbsaveclusters(data)
    % Save modified clusters to database
      data.dbopen();
      % Prune
      % No need to prune, just skip these when updating
      %data.pruneclusters('mincnt',2,'minmembers',2);  % Get rid of clusters with 1 member, 1 cnt

      % Check if any flagged members are already in clusters
      newmembers=[data.clusters.members];
      newmembers=-newmembers(newmembers<0);
      if isempty(newmembers)
        fprintf('No new cluster members to save\n');
        return;
      end
      fprintf('Saving up to %d new clusters and %d new members (but skipping singleton unnamed clusters)...',sum([data.clusters.cluster]<0),length(newmembers));
      nclusters=0;
      nmembers=0;
      
      cmd=sprintf('select naseq from ngs.clustermembers where naseq in (%s0)',sprintf('%d,',newmembers));
      already=mysql(cmd);
      if length(already)>0
        fprintf('Already have %d/%d potential new members in database - skipping them\n', length(already), length(newmembers));
      end
      
      for i=1:length(data.clusters)
        c=data.clusters(i);
        if c.cluster<0
          % New cluster, add to DB
          if length(c.members)==1 && isempty(c.name)
            %fprintf('-');   % Skip it
            continue;
          else
            cmd=sprintf('insert into ngs.clusters(root) values(%d)',c.root);
            xx=mysql(cmd);
            c.cluster=mysql('select LAST_INSERT_ID()');
            %fprintf('%s->%d\n',cmd,c.cluster);
            c.members=-abs(c.members);   % All need to be added
            fprintf('+');
            nclusters=nclusters+1;
          end
        end
        if any(c.members<0)
          % New members (break into max 1000 entries per SQL stmt)
          for j1=1:1000:length(c.members)
            cmd='insert into ngs.clustermembers(cluster,naseq,cnt,dist) values';
            added=false;
            for j=j1:min(length(c.members),j1+999)
              if c.members(j)<0 && ~ismember(-c.members(j),already)
                cmd=[cmd,sprintf('(%d,%d,%d,%d),\n',c.cluster,-c.members(j),c.mcnt(j),c.dist(j))];
                c.members(j)=-c.members(j);
                nmembers=nmembers+1;
                added=true;
              end
            end
            if added
              cmd=cmd(1:end-2);
              xx=mysql(cmd);
              fprintf('M');
            end
          end
        end
        data.clusters(i)=c;
      end
      fprintf('done\nAdded %d clusters and %d members to database\n',nclusters,nmembers);
    end
    
    function [badids,badseqs]=verifyclusters(obj)
    % Check all clusters
      fprintf('Verifying %d clusters with %d members...',length(obj.clusters),length([obj.clusters.members]));
      badids=[];badseqs=[];
      % Preload all the root seqs at once
      preload=setdiff(union([obj.clusters.root],[obj.clusters.members]),[obj.naseq]);
      fprintf('Preloading %d seqs...',length(preload));
      obj.getorloadseq(preload);
      fprintf('done\n');
      
      fprintf('Verifying %d clusters...',length(obj.clusters));
      for i=1:length(obj.clusters)
        if mod(i,100)==0
          fprintf('%d...',i);
        end
        c=obj.clusters(i);
        for j=1:length(c.members)
          if abs(c.members(j))==abs(c.root)
            assert(c.dist(j)==0);
          else
            [~,~,dist]=obj.compare(abs(c.root),abs(c.members(j)),'silent',true);
            if dist>c.dist(j)   % Sometimes insertion/deletions reduce dist by 1 compared to substitutions
              fprintf('Cluster(%d) ID=%d has dist between %d and %d as %d, but measured is %d\n', i, c.cluster, c.root, c.members(j), c.dist(j), dist);
              badids=[badids,c.cluster];
              badseqs=[badseqs,c.members(j)];
            end
          end
        end
      end
      badids=unique(badids);
    end
    
    function dbsplitcluster(obj,naseq)
    % Split a cluster by moving naseq to its own cluster and then reassign all members
      cluster=mysql(sprintf('select cluster from ngs.clustermembers where naseq=%d',naseq));
      if isempty(cluster)
        fprintf('Naseq %d not found in any cluster\n',naseq)
        return;
      end
      isroot=mysql(sprintf('select count(*) from ngs.clusters where root=%d',naseq));
      if isroot
        fprintf('Naseq %d is already the root of cluster %d\n', naseq, cluster);
        return;
      end
      
      [others,olddist]=mysql(sprintf('select naseq,dist from ngs.clustermembers where cluster=%d and (naseq=%d or dist>1)',cluster,naseq));
      fprintf('Have %d other seqs to possibly reassign\n', length(others));
      seq=obj.getorloadseq(naseq);
      tomove=[];
      newdist=[];
      for i=1:length(others)
        [~,~,d]=obj.compare(naseq,others(i),'silent',true,'exact',true,'maxchanges',99);
        if ~isempty(d) && d<olddist(i)
          fprintf('%d is closer to %d (%d) than to old root (%d)\n', others(i), naseq, d, olddist(i));
          tomove(end+1)=others(i);
          newdist(end+1)=d;
        end
      end
      if ~isempty(tomove)
        cmd=sprintf('insert into ngs.clusters(root) values(%d)',naseq)
        mysql(cmd);
        newcluster=mysql('select LAST_INSERT_ID()');
        for k=1:length(tomove)
          cmd=sprintf('update ngs.clustermembers set cluster=%d,dist=%d where naseq=%d',newcluster,newdist(k),tomove(k));
          mysql(cmd);
        end
        fprintf('Need to reload clusters...\n');
      end
    end
    
    function dbdropcluster(obj,cluster,varargin)
    % Drop the given cluster and reassign its members
      defaults=struct('dropnamed',0,'reassign',true);
      args=processargs(defaults,varargin);

      obj.dbopen();
      others=[];
      for i=1:length(cluster)
        c=cluster(i);
        sel=find([obj.clusters.cluster]==c);
        if isempty(sel)
          error('Cluster %d not loaded\n',c);
        end
        if ~isempty(obj.clusters(sel).name) 
          if args.dropnamed
            fprintf('About to drop named cluster %s...pause(5)...',obj.clusters(sel).name);
            pause(5);
            fprintf('done\n');
          else
            error('Cannot delete named cluster %s (%d)\n', obj.clusters(sel).name, c);
          end
        end
        [o,olddist]=mysql(sprintf('select naseq,dist from ngs.clustermembers where cluster=%d',c));
        others=[others;o];
        cmd=sprintf('delete from ngs.clusters where cluster=%d',c);
        mysql(cmd);
      end

      fprintf('Have %d other seqs to possibly reassign\n', length(others));
      obj.dbloadclusters();
      if args.reassign
        obj.clusterseqs(others,ones(size(others)),'mincnt',1,'minlength',5);
      end
      obj.dbsaveclusters();
    end
    
    function dbmovetonamedclusters(obj)
    % Move any cluster members to named clusters if possible
      obj.dbopen();
      for i=1:length(obj.clusters)
        c=obj.clusters(i);
        if isempty(c.name)
          continue;
        end
        if ~ismember(c.root,c.members)
          % Not loaded
          continue;
        end
        [neighbors,dists]=getcloseseqs(obj,c.root,'maxdist',5);
        mods=[];
        for k=1:length(neighbors)
          ind=obj.getclusterindex(neighbors(k));
          cl=obj.clusters(ind);
          olddist=cl.dist(cl.members==neighbors(k));
          if isempty(cl.name) || dists(k)<olddist
            fprintf('Can move %d from cluster %d (dist %d) to cluster %s (%d) with dist %d\n', neighbors(k), cl.cluster, olddist, c.name, c.cluster,dists(k));
            mods=[mods,neighbors(k)];
          end
        end
        if length(mods)>0
          cmd=sprintf('update ngs.clustermembers set cluster=%d where naseq in (%s-1);',c.cluster,sprintf('%d,',mods))
          nc=mysql(cmd);
          fprintf('Updated %d rows\n',nc);
          assert(nc==length(mods));
        end
      end
      obj.dbloadclusters();
    end
    
    function pruneclusters(obj,varargin)
      % Only keep ones with cnt>=mincnt, members>=minmembers, or named, or from DB
      defaults=struct('mincnt',10,'minmembers',2);
      args=processargs(defaults,varargin);
      if isempty(obj.clusters)
        return;
      end
      cnt=arrayfun(@(z) sum(z.mcnt), obj.clusters);
      named=arrayfun(@(z) ~isempty(z.name), obj.clusters);
      nmembers=arrayfun(@(z) length(z.members), obj.clusters);

      sel=(cnt>=args.mincnt & nmembers>=args.minmembers) | named | [obj.clusters.cluster]>0;
      if any(~sel)
        fprintf('Removing %d/%d unnamed clusters with less than %d reads or less than %d members\n', sum(~sel), length(sel), args.mincnt, args.minmembers);
        obj.clusters=obj.clusters(sel);
      end
    end

    function consensus=getconsensus(obj,c,varargin)
    % Get consensus seqs for clusters
      defaults=struct('force',false,'usefrac',.75,'debug',false);
      args=processargs(defaults,varargin);
      if length(c.members)<3
        % Not enough members for multialign, use most frequent one
        [~,sel]=max([c.mcnt]);
        consensus=obj.getseq(c.members(sel));
        consensus=consensus{1};
      else
        cseq=obj.getorloadseq(abs(c.members));
        mcnt=min(c.mcnt,1);
        mcnt(1)=mcnt(1)+1;   % Root gets bonus
        [mcnt,ord]=sort(mcnt,'desc');
        cseq=cseq(ord);
        keep=find(cumsum(mcnt)>=args.usefrac*sum(mcnt),1);   % Keep only the main ones to make the multialign faster
        if isempty(keep)
          keep=length(mcnt);
        end
        keep=max(min(keep,10),3);   % Need >=3 for multialign
        cseq=cseq(1:keep);mcnt=mcnt(1:keep);
        ma=multialign(cseq);
        if args.debug
          ma
          mcnt
        end
        nts='ACGT-';
        ntcnt=[];
        for j=1:length(nts)
          ntcnt(j,:)=mcnt*(ma==nts(j));
        end
        [~,maxind]=max(ntcnt);
        consensus=nts(maxind(maxind<=4));
        if args.debug
          fprintf('Consensus: %s\n', conseq);
        end
      end
    end
    
    function listconsensusdiffs(obj)
    % List clusters with consensus different from root
      for i=1:length(obj.clusters)
        c=obj.clusters(i);
        consensus=obj.getconsensus(c);
        if ~strcmp(c.rootseq,consensus)
          fprintf('\nCluster %d %s:\n', c.cluster,c.name);
          [d,al]=obj.seqcompare(c.rootseq,consensus);
          fprintf('%s Root\n%s\n%s Consensus\n', al(1,:), al(2,:),al(3,:));
        end
      end
    end
    
    function dist=getclusterdists(obj,varargin)
    % Find distance between cluster roots up to and including maxdist
    % Negative values indicate the minimum distance (e.g. -6 means at least 6)
      defaults=struct('maxdist',10,'onlynamed',false,'onlycluster',[]);
      args=processargs(defaults,varargin);

      nclust=length(obj.clusters);
      dist=zeros(nclust,nclust,'int8');
      dist(:)=-1;
      if ~isempty(obj.clusterdists)
        [~,ia,ib]=intersect([obj.clusters.cluster],obj.clusterdists.cluster);
        for i=1:length(ia)
          dist(ia(i),ia)=obj.clusterdists.dist(ib(i),ib);
        end
        fprintf('Using cached cluster dists for %d/%d\n',length(ia),nclust);
      end
      fprintf('Computing pairwise distances between %d clusters up to maxdist=%d...',nclust,args.maxdist);

      nc=0;
      for i=1:nclust
        if mod(i,200)==0
          fprintf('%d...',i);
        end
        if isempty(args.onlycluster) || ismember(obj.clusters(i).cluster,args.onlycluster)
          jvals=i+1:nclust;
        else
          jvals=find(ismember([obj.clusters.cluster],args.onlycluster));
          jvals=jvals(jvals>i);
        end
        if all(dist(i,jvals)>0 | -dist(i,jvals) > args.maxdist)
          continue;   % Already done
        end
        for j=jvals
          if dist(i,j)>0 || -dist(i,j) > args.maxdist
            continue;  % Already computed
          end
          if args.onlynamed && isempty(obj.clusters(i).name) && isempty(obj.clusters(j).name)
            continue;
          end
          
          lendiff=abs(length(obj.clusters(i).rootseq)-length(obj.clusters(j).rootseq));
          if  lendiff> args.maxdist
            dist(i,j)=-lendiff;
          else
            dist(i,j)=obj.seqcompare(obj.clusters(i).rootseq,obj.clusters(j).rootseq);
            nc=nc+1;
          end
        end
        dist(i,1:i-1)=dist(1:i-1,i)';  % Transpose
        dist(i,i)=0;
      end
      fprintf('done\nUsed %d seq compares\n',nc);
      obj.clusterdists=struct('cluster',[obj.clusters.cluster],'dist',dist);
    end
    
    function res=listsimilarclusters(obj,varargin)
    % Find clusters whose root differ by <= maxdist
      defaults=struct('maxdist',5,'onlynamed',false,'debug',false,'unnamed',false);
      args=processargs(defaults,varargin);

      fprintf('Comparing %d clusters...',length(obj.clusters));
      dist=obj.getclusterdists('maxdist',args.maxdist,'onlynamed',args.onlynamed);
      res=[];
      for i=1:length(obj.clusters)
        for j=i+1:length(obj.clusters)
          if abs(dist(i,j))<=args.maxdist
            if args.onlynamed && isempty(obj.clusters(i).name) && isempty(obj.clusters(j).name)
              continue;
            end
            if args.unnamed && ~isempty(obj.clusters(i).name) && ~isempty(obj.clusters(j).name)
              continue;
            end
            fprintf('Clusters %d (%s,N=%d) and %d (%s,N=%d) have a distance of %d\n', obj.clusters(i).cluster,obj.clusters(i).name, length(obj.clusters(i).members), obj.clusters(j).cluster,obj.clusters(j).name,length(obj.clusters(j).members),dist(i,j));
            if args.debug
              obj.compare(obj.clusters(i).root,obj.clusters(j).root,'maxchanges',dist(i,j));
              %            fprintf(' %s\n %s\n %s\n',al(1,:),al(2,:),al(3,:));
            end
            res(end+1,:)=[obj.clusters(i).cluster,obj.clusters(j).cluster,double(dist(i,j))];
          end
        end
      end
    end
    
    function dbupdateprimermutants(obj)
    % Find any named sequences that are primer mutants and tag them
      for naseq=[obj.seqnames.naseq]
        seq=obj.getseq(naseq);
        seq=seq{1};
        s=regexp(seq,'^GCTG.*GGACGAAACAGC$');
        if isempty(s)
          fprintf('Seq %d: %s - primer mutant\n',naseq,seq);
          % Figure out what this is likely corrected to:
          corrected=seq;
          if strncmp(corrected,'CTG',3)
            corrected=['G',corrected];
          end
          if sum(corrected(1:4)~='GCTG')==1
            corrected(1:4)='GCTG';
          end
          if sum(corrected(end-11:end)~='GGACGAAACAGC')<=2
            corrected(end-11:end)='GGACGAAACAGC';
          end
          if strcmp(seq,corrected)
            fprintf('Unable to figure out corrected seq\n');
            obj.dbaddtag(naseq,'primermutant','?',true);
          else
            other=obj.findseq(corrected);
            fprintf('After PCR, gives: %s (%d)\n', corrected,other);
            obj.dbaddtag(naseq,'primermutant',other,true);
          end
          fprintf('\n');
        end
      end
      obj.dbloadtags();
    end

    function [reassign,mcnt]=reassignclusters(obj,varargin)
      defaults=struct('maxdist',5,'cluster',[],'movefromnamed',false);
      args=processargs(defaults,varargin);

      fprintf('Comparing %d clusters...',length(obj.clusters));
      dist=obj.getclusterdists('maxdist',args.maxdist,'onlycluster',args.cluster);
      reassign=[]; 
      mcnt=[];   % Prior count of reads of this member
      if isempty(args.cluster)
        cindices=1:length(obj.clusters)
      else
        cindices=find(ismember([obj.clusters.cluster],args.cluster));
      end
      for ii=1:length(cindices)
        i=cindices(ii);
        if mod(ii-1,100)==0
          fprintf('%d...',ii);
        end
        if ~isempty(obj.clusters(i).name) && ~args.movefromnamed
          % Don't move away from from a named cluster
          continue;
        end
        obj.getorloadseq(abs(obj.clusters(i).members));  % Discard, but loads into otherseqs
        for j=1:length(obj.clusters)
          if i==j
            continue;
          end
          if abs(dist(i,j))<=5*2-1
            for k=1:length(obj.clusters(i).members)
              d1=obj.clusters(i).dist(k);
              if 2*d1 > abs(dist(i,j))
                % Possible that seq could be closer to cluster(j)
                naseq=abs(obj.clusters(i).members(k));
                seq=obj.getorloadseq(naseq);
                seq=seq{1};
                d2=obj.seqcompare(obj.clusters(j).rootseq,seq);
                if d2<d1
                  fprintf('Member %d is closer to cluster %d (d=%d) than cluster %d (d=%d)\n', naseq, j, d2, i, d1);
                  reassign(end+1)=naseq;
                  mcnt(end+1)=obj.clusters(i).mcnt(k);
                end
              end
            end
          end
        end
      end
      fprintf('Have %d sequences that should be reassigned\n', length(reassign));
      if length(reassign)>0
        % Remove these sequences from their current clusters
        obj.dbclusterremovemember(reassign);
        obj.clusterseqs(reassign,mcnt);
        obj.dbsaveclusters();
        fprintf('Removed %d sequences from their current clusters and reassigned\n',length(reassign));
      end
    end

    function sam=bowtiealign(obj,ref,naseqs)
    % Align sequences to a reference using BowTie2
    % Note: shapemapper2 command line looks like the following (with defaults after):
    % --wrapper basic-0 
    % -p 4  (num threads, default 1)
    % --local (default is --end-to-end )
    % --sensitive-local  (default in --local mode,  -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 )
    % --mp 3,1 (max/min mismatch penalties, default is 6,2)
    % --rdg 5,1 (read gap penalties [open,extend], default is 5,3 )
    % --rfg 5,1 (reference gap penalties, default is 5,3)
    % --dpad 30 (padding of dynamic programming problems, default is 15)
    % --maxins 800 (fragment length for paired end alignments, default is 500)
    % --ignore-quals (ignore-qualities)
    % --tab6 (file format)
    % bowtie2 defaults are (from matlab alignoptions.getBowtie2Command() )
    % --np 1 --phred33 --ma 0 --n-ceil L,0,0.15 --score-min L,-0.6,-0.6 --mp 6,2 --end-to-end --gbar 4 -R 2 -D 15 -N 0 -p 1 --dpad 15 --rdg 5,3 --rfg 5,3 --seed 0 -i S,1,1.15 -L 20 -s 0 -3 0 -5 0'

    % TODO: if the read is shorter than the reference, the alignment (and its score) doesn't consider the unmatched reference.  In our case, it might be preferable to enforce that the 5' end of the cDNA starts with the reference since the cDNA should start there (since the RT primer was already aligned to the read in prior processing), though the 3' end may be truncated if the RT aborted early.
    % Temp folder
      tmp=tempname('/tmp');
      mkdir(tmp);
      cd(tmp);
      % Build index
      refseq=obj.getorloadseq(ref); refseq=refseq{1};
      refname=sprintf('N%d',ref);
      fastawrite('ref.fa',refname,refseq);
      if bowtie2build('ref.fa','ref_index')~=0
        error('bowtie2build failed');
      end
      % Save seqs as fastq
      seqs=obj.getorloadseq(naseqs);
      ise=cellfun(@(z) isempty(z), seqs);
      if any(ise)
        fprintf('bowtiealign: empty sequences omitted\n');
        naseqs=naseqs(~ise);
        seqs=seqs(~ise);
      end
      names=arrayfun(@(z) sprintf('N%d',z), naseqs,'Unif',false);
      qual=cellfun(@(z) repmat('A'+30,1,length(z)),seqs,'Unif',false);
      fastqwrite('targets.fq',struct('Header',names,'Sequence',seqs,'Quality',qual));
      % Align
      alignOpt = Bowtie2AlignOptions;
      alignOpt.IgnoreQuality=true;
      alignOpt.ReadGapCosts=[5,1];
      alignOpt.RefGapCosts=[5,1];
      alignOpt.MismatchPenalty=[3,1];
      alignOpt.PadPositions=30;   % Doesn't make much difference from 15, but 0 results in unwanted gaps at ends
      alignOpt.SeedLength=22;
      alignOpt.NumThreads=4;
      % Different from ShapeMapper2
      alignOpt.Mode='EndToEnd';  % ShapeMapper uses local, but here we want full matches
                                 % alignOpt.NoGapPositions=1;  % Don't treat start/end any differently for gaps (since there was extra sequence around before)
                                 %arglist=alignOpt.getBowtie2Command('IncludeAll',true);
      status=bowtie2('ref_index','targets.fq','','align.sam',alignOpt);
      if status~=0
        error('bowtie2 failed');
      end
      % Read SAM file
      sam=samread('align.sam');
      verify=arrayfun(@(z) str2double(z.QueryName(2:end)),sam);
      [~,~,ib]=intersect(naseqs,verify,'stable');
      sam=sam(ib);
      verify=arrayfun(@(z) str2double(z.QueryName(2:end)),sam);
      assert(all(verify==naseqs));
      for i=1:length(sam)
        sam(i).Tags.naseq=verify(i);
        sam(i).Tags.ref=ref;
      end
      rmdir(tmp,'s');
    end

    function m=sam2muts(obj,sam)
    % Convert a sam record into a string of mutations
      if strcmp(sam.ReferenceName,'*')
        % Unaligned
        m={};
        return;
      end
      seq=sam.Sequence;
      cig=sam.CigarString;
      refnaseq=str2double(sam.ReferenceName(2:end));
      refseq=obj.getorloadseq(refnaseq); refseq=refseq{1};

      qpos=1; rpos=int32(sam.Position);
      s=repmat('-',1,rpos-1);
      r=refseq(1:rpos-1);
      while ~isempty(cig)
        i=find(cig<'0' | cig>'9',1);
        n=str2double(cig(1:i-1));
        c=cig(i);
        cig=cig(i+1:end);
        if c=='M'
          s=[s,seq(qpos:qpos+n-1)];
          r=[r,refseq(rpos:rpos+n-1)];
          qpos=qpos+n;
          rpos=rpos+n;
        elseif c=='I'||c=='S'
          s=[s,seq(qpos:qpos+n-1)];
          r=[r,repmat('-',1,n)];
          qpos=qpos+n;
        elseif c=='D'
          s=[s,repmat('-',1,n)];
          r=[r,refseq(rpos:rpos+n-1)];
          rpos=rpos+n;
        else
          assert(false);
        end
      end
      s=[s,repmat('-',1,length(refseq)-rpos+1)];
      r=[r,refseq(rpos:end)];
      assert(length(r)==length(s));
      al(1,:)=s;
      al(3,:)=r;
      al(2,:)=blanks(length(r));
      al(2,r==s)='|';
      fprintf('%s\n%s\n%s\n',al(1,:),al(2,:),al(3,:));
      if isfield(sam.Tags,'adducts')
        astr=al(3,:);
        astr(astr~='-')=sam.Tags.adducts;
        fprintf('%s\n',astr);
      end
      if isfield(sam.Tags,'coverage')
        cstr=al(3,:);
        cstr(cstr~='-')=sam.Tags.coverage;
        fprintf('%s\n',cstr);
      end
      
      changes=find(s~=r);
      m=cell(1,length(changes));
      for ii=1:length(changes)
        i=changes(ii);
        rpos=sum(r(1:i)~='-');
        m{ii}=sprintf('%c%d%c',r(i),rpos,s(i));
      end
    end

  end % methods

end % classdef
