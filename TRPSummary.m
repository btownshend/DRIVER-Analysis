% Class to encapsulate summary data of TRP experiments over multiple NGS runs
% TODO: Check out https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/PASS/Confidence_Intervals_for_the_Odds_Ratio_of_Two_Proportions.pdf  for confidence intervals
classdef TRPSummary < handle
  properties
    run;	% NGS run in which this summary was built
    naseq;	% Sequences with data
    names;	% Names for each sequence
    cleavage;	% Cleavage data for the list of seq at each condition
    stdLogRatio; % Expected std(log(ratio))
    DFLogRatio;  % D.F. for above
    metric;      % metric used for fold change computations
    cirange;     % C.I. range for fold change CIs
    tree;	 % structure for hierarchical clustering
    refratio;	 % Ratio for references
    rename;	 % Map to rename targets
    mergenaseqs; % List of naseqs used in each cluster after merging.
  end
  
  
  methods
    function summary=TRPSummary(run,naseqs)
      summary.run=run;
      summary.naseq=unique(naseqs(:));
      summary.names=arrayfun(@(z) '',summary.naseq,'UniformOutput',false);
      summary.rename=containers.Map();
    end
    
    function save(summary,file)
    % Save data in a .mat file
    % Use HDF (v7.3) since the file is usually >2GB
      if isempty(summary.mergenaseqs)
        vname=sprintf('summary%s',strrep(summary.run,'NGS',''));
      else
        vname=sprintf('merged%s',strrep(summary.run,'NGS',''));
      end
      if nargin<2
        file=sprintf('%s.mat',vname);
      end
      eval(sprintf('%s=summary;',vname));
      fprintf('Saving %s in %s...',vname,file);
      tic
      save(file,'-v7.3',vname);
      elapsed=toc;
      d=dir(file);
      fprintf('done. Took %.1f seconds - %.2f Gbytes\n', elapsed, d.bytes/1e9);
    end

    function setnames(obj,ngs)
      for i=1:length(obj.naseq)
        if isempty(obj.names{i})
          obj.names{i}=ngs.seqs.getname(obj.naseq(i));
        end
        if isempty(obj.names{i})
          % No name, use cluster name if any
          obj.names{i}=ngs.seqs.getcluster(obj.naseq(i));
        end
        if isempty(obj.names{i})
          obj.names{i}=sprintf('New-%04d',mod(obj.naseq(i),10000));
        end
      end
    end
    
    function renametarget(obj,ngsname,summaryname)
      obj.rename(ngsname)=summaryname;
    end
    
    function addsummary(obj,ngs,trpexpt,varargin)
    % Save trpexpt data from NGS into obj, updating ngs.targets and ngs.trpruns as needed
      defaults=struct('mincnt',20,'concfuzz',0.05);
      args=processargs(defaults,varargin);

      % Add any missing names
      for i=1:length(obj.names)
        if isempty(obj.names{i})
          obj.names{i}=ngs.seqs.getname(obj.naseq(i));
        end
        if isempty(obj.names{i})
          % No name, use cluster name if any
          obj.names{i}=ngs.seqs.getcluster(obj.naseq(i));
        end
      end
      
      ngs.trpanalyze('mincnt',args.mincnt);
      sel=[ngs.trpexpts.trpexpt]==trpexpt;
      if ~any(sel)
        error('TRP expt %d not found',trpexpt);
      end
      
      t=ngs.trpexpts(sel);
      if isnan(t.robotrun)
        error('No robotrun for %s',t.descr);
      end
      if isempty(t.mixture) || isnan(t.mixture)
        error('No mixture for %s (target=%s)',t.descr,t.target)
      end
      target=t.target;
      if isKey(obj.rename,target)
        target=obj.rename(target);
      end
      ind=[];
      for i=1:length(obj.cleavage)
        c=obj.cleavage(i);
        if c.mixture==t.mixture && (isempty(t.targetconc) ||  abs((c.conc-t.targetconc)/c.conc)<args.concfuzz)
            % Matched existing
            ind=i;
            break;
        end
      end
      if isempty(ind)
        obj.cleavage(end+1).target=target;
        ind=length(obj.cleavage);
        obj.cleavage(ind).mixture=t.mixture;
        if isempty(t.targetconc)
          obj.cleavage(ind).conc=0;
        else
          obj.cleavage(ind).conc=t.targetconc;
        end
        obj.cleavage(ind).ratioboot={};
        obj.cleavage(ind).rawratio=[];
        obj.cleavage(ind).rawratioci90=[];
        obj.cleavage(ind).cnt=[];
        obj.cleavage(ind).runs={};
        obj.cleavage(ind).trpexpt=[];
        obj.cleavage(ind).negatives={};
      end
      
      % Add data from this experiment
      a=t.analysis;
      [~,ia,ib]=intersect(a.naseq,obj.naseq);   % Line up measurements
      ratio=nan(length(obj.naseq),1);
      ratioboot=nan(length(obj.naseq),size(a.ratioboot,2));
      r90=nan(length(obj.naseq),2);
      n=zeros(length(obj.naseq),1);
      cnt=zeros(length(obj.naseq),3);
      
      ratio(ib)=a.ratio(ia);
      ratioboot(ib,:)=a.ratioboot(ia,:);
      r90(ib,:)=prctile(a.ratioboot(ia,:),[5,95],2);
      n(ib)=sum(a.cnt(ia,2:3),2);
      cnt(ib,:)=a.cnt(ia,:);
      % Remove ones with too few counts
      ratio(n<args.mincnt)=nan;
      r90(n<args.mincnt,:)=nan;
      
      pos=length(obj.cleavage(ind).runs)+1;
      obj.cleavage(ind).ratioboot{pos}=single(ratioboot);
      obj.cleavage(ind).rawratio(:,pos)=ratio;
      obj.cleavage(ind).rawratioci90(:,:,pos)=r90;
      assert(all(all(isfinite(cnt(:,2:3)))));
      obj.cleavage(ind).cnt(:,:,pos)=cnt;
      ucode=strrep(ngs.subruns.codes{ngs.subruns.subrun==t.uncleavedsr},'/X','');
      ccode=strrep(ngs.subruns.codes{ngs.subruns.subrun==t.cleavedsr},'/X','');
      obj.cleavage(ind).runs{pos}=sprintf('%s;%d;%s->%s',ngs.run,t.robotrun,ucode,ccode);
      obj.cleavage(ind).trpexpt(pos)=trpexpt;
      obj.cleavage(ind).negatives{pos}=[];
    end

    function checkforduplicates(summary)
    % Check for technical replicates -- the same cleaveseq that was barcoded and included multiple times
      for i=1:length(summary.cleavage)
        c=summary.cleavage(i);
        robotrun=[];
        for j=1:length(c.runs)
          s=strsplit(c.runs{j},';');
          robotrun(j)=str2double(s{2});
        end
        for j=1:length(robotrun)
          for k=j+1:length(robotrun)
            if robotrun(j)==robotrun(k)
              fprintf('cleaveage(%d) - %s - has %s and %s\n', i, c.target,c.runs{j}, c.runs{k});
            end
          end
        end
      end
    end
    
    function setscaling(summary,refind,varargin)
    % Adjust the scaling of all summary data to match summary.cleavage(refind)
    % Setup base average
      defaults=struct('nonswitching',[],'prctile',10);
      args=processargs(defaults,varargin);

      c=summary.cleavage;
      % First rescale each of the runs within each target entry 
      for i=1:length(c)
        meanratio=exp(nanmean(log(c(i).rawratio),2));
        for j=1:length(c(i).runs)
          c(i).scaling(j)=nanmedian(meanratio./c(i).rawratio(:,j));
          c(i).ratiosc(:,j)=c(i).rawratio(:,j)*c(i).scaling(j);
        end
        c(i).ratioavg=exp(nanmean(log(c(i).ratiosc),2));
      end
      
      % Now rescale those to the reference target (refind)
      refratio=c(refind).ratioavg;
      
      if isempty(args.nonswitching)
        fprintf('Assuming %d%% of all sequences are non-switching\n',args.prctile);
        nonsw=true(size(refratio,1),1);
      else
        nonsw=ismember(summary.naseq,args.nonswitching)&isfinite(refratio);
        fprintf('Using %d non-switching sequences to set scaling\n',sum(nonsw));
      end
        
      for i=1:length(c)
        if i==refind
          continue;
        end
        if isempty(args.nonswitching)
          s2=exp(prctile(log(refratio./c(i).ratioavg),args.prctile));
        else
          s2=exp(nanmedian(log(refratio(nonsw)./c(i).ratioavg(nonsw))));
        end
        %fprintf('i=%d,s2=%.2f\n',i,s2);
        c(i).scaling=c(i).scaling*s2;
        c(i).ratiosc=c(i).ratiosc*s2;
        c(i).ratioavg=c(i).ratioavg*s2;
      end
      summary.cleavage=c;
      summary.refratio=refratio;
    end
    
    function computefolds_old(summary,varargin)
    % Compute fold changes between +target and -target conditions
      defaults=struct('mincnt',10,'ci',normcdf([-1,1])*100,'metric','ratio_of_ratio');
      args=processargs(defaults,varargin);
      summary.metric=args.metric;
      summary.cirange=args.ci;

      r0=summary.cleavage(1).ratiosc;
      obj=NGSVitro('tmp');
      for i=1:length(summary.cleavage)
        r=summary.cleavage(i).ratiosc;
        summary.cleavage(i).ofold=nan(size(r,1),1);
        summary.cleavage(i).ofoldci=nan(size(r,1),2);
        for j=1:size(r,1)
          r0j=r0(j,isfinite(r0(j,:)))';
          rj=r(j,isfinite(r(j,:)))';
          if ~isempty(r0j) && ~isempty(rj)
            % Generate bootstrap samples of the means (sigma=population sigma/sqrt(n))
            nboot=1000;
            if length(r0j)>2
              r0jStdOfMean=std(log(r0j))/sqrt(length(r0j));
              r0jDOF=length(r0j)-1;
            else
              fprintf('Only %d samples of -target for %s: using overall std\n', length(r0j), summary.names{j});
              r0jStdOfMean=summary.stdLogRatio/sqrt(length(r0j));
              r0jDOF=summary.DFLogRatio;
            end
            
            if length(rj)>2
              rjStdOfMean=std(log(rj))/sqrt(length(rj));
              rjDOF=length(rj)-1;
            else
              % Use r0 stats for sigma
              rjStdOfMean=r0jStdOfMean;
              rjDOF=r0jDOF;
            end
            r0jboot=max(1,min(20,exp(trnd(r0jDOF,nboot,1)*r0jStdOfMean+mean(log(r0j)))));
            rjboot=max(1,min(20,exp(trnd(rjDOF,nboot,1)*rjStdOfMean+mean(log(rj)))));

            [meanlogfold,stdlogfold,foldsamps]=obj.bootfold(rjboot,r0jboot,'metric',args.metric);
            if ~any(isfinite(meanlogfold))
              keyboard;
            end
            summary.cleavage(i).ofold(j)=exp(meanlogfold);
            summary.cleavage(i).ologfoldsem(j)=0.5*sqrt(rjStdOfMean^2+r0jStdOfMean^2);   % Scale by 0.5 since ratio of ratio has sqrt(r1/r2)
            summary.cleavage(i).ofoldci(j,:)=prctile(foldsamps,args.ci);
            summary.cleavage(i).ofoldsamps(j,:)=foldsamps;
          end
        end
      end
    end
    
    function showswitch(summary,naseq,varargin)
    % Dump all info about a particular switch
      defaults=struct('minfold',1.5,'extra',[]);
      args=processargs(defaults,varargin);

      if ischar(naseq)
        ind=find(strcmp(summary.names,naseq));
        if length(ind)>1
          fprintf('%s is ambiguous, could be any of [%s]\n', naseq, sprintf('%d ',summary.naseq(ind)));
          return;
        end
      else
        ind=find(summary.naseq==naseq);
      end
      if isempty(ind)
        error('Not found\n');
      end
      fprintf('Summary for %s (%d) maxfold=%.2f, listing targets with fold>=%.2f\n', summary.names{ind}, summary.naseq(ind),nanmax(arrayfun(@(z) z.fold(ind),summary.cleavage)),args.minfold);
      for i=1:length(summary.cleavage)
        c=summary.cleavage(i);
        if c.foldci(ind,1) > 2
          flag='**';
        elseif c.foldci(ind,1) > 1.5
          flag='*';
        else
          flag='';
        end
        if i==1 || c.fold(ind)>=args.minfold
          target=c.target;
          fprintf('%5d %-19.19s %5s: %5.2f  [%5.2f,%5.2f] %3.3s',c.mixture,target,concfmt(c.conc), c.fold(ind), c.foldci(ind,:),flag);
          if ~isempty(args.extra) && isKey(args.extra,target)
            fprintf(' %s',args.extra(target));
          end
          fprintf('\n');
        end
      end
    end
    
    function computefolds(summary,varargin)
    % Compute fold changes between +target and -target conditions
    % Match up +ves and -ves from same run, comparing individual subruns
    % Bootstrap those results (random reselecting) to give overall fold bootstrap samples
      defaults=struct('mincnt',10,'ci',normcdf([-1,1])*100,'metric','ratio_of_ratio','nboot',200,'useallnegatives',false);
      args=processargs(defaults,varargin);
      summary.metric=args.metric;
      summary.cirange=args.ci;

      c0=summary.cleavage(1);
      obj=NGSVitro('tmp');
      for i=1:length(summary.cleavage)
        c=summary.cleavage(i);
        r=summary.cleavage(i).ratiosc;
        foldsampsbyrun=single([]);
        for k=1:length(c.trpexpt)
          if args.useallnegatives
            summary.cleavage(i).negatives{k}=1:length(c0.trpexpt);
            c=summary.cleavage(i);
          elseif isempty(c.negatives{k})
            sel0=[];   % runs within c0 to use as -target
            for j=1:length(c0.trpexpt)
              if strcmp(c0.runs{j},c.runs{k})
                sel0(end+1)=j;
              end
            end
            if isempty(sel0)
              fprintf('No negatives in run %s, compare with %s.%d\n', c.runs{k}, c0.runs{end}, c0.trpexpt(end));
              sel0=length(c.runs);
            end
            summary.cleavage(i).negatives{k}=sel0;
            c=summary.cleavage(i);
          end
          assert(~isempty(c.negatives{k}));
          % Average over the negatives on a per-bootsample basis
          scratioboot=arrayfun(@(z) c0.ratioboot{z}*c0.scaling(z),c.negatives{k},'Unif',false);
          negboot=exp(nanmean(log(cat(3,scratioboot{c.negatives{k}})),3));
          fk=[];
          for m=1:length(summary.naseq)
            [~,~,fk(m,:)]=obj.bootfold(c.ratioboot{k}(m,:)*c.scaling(k),negboot(m,:),'metric',args.metric,'nboot',args.nboot);
          end
          
          assert(size(fk,2)==args.nboot);
          foldsampsbyrun(:,:,k)=fk;
        end
        lfoldsampsbyrun=log(foldsampsbyrun);
        summary.cleavage(i).foldbyrun=exp(squeeze(nanmean(lfoldsampsbyrun,2)));
        summary.cleavage(i).foldcibyrun=squeeze(prctile(foldsampsbyrun,args.ci,2));
        if size(foldsampsbyrun,3)==1
          lfoldsamps=log(foldsampsbyrun);
        else
          summary.cleavage(i).foldsampsbyrun=foldsampsbyrun;
          lfoldsamps=single([]);
          for k=1:size(lfoldsampsbyrun,1)
            f=squeeze(lfoldsampsbyrun(k,:,:));
            lfoldsamps(k,:)=nanmean(reshape(randsample(f(:),length(f(:))),size(f)),2);
          end
        end
        summary.cleavage(i).foldsamps=exp(lfoldsamps);
        assert(strcmp(class(summary.cleavage(i).foldsamps),'single'));
        summary.cleavage(i).fold=exp(nanmean(lfoldsamps,2));
        summary.cleavage(i).foldci=exp(prctile(lfoldsamps,args.ci,2));
      end
    end
    
    function comparefold(summary,ind1,ind2,varargin)
      defaults=struct('alpha',.05,'newfig',true,'list',true,'minfold',0);
      args=processargs(defaults,varargin);

      if args.newfig
        setfig('comparefold');clf;
      end
      c1=summary.cleavage(ind1);
      c2=summary.cleavage(ind2);
      sel=c1.fold>=args.minfold | c2.fold>=args.minfold;
      f1=c1.fold(sel);f1ci=c1.foldci(sel,:);f1s=c1.foldsamps(sel,:);
      f2=c2.fold(sel);f2ci=c2.foldci(sel,:);f2s=c2.foldsamps(sel,:);
      loglog(f1,f2,'ob');
      hold on;
      for i=1:length(f1)
        plot(f1ci(i,:),f2(i)*[1,1],'-g');
        plot(f1(i)*[1,1],f2ci(i,:),'-g');
      end
      loglog(f1,f2,'ob');
      axis equal
      axis tight
      ax=axis;
      plot(ax(1:2),ax(1:2),'r:');
      axis tight
      logticks(1,1);
      xlabel(sprintf('%s fold',c1.target));
      ylabel(sprintf('%s fold',c2.target));
      title(sprintf('%s vs %s',c1.target,c1.target));
      cmp=[];
      for i=1:length(f1)
        if ~all(isfinite([f1ci(i,:),f2ci(i,:)]))
          cmp(i)=0.5;
          continue;
        end
        sel1=isfinite(f1s(i,:)); sel2=isfinite(f2s(i,:));
        cmp(i)=nanmean(nanmean(f1s(i,sel1)>f2s(i,sel2)'));
      end
      [~,ord]=sort(cmp);
      for ii=1:length(cmp)
        i=ord(ii);
        if cmp(i)<=args.alpha || cmp(i)>=1-args.alpha
          if args.list
            fprintf('%2d %9.9s  [%3.1f,%3.1f] vs [%3.1f,%3.1f] p=%6.4f\n', i, summary.names{i}, f1ci(i,:), f2ci(i,:), cmp(i));
          end
          plot(f1(i),f2(i),'or');
        end
      end
    end
    
    function sort(summary)
    % Sort by target, then concentration increasing
      c=summary.cleavage;
      keys={};
      for j=1:length(c)
        if isempty(c(j).target)
          keys{j}='';
        else
          keys{j}=sprintf('%s@%.12f',c(j).target,c(j).conc);
        end
      end
      [~,ord]=sort(keys);
      summary.cleavage=c(ord);
    end
    
    
    function plotcompare(summary,refind)
      setfig('TRPSummaryCompare');clf;
      c=summary.cleavage;
      nrows=ceil(sqrt(length(c))*3/4);
      ncols=ceil((length(c))/nrows);
      cur=1;
      for i=1:length(c)
        subplot(nrows,ncols,cur);
        cur=cur+1;
        loglog(c(refind).ratioavg,c(i).ratiosc,'.');
        cleaveticks(1,1);
        ax=axis;
        hold on;
        plot(ax(1:2),ax(1:2),':');
        if isempty(c(refind).target)
          xlabel('No target');
        else
          xlabel(sprintf('%s@%s',c(refind).target,concfmt(c(refind).conc)));
        end
        if isempty(c(i).target)
          title('No target');
        else
          title(sprintf('%s@%s',c(i).target,concfmt(c(i).conc)));
        end
        if length(c(i).runs)>1
          legend(c(i).runs,'Location','NorthWest');
        end
        axis tight
      end
      pause(0.01);
    end
    
    function calcstats(summary)
    % Compute statistics of repeated measurements
      c=summary.cleavage(1);
      r=c.ratiosc;   % Scaled ratio for -target
      n=sum(isfinite(r),2);
      r=r(n>=min(5,max(n)),:);
      lr=log10(r);
      s=nanstd(lr,0,2);
      m=nanmean(lr,2);
      setfig('summarystats');clf;
      if exist('tiledlayout')
        tiledlayout('flow');
        nexttile;
      end
      [m,ord]=sort(m);
      s=s(ord);
      loglog(10.^m,10.^s,'o');
      xlabel('Mean Ratio');
      ylabel('Multiplicative Std');
      title(sprintf('Variation as a function of ratio (mean std=%.2f)',10.^(mean(s))));
      ax=axis;ax(3)=1;axis(ax);
      logticks(1,0);
      summary.stdLogRatio=median(s);
      summary.DFLogRatio=round(mean(n))-1;

      if exist('tiledlayout')
        nexttile;
      end
      s=[];
      for i=1:length(c.runs)
        s(i)=sqrt(nanmean(log10(c.ratiosc(:,i)./c.ratioavg).^2));
        lbls{i}=sprintf('%s.%d',c.runs{i},c.trpexpt(i));
      end
      bar(10.^s);
      set(gca,'YScale','log');
      ax=axis;ax(3)=1;axis(ax);
      ylabel('Multiplicative Std');
      set(gca,'XTick',1:length(lbls));
      set(gca,'XTickLabels',lbls);
      set(gca,'XTickLabelRotation',90);
      title('Variation as a function of run');
    end
    
    function d=getdesc(summary,i,withconc,usecompound)
    % Get description of conditions for cleavage(i)
      if nargin<3
        withconc=true;
      end
      if nargin<4
        usecompound=true;
      end
      c=summary.cleavage(i);
      if isempty(c.target)
        d='-target';
      else
        if length(c.contains)==1 && usecompound
          co=Compounds.instance().get(c.contains);
          nm=co.name;
          m=Mixtures.instance().get(c.mixture);
          if length(m.name)==4
            % Special name; e.g. 8740
            nm=sprintf('%s(%s)',nm,m.name);
          end
        else
          m=Mixtures.instance().get(c.mixture);
          nm=m.name;
        end
        if withconc
          d=sprintf('%s@%s',nm,concfmt(c.conc,[],1));
        else
          d=nm;
        end
      end
    end
    
    function getcoverage(summary)
      for i=1:length(summary.cleavage)
        n=sum(isfinite(summary.cleavage(i).rawratio),2);
        fprintf('%25.25s: %.1f\n', summary.getdesc(i), median(n(n>0)));
      end
    end
    
    function x=summarytable(summary)
    % Create a table based on summary data
      x=table();
      x.naseq=int32(summary.naseq);
      x.Properties.VariableDescriptions{1}='ID#';
      c=summary.cleavage;
      for i=1:length(c)
        r=c(i).ratiosc
        if size(c(i).ratiosc,2)>=2
          % If more than 1, return individual points, mean, mean-sem, mean+sem
          r(:,end+1)=c(i).ratioavg;
          rstd=nanstd(log(c(i).ratiosc),[],2);
          rsem=rstd/size(c(i).ratiosc,2);   % FIXME: This is not right since using the sample std is a student t-distribution, not normal
          r(:,end+1)=exp(log(c(i).ratioavg)-rsem);
          r(:,end+1)=exp(log(c(i).ratioavg)+rsem);
        end
        cl=r./(1+r)*100;
        x.(sprintf('Run%d',i))=round(cl,0);
        if isempty(c(i).target)
          descr=sprintf('-target');
        else
          descr=strrep(sprintf('+%s@%s',c(i).target,concfmt(c(i).conc)),' ','');
        end
        x.Properties.VariableDescriptions{i+1}=descr;
        vname=strrep(strrep(strrep(strrep(descr,'-','_'),'.','_'),'+',''),'@','_');
        if vname(1)=='_'
          vname=['no',vname];
        end
        x.Properties.VariableNames{i+1}=vname;
      end
    end

    function [t,ord]=targetsort(summary,t)
    % Sort target names in a meaningful way
      keys=t;
      for i=1:length(keys)
        k=keys{i};
        % Change each numeric string into a 3-digit leading 0 string
        isnum=[false,k>='0' & k<='9',false];
        startnum=find(isnum(2:end) & ~isnum(1:end-1));
        endnum=find(~isnum(2:end) & isnum(1:end-1))-1;
        assert(length(startnum)==length(endnum));
        for j=length(startnum):-1:1
          v=str2num(k(startnum(j):endnum(j)));
          k=[k(1:startnum(j)-1),sprintf('%04d',v),k(endnum(j)+1:end)];
        end
        if any(k=='-')
          k=['-',k];
        end
        keys{i}=k;
      end
      [~,ord]=sort(keys);
      t=t(ord);
    end
    
    function x=plotgradsummary(summary,obj,varargin)
    % Generate a pseudocolor plot with sequences vs expts (summary conditions); switching ratio is shown in color
    % The first expt is assumed to be the -target condition
      defaults=struct('naseqs',[],'expts',[],'onlyvalid',true,'showids',false,'shownaseqs',true,'metric','ratio_of_uncleavage','aptprefixes',[],'ignorenaseqs',[],'minswitching',2,'sort','name','onlyconc',[],'alltargets',false,'maxhits',10000,'mintargets',[],'divideexpts',false,'showexptnames',true);
      args=processargs(defaults,varargin);

      if ~isempty(args.naseqs)
        naseqs=args.naseqs;
      else
        % Use all analyzed sequences
        naseqs=summary.naseq;
      end
      
      if isempty(args.expts)
        args.expts=find(~strcmp({summary.cleavage.target},''));
      elseif ischar(args.expts)
        args.expts=find(summary.exptsel(args.expts));
      elseif islogical(args.expts)
        args.expts=find(args.expts);
      end
      if ~isempty(args.onlyconc)
        args.expts=args.expts(ismember([summary.cleavage(args.expts).conc],args.onlyconc));
      end
      
      if strcmp(args.sort,'tree')
        if isempty(summary.tree) || isempty(summary.tree.useseqs)
          error('Attempt to sort by tree, but tree not built');
        end
        naseqs=summary.naseq(summary.tree.useseqs);
        naseqs=naseqs(summary.tree.leaforder);
        [~,ord]=sort(summary.names);
        n2=summary.naseq(ord);
        naseqs=union(naseqs,n2,'stable');   % Append any other naseqs, sorted by name
      elseif strcmp(args.sort,'name')
        [~,ord]=sort(summary.names);
        n2=summary.naseq(ord);
        naseqs=intersect(n2,naseqs,'stable');
      elseif strcmp(args.sort,'none') || strcmp(args.sort,'switching')
        ;  % Whatever order naseqs are in
      elseif strcmp(args.sort,'target');
        ;  % Handled later
      else
        error('Unknown sort key: %s (expected "name" or "tree")',args.sort);
      end

      if args.onlyvalid
        ignoreseqs=[obj.seqs.tags(ismember({obj.seqs.tags.name},{'ref','amplicon','invalid','grz'})).naseq];
        ignoreseqs=union(ignoreseqs,[54029,1613602,202280]);  % Remove sTRSV, sTRSVCtl
        nremove=length(intersect(ignoreseqs,naseqs));
        if nremove>0
          fprintf('Removing %d/%d invalid (or sTRSV) sequences\n', nremove,length(naseqs));
          naseqs=setdiff(naseqs,ignoreseqs,'stable');
        end
      end
      if ~isempty(args.ignorenaseqs)
        naseqs=setdiff(naseqs,args.ignorenaseqs,'stable');
      end
      % Only use ones in summary
      naseqs=intersect(naseqs,summary.naseq,'stable');
      
      [~,targetord]=summary.targetsort({summary.cleavage(args.expts).target});
      args.expts=args.expts(targetord);
      
      % Remove everything after '-' from names
      names={};
      for i=1:length(naseqs)
        names{i}=summary.names{summary.naseq==naseqs(i)};
      end
      
      aptprefixes=names;
      ids=names;
      for i=1:length(aptprefixes)
        t=aptprefixes{i};
        lastpos=find(t=='-',1);
        if ~isempty(lastpos)
          aptprefixes{i}=t(1:lastpos-1);
          ids{i}=t(lastpos+1:end);
        else
          aptprefixes{i}='?';
          ids{i}=t;
        end
      end
      if ~strcmp(args.sort,'name')
        ids=names;
      end
      if ~isempty(args.aptprefixes)
        sel=ismember(aptprefixes,args.aptprefixes);
        aptprefixes=aptprefixes(sel);
        fprintf('Keeping %d/%d seqs that are in aptprefixes list\n', sum(sel), length(sel));
        names=names(sel);
        ids=ids(sel);
        naseqs=naseqs(sel);
      end

      % Collect fold change data
      switching=nan(length(args.expts),length(naseqs));
      for i=1:length(args.expts)
        t=summary.cleavage(args.expts(i));
        for j=1:length(naseqs)
          sel=summary.naseq==naseqs(j);
          if sum(sel)==1
            switching(i,j)=t.fold(sel);
          end
        end
      end
      
      if ~args.alltargets
        keepexpts=max(switching,[],2)>=args.minswitching;
        fprintf('Keeping %d/%d experiments that have switching >= %.1f\n',sum(keepexpts),length(keepexpts),args.minswitching);
        switching=switching(keepexpts,:);
        args.expts=args.expts(keepexpts);
      end
      
      ord=1:length(naseqs);
      if strcmp(args.sort,'switching')
        [~,ord]=sort(max(switching,[],1),'desc');
      elseif strcmp(args.sort,'target')
        % Sort by target sensitivity
        maxtgt=zeros(size(naseqs));
        stmp=switching; stmp(isnan(stmp))=0;
        for j=1:length(naseqs)
          [~,ord]=sort(stmp(:,j),'desc');
          for k=1:min(3,length(ord))
            maxtgt(j)=maxtgt(j)*100+ord(k);
          end
        end
        [~,ord]=sort(maxtgt,'desc');
      end

      naseqs=naseqs(ord);
      names=names(ord);
      aptprefixes=aptprefixes(ord);
      ids=ids(ord);
      switching=switching(:,ord);
      
      if isempty(args.mintargets)
        args.mintargets=length(args.expts)/2;
      end

      seqsel=sum(isfinite(switching))>=args.mintargets & nanmax(switching)>args.minswitching & sum(switching>args.minswitching)<=args.maxhits;
      if args.maxhits<=length(args.expts)
        fprintf('Keeping %d/%d sequences that have data on at least %.0f targets and show at least %.1f switching and have <= %d hits\n',sum(seqsel), length(seqsel), args.mintargets, args.minswitching,args.maxhits);
      else
        fprintf('Keeping %d/%d sequences that have data on at least %.0f targets and show at least %.1f switching\n',sum(seqsel), length(seqsel), args.mintargets, args.minswitching);
      end
      switching=switching(:,seqsel);
      naseqs=naseqs(seqsel);
      aptprefixes=aptprefixes(seqsel);
      ids=ids(seqsel);
      names=names(seqsel);
      data=switching; data(end+1,:)=nan; data(:,end+1)=nan;

      setfig(sprintf('gradsummary %d',length(summary.naseq)));clf;
      pcolor(log(data));
      cax=caxis;
      cax(1)=0;
      cax(2)=min(log(5),cax(2));
      caxis(cax);
      h=colorbar;
      set(get(h,'Label'),'string','Switching')
      set(h, 'Ticks',log([1,2,3,4,5,10,20,30,50,100]));
      set(h,'TickLabels',arrayfun(@(z) sprintf('%.0f',exp(z)), get(h,'Ticks'),'UniformOutput',false)); 
      exptnames={};
      for i=1:length(args.expts)
        t=summary.cleavage(args.expts(i));
        if isempty(t.target)
          exptnames{end+1}='no target';
        elseif length(args.onlyconc)==1
          exptnames{end+1}=summary.getdesc(args.expts(i),false);
        else
          exptnames{end+1}=summary.getdesc(args.expts(i),true);
        end
      end
      if args.showexptnames
        set(gca,'YTick',(1:length(args.expts))+0.5);
        set(get(gca,'YAxis'),'FontSize',8);
        set(gca,'YTickLabel',exptnames);
        set(gca,'TickLabelInterpreter','None');
        set(get(gca,'YAxis'),'FontSize',max(6,min(10,800/length(exptnames))));
      end
      
      shading flat;
      ax=axis;
      hold on;
      if strcmp(args.sort,'name')
        lasti=1;
        ticks=[]; ticklabel={};
        for i=2:length(naseqs)
          if ~strcmp(aptprefixes{i},aptprefixes{i-1})
            plot(i*[1,1],ax(3:4),'k-','LineWidth',1.5);
            ticks(end+1)=(lasti+i)/2;
            ticklabel{end+1}=aptprefixes{i-1};
            lasti=i;
          end
        end
        ticks(end+1)=(lasti+length(naseqs)+1)/2;
        ticklabel{end+1}=aptprefixes{end};
        set(gca,'XTick',ticks);
        set(gca,'XTickLabel',ticklabel);
        set(gca,'XTickLabelRotation',45);
      end

      if args.showids
        % Add IDs along top
        for i=1:length(ids)
          text(i+0.5,length(args.expts)+1.1,ids{i},'Rotation',90,'HorizontalAlignment','Left','VerticalAlignment','middle','FontSize',max(6,min(10,1000/length(ids))));
        end
      elseif args.shownaseqs
        for i=1:length(ids)
          text(i+0.5,length(args.expts)+1.1,sprintf('%d',naseqs(i)),'Rotation',90,'HorizontalAlignment','Left','VerticalAlignment','middle','FontSize',max(6,min(10,1000/length(ids))));
        end
      end
        
      if args.divideexpts
        etmp=cellfun(@(z) z(z<'0' | z>'9'),exptnames,'UniformOutput',false);  % Remove any numbers for expt names
        etmp=cellfun(@(z) regexprep(z,'^R[A-H].*','R'),etmp,'UniformOutput',false);  % Remove any row labels
        for i=2:length(etmp)
          if ~strncmp(etmp{i},etmp{i-1},4)
            plot(ax(1:2),i*[1,1],'k-','LineWidth',1.5);
          end
        end
      end
      xlabel(' ');  % Pushes up bottom margin
      if length(args.onlyconc)==1
        ylabel(sprintf('Ligand (@%.1f \\mu{}M)',args.onlyconc*1e6),'Interpreter','tex');
      else
        ylabel('Ligand @ Conc (\mu{}M)','Interpreter','tex');
      end
      xlabel('Sensor');
      x=struct();
      x.naseq=naseqs;
      x.switching=switching;
      x.exptnames=exptnames;
      x.details=summary.cleavage(args.expts);
      x.allnaseqs=summary.naseq;
    end

    function expect=runmodel(summary,x,conc)
      expect=x(1)+x(2).*conc./(conc+x(3));
    end
    
    function [kd,rmin,rmax]=kdfit(summary,conc,ratio)
      lratio=log(ratio);
      x0=[max(lratio),-(max(lratio)-min(lratio)),mean(conc)];
      options=optimset('display','notify','TolX',1e-8,'TolFun',1e-8);
      x=fminsearch(@(x) sum((summary.runmodel(x,conc)-lratio).^2),x0, options);
      kd=x(3);
      rmax=x(1);
      rmin=x(1)+x(2);
      for i=1:length(conc)
        fprintf('%.1f@%.2g ',ratio(i),conc(i));
      end
      fprintf(' -> ');
      fprintf('kd=%.1g, rmin=%.2f, rmax=%.2f\n', kd, rmin, rmax);
      clf;
      loglog(conc*1e6,ratio); hold on;
      ctmp=logspace(log10(min(conc(conc>0))/10),log10(max(conc)*10));
      plot(ctmp*1e6,exp(summary.runmodel(x,ctmp)),':');
      ax=axis;
      plot(kd*[1,1]*1e6,ax(3:4),'r:');
      cleaveticks(0,1);
      xlabel('Conc (uM)');
      ylabel('Cleavage');
      figure(gcf);
    end
    
    function [kd,subweights]=fit(summary,naseq,varargin)
    % Build model for given seq
      defaults=struct('total','HTSLib','subparts',{{'C','R','P'}});
      args=processargs(defaults,varargin);
      targets={summary.cleavage.target};
      targetforms=regexprep(targets,'[0-9]','');
      total=strncmp(targets,args.total,length(args.total))';
      total(1)=true;   % Include no-target
      subparts=[];
      for i=1:length(args.subparts)
        subparts(:,i)=strcmp(targetforms,args.subparts{i});
      end
      naseqsel=summary.naseq==naseq;
      conc=[summary.cleavage(total).conc];
      ratio=arrayfun(@(z) z.ratioavg(naseqsel),summary.cleavage(total));
      [kd,fmax]=summary.kdfit(conc,ratio);
    end

   function [x,ratio,ratioci,fold,foldci]=plotgradient(summary,obj,target,varargin)
   % For each of the naseqs given, produce a bar graph across the conditions for which we have data
   % Pass in an NGS object (obj) for accessing seqs 
   % trpexpts are indices into summary.cleavage() 
      defaults=struct('naseqs',[],'mincnt',10,'showall',false,'units','M','minratio',2,'forcenamed',[],...
                      'oneframe',0,'longlegend',0,'expts',[],'title','','newfig',true);
      args=processargs(defaults,varargin);
      %obj.trpanalyze('mincnt',args.mincnt);

      naseqs=summary.naseq;
      if ~isempty(args.forcenamed)
        % Add any naseq that starts with forcenamed
        sel=strncmp({obj.seqs.seqnames.name},args.forcenamed,length(args.forcenamed));
        forced=[obj.seqs.seqnames(sel).naseq];
      else
        forced=[];
      end

      % Locate relevant entries in summary
      selminus=find(strcmp({summary.cleavage.target},''));
      if isempty(args.expts)
        args.expts=find(strcmp({summary.cleavage.target},target));
      end
      % Make sure negative is not in expts list
      args.expts=setdiff(args.expts,selminus);
      if isempty(args.expts)
        error('No data available with target %s',target);
      end
      
      % Verify that there is some overlap in the contains set of all expts
      allcontain=summary.cleavage(args.expts(1)).contains;
      for i=2:length(args.expts)
        allcontain=intersect(allcontain,summary.cleavage(args.expts(i)).contains);
      end
      if isempty(allcontain)
        error('No overlap of compounds in set of expts\n');
      end
      compounds=Compounds.instance();
      fprintf('Expts all contain: [%s]\n', strjoin(arrayfun(@(z) compounds.get(z).name,allcontain,'Unif',false),','));
      % Order by increase conc
      [~,ord]=sort([summary.cleavage(args.expts).conc]);

      % Keep -target condition first in list
      sel=[selminus,args.expts(ord)];
      
      % Collect data for plot in ratio(naseq,expt),ratioci(naseq,expt,2), fold(naseq,expt), foldci(naseq,expt,2)
      ratio=nan(length(naseqs), length(sel));
      ratioci=nan(length(naseqs), length(sel),2);
      fold=nan(length(naseqs), length(sel));
      foldci=nan(length(naseqs), length(sel),2);
      cnt=nan(length(naseqs), length(sel));
      fold(:,1)=1.0;
      foldci(:,1,:)=1.0;

      r0=summary.cleavage(selminus).ratiosc;

      for ii=1:length(sel)
        i=sel(ii);
        %r=summary.cleavage(i).ratiosc;
        ratio(:,ii)=summary.cleavage(i).ratioavg;
        %ratioci(:,ii,:)=prctile(summary.cleavage(i).ratiosc,summary.cirange,2); % TODO: This is WRONG; needs to be done over the scaled bootstrap samples
        allboot=[];
        for j=1:length(summary.cleavage(i).ratioboot)
          allboot=[allboot,summary.cleavage(i).ratioboot{j}*summary.cleavage(i).scaling(j)];
        end
        ratioci(:,ii,:)=prctile(allboot,summary.cirange,2); 
        
        fold(:,ii)=summary.cleavage(i).fold;
        foldci(:,ii,:)=summary.cleavage(i).foldci;
        cnt(:,ii)=sum(sum(summary.cleavage(i).cnt(:,2:3,:),3),2);
        % ex=1;fprintf('sel(%d)=%d,naseq %d:ratio=%f[%f,%f],fold=%f[%f,%f],cnt=%d\n',ii,i,ex,ratio(ex,ii),ratioci(ex,ii,:),fold(ex,ii),foldci(ex,ii,:),cnt(ex,ii));
      end
      
      tgtconc=[summary.cleavage(sel).conc];  % Conc in nM
      if length(unique(tgtconc))<length(tgtconc)
        fprintf('Have multiple expts with the same concentration, averaging.\n');
        utgtconc=unique(tgtconc);
        for j=1:length(utgtconc)
          sel=tgtconc==utgtconc(j);
          uratio(:,j)=nanmean(ratio(:,sel),2);
          uratioci(:,j,:)=nanmean(ratioci(:,sel,:),2);
          ufold(:,j)=nanmean(fold(:,sel),2);
          ufoldci(:,j,:)=nanmean(foldci(:,sel,:),2);
          ucnt(:,j)=nansum(cnt(:,sel),2);
        end
        ratio=uratio;
        ratioci=uratioci;
        fold=ufold;
        foldci=ufoldci;
        cnt=ucnt;
        tgtconc=utgtconc;
      end
      if strcmp(args.units,'\muM')
        tgtconc=tgtconc*1e6;
      elseif strcmp(args.units,'nM')
        tgtconc=tgtconc*1e9;
      elseif ~strcmp(args.units,'M')
        error('Bad units: %s\n', args.units);
      end
      tgtconc(isnan(tgtconc))=0;
      
      % Compute max fold change for each naseq
      fmax=[];
      for i=1:length(naseqs)
        ftmp=fold(i,isfinite(fold(i,:)));
        if isempty(ftmp)
          fmax(i)=nan;
        else
          fmax(i)=ftmp(end);
        end
      end
      
      if ~args.showall
        if ~isempty(args.naseqs)
          keep=ismember(naseqs,args.naseqs);
        else
          % Only keep sequences that show switching
          keep=false(size(naseqs));
          for i=1:length(naseqs)
            keep(i)=fmax(i)>=args.minratio && sum(isfinite(ratio(i,:)))>1;
            keep(i)=keep(i)&min(cnt(i,:))>=args.mincnt;
            if ~keep(i) && ismember(naseqs(i),forced) && sum(isfinite(ratio(i,:)))>1
              fprintf('Displaying %d (%s) even though it is switching with ratio %.1f < %.1f\n', naseqs(i), obj.seqs.getname(naseqs(i)),max(fold(i,:)),args.minratio);
              keep(i)=true;
            end
          end
        end
        fprintf('Keeping %d/%d sequences\n', sum(keep), length(keep));
        if ~any(keep)
          error('No sequences kept\n');
        end
        naseqs=naseqs(keep);
        ratio=ratio(keep,:);
        ratioci=ratioci(keep,:,:);
        fold=fold(keep,:);
        foldci=foldci(keep,:,:);
        fmax=fmax(keep);
      end
      
      % Title
      if isempty(args.title)
        ti=['Cleavage Gradient: ',target];
      else
        ti=args.title;
      end
      tgtconcnz=tgtconc;
      tgtconcnz(tgtconc==0)=min(tgtconc(tgtconc>0))/10;
      jitter=(max(tgtconcnz)/min(tgtconcnz))^(1/50);  % Jitter by 1/50 of plot width
      if args.newfig
        setfig(ti);clf;
      end
      for sp=1:2
        if args.oneframe==0
          subplot(1,2,sp);
        elseif args.oneframe~=sp
          continue;
        end
        if sp==1
          val=ratio;
          ci=ratioci;
        else
          val=fold;
          ci=foldci;
        end
        
        h=loglog(tgtconcnz,val','o-');
        %set(gca,'XTick',tgtconcnz);
        %set(gca,'XTickLabel',arrayfun(@(z) regexprep(regexprep(sprintf('%.3f',z),'0*$',''),'\.$',''),tgtconc,'UniformOutput',false));
        hold on;
        % Connect lines with breaks due to NaNs
        for i=1:length(naseqs)
          for j=2:length(tgtconc)-1
            if isnan(val(i,j)) & ~isnan(val(i,j-1)) & ~isnan(val(i,j+1))
              plot(tgtconcnz([j-1,j+1]),val(i,[j-1,j+1]),'-','Color',get(h(i),'Color'));
            end
          end
        end
        
        % Add error bars
        for i=1:length(naseqs)
          if length(naseqs)==1
            dodge=1;   % No jitter if only 1 plot
          else
            dodge=exp((((i-1)/(length(naseqs)-1))*2-1)*log(jitter));
          end
          for j=1:length(tgtconc)
            plot(tgtconcnz(j)*[1,1]*dodge,squeeze(ci(i,j,:)),'r-','Color',get(h(i),'Color'));
            plot(tgtconcnz(j)*[.95,1.05]*dodge,ci(i,j,1)*[1,1],'r-','Color',get(h(i),'Color'));
            plot(tgtconcnz(j)*[.95,1.05]*dodge,ci(i,j,2)*[1,1],'r-','Color',get(h(i),'Color'));
          end
        end
        xlabel(sprintf('[%s] (%s)',target, args.units));
        if sp==2 || args.oneframe~=0
          % Add legend
          for i=1:length(naseqs)
            nm=obj.seqs.getname(naseqs(i));
            conc2=nanmin(tgtconcnz(val(i,:)>=2))*1e9;
            if args.longlegend
              leg{i}=sprintf('%d %3.1f %.0fnM %s',naseqs(i), fmax(i), conc2, obj.seqs.getlabels(naseqs(i)));
            elseif ~isempty(nm)
              leg{i}=sprintf('%s %3.1f %.0fnM',nm, fmax(i), conc2);
            else
              leg{i}=sprintf('%d %3.1f %.0fnM',naseqs(i), fmax(i), conc2);
            end
          end
          [~,ord]=sort(fmax,'descend');
          ord=ord(1:min(end,20));
          if args.oneframe
            legend(h(ord),leg(ord),'Location','EastOutside');
          else
            legend(h(ord),leg(ord),'Location','Best');
          end
        end
        % Give a little more room at sides
        c=axis;
        c(1)=min(tgtconcnz)/jitter^2;
        c(2)=max(tgtconcnz)*jitter^2;
        if sp==1
          c(4)=max([c(4),ci(:)']);
        else
          c(3)=0.2;
          c(4)=max([ci(:)',30]);
        end
        axis(c);
        if sp==1
          ylabel('Fraction Cleaved');
          cleaveticks(0,1);
        else
          ylabel('Fold Change of Cleavage','Interpreter','none');
          logticks(0,1);
        end
      end
      if args.oneframe==0
        sti=suptitle(ti);
      else
        sti=title(ti);
      end
      set(sti,'Interpreter','none');
      x=[];
      for i=1:length(naseqs)
        for j=2:size(fold,2)
          x=[x,struct('name',obj.seqs.getname(naseqs(i)),'naseq',int64(naseqs(i)),'target',target,'conc',tgtconc(j),'fold',round(fold(i,j),2),'foldcilow',round(foldci(i,j,1),2),'foldcihigh',round(foldci(i,j,2),2))];
        end
      end
    end
    
    function [data,concs,runs,seqsel]=plotrundependence(summary,obj,target,varargin)
    % Check dependency against run by averaging across sequences that are sensitive to the target
    % TODO- could instead use statistical tests based on cnt to see if cleavage is different for this target across runs
      defaults=struct('minfold',2,'fulldata',true);
      args=processargs(defaults,varargin);

      exptsel=[1,find(strcmp({summary.cleavage.target},target))];
      concs=[summary.cleavage(exptsel).conc];
      [concs,ord]=sort(concs);
      exptsel=exptsel(ord);
      runs={};
      for i=2:length(exptsel)  % Only runs with target  
        runs=union(runs,summary.cleavage(exptsel(i)).runs);
      end
      runs=sort(runs);
      
      seqsel=find(summary.cleavage(exptsel(end)).fold>=args.minfold);
      fprintf('Averaging across %d seqs\n', length(seqsel));
      if length(seqsel)==0
        error('No sequences with fold change > %f for target %s', args.minfold, target);
      end

      data=nan(length(exptsel),length(runs),length(seqsel));
      for i=1:length(exptsel)
        c=summary.cleavage(exptsel(i));
        for k=1:length(runs)
          runsel=find(strcmp(c.runs,runs{k}));
          if ~isempty(runsel)
            data(i,k,:)=mean(c.ratiosc(seqsel,runsel),2);
          end
        end
      end
      allfinite=all(any(isfinite(data)));
      if args.fulldata
        fprintf('Keeping %d/%d seqs with full data\n', sum(allfinite), length(allfinite));
        seqsel=seqsel(allfinite);
        data=data(:,:,allfinite);
      end
      for i=1:length(seqsel)
        r=nanmean(data(:,:,i),2);
        c=r./(r+1);
        %fprintf('%d: [%s] %s \n', summary.naseq(seqsel(i)),sprintf('%2.0f ',c*100),obj.seqs.getlabels(summary.naseq(seqsel(i))));
      end
      setfig(['rundependence - ',target]);clf;
      for i=1:length(runs)
        v=nanmean(data(:,i,:),3);
        sel=isfinite(v);
        loglog(max(concs(sel),concs(2)/10)*1e6,v(sel),'o-');
        hold on;
      end
      xlabel('Concentration (nM)');
      legend(runs);
      cleaveticks(0,1);
      logticks(1,0);
      xticks=get(gca,'XTick');
      xticklabels=get(gca,'XTickLabel');
      keepticks=xticks>concs(2)/10*1e6*1.1;
      xticks=[concs(2)/10*1e6,xticks(keepticks)];
      xticklabels={'0',xticklabels{keepticks}};
      set(gca,'XTick',xticks);
      set(gca,'XTickLabel',xticklabels);

      ylabel('Cleavage');
      title(sprintf('Average Cleavage for %s over %d Sensors',target, length(seqsel)));
    end
    
    function dendrogram(summary,varargin)
      defaults=struct();
      args=processargs(defaults,varargin);

      ti=sprintf('Dendrogram(%s)',summary.tree.method);
      setfig(ti);clf;
      nleafs=sum(summary.tree.linkages(:,3)>0)+1;   % Collapse 0-distance leaves
      dendrogram(summary.tree.linkages,nleafs,'Labels',summary.tree.labels,'Reorder',summary.tree.leaforder);
      set(gca,'XTickLabelRotation',90)
      title(ti);
    end
    
    function sel=exptsel(summary,expts,conc)
      if strcmp(expts,'all')
        sel=true(length(summary.cleavage),1);
      elseif strncmp(expts,'singles',length(expts))
        sel=arrayfun(@(z) length(z.contains)==1,summary.cleavage);
      elseif strcmp(expts,'V256')
        sel=arrayfun(@(z) length(z.contains)==256,summary.cleavage);
      elseif strcmp(expts,'CDIQ')
        sel=arrayfun(@(z) strncmp(z.target,'CDIQ',4),summary.cleavage);
      else
        error('Expected "all","singles","CDIQ", or "V256" for expts arg');
      end
      if nargin>=3 && ~isempty(conc)
        sel=sel & ismember([summary.cleavage.conc],conc);
      end
      if ~any(sel)
        fprintf('exptsel(%s) had no results\n', expts);
      end
    end
    
    function clusterbytarget(summary,varargin)
    % Cluster the results by pattern of target activation
      defaults=struct('method','average','expts','all','minfold',1.5);
      args=processargs(defaults,varargin);

      useseqs=isfinite(summary.cleavage(1).ratioavg);   % Only ones that we have at least -target data for
      useseqs=useseqs & summary.getminimumfold()>args.minfold;   % And are certain to have at least minfold fold for some condition (using lower bound of CI)
      fprintf('Running clusterbytarget using %d/%d seqs that have at least one condition with at least %.1f fold\n', sum(useseqs), length(useseqs),args.minfold);
      x=summary.compareseqs(summary.naseq(useseqs),'showdiffs',false,'expts',args.expts);
      dist=squareform(x.dpmax);
      linkages=linkage(dist,args.method);  
      leaforder=optimalleaforder(linkages,dist,'criteria','group');
      labels=arrayfun(@(z) sprintf('%d',z),find(useseqs),'Unif',false);
      summary.tree=struct('method',args.method,'dist',dist,'compare',x,'linkages',linkages,'leaforder',leaforder,'labels',{labels},'useseqs',useseqs,'useconditions',x.expts);
    end
    
    function x=compareseqs(summary,naseqs,varargin)
    % Compare seqs across all conditions
      defaults=struct('thresh',1,'showdiffs',true,'expts','all','details',false);
      args=processargs(defaults,varargin);

      sel=summary.exptsel(args.expts);
      
      fold=cat(2,summary.cleavage.fold);  % fold(i,k) - seq i, condition k
      foldci=cat(3,summary.cleavage.foldci);   % foldci(i,j,k) - seq i, bound j (1 or 2), condition k
      sel=find(sel);
      fold=fold(:,sel);
      foldci=foldci(:,:,sel);
      
      ind=[];
      for i=1:length(naseqs)
        ind(i)=find(naseqs(i)==summary.naseq);
      end
      fold=fold(ind,:);
      foldci=foldci(ind,:,:);
      lfold=log(fold);
      lfoldci=log(foldci);
      sigma1=squeeze(lfoldci(:,1,:))-lfold;
      sigma2=squeeze(lfoldci(:,2,:))-lfold;
      fprintf('Comparing %d seqs...', length(naseqs));
      sigmadiff1=zeros(length(naseqs),length(naseqs),size(lfold,2),'single');
      sigmadiff2=zeros(length(naseqs),length(naseqs),size(lfold,2),'single');
      diff=zeros(length(naseqs),length(naseqs),size(lfold,2),'single');

      for i=1:length(naseqs)
        if mod(i,100)==0
          fprintf('%d...',i);
        end
        for j=i+1:length(naseqs)
          diff(i,j,:)=squeeze(lfold(i,:)-lfold(j,:));
          sigmadiff1(i,j,:)=sqrt(1/2*(sigma1(i,:).^2+sigma2(j,:).^2));
          sigmadiff2(i,j,:)=sqrt(1/2*(sigma2(i,:).^2+sigma1(j,:).^2));
        end
      end
      fprintf('done\n');
      sigmadiff=sigmadiff1;
      % if diff<0 then fold(1)<fold(2), so should use sigma2(1) and sigma1(2), which is sigmadiff2
      sigmadiff(diff<0)=sigmadiff2(diff<0);
      % Fold around diagonal
      diff=diff+permute(diff,[2,1,3]);
      sigmadiff=sigmadiff+permute(sigmadiff,[2,1,3]);
      dprime=diff./sigmadiff;
      for i=1:size(dprime,1)
        dprime(i,i,:)=0;  % Force diagonals to zero (may be nan if no data)
      end
      
      if args.showdiffs
        for i=1:length(naseqs)
          for j=i+1:length(naseqs)
            first = true;
            for k=1:size(dprime,3)
              if abs(dprime(i,j,k))>args.thresh
                if first
                  fprintf('%d(%d) vs %d(%d)\n', naseqs(i),ind(i),naseqs(j),ind(j));
                  first=false;
                end
                fprintf('%-18.18s: %.1f [%.1f,%.1f] vs %.1f [%.1f, %.1f];  d''=%5.1f %s\n', summary.getdesc(sel(k)), fold(i,k),foldci(i,:,k),fold(j,k),foldci(j,:,k),dprime(i,j,k),repmat('*',1,round(abs(2*dprime(i,j,k)))));
              end
            end
          end
        end
      end
      dpmax=nanmax(abs(dprime),[],3);
      nmissing=sum(isnan(dpmax(:)));
      if nmissing>0
        fprintf('Have %d seq comparison that have no overlap in the experiments (over %d naseqs)\n', nmissing, sum(any(isnan(dpmax))));
        [si,sj]=ind2sub(size(dpmax),find(isnan(dpmax(:)),1));
        fprintf('e.g. %d vs %d\n', naseqs([si,sj]));
        dpmax(isnan(dpmax))=0;
      end
      assert(all(isfinite(dpmax(:))));
      if nargout>0
        x=struct('naseqs',naseqs,'ind',ind,'expts',find(sel),'dpmax',dpmax);
        if args.details
          x.dprime=dprime;   % Only save if requested since it can be very large
        end
      end
    end
    
    function targetmdscaling(summary,varargin)
      % Cluster the targets by pattern of fold change as a function of sequence
      defaults=struct('minfold',1.2,'mindata',10,'thresh',1.5);
      args=processargs(defaults,varargin);

      fold=[];
      targets={};
      cind=[];
      for i=1:length(summary.cleavage)
        c=summary.cleavage(i);
        issingle=~isempty(regexp(c.target,'^[0-9][0-9]*[A-H][0-9][0-9]*$'));
        if ~issingle
          continue;
        end
        for k=1:length(c.trpexpt)
          fold(:,end+1)=summary.cleavage(1).ratioavg./c.ratiosc(:,k);
          if length(c.trpexpt)==1
            targets{end+1}=c.target;
          else
            targets{end+1}=sprintf('%s.%d',c.target,k);
          end
          cind(end+1)=i;
        end
      end
      
      %fold(abs(fold)<args.thresh)=1.0;
      nseqs=sum(isfinite(fold),1);  % Number of sequences with data for condition i
      nconditions=sum(isfinite(fold),2);  % Number of conditions for seq i

      useconditions=nseqs>=args.mindata;   % Only conditions for which we have data for most sequences
      useconditions=useconditions & nanmax(fold,[],1)>=args.minfold;
      
      fold=fold(:,useconditions);
      targets=targets(useconditions);
      cind=cind(useconditions);
      fprintf('Using %d/%d conditions\n', sum(useconditions), length(useconditions));

      useseqs=nconditions>=args.mindata;  % Use only these sequences for analysis
      useseqs=useseqs & nanmax(fold,[],2)>=args.minfold;
      fprintf('Using %d/%d seqs\n', sum(useseqs), length(useseqs));
      
      % Set switching for unused seqs, conditions to 1
      fold=fold(useseqs,:);

      dist=zeros(length(targets));
      for i=1:length(targets)
        for j=i+1:length(targets)
          c=corrcoef(log(fold(:,[i,j])),'rows','pair');
          dist(i,j)=1-c(1,2);
        end
      end
      dist=dist+dist';
      %dist=pdist(log(fold'),'correlation');
      y=mdscale(dist,3,'criterion','metricstress');

      ti='MDScale Targets'
      setfig(ti);clf;
      plot3(y(:,1),y(:,2),y(:,3),'ro'); 
      axis off;
      hold on;
      text(y(:,1),y(:,2),y(:,3),targets);
      for i=2:length(cind)
        if cind(i)==cind(i-1)
          plot3(y([i-1,i],1),y([i-1,i],2),y([i-1,i],3),'r');
        end
      end
      
      % Check distance for nearby targets
      d3=zeros(size(dist));
      for i=1:size(dist,1)
        for j=i+1:size(dist,2)
          d3(i,j)=norm(y(i,:)-y(j,:));
          if dist(i,j)<0.06 || cind(i)==cind(j)
            fprintf('%6.6s(%2d) - %6.6s(%2d) - d=%4.2f, d3=%4.2f\n', targets{i}, i, targets{j}, j, dist(i,j), norm(y(i,:)-y(j,:)));
          end
        end
      end
      d3=d3+d3';
      d3=d3+eye(size(d3,1));
      setfig('dcompare');clf;
      plot(dist(:),d3(:),'.');
      xlabel('1-correlation');
      ylabel('d3');
      
      title(ti);
    end


    function checkconsistency(summary,varargin)
    % Check consistency of summaries
    % Print any switches with ratiosc's over different trpexpts that differ by a factor >= thresh
      defaults=struct('thresh',2,'nsref',false,'target',[],'mincnt',100);
      args=processargs(defaults,varargin);
      allruns= unique(horzcat(summary.cleavage.runs));
      allh=nan(length(allruns),1);
      colors=nan(length(allruns),3);
      if args.nsref
        seqsel=strncmp(summary.names,'NSRef',5);
      else
        seqsel=true(size(summary.names));
      end
      setfig('Consistency');clf;
      for i=1:length(summary.cleavage)
        s=summary.cleavage(i);
        ratiosc=s.ratiosc;
        ratiosc(squeeze(sum(s.cnt(:,2:3,:),2))<args.mincnt)=nan;
        if ~isempty(args.target) && ~strcmp(s.target,args.target)
          continue;
        end
        if length(s.trpexpt)==1
          ;%fprintf('Only 1 experiment for %s@%s\n', s.target,concfmt(s.conc));
        else
          c=ratiosc./(ratiosc+1);
          ratioavg=nan(size(ratiosc));
          for j=1:size(ratioavg,2)
            rtmp=ratiosc;
            rtmp(:,j)=nan;  % Average over all other measurements
            ratioavg(:,j)=exp(nanmean(log(rtmp),2));
          end

          cm=s.ratioavg./(s.ratioavg+1);
          h=loglog(ratioavg(seqsel,:),ratiosc(seqsel,:),'.');
          hold on;
          for j=1:length(h)
            sel=strcmp(s.runs{j},allruns);
            if isnan(colors(sel,1))
              allh(sel)=h(j);
              colors(sel,:)=get(h(j),'Color');
            else
              set(h(j),'HandleVisibility','off');
            end
            set(h(j),'Color',colors(sel,:));
          end
          ratio=[];
          for j=1:size(ratiosc,1)
            ratio(j)=max(max(ratiosc(j,:)./ratioavg(j,:),ratioavg(j,:)./ratiosc(j,:)));
          end
          outliers=find(ratio>args.thresh);
          if ~isempty(outliers)
            fprintf('Have %d outliers for %s@%g - max ratio=%.1f\n', length(outliers), s.target,s.conc,max(ratio));
            x=table();
            x.name=summary.names(outliers);
            x.naseq=int32(summary.naseq(outliers));
            for m=1:size(ratiosc,2)
              cname=sprintf('%s_%d',s.runs{m},s.trpexpt(m));
              cname=strrep(cname,';','_');
              cname=strrep(cname,'>','_');
              cname=strrep(cname,'-','_');
              x.(cname)=round(ratiosc(outliers,m),2);
            end
            fprintf(' %-10.10s %9.9s %5.5s','Name','naseq','max');
            for m=1:size(ratiosc,2)
              fprintf('%20.20s',sprintf('%s',s.runs{m}));
            end
            fprintf('%10.10s\n',' Mean');
            fprintf(' %-10.10s %9.9s %5.5s','','','ratio');
            for m=1:size(ratiosc,2)
              fprintf(' %19.19s',sprintf('%d',s.trpexpt(m)));
            end
            fprintf('%10.10s\n','');
            for kk=1:length(outliers)
              k=outliers(kk);
              fprintf(' %-10.10s %9d %5.2f ',summary.names{k},summary.naseq(k),ratio(k));
              for m=1:size(ratiosc,2)
                fprintf('%19.2f',c(k,m));
                if ratioavg(k,m)/ratiosc(k,m)>args.thresh
                  fprintf('<');
                elseif ratiosc(k,m)/ratioavg(k,m)>args.thresh
                  fprintf('>');
                else
                  fprintf(' ');
                end
              end
              fprintf('%9.2f\n',cm(k));
            end
          end
        end
      end
      cleaveticks(1,1);
      ax=axis;
      plot(ax(1:2),ax(1:2),'r:');
      plot(ax(1:2),ax(1:2)*args.thresh,'r:');
      plot(ax(1:2),ax(1:2)/args.thresh,'r:');
      legend(allh(isfinite(allh)),allruns(isfinite(allh)),'Location','best');
      xlabel('Mean of other measurements');
      ylabel('Specific values');
      title(sprintf('Consistency - mincnt=%d, thresh=%.1f',args.mincnt,args.thresh));
    end

    function res=significantswitches(summary,varargin)
    % Find switch/target combinations that have signficant switching
    % Return issig - 0 not a switch, 1 is a switch, nan - uncertain
      defaults=struct('minfold',2,'alpha',.05);
      args=processargs(defaults,varargin);
      issig=nan(length(summary.cleavage),length(summary.naseq));
      for i=1:length(summary.cleavage)
        c=summary.cleavage(i);
        for j=1:size(c.foldsamps,1)
          sel=isfinite(c.foldsamps(j,:));
          if sum(sel)>0
            issig(i,j)=mean(c.foldsamps(j,sel)>=args.minfold);
          end
        end
        %issig(i,isfinite(c.foldci(:,1)))=0.5;
        %sig=find(c.foldci(:,1)>=args.minfold);
        %issig(i,sig)=1;
        %issig(i,c.foldci(:,2)<args.minfold)=0;
        sig=find(issig(i,:)>(1-args.alpha));
        if ~isempty(sig)
          fprintf('%s@%-s:  ',c.target, concfmt(c.conc));
          for k=1:length(sig)
            fprintf('%s(%.1f) ',summary.names{sig(k)},c.fold(sig(k)));
          end
          fprintf('\n');
        end
      end
      
      issingle=cellfun(@(z) ~isempty(z), regexp({summary.cleavage.target},'^[0-9][0-9]*[A-H][0-9][0-9]*$'));
      targets={summary.cleavage(issingle).target};
      if isempty(targets)
        fprintf('No single-target data\n');
        return;
      end
      issig=issig(issingle,:);
      keepseq=nansum(issig>=1-args.alpha)>0;  % At least one target for which it switches
      issig=issig(:,keepseq);
      seqs=summary.names(keepseq);
      keeptarget=nansum(issig'>=1-args.alpha)>0; % At least one seq that switches for this target
      issig=issig(keeptarget,:);   
      targets=targets(keeptarget);
      
      % Determine equivalent targets
      dsig=nan(size(issig));
      dsigalpha=args.alpha;
      dsig(issig>(1-dsigalpha))=1;
      dsig(issig<dsigalpha)=-1;
      dsig(isnan(issig))=nan;
      tdist=zeros(length(targets));
      for i=1:length(targets)
        for j=i+1:length(targets)
          sel=isfinite(dsig(i,:)) & dsig(i,:)~=0 & isfinite(dsig(j,:)) & dsig(j,:)~=0;
          if all(dsig(i,sel)==dsig(j,sel)) && sum(sel)>0
            fprintf('%s and %s are equivalent over %d active and %d inactive sequences\n', targets{i}, targets{j}, sum(dsig(i,sel)==1), sum(dsig(i,sel)==-1));
          end
          tdist(i,j)=sum(dsig(i,sel)~=dsig(j,sel));
        end
      end
      tdist=tdist+tdist';
      tord=optimalleaforder(linkage(tdist),tdist);
      targets=targets(tord);
      issig=issig(tord,:);
      dsig=dsig(tord,:);
      tdist=tdist(tord,tord);
      
      % Determine equivalent seqs
      sdist=zeros(length(seqs));
      for i=1:length(seqs)
        for j=i+1:length(seqs)
          sel=isfinite(dsig(:,i)) & dsig(:,i)~=0 & isfinite(dsig(:,j)) & dsig(:,j)~=0;
          if all(dsig(sel,i)==dsig(sel,j)) && sum(sel)>0
            fprintf('%s and %s are equivalent over %d active and %d inactive targets\n', seqs{i}, seqs{j}, sum(dsig(sel,i)==1), sum(dsig(sel,i)==-1));
          end
          sdist(i,j)=sum(dsig(sel,i)~=dsig(sel,j));
        end
      end
      sdist=sdist+sdist';
      sord=optimalleaforder(linkage(sdist),sdist);
      seqs=seqs(sord);
      issig=issig(:,sord);
      sdist=sdist(sord,sord);
      
      res=struct('minfold',args.minfold,'seqs',{seqs},'targets',{targets},'issig',issig,'tdist',tdist,'sdist',sdist);
      setfig('sigswitches');clf;
      data=log10(issig./(1-issig));
      data(end+1,:)=nan;data(:,end+1)=nan;
      pcolor(data);
      caxis([-2,2]);
      shading flat;
      hold on;
      ax=axis;
      % Draw dividing line between distinct conditions
      groupstart=1;
      for i=2:size(sdist,1)
        if any(sdist(i,groupstart:i-1)~=0)
          plot(i*[1,1],ax(3:4),'k','LineWidth',3);
          groupstart=i;
        end
      end
      prev=1;
      for i=2:size(tdist,1)
        if tdist(i,prev)~=0
          plot(ax(1:2),i*[1,1],'k','LineWidth',3);
          prev=i;
        end
      end

      set(gca,'XTick',(1:length(seqs))+0.5);
      set(gca,'XTickLabel',seqs);
      set(gca,'XTickLabelRotation',90);
      set(gca,'YTick',(1:length(targets))+0.5);
      set(gca,'YTickLabel',targets);
      h=colorbar;
      set(get(h,'Label'),'string',sprintf('P[fold > %.1f]',args.minfold))
      p=[0.01,.05,0.1,0.2,0.5,.8,.9,.95,.99];
      set(h, 'Ticks',log10(p./(1-p)));
      set(h,'TickLabels',arrayfun(@(z) sprintf('%.2f',z), p,'UniformOutput',false)); 
      title(sprintf('Seq-Targets with >%.1f fold switching (alpha=%.3f)',args.minfold,args.alpha));
    end
    
    function plotrunvariation(summary,ind,varargin)
      defaults=struct('conc',2e-6);
      args=processargs(defaults,varargin);
      if ischar(ind)
        ind={ind};
      end
      if iscell(ind)
        ind=find(ismember({summary.cleavage.target},ind) & [summary.cleavage.conc]==args.conc);
      end
      ind=[1,setdiff(ind,1)];
      c=summary.cleavage(ind);
      setfig('runvariation');clf;
      [~,ord]=sort(c(1).ratioavg);
      t=tiledlayout(length(ind),1,'TileSpacing','compact');
      ax=[];
      for i=1:length(c)
        nexttile();
        boxplot(c(i).ratiosc(ord,:)');
        hold on;
        plot(c(1).ratioavg(ord),'-g');
        set(gca,'YScale','log');
        cleaveticks(0,1);
        set(gca,'XTickLabel',summary.names(ord));
        set(gca,'XTickLabelRotation',90);
        if isempty(c(i).target)
          ylabel('-target');
        else
          ylabel(c(i).target);
        end
        ax(i)=gca;
      end
      linkaxes(ax);
      xlabel(t,'Sensor');
      ylabel(t,'Cleavage');
      if length(ind)==1
        title(sprintf('Target: %s, Num runs: %d', c.target,size(c.ratiosc,2)));
      else
        title(t,'Run Variations');
        %legend({summary.cleavage(ind).target});
      end
    end

    function [p,foldci]=foldprob(summary,fold,varargin)
    % Compute matrix of probability that fold change for each target,aptamer combination is >= fold
    % Use row, column, plate, diag as bounds
    % p is a high bound on the prob that each sensor is a switch with >= fold change
      defaults=struct('minincprob',0.9,'ci',[5,95],'minhitprob',0.5,'basengs',[]);
      args=processargs(defaults,varargin);
      p=nan(8,10,12,length(summary.names));
      foldci=nan(8,10,12,length(summary.names),2);
      foldci(:,:,:,:,1)=1;   % Lower bound of CI is initial 1.0 fold
      for plate=1:12
        pname{plate}=sprintf('P%02d',(plate-1)*10+1);
        ind=find(strcmp(pname{plate},{summary.cleavage.target})&[summary.cleavage.conc]==2e-6);
        assert(length(ind)==1);
        pfold(plate,:,:)=summary.cleavage(ind).foldsamps;
      end
      for col=1:10
        cname{col}=sprintf('C%d',col+1);
        ind=find(strcmp(cname{col},{summary.cleavage.target})&[summary.cleavage.conc]==2e-6);
        assert(length(ind)==1);
        cfold(col,:,:)=summary.cleavage(ind).foldsamps;
      end
      for row=1:8
        rname{row}=sprintf('R%c',char(row-1+'A'));
        ind=find(strcmp(rname{row},{summary.cleavage.target})&[summary.cleavage.conc]==2e-6);
        assert(length(ind)==1);
        rfold(row,:,:)=summary.cleavage(ind).foldsamps;
      end
      for dprdiag=1:8
        dprname{dprdiag}=sprintf('D%d',dprdiag);
        ind=find(strcmp(dprname{dprdiag},{summary.cleavage.target})&[summary.cleavage.conc]==2e-6);
        assert(length(ind)==1);
        dprfold(dprdiag,:,:)=summary.cleavage(ind).foldsamps;
      end
      for dpcdiag=1:10
        dpcname{dpcdiag}=sprintf('DPC%d',dpcdiag);
        ind=find(strcmp(dpcname{dpcdiag},{summary.cleavage.target})&[summary.cleavage.conc]==2e-6);
        assert(length(ind)==1);
        dpcfold(dpcdiag,:,:)=summary.cleavage(ind).foldsamps;
      end
      for row=1:8
        for col=1:10
          for plate=1:12
            dprdiag=mod(row-plate,8)+1;
            dpcdiag=mod(col-plate,10)+1;
            f=cat(1,pfold(plate,:,:),cfold(col,:,:),rfold(row,:,:),dprfold(dprdiag,:,:),dpcfold(dpcdiag,:,:));
            minf=squeeze(nanmin(f));   % Find min of the 4 measures for each boot,each aptamer
            p(row,col,plate,:)=mean(minf>=fold,2);   % prob >= fold
            foldci(row,col,plate,:,2)=prctile(minf,args.ci(2),2);
            mn=mean(minf,2);   % To check for nans
            p(row,col,plate,isnan(mn))=nan;
            foldci(row,col,plate,isnan(mn),2)=nan;
          end
        end
      end
      fprintf('Possible hits (with p>=%.1f):\n',args.minhitprob);
      for apt=1:length(summary.names)
        missing=[pname(isnan(pfold(:,apt,1))),rname(isnan(rfold(:,apt,1))),cname(isnan(cfold(:,apt,1))),dprname(isnan(dprfold(:,apt,1))),dpcname(isnan(dpcfold(:,apt,1)))];
        fprintf('%s: (%s)',summary.names{apt},strjoin(missing,','));
        [r,j,v]=find(p(:,:,:,apt)>args.minhitprob);
        [c,pl]=ind2sub([10,12],j);
        if length(r)>30
          fprintf(' Have %d possible hits', length(r));
        else
          for i=1:length(r)
            wellname=sprintf('%02d%c%d',(pl(i)-1)*10+1,'A'+r(i)-1,c(i)+1);
            fprintf(' %s',wellname);
            wind=find(strcmp({summary.cleavage.target},wellname));
            if length(wind)==1
              fprintf('(%.2f)',summary.cleavage(wind).fold(apt));
            end
          end
        end
        fprintf('\n');
      end

      % List any inconsistencies between plate,row,col,dprdiag,dpcdiag measurements and single target
      fprintf('\nFlagging aptamers that are inconsistent with prob >= %.2f:\n', args.minincprob);
      aptlist=[];
      exptlist=find(arrayfun(@(z) z.conc==2e-6 && (length(z.target)<1 || z.target(1)=='R' || z.target(1)=='C' || z.target(1)=='D' || z.target(1)=='P'),summary.cleavage));
      for row=1:8
        for col=1:10
          for plate=1:12
            wellname=sprintf('%02d%c%d',(plate-1)*10+1,'A'+row-1,col+1);
            ind=find(strcmp({summary.cleavage.target},wellname) & [summary.cleavage.conc]==2e-6);
            if ~isempty(ind)
              c=summary.cleavage(ind);
              % Override foldci
              foldci(row,col,plate,:,:)=prctile(c.foldsamps,args.ci,2);   % Ignores mux data
              dprdiag=mod(row-plate,8)+1;
              dpcdiag=mod(col-plate,10)+1;
              p1=squeeze(p(row,col,plate,:));   % Prob(fold>2) < p1
              p2=mean(c.foldsamps>fold,2);   % Prob(fold>2) = p2
              pinc=(1-p1).*p2;   % Probability that mux says its not a switch and specific says it is
              if any(pinc>args.minincprob)
                % Inconsistent
                badapt=find(pinc>args.minincprob);
                for kk=1:length(badapt)
                  k=badapt(kk);
                  fprintf('%s(%d) %s(%d): p(inc)=%.2f:  mux says p(fold>%.1f) < %.2f (P%02d=%.2f,R%c=%.2f,C%d=%.2f,DPR%d=%.2f,DPC%d=%.2f), but target measurement gives p(fold>%.1f) = %.2f\n',...
                          summary.names{k},k, c.target, ind, pinc(k), fold, p1(k), ...
                          (plate-1)*10+1,mean(pfold(plate,k,:)>fold), row+'A'-1,mean(rfold(row,k,:)>fold), col+1, mean(cfold(col,k,:)>fold), dprdiag, mean(dprfold(dprdiag,k,:)>fold), dpcdiag,  mean(dpcfold(dpcdiag,k,:)>fold), fold, p2(k));
                end
                aptlist=[aptlist;badapt];
                exptlist=[exptlist,ind];
                fprintf('\n');
              else
                ; % fprintf('%s ok (max(pinc)=%.2f)\n', c.target, max(pinc));
              end
            end
          end
        end
      end
      if ~isempty(args.basengs)
        summary.plotgradsummary(args.basengs,'metric','ratio_of_ratio','minswitching',1.0,'expts',exptlist,'naseqs',summary.naseq(aptlist));
      end
      p1=foldci(:,:,:,:,1);
      p2=foldci(:,:,:,:,2);
      fprintf('Using a [%.0f,%.0f] CI, %.1f%% are < %.1f fold, %.1f%% are > %.1f fold, and %.1f%% are indeterminate\n', args.ci, 100*mean(p2(:)<fold), fold, 100*mean(p1(:)>fold), fold, (1-mean(p2(:)<fold)-mean(p1(:)>fold))*100);
      p1=any(foldci(:,:,:,:,1)>fold,4);
      p2=all(foldci(:,:,:,:,2)<fold,4);
      fprintf('Using a [%.0f,%.0f] CI, %d targets are < %.1f fold, %d are > %.1f fold, and %d are indeterminate\n', args.ci, sum(p2(:)), fold, sum(p1(:)), fold, 960-sum(p2(:))-sum(p1(:)));
    end
    
    function plotcoverage(summary)
    % Show coverage of each target in terms of number of aptamers with 0,1,2,>3 measurements
      cov=[];
      for i=1:length(summary.cleavage)
        cov(i,:)=hist(sum(isfinite(summary.cleavage(i).foldbyrun),2),0:3);
      end
      setfig('lowcoverage');
      %bar(sum(cov(:,3:end),2));
      bar(cov,'stacked');
      set(gca,'XTick',1:length(summary.cleavage));
      ticks={summary.cleavage.target};
      for i=1:length(ticks)
        if summary.cleavage(i).conc~=2e-6
          ticks{i}=sprintf('%s@%s',ticks{i},concfmt(summary.cleavage(i).conc));
        end
      end
      set(gca,'XTickLabel',ticks);
      set(gca,'XTickLabelRotation',90);
      ylabel('Num Aptamers');
      legend('0','1','2','>2','location','eastoutside');
    end
    
    function dumphits(summary,varargin)
    % Create a file for running ML tests
    % Output each of the tested targets with one column per aptamer
    % Use only aptamers with at least minfold switching for some target with given certainty
    % If compounds is given, then use to find molecules which are not in massspec and thus should have nan for whether they are hits
      defaults=struct('minfold',1.3,'certainty',0.9,'minhitsperaptamer',2,'filename','hits.csv','compounds',[]);
      args=processargs(defaults,varargin);

      [p,foldci]=summary.foldprob(args.minfold,'ci',[1-args.certainty,args.certainty]*100);
      ishit=nan(size(p));
      ishit(foldci(:,:,:,:,1)>=args.minfold)=1;
      ishit(foldci(:,:,:,:,2)<=args.minfold)=0;
      fprintf('Of %d conditions, have %d hits, %d misses, %d indeterminate\n', length(ishit(:)),sum(ishit(:)==1),sum(ishit(:)==0),sum(isnan(ishit(:))));
      activeaptamers=nansum(reshape(ishit,960,length(summary.naseq)))>=args.minhitsperaptamer;
      fprintf('%d/%d aptamers have at least %d definite hits\n', sum(activeaptamers), length(activeaptamers),args.minhitsperaptamer);

      if ~isempty(args.compounds)
        verified=~isnan(args.compounds.meantime);
        fprintf('Mass spec verified presence of %d/%d compounds\n',sum(verified),length(verified));
      end
      
      fd=fopen(args.filename,'w');
      fprintf(fd,'Target');
      fprintf(fd,',%s',summary.names{activeaptamers});
      fprintf(fd,'\n');

      nancnt=0;
      for plate=1:size(ishit,3)
        for row=1:size(ishit,1)
          for col=1:size(ishit,2)
            tname=sprintf('%d%c%02d',(plate-1)*10+1,row+'A'-1,col+1);
            if ~isempty(args.compounds)
              cind=strcmp(args.compounds.names,tname);
              if sum(cind)~=1
                error('Compounds does not include entry for %s', tname);
              end
              verified=~isnan(args.compounds.meantime(cind));
              if ~verified
                if any(ishit(row,col,plate,activeaptamers))
                  fprintf('Have %d hit(s) for %s even though mass spec has not verified presence -- assuming it really was present\n', nansum(ishit(row,col,plate,activeaptamers)), tname);
                else
                  ishit(row,col,plate,:)=nan;   % Unknown state for these
                end
              end
            end
            fprintf(fd,'%s',tname);
            fprintf(fd,',%.0f',ishit(row,col,plate,activeaptamers));
            fprintf(fd,'\n');
            nhits=nansum(ishit(row,col,plate,activeaptamers));
            if nhits>1
              fprintf('%s: %d hits\n', tname, nhits);
            end
          end
        end
      end
      fclose(fd);

      return;
      

      targets=false(length(summary.cleavage),1);
      activeaptamers=false(length(summary.naseq),1);
      for i=1:length(summary.cleavage)
        c=summary.cleavage(i);
        split=find(c.target>='A');
        if length(split)==1 && split>1
          targets(i)=true;
          activeaptamers=activeaptamers | c.foldci(:,1)>=args.minfold;
        end
      end
      fprintf('Have %d targets and %d aptamers with at least %.1f fold switching\n', sum(targets), sum(activeaptamers),args.minfold);

      tnames={summary.cleavage(targets).target};
      % Remove leading 0 from target names
      for i=1:length(tnames)
        if tnames{i}(1)=='0'
          tnames{i}=tnames{i}(2:end);
        end
      end

      fold=nan(sum(targets),sum(activeaptamers));
      ishit=nan(sum(targets),sum(activeaptamers));
      csel=summary.cleavage(targets);
      for i=1:sum(targets)
        c=csel(i);
        fold(i,:)=c.fold(activeaptamers);
        phigh=mean(c.foldsamps(activeaptamers,:)>args.minfold,2);
        ishit(i,phigh>args.certainty)=1;
        ishit(i,phigh<1-args.certainty)=0;
      end

      filenames={'fold.csv','hits.csv'};
      for of=1:2
        fd=fopen(filenames{of},'w');
        fprintf(fd,'Target');
        fprintf(fd,',%s',summary.names{activeaptamers});
        fprintf(fd,'\n');

        nancnt=0;
        for i=1:size(fold,1)
          fprintf(fd,'%s',tnames{i});
          if of==1
            fprintf(fd,',%.2f',fold(i,:));
          else
            fprintf(fd,',%.0f',ishit(i,:));
          end
          fprintf(fd,'\n');
        end
        nancnt=sum(isnan(fold(:)));
        ntotal=length(fold(:));
        fprintf('Have %d NaNs over %d measurements (%.2f%%)\n', nancnt, ntotal,nancnt/ntotal*100);
        fclose(fd);
      end
    end


    function testEqualCleavage(summary,varargin)
    % For each experiment, run pairwise test between each sequence for significant difference in cleavage
    % Store p-value for test in summary.cleavage(i).psameclv(m,n) 
    % If p-value< alpha (.05, for example), then reject H0: cleavages are equal
      defaults=struct('mincnt',100,'force',false);
      args=processargs(defaults,varargin);

      if isfield(summary.cleavage,'psameclv') && ~isempty(summary.cleavage(1).psameclv)
        if nargin<2 || ~args.force
          return;   % Already computed
        end
      end
      fprintf('Computing alpha for hypothesis of equal cleavage across %d conditions...',length(summary.cleavage));
      for i=1:length(summary.cleavage)
        fprintf('%d...',i);
        c=summary.cleavage(i);
        psameclv={};
        for j=1:length(c.runs)
          if length(c.runs)>1
            fprintf('(%c)','a'+j-1);
          end
          cnt=c.cnt(:,:,j);
          psameclv{j}=zeros(length(summary.naseq),length(summary.naseq));
          for m=1:length(summary.naseq)
            or=OddsRatio(cnt(m,2),cnt(m,3),cnt(m+1:end,2),cnt(m+1:end,3));
            psameclv{j}(m,m+1:end)=or.getAlphaNorm();
          end
          psameclv{j}=eye(length(summary.naseq))+psameclv{j}+psameclv{j}';   % Make symmetric
          % Use mincnt to filter out and set psameclv to nan when <mincnt on either side of comparison
          lowcnt=sum(cnt(:,2:3),2)<args.mincnt;
          psameclv{j}(lowcnt,:)=nan;
          psameclv{j}(:,lowcnt)=nan;
        end
        summary.cleavage(i).psameclv=psameclv;
      end
      fprintf('done\n');
    end
    
    function [ndiffer,ntest,diffbyexpt]=testAllEqualCleavage(summary,varargin)
    % Test for significantly different cleavage in at least one experiment using summary.cleavage(*).psameclv
    % Also, return number of tests used (differ will be nan for conditions with <mintest tests)
    % If multiple runs exist in any experiment, only use the one with the most total reads for these sequences
    % Alpha will be adjusted for FWER (Holms-Bonferroni correction)
      defaults=struct('naseq',[],'expt',[],'alpha',0.05,'mintest',[]);
      args=processargs(defaults,varargin);
      if isempty(args.expt)
        args.expt=1:length(summary.cleavage);
      end
      if isempty(args.naseq)
        args.naseq=summary.naseq;
      end
      if isempty(args.mintest)
        args.mintest=ceil(length(args.expt)/2);  
      end
      summary.testEqualCleavage();
      naseqsel=arrayfun(@(z) find(summary.naseq==z),args.naseq);
      if length(naseqsel) ~= length(args.naseq)
        error('Only found %d/%d of the requested naseqs', length(naseqsel),length(args.naseq));
      end
      ndiffer=zeros(length(naseqsel));
      % Count the number of tests with valid data for each condition
      ntest=zeros(length(naseqsel));
      rsel=[];
      for i=args.expt
        c=summary.cleavage(i);
        if length(c.runs)==1
          rsel(i)=1;
        else
          [~,rsel(i)]=max(sum(sum(c.cnt(naseqsel,2:3,:))));
        end
        ntest=ntest+isfinite(c.psameclv{rsel(i)}(naseqsel,naseqsel));
      end
      p1=cat(2,summary.cleavage(args.expt).psameclv);
      p2=cat(3,p1{:});   % Now have all of the alphas as an nseq*nseq*nexpt matrix
      p2=p2(naseqsel,naseqsel,:);
      [p2,ord]=sort(p2,3);   % Sort by increasing alpha
      ntest=sum(isfinite(p2),3);
      nremain=ntest;
      reject=false(size(p2));

      for i=1:size(p2,3)
        reject(:,:,i)=p2(:,:,i)<args.alpha./nremain;
        if i>1
          reject(:,:,i)=reject(:,:,i)|reject(:,:,i-1);   % Once we hit a failure, all others are failures
        end
        nremain=nremain-isfinite(p2(:,:,i));
      end
      assert(all(nremain(:)==0));
      ndiffer=sum(reject,3);
      
      diffbyexpt=false(length(naseqsel),length(naseqsel),length(args.expt));
      for i=args.expt
        c=summary.cleavage(i);
        diffbyexpt(:,:,i)=(c.psameclv{rsel(i)}(naseqsel,naseqsel)<args.alpha./ntest); % Make sure nans don't affect; apply Bonferonni correction
      end
      ndiffer(ntest<args.mintest)=nan;
      % Clear diagonal
      ndiffer(1:size(ntest,1)+1:end)=0;
    end
    

    function c=copy(summary)
      c=TRPSummary(summary.run,summary.naseq);
      fn=fieldnames(summary);
      for i=1:length(fn)
        c.(fn{i})=summary.(fn{i});
      end
    end
    
    function minfold=getminimumfold(summary)
    % Find minimum fold change for each sequence (i.e. max(foldci(:,1)) )
      minfold=summary.cleavage(1).foldci(:,1);
      for i=2:length(summary.cleavage)
        minfold=max(minfold,summary.cleavage(i).foldci(:,1));
      end
    end
    
    function [refs,ratio,fold]=chooserefs(summary,varargin)
      defaults=struct('num',15,'minfold',1.2,'mincoverage',0.8);
      args=processargs(defaults,varargin);

      cnt=arrayfun(@(z) sum(sum(z.cnt(:,2:3,:),2),3),summary.cleavage,'Unif',false);
      cnt=horzcat(cnt{:});
      coverage=mean(cnt>0,2);
      minfold=summary.getminimumfold();
      sel=find(coverage>=args.mincoverage & minfold<=args.minfold);
      ratio=summary.cleavage(1).ratioavg(sel);
      [ratio,ord]=sort(ratio);
      sel=sel(ord);
      fprintf('Initially have %d potential refs with fold<=%.1f\n', length(sel), args.minfold);
      
      while length(sel)>args.num
        % Eliminate one with smallest ratio gap between its neighbors (but not end ones)
        rratio=ratio(3:end)./ratio(1:end-2);
        [~,drop]=min(rratio);
        drop=drop+1;
        fprintf('Drop [%.2f,%.2f,%.2f] with rratio=%.2f\n', ratio(drop-2:drop), rratio(drop-1));
        sel=sel([1:drop-1,drop+1:end]);
        ratio=summary.cleavage(1).ratioavg(sel);
      end
      fold=minfold(sel);
      refs=summary.naseq(sel);
    end
        
    function removenonswitching(summary,varargin)
    % Remove any sequences for which all the foldci(:,1)<minfold
      defaults=struct('minfold',1.2,'refs',[]);
      args=processargs(defaults,varargin);
      minfold=summary.getminimumfold();
      keep=minfold>=args.minfold;
      keep=keep | ismember(summary.naseq,args.refs);
      fprintf('Removing %d/%d sequences that have fold<%.1f\n', sum(~keep), length(keep), args.minfold);
      if all(keep)
        return;
      end
      for i=1:length(summary.cleavage)
        c=summary.cleavage(i);
        c.rawratio=c.rawratio(keep,:);
        c.rawratioci90=c.rawratioci90(keep,:,:);
        c.cnt=c.cnt(keep,:,:);
        c.ratiosc=c.ratiosc(keep,:);
        c.ratioavg=c.ratioavg(keep);
        c.foldsampsbyrun=c.foldsampsbyrun(keep,:,:);
        c.foldbyrun=c.foldbyrun(keep,:);
        c.foldcibyrun=c.foldcibyrun(keep,:,:);
        c.foldsamps=c.foldsamps(keep,:);
        c.fold=c.fold(keep);
        c.foldci=c.foldci(keep,:);
        summary.cleavage(i)=c;
      end
      summary.naseq=summary.naseq(keep);
      summary.names=summary.names(keep);
      summary.refratio=summary.refratio(keep);
      if isfield(summary,'mergenaseqs')
        summary.mergenaseqs=summary.mergenaseqs(keep);
      end
    end
    
    function merged=mergeSimilarSeqs(summary,ngs,varargin)
    % Group together sequences with the same cleavage pattern (ndiffer <= cutoff)
      defaults=struct('usecleavage',false,'cutoff',[],'alpha',.05,'dump',false,'refs',[],'seqclusters',true);
      args=processargs(defaults,varargin);

      if args.usecleavage
        if isempty(args.cutoff)
          args.cutoff=0.99;
        end
        [ndiffer,ntest]=summary.testAllEqualCleavage('alpha',args.alpha);
        ndiffer(isnan(ndiffer))=max(ndiffer(:));
        z=linkage(squareform(ndiffer),'average');
      else
        z=summary.tree.linkages;
        if isempty(args.cutoff)
          args.cutoff=2;
        end
      end
      cl=zeros(size(summary.naseq));
      cl(summary.tree.useseqs)=cluster(z,'cutoff',args.cutoff,'criterion','distance');
      % Put each ref into its own cluster so we keep them all even if they weren't showing any fold change
      fprintf('Putting %d refs into their own clusters\n', length(args.refs));
      for i=1:length(args.refs)
        ind=summary.naseq==args.refs(i);
        assert(sum(ind)==1);
        cl(ind)=max(cl)+1;
      end
      
      ucl=unique(cl(cl>0));
      fprintf('Formed %d clusters from %d seqs with cutoff=%f\n', length(ucl), sum(cl>0), args.cutoff);
      if args.seqclusters
        % Split these clusters based on sequence clustering as well
        for i=1:length(ucl)
          sel=cl==ucl(i);
          seqclusts=ngs.seqs.getclusterindex(summary.naseq(sel));
          usc=unique(seqclusts);
          if length(usc)>1
            cl(sel)=cl(sel)+seqclusts/max(seqclusts);
          end
        end
        ucl=unique(cl(cl>0));
        fprintf('Split into %d clusters based on sequence dissimilarity\n', length(ucl));
      end

      % Construct a new summary structure with the sequences grouped
      % Keep only one sequence from each group and note the others in merged.mergenaseqs
      merged=summary.copy();
      merged.tree=[];
      
      % Clustering
      merged.mergenaseqs={};
      sel=false(size(summary.naseq));
      totalreads=zeros(size(summary.naseq));
      for i=1:length(summary.cleavage)
        totalreads=totalreads+sum(sum(summary.cleavage(i).cnt(:,2:3,:),3),2);
      end
      fsel=[];
      for i=1:length(ucl)
        if isempty(summary.mergenaseqs)
          merged.mergenaseqs{i}=summary.naseq(cl==ucl(i));
        else
          naseqs=summary.naseq(cl==ucl(i));
          msel=find(cellfun(@(z) ~isempty(intersect(z,naseqs)),summary.mergenaseqs));
          merged.mergenaseqs{i}=vertcat(summary.mergenaseqs{msel});
          assert(all(ismember(naseqs,merged.mergenaseqs{i})));
        end
        
        rds=max(totalreads(cl==ucl(i)));
        % Keep the one in the group with the most reads
        ind=find(cl==ucl(i) & totalreads==rds,1);
        assert(length(ind)==1);
        assert(~ismember(ind,fsel));
        fsel(i)=ind;
      end

      assert(length(fsel)==length(ucl));
      % filter out only sel sequences
      merged.naseq=summary.naseq(fsel);
      merged.names=summary.names(fsel);
      merged.refratio=summary.refratio(fsel);
      
      for i=1:length(merged.cleavage)
        c=merged.cleavage(i);
        s=summary.cleavage(i);
        c.ratioboot=cell(size(s.ratioboot));
        c.cnt=nan(length(ucl),size(s.cnt,2),size(s.cnt,3));
        c.rawratio=nan(length(ucl),size(s.rawratio,2));
        c.rawratioci90=nan(length(ucl),2,size(s.rawratioci90,3));
        for j=1:length(ucl)
          sel=find(cl==ucl(j));
          for k=1:length(c.trpexpt)
            c.ratioboot{k}(j,:)=exp(nanmean(log(s.ratioboot{k}(sel,:)),1));
            c.cnt(j,:,k)=sum(s.cnt(sel,:,k),1);
            c.rawratio(j,k)=exp(nanmean(log(s.rawratio(sel,k)),1));
            c.rawratioci90(j,:,k)=prctile(c.ratioboot{k}(j,:),[5,95]);
          end
        end
        % Need to reset scaling, etc.
        c.scaling=[];
        c.ratiosc=[];
        c.ratioavg=[];
        c.foldsampsbyrun=single([]);
        c.foldbyrun=single([]);
        c.foldcibyrun=single([]);
        c.foldsamps=single([]);
        c.fold=single([]);
        c.foldci=single([]);
        merged.cleavage(i)=c;
      end

      % Remove invalid psameclv
      if isfield(merged.cleavage,'psameclv')
        merged.cleavage=rmfield(merged.cleavage,'psameclv');
      end

      fprintf('Need to rerun setscaling() and computefolds()\n');
      if args.dump
        % Dump cluster information
        for i=1:length(merged.mergenaseqs)
          naseq=merged.mergenaseqs{i};
          if length(naseq)<2
            continue;
          end
          fprintf('Cl%d: N=%d\n',ucl(i),length(naseq));
          d={}; seqcl=[];
          for j=1:length(naseq)
            [d{j},tmp]=ngs.seqs.getcluster(naseq(j));
            seqcl(j)=tmp.cluster;
          end
          useqcl=unique(seqcl);
          for j=1:length(useqcl)
            fprintf('\t[%s]\n',strjoin(d(seqcl==useqcl(j)),','));
          end
        end
      end

    end
    
    
    function runCleavageTests(summary)
    % Run the suite of cleavage tests (i.e. not fold change, just cleavage)
      [ndiffer,ntest,diffbyexpt]=summary.testAllEqualCleavage();
      nexpt=length(summary.cleavage);
      totaldiff=[];
      for n=0:nexpt
        totaldiff(n+1)=sum(ndiffer(:)==n);
      end
      for n=[0:2,nexpt-2:nexpt]
        fprintf('ndiff=%2d: %.2f%%\n', n, totaldiff(n+1)/sum(totaldiff)*100);
      end
      setfig('ndiff');clf;
      bar(0:nexpt,totaldiff/sum(totaldiff)*100);
      set(gca,'YScale','log');
      ax=axis;ax(3)=ax(3)/1.1; ax(4)=ax(4)*1.1; axis(ax);  % R2020b seems to set axis bounds to exact min/max
      xlabel('Number of experiments that show significant differences in cleavage');
      ylabel('Percentage of sequences');
      
      setfig('Outlier Experiments');clf;
      sel=ntest==max(ntest(:)) & ndiffer==1;   % Only differ in one experiment; maybe that's an outlier?
      for i=1:size(diffbyexpt,3)
        dtmp=diffbyexpt(:,:,i);
        n(i)=sum(dtmp(sel));
      end
      barh(n);
      ylabel('Experiment');
      xlabel('Number of sequence pairs that have different cleavage only in this expt');
      enames={};
      for i=1:length(summary.cleavage)
        c=summary.cleavage(i);
        if isempty(c.target)
          enames{i}='DMSO';
        elseif c.conc==2e-6 || c.conc==0
          enames{i}=c.target;
        else
          enames{i}=sprintf('%s@%s',c.target,concfmt(c.conc));
        end
        enames{i}=sprintf('%s (%d)',enames{i},i);
      end
      
      set(gca,'YTick',1:length(enames));
      set(gca,'YTickLabel',enames);
    end
    
    function plotcleavagebyexpt(summary,naseq,varargin)
      defaults=struct('expts','all','conc',[],'minfold',[],'sort',true);
      args=processargs(defaults,varargin);

      sel=summary.exptsel(args.expts,args.conc);

      if ~any(sel)
        error('No experiment conditions selected');
      end
      sel(1)=true;  % Always include -target
      ind=find(summary.naseq==naseq);
      assert(length(ind)==1);
      ti=sprintf('Cleavage by expt for %d',summary.naseq(ind));
      if ~strcmp(args.expts,'all')
        ti=[ti,', ',args.expts];
      end
      if ~isempty(args.conc)
        ti=sprintf('%s,%s',ti,concfmt(args.conc));
      end
      if ~isempty(args.minfold)
        ti=sprintf('%s, fold>=%.1f',ti,args.minfold);
      end
      setfig(ti);clf;
      if false
        rlow=[];rhigh=[];c=[];
        for i=1:length(summary.cleavage)
          cnt=nansum(summary.cleavage(i).cnt(ind,2:3,:),3);
          [phat,pci]=binofit(cnt(1),sum(cnt));
          rlow(i)=pci(1);
          rhigh(i)=pci(2);
          c(i)=phat;
        end
      else
        rlow=arrayfun(@(z) nanmean(z.scaling.*squeeze(z.rawratioci90(ind,1,:))'), summary.cleavage(sel));
        rhigh=arrayfun(@(z) nanmean(z.scaling.*squeeze(z.rawratioci90(ind,2,:))'), summary.cleavage(sel));
        ratio=arrayfun(@(z) z.ratioavg(ind), summary.cleavage(sel));
      end
      other=[];
      
      targets=arrayfun(@(z) summary.getdesc(z,isempty(args.conc)), find(sel),'Unif',false);
      if ~isempty(args.minfold)
        keep=rlow<=ratio(1)/(args.minfold^2);
        fprintf('Showing %d/%d conditions that may have fold>=%.1f\n', sum(keep),length(keep), args.minfold);
        keep(1)=true;
        rlow=rlow(keep);
        rhigh=rhigh(keep);
        other=ratio(~keep);
        ratio=ratio(keep);
        targets=targets(keep);
      end
      if args.sort
        [~,ord]=sort(ratio(2:end),'desc');
        ord=[1,ord+1];
        ratio=ratio(ord);
        rhigh=rhigh(ord);
        rlow=rlow(ord);
        targets=targets(ord);
      end
      errorbar(ratio,1:length(ratio),0*ratio,0*ratio,ratio-rlow,rhigh-ratio,'o');
      hold on;
      if length(other)>0
        plot(other,other*0+length(ratio)+1,'.r');   % Other points
        targets{end+1}=sprintf('Other(%d)',length(other));
      end
      cleaveticks(1,0);
      ax=axis;
      ax(3)=0.5;ax(4)=length(ratio)+1.5;
      axis(ax);
      plot(ratio(1)*[1,1],ax(3:4),'-r');
      for fold=2:10
        plot(ratio(1)/(fold^2)*[1,1],ax(3:4),'-m');
      end
      set(gca,'XScale','log');
      xlabel('Cleavage');
      set(gca,'YTick',1:length(targets));
      set(gca,'YTickLabel',targets);
      set(get(gca,'YAxis'),'FontSize',min(10,800/length(ratio)));
      set(gca,'YGrid','on');
      pos=get(gcf,'Position');
      minheight=min(800,length(ratio)*11+50);
      pos(4)=minheight;
      set(gcf,'Position',pos);
      title(ti);
    end

    function comparecleavage(summary,naseqs)
    % Plot a comparison of cleavage of naseqs(1) vs naseqs(2:end) over all the experiments
      ti=sprintf('comparecleavage %s',sprintf('%d ',naseqs));
      setfig(ti);clf;
      assert(length(naseqs)>=2);
      sel=[];
      for j=1:length(naseqs)
        ind=find(naseqs(j)==summary.naseq);
        if isempty(ind)
          error('naseq %d not found in summary',naseqs(j));
        end
        sel(end+1)=ind;
      end
      
      allphat=[]; allpci=[];totalcnt=[];allnphat=[];allnpci=[];
      for i=1:length(summary.cleavage)
        c=summary.cleavage(i);
        for j=1:length(c.runs)
          cnt=c.cnt(sel,2:3,j);
          [phat,pci]=binofit(cnt(:,1),sum(cnt,2));
          [nphat,npci]=binofit(sum(cnt(:,1))-cnt(:,1),sum(sum(cnt,2))-sum(cnt,2));
          allphat(1:length(naseqs),end+1)=phat;
          allpci(1:length(naseqs),1:2,size(allphat,2))=pci;
          allnphat(1:length(naseqs),end+1)=nphat;
          allnpci(1:length(naseqs),1:2,size(allphat,2))=npci;
          totalcnt(1:length(naseqs),end+1)=sum(cnt,2);
        end
      end

      allphat=allphat./(1-allphat);
      allpci=allpci./(1-allpci);
      allnphat=allnphat./(1-allnphat);
      allnpci=allnpci./(1-allnpci);
      
      leg={};h=[];
      ypos=squeeze(allpci(:,2,:))-allphat;
      yneg=allphat-squeeze(allpci(:,1,:));
      xpos=squeeze(allnpci(:,2,:))-allnphat;
      xneg=allnphat-squeeze(allnpci(:,1,:));

      for i=1:length(naseqs)
        h(end+1)=errorbar(allnphat(i,:),allphat(i,:),yneg(i,:),ypos(i,:),xneg(i,:),xpos(i,:),'o');
        hold on;
        leg{end+1}=sprintf('%s med N=%.0f',summary.names{sel(i)},median(totalcnt(i,:)));
      end
      n0=length(summary.cleavage(1).runs);
      plot(allnphat(:,1:n0),allphat(:,1:n0),'xk','markersize',20);   % Mark -target
      xlabel('Combined less specific Seq');
      ylabel('Specific Seq');
      set(gca,'XScale','log');
      set(gca,'YScale','log');
      axis tight;
      lows=allpci(:,1,:); lows=lows(isfinite(lows(:)));
      highs=allpci(:,2,:); highs=highs(isfinite(highs(:)));
      ax=axis;
      plot(ax(1:2),ax(1:2),'r:');
      legend(h,leg,'location','best');
      cleaveticks(1,1);
      axis(ax);
    end

    function cnt=allcnts(summary)
    % Get counts across all conditions
    % Return c(i,3,j) - naseq i, experiment j, 3 count values
    % Sum cnts across for conditions with replicates
      cnt=nan(length(summary.naseq),3,length(summary.cleavage));
      for i=1:length(summary.cleavage)
        cnt(:,:,i)=sum(summary.cleavage(i).cnt,3);
      end
    end
    
    function setcontains(summary,vecs)
    % Setup contains based on vecs
      for i=1:length(summary.cleavage)
        if isempty(summary.cleavage(i).mixture)
          summary.cleavage(i).contains=[];
        else
          mix=Mixtures.instance().get(summary.cleavage(i).mixture);
          assert(length(mix)<=1);
          if length(mix)==1
            summary.cleavage(i).contains=mix.contents;
          end
        end
      end
    end

    function classifyhits(obj,varargin)
    % What targets are distinguishable by the pattern of aptamer activation?
      defaults=struct('minfold',2,'alpha',0.05,'conc',10e-6);
      args=processargs(defaults,varargin);


      % Only single-target tests
      minmaxfold=arrayfun(@(z) max(z.foldci(:,1)),obj.cleavage);
      single=arrayfun(@(z) length(z.contains)==1,obj.cleavage);
      if isempty(args.conc)
        concmatch=true(size(obj.cleavage));
      else
        concmatch=arrayfun(@(z) z.conc==args.conc,obj.cleavage);
      end
      sel=single&concmatch&minmaxfold>=args.minfold;
      orig=find(sel);
      fprintf('Have %d single target experiments, %d of which have at least %.1f fold for some seq\n',sum(single),sum(sel), args.minfold);

      % Retrieve the bootstrap fold samples for each and the compound pks
      fold=cat(2,obj.cleavage(sel).fold);
      compounds=[obj.cleavage(sel).contains];
      conc=[obj.cleavage(sel).conc];
      % For each pair (i,j), check if there are any aptamers such that P(fold(i)>fold(j))>1-args.alpha
      dist=[];
      for i=1:size(fold,2)
        fi=fold(:,i);
        for j=1:size(fold,2)
          fj=fold(:,j);
          sel=isfinite(fi) & isfinite(fj) & (fi>=args.minfold | fj>=args.minfold);
          sc=sum(fj(sel))/sum(fi(sel));
          dvec=log(fj/sqrt(sc))-log(fi*sqrt(sc));
          dist(i,j)=sqrt(mean(dvec(sel).^2));
          if compounds(i)==compounds(j) && i<j
            fprintf('%s vs %s: %.2f\n', obj.getdesc(orig(i)), obj.getdesc(orig(j)), dist(i,j));
          end
        end
      end

      % Dendrogram based on average distance
      link=linkage(squareform(double(dist)),'average');
      setfig('classifyhits-dendrogram');clf;
      dendrogram(link,size(dist,1),'Labels',{Compounds.instance().get(compounds).name});
      set(gca,'XTickLabelRotation',90)

      % Cluster
      clusters=cluster(link,0.0001);
      pklist={};titles={};
      for i=1:max(clusters)
        sel=clusters==i;
        if sum(sel)>1
          %fprintf('%d: %s\n', i, strjoin({Compounds.instance().get(compounds(sel)).name},','));
          pklist{end+1}=compounds(sel);
          titles{end+1}=sprintf('Cluster %d',i);
        end
      end
      fprintf('Have %d groups of targets, %d of which have a single target\n', length(unique(clusters)),length(unique(clusters))-length(pklist));
      % Plot them
      Compounds.instance().savehtml('tgtclusters',pklist,'title',titles);

      % Verify some
      setfig('classifyhits');clf;
      tiledlayout('flow');
      for i=1:length(pklist)
        if length(pklist{i})>=2
          if isempty(args.conc)
            sel=find(arrayfun(@(z) length(z.contains)==1 && ismember(z.contains,pklist{i}),obj.cleavage));
          else
            sel=find(arrayfun(@(z) length(z.contains)==1 && ismember(z.contains,pklist{i}) && z.conc==args.conc,obj.cleavage));
          end
          fprintf('%s: %s\n%s: %s\n', titles{i}, strjoin({obj.cleavage(sel).target},','),blanks(length(titles{i})),strjoin({Compounds.instance().get([obj.cleavage(sel).contains]).name},','));
          nexttile
          obj.comparefold(sel(1),sel(2),'alpha',args.alpha,'newfig',false,'list',false,'minfold',args.minfold);
          title(titles{i});
          axis([0.5,10,0.5,10]);
        end
      end
    end
    
    function plothits(obj,ngs,varargin)
    % Plot the hit compounds
    % Note that there is also TargetID.plothits() with similar output
    % TargetID version plots hits inferred by either singles or vectors.  
    %   Ones with fold>=thresh are shown
    % TRPSummary version is based only on singles.  It includes compounds if the high 
    %   end of their 90% confidence interval are >thresh
      defaults=struct('thresh',2,'maxshow',12,'dirname','./summary','conc',10e-6,'naseq',[]);
      % debug - list of targets to emit debug messages
      args=processargs(defaults,varargin);

      captions={}; ti={}; pk={}; im={};
      exptsel=find(arrayfun(@(z) z.conc==args.conc && length(z.contains)==1,obj.cleavage));
      allcompounds=arrayfun(@(z) z.contains, obj.cleavage(exptsel));
      allfold=horzcat(obj.cleavage(exptsel).fold);
      allfoldci=cat(3,obj.cleavage(exptsel).foldci);
      maxfold=nanmax(allfold,[],2);
      if isempty(args.naseq)
        seqsel=find(maxfold>=args.thresh);
        fprintf('Have %d sequences with at least %.1f fold change over %d experiments\n', length(seqsel),args.thresh,length(exptsel));
      else
        seqsel=find(ismember(obj.naseq,args.naseq));
        fprintf('Have %d sequences over %d experiments\n', length(seqsel),length(exptsel));
      end
      allfold=allfold(seqsel,:);
      allfoldci=allfoldci(seqsel,:,:);
      for ii=1:length(seqsel)
        naseq=obj.naseq(seqsel(ii));
        fold=allfold(ii,:);
        foldci=squeeze(allfoldci(ii,:,:))';
        [fold,ord]=sort(fold,'desc');
        foldci=foldci(ord,:);
        assert(~any(fold'>foldci(:,2)));
        assert(~any(fold'<foldci(:,1)));
        compounds=allcompounds(ord);
        sel=find(foldci(:,2)>=args.thresh);
        if isempty(sel)
          [~,sel]=max(fold);  % Display at least 1
        end
        sel=sel(1:min(end,args.maxshow));
        fprintf('Plotting %d structures with minfold>=%.1f for naseq %d\n', length(sel), args.thresh, naseq);
        notelist={};
        for kk=1:length(sel)
          k=sel(kk);
          note=sprintf('s=%.1f',fold(k));
          notelist{end+1}=note;
        end
        captions{end+1}=notelist;
        seq=ngs.seqs.naseqformat(naseq);
        ti{end+1}=sprintf('%d %s',naseq,seq);
        pk{end+1}=compounds(sel);
        % Plot cleavage
        obj.plotcleavagebyexpt(naseq,'expts','single','minfold',args.thresh);
        set(gcf,'Color',[1,1,1]);  % White background
        im{end+1}=imtrim(getframe(gcf).cdata);
        close(gcf);
      end
      Compounds.instance().savehtml(args.dirname,pk,'captions',captions,'title',ti,'images',im,'maintitle',sprintf('Summary hits with fold >= %.1f (%s)',args.thresh,datestr(now)));
      fprintf('Saved output in %s/index.html\n',args.dirname);
    end

    function dupewells(summary)
    % Identify samples that appear to be the same CleaveSeq sample sequenced multiple times
      for i=1:length(summary.cleavage)
        c=summary.cleavage(i);
        for j=1:length(c.runs)
          sj=strsplit(c.runs{j},';');
          for k=j+1:length(c.runs)
            sk=strsplit(c.runs{k},';');
            if strcmp(sj{2},sk{2})
              fprintf('cleavage(%d) with target %s has robotrun %s twice with cnts [%d,%d]\n', i, c.target, sj{2},sum(sum(c.cnt(:,2:3,[j,k]),1),2));
            end
          end
        end
      end
    end
    
    function removebootstraps(summary)
    % Remove ratioboot from cleavage (to conserve memory -- large summary won't save)
      for i=1:length(summary.cleavage)
        if length(summary.cleavage(i).runs)==1
          summary.cleavage(i).foldsampsbyrun=[];
        end
      end
    end
    
  end % methods
end % TRPSummary
