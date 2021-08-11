classdef TRPExpt < matlab.mixin.Copyable
  properties
    run;	% Database
    trpexpt;	% Index in database to this TRP exp
    descr;	% Description
    cleavedsr;	% Subrun for cleaved
    uncleavedsr;% Subrun for uncleaved
    insr;	% Subrun for input
    bulkcleavage; % Bulk cleavage (qPCR), nan if not known
    target;	% Target added, or '' if no target
    targetconc;	% Conc of target
    concunits;	% Units of targetconc
    duration;	% Duration (min)
    idsr;	% Subrun of ID for source of sequences (or [])
    robotrun;	% Index into robot.runs
    mixture;	% Index into compounds.mixtures
    analysis;	% Analysis 
  end
  
  methods
    function obj=TRPExpt(s)  % Convert from struct to class
      if nargin>0
        fn=fieldnames(s);
        for i=1:length(fn)
          obj.(fn{i})=s.(fn{i});
        end
      end
    end
    
    function w=getwell(obj)
    % Retrieve well name corresponding to this sample
      NGSDatabase.open();
      query=sprintf('SELECT well FROM robot.samples WHERE name=''%s'' AND program=(SELECT program FROM robot.runs WHERE run=%d)',...
                    obj.descr,obj.robotrun);
      res=mysql(query);
      w='';
      if isempty(res)
        fprintf('Unable to location well for robotrun %d, name %s\n', obj.robotrun, obj.descr);
      elseif length(res)>1
        error('Got multiple wells for robotrun %d, name %s\n', obj.robotrun, obj.descr);
      else
        w=res{1};
      end
    end
    
  end
end
