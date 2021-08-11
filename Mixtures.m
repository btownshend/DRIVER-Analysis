% Shadow of database compounds.Mixtures
classdef Mixtures < handle
  properties
    m;    % array of individual mixtures, each with fields: mixture,name,contents
  end

  methods(Access=private)
    function obj=Mixtures()
      obj.load();
    end
  end
  
  methods(Static)
    function obj=instance()
      persistent theInstance
      if isempty(theInstance)
        theInstance=Mixtures();
      end
      obj=theInstance;
    end
  end
  
  methods
    function load(obj,force)
    % Load from database
      NGSDatabase.open();
      [umix,name]=mysql('SELECT m.mixture,m.name FROM compounds.mixtures m');
      obj.m=struct('mixture',num2cell(umix),'name',name);
      [mixture,compound,conc]=mysql('SELECT mixture,compound,concentration FROM compounds.contents');
      for i=1:length(obj.m)
        obj.m(i).contents=compound(mixture==obj.m(i).mixture);
      end
    end
    
    function c=get(obj,pk)
      c=[];
      for i=1:length(pk)
        c=[c,obj.m([obj.m.mixture]==pk(i))];
      end
      if length(c)~=length(pk)
        error('Only found %d/%d requested mixtures', length(c), length(pk));
      end
    end

    function c=find(obj,name)
      c=[];
      if ischar(name)
        name={name};
      end
      for i=1:length(name)
        c=[c,obj.m(strcmp({obj.m.name},name{i}))];
      end
      if length(c)~=length(name)
        error('Only found %d/%d requested mixtures', length(c), length(name));
      end
    end

    function c=findcontains(obj,compoundid)
      sel=arrayfun(@(z) any(compoundid==z.contents),obj.m);
      c=obj.m(sel);
    end

    function c=intersect(obj,names)
      for i=1:length(names)
        mix=obj.find(names{i});
        if i==1
          c=mix.contents;
        else
          c=intersect(c,mix.contents);
        end
      end
    end

    function [c,cnt]=majority(obj,names)
      cnt=[];
      for i=1:length(names)
        mix=obj.find(names{i});
        c=mix.contents;
        if max(c)>length(cnt)
          cnt(max(c))=0;
        end
        cnt(c)=cnt(c)+1;
      end
      c=find(cnt==max(cnt));
      cnt=max(cnt);
    end

  end
end