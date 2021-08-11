% Class for database access
% Singleton (only use static functions)

classdef NGSDatabase < handle
  properties
    host;
    user;
    password;
    schema;   % Current schema 'use'd
    isopen;
  end
  
  properties (Constant)
  end

  methods(Access=private)
    function obj=NGSDatabase()
      fprintf('Creating NGSDatabase instance\n');
      obj.host='XXXXXX';
      obj.user='XXXXXX';
      obj.password='XXXXXX';
      obj.schema='';
      obj.isopen=false;
    end
  end
  
  methods(Static)
    function inst=getinstance()
      persistent theInstance;
      if isempty(theInstance)
        theInstance=NGSDatabase();
      end
      inst=theInstance;
    end
    
    function open(schema)
    % Open a database connection and connect to given schema, close all other connections
      inst=NGSDatabase.getinstance();
      if ~inst.isopen || mysql('status')~=0
        mysql('closeall');
        mysql('open',inst.host,inst.user,inst.password);
        inst.isopen=true;
        inst.schema='';
        inst
      end
      if nargin>=1
        NGSDatabase.use(schema)
      end
    end
    
    function use(schema)
      inst=NGSDatabase.getinstance();
      if ~inst.isopen
        NGSDatabase.open(schema);
      elseif ~strcmp(inst.schema,schema)
        inst.schema=schema;
        mysql('use',inst.schema);
      end
    end
    
    function close()
      inst=NGSDatabase.getinstance();
      mysql('closeall');
      inst.schema='';
      inst.isopen=false;
    end
    
    function enablewrites(password)
    % Enable writes by changing which user access the instbase
      inst=NGSDatabase.getinstance();
      inst.user='XXXXXX';
      inst.password=password;
      mysql('closeall');
      NGSDatabase.open();
    end
    
    function disablewrites()
    % Disable writes by changing which user access the instbase
      inst=NGSDatabase.getinstance();
      inst.user='XXXXXX';
      inst.password='';
      mysql('closeall');
      NGSDatabase.open();
    end

  end

end
