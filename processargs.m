% Process args
% Usage:
%  defaults=struct(name,value,name,value...)
%  args=processargs(defaults,varargin);
function args=processargs(args,va)
i=1;
while i<=length(va)
  if ~isfield(args,va{i})
    s=sprintf('Unknown option: "%s"\nValid options are:',va{i});
    fn=fieldnames(args);
    for i=1:length(fn)
      if ~strncmp(fn{i},'SET',3)
        s=[s,sprintf('%s ',fn{i})];
      end
    end
    error(s);
  end
  args.(['SET',va{i}])=true;   % Flag that it was set explicitly
  if islogical(args.(va{i})) && (i==length(va) || ischar(va{i+1}))
    % Special case of option string without value for boolean, assume true
    args.(va{i})=true;
  else
    args.(va{i})=va{i+1};
    i=i+1;
  end
  i=i+1;
end
