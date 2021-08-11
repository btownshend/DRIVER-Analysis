% Convert a list of integers to a formatted string
function s=intlist(x)
s=sprintf('%d',x(1));
if length(x)>1
  s=['[',s,sprintf(',%d',x(2:end)),']'];
end
end


  