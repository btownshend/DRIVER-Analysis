% Concatenate 2 structures
% Makes the set of fields the union of the two
function s=concatstructs(s1,s2)
s=s1;
if ~isempty(s2)
  f2=fieldnames(s2);
  for i=1:length(s2)
    for j=1:length(f2)
      s(length(s1)+i).(f2{j})=s2(i).(f2{j});
    end
  end
end
