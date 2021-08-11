% Display a table with any unneeded space removed
function packtable(x)
t=evalc('x');
t=strrep(strrep(t,'<strong>',''),'</strong>','');
s=char(strsplit(t,'\n'));  % Character array
s=s(4:end,:);   % Remove header
s(s==''''| s=='{'|s=='}')=' ';
for i=1:size(s,2)
  if i==1
    repress=s(:,i)==' ';
  else
    repress=s(:,i)==' ' & (s(:,i-1)==' ' | s(:,i-1)=='@');
  end
  s(repress,i)='@';
end
keep=any(~ismember(s([1,3:end],:),['''','@','{','}']))|s(2,:)==' ';
keeph=keep;
for i=length(keeph):-1:2
  if ~keeph(i) && s(1,i)~='@'
    % Move back space removal
    for k=i-1:-1:1
      if keeph(k)
        break;
      end
    end
    keeph(k)=0;
    keeph(i)=1;
  end
end
sout=[s(1,keeph);s(2:end,keep)];
sout(:)=strrep(sout(:)','@',' ');
%sout(2,all(sout(3:end,:)==' ',1))=' ';
res=strjoin(cellstr(sout),'\n');
disp(res);
end

