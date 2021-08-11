function im=imtrim(im)
% Trim an image to bounding box
r1=max(range(im),[],3);
r2=max(range(im,2),[],3);
sel1=find(r2~=0,1,'first'):find(r2~=0,1,'last');
sel2=find(r1~=0,1,'first'):find(r1~=0,1,'last');
im=im(sel1,sel2,:);
end
