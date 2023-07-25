function imout=im2color2(im,lut)

imout=zeros(size(im,1),size(im,2),3);
for ch=1:3
    imout(:,:,ch)=sum(im.*reshape(lut(:,ch),[1 1 size(im,3)]),3);
end
end