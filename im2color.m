function imout=im2color(im,lut)

imout=zeros(size(im,1),size(im,2),3);
for ch=1:3
    imout(:,:,ch)=im*lut(ch);
end
end