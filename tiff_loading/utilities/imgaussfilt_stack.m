function I = imgaussfilt_stack(I,sigma)
for i=1:size(I,3)
    for p=1:size(I,4)
        I(:,:,i,p)=imgaussfilt(I(:,:,i,p),sigma);
    end
end