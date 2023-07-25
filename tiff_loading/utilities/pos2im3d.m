function I=pos2im3d(pos)

I=zeros(256,512,21*4);

pos = max(round(pos/2),1);

for i=1:size(pos,1)
    I(pos(i,2),pos(i,1),pos(i,3))=i;
end
% I=flip(I,3);
I = reshape(I,[256 512 4 21]);
I = repmat(squeeze(max(I,[],3)),[1 1 1 4]);
end
