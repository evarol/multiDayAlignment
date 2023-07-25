function [ufield,phifield,coormap]=rigid2ufield(rigidTform,ref)


% [X,Y]=meshgrid(1:size(ref,1),1:size(ref,2));
% X=permute(X,[2 1]);
% Y=permute(Y,[2 1]);
% coors0=[X(:) Y(:)];
% % A=inv(rigidTform);coorsTformed=coors0*A(1:2,1:2) + A(3,[2 1]);
% A=rigidTform;coorsTformed=coors0*A(1:2,1:2) + A(3,[2 1]);
% % coorsTformed=(coors0-rigidTform(3,[2 1]))*rigidTform(1:2,1:2);
% % coorsTformed=(coors0-A(3,1:2)))*rigidTform(1:2,1:2);
%
% ufield=zeros(size(ref));
% phifield=zeros(size(ref));
% coormap=zeros(size(ref));
% for d=1:size(ref,3)
%     ufield(:,:,d)=reshape(coorsTformed(:,d)-coors0(:,d),size(ref(:,:,d)));
%     phifield(:,:,d)=reshape(coorsTformed(:,d),size(ref(:,:,d)));
%     coormap(:,:,d)=reshape(coors0(:,d),size(ref(:,:,d)));
% end
% ufield=ufield(:,:,[2 1]);
% phifield=phifield(:,:,[2 1]);
% coormap=coormap(:,:,[2 1]);

[M, N,~] = size(ref);
[x, y] = meshgrid(1:N, 1:M);
warpedPoints = transformPointsForward(affine2d(rigidTform), [x(:), y(:)]);
% Extracting the warped x and y coordinates
ufield(:,:,1) = reshape(warpedPoints(:, 1)-x(:), M, N);
ufield(:,:,2) = reshape(warpedPoints(:, 2)-y(:), M, N);
phifield(:,:,1) = reshape(warpedPoints(:, 1), M, N);
phifield(:,:,2) = reshape(warpedPoints(:, 2), M, N);
coormap(:,:,1)= reshape(x(:), M, N);
coormap(:,:,2)= reshape(y(:), M, N);





end