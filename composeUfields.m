function ufieldComposed=composeUfields(ufield1,ufield2,coormap,fast)
%Outputs: u_total = u_2(u_1(x))
if nargin<4
    fast=1;
end
vec=@(x)(x(:));
phifield1=coormap+ufield1;
phifield2=coormap+ufield2;
F1=griddedInterpolant(coormap(:,:,1)',coormap(:,:,2)',phifield2(:,:,1)','linear');
F2=griddedInterpolant(coormap(:,:,1)',coormap(:,:,2)',phifield2(:,:,2)','linear');
if fast==1
    grab = @(a,b)(a(b));
    inRange1=and(and(phifield1(:,:,1)<=size(coormap,1),phifield1(:,:,1)>=1),and(phifield1(:,:,2)<=size(coormap,2),phifield1(:,:,1)>=1));
    A=zeros(size(coormap));tmp=A(:,:,1);
    tmp1=F1(grab(phifield1(:,:,1),inRange1==1),grab(phifield1(:,:,2),inRange1==1));
    tmp2=F2(grab(phifield1(:,:,1),inRange1==1),grab(phifield1(:,:,2),inRange1==1));
    tmp(inRange1==1)=tmp1;A(:,:,1)=tmp;
    tmp(inRange1==1)=tmp2;A(:,:,2)=tmp;
else
    A=zeros(size(coormap));
    A(:,:,1)=reshape(F1(vec(phifield1(:,:,1)),vec(phifield1(:,:,2))),size(A,1),size(A,2));
    A(:,:,2)=reshape(F2(vec(phifield1(:,:,1)),vec(phifield1(:,:,2))),size(A,1),size(A,2));
end

ufieldComposed=A-coormap;
% ufieldComposed(isnan(ufieldComposed))=0;
end