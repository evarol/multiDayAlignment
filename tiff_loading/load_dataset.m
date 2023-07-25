% clear all
% clc
% close all
vec=@(x)(x(:));

addpath ./jsonlab/
addpath ./utilities/
addpath(genpath('./Fiji.app'));
javaaddpath('./Fiji.app/mij.jar')

% I=load_tiff('/Users/erdem/Dropbox/LateBoi/sah022/functional_z_stacks/sah022_zstack_final.tiff');

% I=load_tiff('/Users/erdem/Dropbox/LateBoi/sah022/Chubstack/2020-09-10_OligodT_stitched_Ch1-1.tif');
% I=load_tiff('/Users/erdem/Dropbox/LateBoi/byr005_first_o_pool/F1_UG.tiff');
% I=load_tiff('/Users/erdem/Dropbox/LateBoi/sah046/sah046_zstack_final.tiff');

if ~exist('file');
    disp('Select ims files to stitch...');
    [file,path] = uigetfile('*','Select tiff file to load','MultiSelect','on');
end

I=load_tiff([path file]);

% I=load_tiff('/Users/erdem/Downloads/dendritic_zstacks/pre/time_avg_Ch1.tif');
% I11=load_tiff('/Users/erdem/Downloads/dendritic_zstacks/pre/time_avg_Ch1.tif');
% I12=load_tiff('/Users/erdem/Downloads/dendritic_zstacks/pre/time_avg_Ch2.tif');
% I21=load_tiff('/Users/erdem/Downloads/dendritic_zstacks/post/time_avg_Ch1.tif');
% I22=load_tiff('/Users/erdem/Downloads/dendritic_zstacks/post/time_avg_Ch2.tif');
try
Zproj=double(squeeze(max(I,[],3)));%Zproj=Zproj./imgaussfilt(Zproj,100);
Xproj=double(squeeze(max(I,[],2)));%Xproj=Xproj./imgaussfilt(Xproj,100);
Yproj=double(squeeze(max(I,[],1))');%Yproj=Yproj./imgaussfilt(Yproj,100);

A=[log1p(Zproj) 10*ones(size(I,1),1) log1p(Xproj);10*ones(1,size(I,2)+size(I,3)+1);log1p(Yproj) 10*ones(size(I,3),1) mean(log1p(Zproj(:)))*ones(size(I,3),size(I,3))];

figure('units','normalized','outerposition',[0 0 1 1])
imagesc(A);colormap(gray(256));axis equal;axis off
text(size(I,1)+20,size(I,2)+20,'X-Y-Z Projections','Color','w','FontWeight','bold','FontSize',20);
end
% export_fig('aod_projections.png');