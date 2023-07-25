clear all
% clc
close all

%% parameters

brightness=1; %increase to make images easier to see
channel_to_use=2; %for multi channel images, select the channel for alignment e.g. where GCaMP is
is3d=0; %toggle 1 if the images are 3D


%% paths to utilities + functions
addpath ./tiff_loading/utilities
addpath(genpath('./tiff_loading/Fiji.app'));
javaaddpath('./tiff_loading/Fiji.app/mij.jar');
myimfuse = @(x,y)(imfuse(x,y,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
strpart = @(x,y)(x{y});
minmax = @(x)((x-min(x(:)))./max(x(:)-min(x(:))));

disp('Select images to stitch (IMPORTANT: click on them sequentially in increasing order)');
[file,path] = uigetfile('*','MultiSelect','on');

for i=1:length(file)
    filename{i}=[path file{i}];
end

%%
for i=1:length(filename)
    if is3d==0
        try
            I{i}=double(load_tiff(filename{i}));
            I{i}=zscore(I{i}(:,:,channel_to_use),[],'all');
        catch
            I{i}=double(imread(filename{i}));
            I{i}=zscore(I{i}(:,:,channel_to_use),[],'all');
        end
    else
        try
            I{i}=double(load_tiff(filename{i}));
            I{i}=max(zscore(I{i}(:,:,:,channel_to_use),[],'all'),[],3);
        catch
            I{i}=double(imread(filename{i}));
            I{i}=max(zscore(I{i}(:,:,:,channel_to_use),[],'all'),[],3);
        end
    end
end
MIJ.closeAllWindows()
%%
dims=zeros(1,2);
for i=1:length(I)
    dims=max([dims; size(I{i})],[],1);
end

for i=1:length(I)
    I{i}=padarray(I{i},[dims(1)-size(I{i},1) dims(2)-size(I{i},2) 0],'post');
end

%%
for ii=1:length(I)-1
    sliceA=ii;
    sliceB=ii+1;
    
    I1slice=max(I{sliceA},[],3);
    I2slice=max(I{sliceB},[],3);
    I2slice_warped=I2slice;
    
    notGood=1;
    rigid=1;
    affine=0;
    globTform=eye(3);inv_globTform=eye(3);
    D=zeros(size(I2slice,1),size(I2slice,2),2);
    while notGood==1
        close all
        imagesc(brightness*imfuse(I1slice,I2slice_warped,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
        title('Click successive pairs of matching landmarks - Chose RED first then GREEN  - right click when done');
        
        [x,y]=getpts;
        
        x(end)=[];
        y(end)=[];
        
        matches1=[x(mod(1:length(x),2)==1) y(mod(1:length(x),2)==1)];
        matches2=[x(mod(1:length(x),2)==0) y(mod(1:length(x),2)==0)];
        
        min_pts=min(size(matches1,1),size(matches2,1));
        matches1=matches1(1:min_pts,:);
        matches2=matches2(1:min_pts,:);
        
        
        D_previous=D;
        globTform_previous=globTform;
        inv_globTform_previous=inv_globTform;
        D_inverse_previous=rigid2ufield(inv_globTform_previous,zeros(size(I1slice,1),size(I1slice,2),2));
        DX_inverse_previous=D_inverse_previous(:,:,1);
        DY_inverse_previous=D_inverse_previous(:,:,2);
        if ~isempty(matches2)
            matches2_in_original_space=[matches2(:,1)-DX_inverse_previous(sub2ind(size(DX_inverse_previous),round(matches2(:,1)),round(matches2(:,2)))) matches2(:,2)-DY_inverse_previous(sub2ind(size(DY_inverse_previous),round(matches2(:,1)),round(matches2(:,2))))];%proj(addone(matches2)*globTform_previous);
        end
        %% debug stuff
        %         close all
        %         subplot(1,3,1)
        %         imagesc(brightness*imfuse(I1slice,zeros(size(I2slice)),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
        %         hold on
        %         plot(matches1(:,1),matches1(:,2),'go','MarkerSize',20);
        %         subplot(1,3,2)
        %         imagesc(brightness*imfuse(zeros(size(I1slice)),I2slice,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
        %         hold on
        %         plot(matches2_in_original_space(:,1),matches2_in_original_space(:,2),'ro','MarkerSize',20);
        
        %%
        if and(~isempty(matches1),~isempty(matches2))
            if rigid==1
                [R,T]=wahba(matches1,matches2_in_original_space);
                tform=affine2d([R zeros(2,1);T 1]);
                [R,T]=wahba(matches2_in_original_space,matches1);
                inv_tform=affine2d([R zeros(2,1);T 1]);
            elseif affine==1
                tform = fitgeotrans(matches1,matches2_in_original_space,'affine');
                inv_tform = fitgeotrans(matches2_in_original_space,matches1,'affine');
            end
        else
            tform=affine2d(eye(3));
            inv_tform=affine2d(eye(3));
        end
        
        globTform=tform.T;
        inv_globTform=inv_tform.T;
        D=rigid2ufield(globTform,zeros(size(I1slice,1),size(I1slice,2),2));
        I2slice_warped=imwarp(I2slice,D);
        imagesc(brightness*myimfuse(I1slice,I2slice_warped));
        
        if affine==0
            switch input('Start over? (0) Is it good enough? (1)  Switch to affine? (2) Switch to rigid? (3)');
                case 0
                    rigid=1;
                    affine=0;
                    globTform=eye(3);
                    inv_globTform=eye(3);
                    D=rigid2ufield(globTform,zeros(size(I1slice,1),size(I1slice,2),2));
                    I2slice_warped=imwarp(I2slice,D);
                case 1
                    notGood=0;
                case 2
                    rigid=0;
                    affine=1;
                case 3
                    rigid=1;
                    affine=0;
            end
        else
            notGood=0;
        end
        
    end
    
    radius=40;
    sigma=50;
    
    close all
    I2slice_globtform=imwarp(I2slice,rigid2ufield(globTform,zeros(size(I1slice,1),size(I1slice,2),2)));
    imagesc(brightness*myimfuse(I2slice_globtform,I1slice));
    clear I2transformed*
    for i=1:size(I{sliceB},3)
        I2transformed(:,:,i)=imwarp(I{sliceB}(:,:,i),rigid2ufield(globTform,zeros(size(I1slice,1),size(I1slice,2),2)));
    end
    I1=I{sliceA};
    I2=I{sliceB};
    
    
    close all
    [final_pass,vfield_final]=deformableReg(max(I1,[],3),max(I2transformed,[],3),20,40,200,brightness);
    
    for i=1:size(I{sliceB},3)
        I2transformed_warped(:,:,i)=imwarp(I2transformed(:,:,i),vfield_final);
    end
    
    
    vfieldTotal=pointDeformerIterative(max(I2transformed_warped,[],3),max(I1,[],3),sigma,brightness);
    close all
    
    for i=1:size(I{sliceB},3)
        I2transformed_warped_final(:,:,i)=imwarp(I2transformed(:,:,i),vfieldTotal+vfield_final);
    end
    
    %%
    
    figure
    [ax, ~] = tight_subplot(1, 4, [0.01 0.01], 0.1,0.1);
    axes(ax(1));
    imagesc(brightness*myimfuse(max(I2,[],3),max(I1,[],3)));title('Unregistered');axis equal;axis off;text(100,100,'Unregistered','color','w','FontWeight','bold','FontSize',20);
    axes(ax(2))
    imagesc(brightness*myimfuse(max(I2transformed,[],3),max(I1,[],3)));title('Rigid/Affine');axis equal;axis off;text(100,100,'Rigid/affine','color','w','FontWeight','bold','FontSize',20);
    axes(ax(3))
    imagesc(brightness*myimfuse(max(I2transformed_warped,[],3),max(I1,[],3)));title('Deformable Manual');axis equal;axis off;text(100,100,'Deformable Manual','color','w','FontWeight','bold','FontSize',20);
    
    axes(ax(4))
    imagesc(brightness*myimfuse(max(I2transformed_warped_final,[],3),max(I1,[],3)));title('Deformable Auto');axis equal;axis off;text(100,100,'Deformable Auto','color','w','FontWeight','bold','FontSize',20);
    linkaxes([ax(1) ax(2) ax(3) ax(4)]);
    set(gcf,'color','w');
    
    
    %%
    
    a=strsplit(filename{ii+1},{'/','.'});
    apath=fileparts(filename{ii+1});
    b=strsplit(filename{ii},{'/','.'});
    save([apath '/' a{end-1} '_TO_' b{end-1} '_deformation_map.mat'],'vfieldTotal','vfield_final','globTform');
    
end