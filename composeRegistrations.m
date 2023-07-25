clear all
close all
clc

%% parameters

brightness=1; %increase to make images easier to see
channel_to_use=1; %for multi channel images, select the channel for alignment e.g. where GCaMP is
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
for i=(length(file)):-1:2
    tforms{i}=load([path strpart(strsplit(file{i},'.'),1) '_TO_' strpart(strsplit(file{i-1},'.'),1) '_deformation_map.mat']);
end

dims=zeros(1,2);
for i=1:length(I)
    dims=max([dims; size(I{i})],[],1);
end

for i=1:length(I)
    I{i}=padarray(I{i},[dims(1)-size(I{i},1) dims(2)-size(I{i},2) 0],'post');
end

% for i=2:length(tforms)
%     for d=1:2
%         tforms{i}.vfieldTotal=padarray(tforms{i}.vfieldTotal,[dims(1)-size(tforms{i}.vfieldTotal,1) dims(2)-size(tforms{i}.vfieldTotal,2) 0],'post');
%         tforms{i}.vfield_final=padarray(tforms{i}.vfield_final,[dims(1)-size(tforms{i}.vfield_final,1) dims(2)-size(tforms{i}.vfield_final,2) 0],'post');
%     end
% end

% composing rigid + deformable transformations per slice

for i=2:length(tforms)
    [ufield_rigid,~,coormap]=rigid2ufield(tforms{i}.globTform,zeros(size(tforms{i}.vfieldTotal)));
    ufield_deformable=tforms{i}.vfieldTotal+tforms{i}.vfield_final;
    tic;
    total_ufield{i}=composeUfields(ufield_rigid,ufield_deformable,coormap,1);
    toc
end

total_ufield{1}=zeros(size(total_ufield{end}));
backtracked_ufield{1,1}=zeros(size(coormap));
total_rotation{1}=eye(3);
registered{1}=I{1};
rigid_registered{1}=I{1};

composed_field{1}=zeros(size(I{1},1),size(I{1},2),2);
% back tracking transformations for each slice to slice 1
for i=2:length(tforms)
    composed_field{i}=zeros(size(I{i},1),size(I{i},2),2);
    for z=1:size(I{i},3)
        registered{i}(:,:,z)=I{i}(:,:,z);
    end
    backtracked_ufield{i,i-1}=total_ufield{i};
    for z=1:size(I{i},3)
        registered{i}(:,:,z)=imwarp(registered{i}(:,:,z),total_ufield{i});
    end
    composed_field{i}=composeUfields(composed_field{i},total_ufield{i},coormap,1);
    
    total_rotation{i}=eye(3);
    if i>2
        for j=(i-1):-1:2
            tic
            backtracked_ufield{i,j-1}=composeUfields(backtracked_ufield{i,j},total_ufield{j},coormap,1);
            total_rotation{i}=tforms{j}.globTform*total_rotation{i};
            composed_field{i}=composeUfields(composed_field{i},total_ufield{j},coormap,1);
            for z=1:size(I{i},3)
                registered{i}(:,:,z)=imwarp(registered{i}(:,:,z),total_ufield{j});
            end
            toc
            [i j]
        end
    end
    
end

for i=1:length(I)
    registered_composed{i}=imwarp(I{i},composed_field{i});
end


for i=1:length(I)
    for z=1:size(I{i},3)
        rigid_registered{i}(:,:,z)=imwarp(I{i}(:,:,z),affine2d(total_rotation{i}),"OutputView",imref2d(size(I{i}(:,:,1))));
    end
end


%% final pass
for i=1:length(I)
[final_pass{i},vfield_final]=deformableReg(I{1}(:,:,1),registered{i}(:,:,1),40,40,100,brightness);
end

%% final pass
for i=1:length(I)
[final_pass_composed{i},vfield_final]=deformableReg(I{1}(:,:,1),registered_composed{i}(:,:,1),40,40,100,brightness);
end



%% 
%% Pairwise aligned slices
figure(1)
for i=1:length(I)-1
    ax1(i)=subplot(2,4,i);
    imagesc(brightness*myimfuse(I{i}(:,:,1),imwarp(I{i+1}(:,:,1),total_ufield{i+1})));
    title(['Pairwise: Red: Slice ' num2str(i) ' - Green: Slice ' num2str(i+1)]);drawnow
    set(gca,'FontSize',14,'FontWeight','bold');
%     axis equal;axis off
end
set(gcf,'color','w');
linkaxes(ax1)


%% Globally aligned slices
figure(2)
for i=1:length(I)
    ax2(i)=subplot(2,4,i);
    imagesc(brightness*myimfuse(I{1}(:,:,1),registered{i}(:,:,1)));    
    title({'Global alignment:',['Red: ' strpart(strsplit(filename{1},'/'),length(strsplit(filename{1},'/')))],['Green: ' strpart(strsplit(filename{i},'/'),length(strsplit(filename{i},'/')))]},'interpreter','none');drawnow
    set(gca,'FontSize',14,'FontWeight','bold');
%     axis equal;axis off
end
set(gcf,'color','w');
linkaxes(ax2)

%% un registered
figure(3)
for i=1:length(I)-1
    ax3(i)=subplot(2,4,i);
    imagesc(brightness*myimfuse(I{1}(:,:,1),I{i+1}(:,:,1)));
    title({'Unaligned:',['Red: ' strpart(strsplit(filename{1},'/'),length(strsplit(filename{1},'/')))],['Green: ' strpart(strsplit(filename{i},'/'),length(strsplit(filename{i},'/')))]},'interpreter','none');drawnow
    set(gca,'FontSize',14,'FontWeight','bold');
%     axis equal;axis off
end
set(gcf,'color','w');
linkaxes(ax3)

%% Globally aligned slices - final pass
figure(4)
for i=1:length(I)
    ax4(i)=subplot(2,4,i);
    imagesc(brightness*myimfuse(I{1}(:,:,1),final_pass{i}{end}(:,:,1)));    
    title({'Global alignment:',['Red: ' strpart(strsplit(filename{1},'/'),length(strsplit(filename{1},'/')))],['Green: ' strpart(strsplit(filename{i},'/'),length(strsplit(filename{i},'/')))]},'interpreter','none');drawnow
    set(gca,'FontSize',14,'FontWeight','bold');
%     axis equal;axis off
end
set(gcf,'color','w');

%% Globally aligned slices through composition - final pass
figure(5)
for i=1:length(I)
    ax5(i)=subplot(2,4,i);
    imagesc(brightness*myimfuse(I{1}(:,:,1),final_pass_composed{i}{end}(:,:,1)));    
    title({'Global alignment (composed):',['Red: ' strpart(strsplit(filename{1},'/'),length(strsplit(filename{1},'/')))],['Green: ' strpart(strsplit(filename{i},'/'),length(strsplit(filename{i},'/')))]},'interpreter','none');drawnow
    set(gca,'FontSize',14,'FontWeight','bold');
%     axis equal;axis off
end
set(gcf,'color','w');
linkaxes([ax2 ax4])

%% Globally aligned slices - final pass
figure(6)
ha=tight_subplot(2,length(I),[0.01 0.01],0.1,0.1);
for i=1:length(I)
    axes(ha(i))
    imagesc(brightness*myimfuse(I{1}(:,:,1),final_pass{i}{end}(:,:,1)));    
    title({'Alignment overlay:',['Red: ' strpart(strsplit(filename{1},'/'),length(strsplit(filename{1},'/')))],['Green: ' strpart(strsplit(filename{i},'/'),length(strsplit(filename{i},'/')))]},'interpreter','none');drawnow
    axis equal;set(gca,'xtick',[],'ytick',[],'ticklength',[0 0]);
    axes(ha(length(I)+i))
    imagesc(brightness*myimfuse(zeros(size(I{i}(:,:,1))),minmax(final_pass{i}{end}(:,:,1))));    
    xlabel({'Aligned individual day only:',['Green: ' strpart(strsplit(filename{i},'/'),length(strsplit(filename{i},'/')))]},'interpreter','none');drawnow
    set(gca,'FontWeight','bold');
    axis equal;set(gca,'xtick',[],'ytick',[],'ticklength',[0 0]);
end
set(gcf,'color','w');
linkaxes(ha)

%% Globally aligned slices through composition - final pass
figure(7)
ha2=tight_subplot(2,length(I),[0.01 0.01],0.1,0.1);
for i=1:length(I)
    axes(ha2(i))
    imagesc(brightness*myimfuse(zeros(size(I{i}(:,:,1))),minmax(final_pass{i}{end}(:,:,1))));  
    title({'Alignment overlay:',['Red: ' strpart(strsplit(filename{1},'/'),length(strsplit(filename{1},'/')))],['Green: ' strpart(strsplit(filename{i},'/'),length(strsplit(filename{i},'/')))]},'interpreter','none');drawnow
    axis equal;set(gca,'xtick',[],'ytick',[],'ticklength',[0 0]);
    axes(ha2(length(I)+i))
    imagesc(brightness*myimfuse(zeros(size(I{i}(:,:,1))),minmax(final_pass_composed{i}{end}(:,:,1))));
    xlabel({'Aligned individual day only:',['Green: ' strpart(strsplit(filename{i},'/'),length(strsplit(filename{i},'/')))]},'interpreter','none');drawnow
    set(gca,'FontWeight','bold');
    axis equal;set(gca,'xtick',[],'ytick',[],'ticklength',[0 0]);
end
set(gcf,'color','w');
linkaxes(ha2)
%% save data
% 
% for i=1:length(I)
%     disp(['Saving aligned slice ' num2str(i)]);
%     im=registered{i};
%     save([path strpart(strsplit(file{i},'.'),1) '_to_' strpart(strsplit(file{1},'.'),1) '_globally_deformed.mat'],'im');
% end
