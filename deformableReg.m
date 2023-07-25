function [movingApplied_reg,sumVfield]=deformableReg(fixed,moving,numblocks,numiter,sigma,brightness)
fixedApplied=fixed;
movingApplied=moving;
myimfuse = @(x,y)(imfuse(x,y,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
moving_reg{1}=moving;
for t=1:numiter
    tic;
    vfield{t}=repmat(zeros(size(fixed)),[1 1 2]);
    efield{t}=zeros(size(fixed));
    dims=size(fixed);
    
    for i=1:numblocks
        for j=1:numblocks
            block_idx{i,j,1}=(i-1)*floor(dims(1)/numblocks)+1:i*floor(dims(1)/numblocks);
            block_idx{i,j,2}=(j-1)*floor(dims(2)/numblocks)+1:j*floor(dims(2)/numblocks);
            fixed_block{i,j}=fft2(fixed(block_idx{i,j,1},block_idx{i,j,2}));
            moving_block{i,j}=fft2(moving_reg{t}(block_idx{i,j,1},block_idx{i,j,2}));
        end
    end
    
    for i=1:numblocks
        for j=1:numblocks
            [output, Greg] = dftregistration(moving_block{i,j},fixed_block{i,j},10);
            if isnan(output(1))
                output(1)=1;
            end
            
            vfield{t}(block_idx{i,j,1},block_idx{i,j,2},1)=output(4)*(1-output(1));
            vfield{t}(block_idx{i,j,1},block_idx{i,j,2},2)=output(3)*(1-output(1));
            efield{t}(block_idx{i,j,1},block_idx{i,j,2})=(1-output(1));
        end
    end
    
    for i=1:2
        vfield{t}(:,:,i)=imgaussfilt(vfield{t}(:,:,i),sigma);
    end
    
    sumVfield=zeros(size(vfield{1}));
    for tt=1:t
        sumVfield=sumVfield+vfield{tt};
    end
    %     for i=1:2
    %         sumVfield(:,:,i)=imgaussfilt(sumVfield(:,:,i),sigma);
    %     end
    moving_reg{t+1}=imwarp(moving,sumVfield);
    
    for i=1:size(movingApplied,3)
        movingApplied_reg{t}(:,:,i)=imwarp(movingApplied(:,:,i),sumVfield);
    end
    obj(t)=corr(moving_reg{t+1}(:),fixed(:));
    subplot(2,4,1)
    plot(obj,'LineWidth',2);
    title('Corr');
    xlabel('Iter');
    ylabel('Corr');
    set(gca,'FontSize',20,'FontWeight','bold');drawnow
    subplot(2,4,[2:4 6:8])
    imagesc(brightness*myimfuse(fixed,moving_reg{t+1}));axis off;set(gcf,'color','w');drawnow
    subplot(2,4,5)
    imagesc(sqrt(sum(sumVfield.^2,3)));colorbar;axis off;drawnow
    title('Displacement magnitutde');
    set(gca,'FontSize',20,'FontWeight','bold');drawnow
    toc
end





end