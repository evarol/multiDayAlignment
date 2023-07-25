function vfieldTotal=pointDeformerIterative(moving,fixed,sigma,brightness)
close all
% sigma=100;


myimfuse = @(x,y)(imfuse(x,y,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));

% I=imread('cameraman.tif');axis equal;axis off
% imagesc(I);
notDone=1;

moving_warped=moving;
vfieldTotal=repmat(zeros(size(moving)),[1 1 2]);
while notDone==1
    imagesc(brightness*myimfuse(moving_warped,fixed));axis equal;axis off
    title('Click successive pairs of matching landmarks - Chose RED first then GREEN  - right click when done');
    
    [x,y]=getpts;
    x(end)=[];
    y(end)=[];
    
    fixed_matches=round([x(mod(1:length(x),2)==1) y(mod(1:length(x),2)==1)]);
    moving_matches=round([x(mod(1:length(x),2)==0) y(mod(1:length(x),2)==0)]);
    
    
    vfield=repmat(zeros(size(moving)),[1 1 2]);
    
    if ~isempty(moving_matches)
    for i=1:size(moving_matches,1)
        vfield(moving_matches(i,2),moving_matches(i,1),2)=fixed_matches(i,2)-moving_matches(i,2);
        vfield(moving_matches(i,2),moving_matches(i,1),1)=fixed_matches(i,1)-moving_matches(i,1);
    end
    for i=1:size(vfield,3)
        vfield(:,:,i)=imgaussfilt(vfield(:,:,i),sigma,'FilterDomain','auto','FilterSize',4*ceil(2*sigma)+1);
    end
    
    
    
    
    clear A b step
    for d=1:2
        for i=1:size(moving_matches,1)
            A(i,d)=vfield(moving_matches(i,2),moving_matches(i,1),d);
            b(i,d)=fixed_matches(i,d)-moving_matches(i,d);
        end
        step(1,1,d)=linsolve(A(:,d),b(:,d));
    end
    
    step(or(isnan(step),isinf(step)))=0;
    
    
    vfieldTotal=vfieldTotal+vfield.*step;
    end
    moving_warped=imwarp(moving,vfieldTotal);
    moving_matches_moved=zeros(size(moving_matches));
    if ~isempty(moving_matches)
    for d=1:2
        for i=1:size(moving_matches,1)
            moving_matches_moved(i,d)=moving_matches(i,d)+vfield(moving_matches(i,2),moving_matches(i,1),d)*step(d);
        end
    end
    end
    
    imagesc(brightness*myimfuse(moving_warped,fixed));axis equal;axis off
    done=input('Done?(1) Add more?(0), Start Over(2)');
    if done==1
        notDone=0;
    elseif done==0
    elseif done==2
        vfieldTotal=zeros(size(vfieldTotal));
        moving_warped=moving;
    end
end
end