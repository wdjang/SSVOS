function [fore_flow, back_flow] = Estimate_Optical_Flow(frame_list, flow_dir)
    
    if exist(flow_dir,'dir') == 0
        mkdir(flow_dir);
    end
    
    %% set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
    alpha = 0.012;
    ratio = 0.75;
    minWidth = 20;
    nOuterFPIterations = 7;
    nInnerFPIterations = 1;
    nSORIterations = 30;

    para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];
    
    if exist(fullfile(flow_dir,'fore_flow_new.mat'),'file')
        load(fullfile(flow_dir,'fore_flow_new.mat'));
        load(fullfile(flow_dir,'back_flow_new.mat'));
    else
        
        %% Forward optical flow
        fore_flow = cell(length(frame_list),1);

        for img_id = 1:length(frame_list)-1
            im1 = frame_list{img_id};
            im2 = frame_list{img_id+1};

%             [vx, vy, ~] = Coarse2FineTwoFrames(im1,im2,para);
            uv = EPPM(im1,im2);
            vx = uv(:,:,1);
            vy = uv(:,:,2);

            fore_flow{img_id}.x = vx;
            fore_flow{img_id}.y = vy;

%             % visualize flow field
%             clear flow;
%             flow(:,:,1) = xflow_list{img_id};
%             flow(:,:,2) = yflow_list{img_id};
%             imflow = flowToColor(flow);
%     
%             figure;imshow(imflow);
        end
        
        fore_flow{length(frame_list)}.x = zeros(size(vx));
        fore_flow{length(frame_list)}.y = zeros(size(vy));
        
        %% Backward optical flow
        back_flow = cell(length(frame_list),1);

        for img_id = 2:length(frame_list)
            im1 = frame_list{img_id};
            im2 = frame_list{img_id-1};

            [vx, vy, ~] = Coarse2FineTwoFrames(im1,im2,para);

            back_flow{img_id}.x = vx;
            back_flow{img_id}.y = vy;

%             % visualize flow field
%             clear flow;
%             flow(:,:,1) = xflow_list{img_id};
%             flow(:,:,2) = yflow_list{img_id};
%             imflow = flowToColor(flow);
%     
%             figure;imshow(imflow);
        end
        
        back_flow{1}.x = zeros(size(vx));
        back_flow{1}.y = zeros(size(vy));
        
        %% Save optical flows
        save(fullfile(flow_dir,'fore_flow_new.mat'),'fore_flow');
        save(fullfile(flow_dir,'back_flow_new.mat'),'back_flow');
        
    end
    
end
