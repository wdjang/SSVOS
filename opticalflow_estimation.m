function flow_out = opticalflow_estimation(frame_img, prev_frame_img, frame_id, frame_interval, result_path, flow_param)

    flow_path = fullfile(result_path, 'data', 'opticalflow/');
    if ~exist(flow_path,'dir')
        mkdir(flow_path);
    end

    para = [flow_param.alpha,flow_param.ratio,flow_param.minWidth,flow_param.nOuterFPIterations,...
        flow_param.nInnerFPIterations,flow_param.nSORIterations];

    [fvx, fvy, ~] = Coarse2FineTwoFrames(frame_img,prev_frame_img,para);

    flow_out.fvy = fvy;
    flow_out.fvx = fvx;

    save(fullfile(flow_path, sprintf('flow_%04d_%d.mat', frame_id, frame_interval)),'flow_out');
    
end