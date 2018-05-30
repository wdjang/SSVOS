function segment_track = SSVOS(seq_path, ann_path, result_path, param_list)

% List of frames
frame_list = dir(fullfile(seq_path,'*.jpg'));
if isempty(frame_list)
    frame_list = dir(fullfile(seq_path,'*.bmp'));
    if isempty(frame_list)
        frame_list = dir(fullfile(seq_path,'*.png'));
    end
end
% User annotation at the first frame
ann_list = dir(fullfile(ann_path,'*.jpg'));
if isempty(ann_list)
    ann_list = dir(fullfile(ann_path,'*.bmp'));
    if isempty(ann_list)
        ann_list = dir(fullfile(ann_path,'*.png'));
    end
end
% Make result directory
if ~exist(result_path,'dir')
    mkdir(result_path);
end

% Make resultant segmentation track
segment_track = cell(length(frame_list),1);

%% At the first frame
frame_id = 1;
fprintf('Processing frame %d\n',frame_id);
prev_frame_img = imread(fullfile(seq_path,frame_list(frame_id).name)); % Load the first frame
prev_sp_map = superpixel_generation(prev_frame_img, frame_id, result_path, param_list.slic); % Generate SLIC superpixels
prev_num_sp = max(prev_sp_map(:)); % Number of superpixels
[h_size, w_size] = size(prev_sp_map); % Height and Width
hw_size = h_size*w_size; % Number of pixels

prev_result_map = im2double(imread(fullfile(ann_path,ann_list(frame_id).name))); % Load a user annotation at the first frame
prev_result_map = imresize(prev_result_map,[h_size,w_size]);
prev_fgprob = zeros(prev_num_sp,1);
for sp_id = 1:prev_num_sp
    prev_fgprob(sp_id) = mean(prev_result_map(prev_sp_map==sp_id));
end
prev_result_vec = double(prev_fgprob > 0.5);
imwrite(prev_result_map,fullfile(result_path,sprintf('%04d.png',frame_id)));
segment_track{frame_id} = prev_result_map;

rgb_p = reshape(im2double(prev_frame_img),hw_size,3);
rgb_sp = zeros(prev_num_sp,3);
for sp_id = 1:prev_num_sp
    find_list = (prev_sp_map==sp_id);
    rgb_sp(sp_id,:) = mean(rgb_p(find_list,:),1);
end
prev_sp_feat = colorspace('Lab<-', rgb_sp)';

%% From the second to the last frames
for frame_id = 2:length(frame_list)
    fprintf('Processing frame %d\n',frame_id);
    frame_img = imread(fullfile(seq_path, frame_list(frame_id).name)); % Load frame
    sp_map = superpixel_generation(frame_img, frame_id, result_path, param_list.slic); % Generate superpixel
    num_sp = max(sp_map(:)); % Number of superpixels
    
    y_list = repmat([1:h_size]',1,w_size);
    x_list = repmat(1:w_size,h_size,1);
    
    flow = opticalflow_estimation(frame_img, prev_frame_img, frame_id, 1, result_path, param_list.flow);
    y_bwarp = max(min(round(y_list + flow.fvy),h_size),1);
    x_bwarp = max(min(round(x_list + flow.fvx),w_size),1);
    bwarp_list = sub2ind([h_size,w_size], y_bwarp, x_bwarp);
    
    if frame_id > 2
        flow = opticalflow_estimation(frame_img, third_frame_img, frame_id, 2, result_path, param_list.flow);
        y_bwarp = max(min(round(y_list + flow.fvy),h_size),1);
        x_bwarp = max(min(round(x_list + flow.fvx),w_size),1);
        bwarp_list2 = sub2ind([h_size,w_size], y_bwarp, x_bwarp);
    end
    
    conn_edge_mat = k_ring_graph_construction(sp_map, param_list.fg_k);
    four_edge_mat = k_ring_graph_construction(sp_map, param_list.bg_k);

    rgb_p = reshape(im2double(frame_img),hw_size,3);
    rgb_sp = zeros(num_sp,3);
    for sp_id = 1:num_sp
        find_list = (sp_map==sp_id);
        rgb_sp(sp_id,:) = mean(rgb_p(find_list,:),1);
    end
    sp_feat = colorspace('Lab<-', rgb_sp)';

    color_dist = dist(sp_feat);
    color_dist = color_dist / max(color_dist(:));
    color_aff = exp(-color_dist.^2/param_list.sigsqr);

    flow_feat = zeros(2,num_sp);
    for sp_id = 1:num_sp
        flow_feat(1,sp_id) = mean(flow.fvy(sp_map==sp_id));
        flow_feat(2,sp_id) = mean(flow.fvx(sp_map==sp_id));
    end
    flow_dist = dist(flow_feat);
    flow_dist = flow_dist / max(flow_dist(:));
    
    bwarp_mat = zeros(num_sp,prev_num_sp);
    for sp_id = 1:num_sp
        warp_hist = hist(prev_sp_map(bwarp_list(sp_map==sp_id)),1:prev_num_sp);
        bwarp_mat(sp_id,:) = warp_hist/sum(warp_hist);
    end
    if frame_id > 2
        bwarp_mat2 = zeros(num_sp,third_num_sp);
        for sp_id = 1:num_sp
            warp_hist = hist(third_sp_map(bwarp_list2(sp_map==sp_id)),1:third_num_sp);
            bwarp_mat2(sp_id,:) = warp_hist/sum(warp_hist);
        end
    end
    
    inter_color_dist = pdist2(sp_feat',prev_sp_feat');
    inter_color_dist = inter_color_dist / max(inter_color_dist(:));
    inter_color_aff = exp(-inter_color_dist.^2/param_list.sigsqr);
    inter_color_trans = inter_color_aff.*bwarp_mat;
    
    fore_infer_r = bwarp_mat*prev_result_vec;
    fore_infer_r = fore_infer_r / sum(fore_infer_r); % Foreground guidance distributions
    back_infer_r = bwarp_mat*(1-prev_result_vec);
    back_infer_r = back_infer_r / sum(back_infer_r); % Background guidance distributions
    if frame_id > 2
        fore_infer_r = 0.5*bwarp_mat*prev_result_vec + 0.5*bwarp_mat2*third_result_vec;
        fore_infer_r = fore_infer_r / sum(fore_infer_r); % Foreground guidance distributions
        back_infer_r = 0.5*bwarp_mat*(1-prev_result_vec) + 0.5*bwarp_mat2*(1-third_result_vec);
        back_infer_r = back_infer_r / sum(back_infer_r); % Background guidance distributions
    end

    back_color_prob = inter_color_trans*(1-prev_result_vec);
    back_color_prob = back_color_prob / sum(back_color_prob);
    
    % Initial probability estimation
    [fore_init_p, back_init_p] = initial_probability_estimation(color_dist, flow_dist, ...
        sp_map, four_edge_mat, back_color_prob, param_list);
    
    % Inference restart probability distributions
    fore_infer_r = fore_infer_r.*fore_init_p;
    fore_infer_r = fore_infer_r / sum(fore_infer_r);
    back_infer_r = back_infer_r.*back_init_p;
    back_infer_r = back_infer_r / sum(back_infer_r);

    % Build color affinity matrices
    conn_color_aff = color_aff.*conn_edge_mat;
    conn_color_aff = conn_color_aff.*(1-eye(num_sp)) + param_list.self_aff*eye(num_sp);
    four_color_aff = color_aff.*four_edge_mat;
    four_color_aff = four_color_aff.*(1-eye(num_sp)) + param_list.self_aff*eye(num_sp);

    fore_color_trans = conn_color_aff*diag(1./(sum(conn_color_aff)+eps));
    back_color_trans = four_color_aff*diag(1./(sum(four_color_aff)+eps));

    % Simulation of multiple random walkers
    fore_p = fore_init_p;
    fore_p = fore_p / sum(fore_p);
    back_p = back_init_p;
    back_p = back_p / sum(back_p);
    
    opt_fore_p = fore_p;
    opt_back_p = back_p;
    prev_fore_r = fore_p;
    prev_back_r = back_p;
        
    iter_cnt = 0;
    while 1
        fore_inter_r = (conn_color_aff*(fore_p/max(fore_p)))./(conn_color_aff*(back_p/max(back_p))+conn_color_aff*(fore_p/max(fore_p)));
        fore_inter_r = fore_inter_r / sum(fore_inter_r);
        back_inter_r = (conn_color_aff*(back_p/max(back_p)))./(conn_color_aff*(fore_p/max(fore_p))+conn_color_aff*(back_p/max(back_p)));
        back_inter_r = back_inter_r / sum(back_inter_r);
        
        fore_inter_r = (1-(param_list.cool_factor^iter_cnt))*prev_fore_r + (param_list.cool_factor^iter_cnt)*fore_inter_r;
        back_inter_r = (1-(param_list.cool_factor^iter_cnt))*prev_back_r + (param_list.cool_factor^iter_cnt)*back_inter_r;
        
        fore_restart = (1-param_list.delta)*fore_inter_r + param_list.delta*fore_infer_r;
        back_restart = (1-param_list.delta)*back_inter_r + param_list.delta*back_infer_r;
        
        fore_p = param_list.epsilon*((eye(num_sp)-(1-param_list.epsilon)*fore_color_trans)\fore_restart);
        back_p = param_list.epsilon*((eye(num_sp)-(1-param_list.epsilon)*back_color_trans)\back_restart);
        p_diff = norm(fore_p-opt_fore_p) + norm(back_p-opt_back_p);
        if p_diff < 1e-4
            break;
        end
        opt_fore_p = fore_p;
        opt_back_p = back_p;
        prev_fore_r = fore_inter_r;
        prev_back_r = back_inter_r;
        iter_cnt = iter_cnt + 1;
    end

    % Foreground / background determination
    p_f = 1/max(opt_fore_p);
    p_b = 1/max(opt_back_p);
    fore_prior = p_f / (p_f + p_b);
    back_prior = p_b / (p_f + p_b);
    fore_posterior = (opt_fore_p*fore_prior)./(opt_fore_p*fore_prior+opt_back_p*back_prior);
    back_posterior = (opt_back_p*back_prior)./(opt_fore_p*fore_prior+opt_back_p*back_prior);
    result_vec = fore_posterior > back_posterior;
    
    result_map = zeros(h_size,w_size);
    for sp_id = 1:num_sp
        result_map(sp_map==sp_id) = result_vec(sp_id);
    end
    
    % Single object selection
    cc = bwconncomp(result_map);
    if cc.NumObjects > 1
        c_num = [];
        for c_id = 1:cc.NumObjects
            c_num = [c_num; length(cc.PixelIdxList{c_id})];
        end
        [~, c_maxi] = max(c_num);
        ccresult_map = zeros(h_size,w_size);
        ccresult_map(cc.PixelIdxList{c_maxi}) = 1;
    else
        ccresult_map = result_map;
    end
    result_vec = zeros(num_sp,1);
    for sp_id = 1:num_sp
        result_vec(sp_id) = mean(ccresult_map(sp_map==sp_id));
    end
    result_vec = double(result_vec > 0.5);
    ccresult_map = ccresult_map/max(ccresult_map(:));
    imwrite(ccresult_map,fullfile(result_path,sprintf('%04d.png',frame_id)));
    segment_track{frame_id} = ccresult_map;

    % Save previous results
    third_sp_map = prev_sp_map;
    third_num_sp = prev_num_sp;
    third_result_vec = prev_result_vec;
    third_frame_img = prev_frame_img;

    prev_sp_map = sp_map;
    prev_num_sp = num_sp;
    prev_sp_feat = sp_feat;
    prev_frame_img = frame_img;
    prev_result_vec = result_vec;

end

end