function [fore_init_p, back_init_p] = initial_probability_estimation(color_dist, flow_dist, ...
    sp_map, edge_mat, back_guide, param_list)

num_sp = max(sp_map(:));

bg_color_dist = color_dist;
bg_color_dist(logical(1-edge_mat)) = Inf;
bg_color_aff = exp(-bg_color_dist.^2/param_list.sigsqr);
bgcolor_trans = bg_color_aff * diag(1./(sum(bg_color_aff)+eps));

bg_flow_dist = flow_dist;
bg_flow_dist(logical(1-edge_mat)) = Inf;
bg_motion_aff = exp(-bg_flow_dist.^2/param_list.sigsqr);
bg_motion_trans = bg_motion_aff * diag(1./(sum(bg_motion_aff)+eps));

fg_color_dist = color_dist;
fg_color_dist(logical(1-edge_mat)) = Inf;
fg_color_aff = exp(-fg_color_dist.^2/param_list.sigsqr);
fg_color_trans = fg_color_aff * diag(1./(sum(fg_color_aff)+eps));

fg_flow_dist = flow_dist;
fg_flow_dist(logical(1-edge_mat)) = Inf;
fg_motion_aff = exp(-fg_flow_dist.^2/param_list.sigsqr);
fg_motion_trans = fg_motion_aff * diag(1./(sum(fg_motion_aff)+eps));

cvx_begin quiet
    cvx_solver gurobi;
    variable back_init_p(num_sp,1);
    minimize( param_list.lambda_c*norm(bgcolor_trans*back_init_p-back_init_p) ...
        + param_list.lambda_m*norm(bg_motion_trans*back_init_p-back_init_p) ...
        + param_list.lambda_g*norm(back_init_p - back_guide) );
    subject to
        0 <= back_init_p <= 1;
        sum(back_init_p) == 1;
cvx_end

fore_guide = exp(-back_init_p/param_list.conv);
fore_guide = fore_guide / sum(fore_guide);

cvx_begin quiet
    cvx_solver gurobi;
    variable fore_init_p(num_sp,1);
    minimize( param_list.lambda_c*norm(fg_color_trans*fore_init_p-fore_init_p) ...
        + param_list.lambda_m*norm(fg_motion_trans*fore_init_p-fore_init_p) ...
        + param_list.lambda_g*norm(fore_init_p - fore_guide) );
    subject to
        0 <= fore_init_p <= 1;
        sum(fore_init_p) == 1;
cvx_end

fore_init_p = fore_init_p / sum(fore_init_p);
back_init_p = back_init_p / sum(back_init_p);

end