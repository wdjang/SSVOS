function param_list = set_param

% SLIC parameters
param_list.slic.area = 1000;
param_list.slic.edge = 20;

% K-ring connectivity
param_list.fg_k = 1;
param_list.bg_k = 4;

% Self affinity
param_list.self_aff = 1;

% Sigma square
param_list.sigsqr = 0.01;

% Conversion parameter from background to foreground
param_list.conv = 0.0002;

% Restart probability
param_list.epsilon = 0.4;

% Restart blending ratio
param_list.delta = 0.7;

% Cooling factor
param_list.cool_factor = 0.995;

% Color, motion, and guidance weights
param_list.lambda_c = 1;
param_list.lambda_m = 2;
param_list.lambda_g = 1;

% Optical flow parameters
param_list.flow.alpha = 0.012;
param_list.flow.ratio = 0.75;
param_list.flow.minWidth = 20;
param_list.flow.nOuterFPIterations = 7;
param_list.flow.nInnerFPIterations = 1;
param_list.flow.nSORIterations = 30;

end