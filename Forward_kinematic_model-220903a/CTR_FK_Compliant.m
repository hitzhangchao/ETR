%% 3管CTR正运动学（典型CTR） - 基于扭转柔性模型
% 输入各管关节变量，输出各管末端位姿和整体形状
% by Zhang Chao
% Date：2022/9/2

% 首先运行ETR_params.m，加载物理参数

% 全局变量
global theta_11 theta_12 theta_13

%% 运动学输入
arm_1_in = [0,40E-3,0,65E-3,0,85E-3];   % arm 1 运动学输入 - 1级rotation和translation；2级rotation和translation；3级rotation和translation

% arm 1的关节变量
theta_11 = deg2rad(arm_1_in(1));        % sheath级tube的rotation
transl_11 = arm_1_in(2);
theta_12 = deg2rad(arm_1_in(3));
transl_12 = arm_1_in(4);
theta_13 = deg2rad(arm_1_in(5));
transl_13 = arm_1_in(6);


%% CTR分段（segmenting）





