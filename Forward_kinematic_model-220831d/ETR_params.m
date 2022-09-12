%% 3管ETR物理参数 - 首先运行该程序
% by Zhang Chao
% Date：2022/7/23

clear;clc;close all;

% 全局变量
global poisson_rate uy_star d

% 材料参数
E = 60*10^9;                % NiTi材料弹性模量，Pa
poisson_rate = 0.3;         % NiTi材料泊松比
G = E/(2*(1+poisson_rate)); % NiTi材料剪切模量，Pa

% sheath级物理参数(3根管是相同的)
phi_1 = 0;                  % 三管离心分布角分别为0°,120°,240°
phi_2 = 2/3*pi;
phi_3 = 4/3*pi;
OD_1 = 2.7E-3;              % 最外管tube外径
ID_1 = 2.0E-3;              % 最外管tube外径
r_1 = 130E-3;               % tube预弯曲半径
s_1 = 50E-3;                % tube伸出弧长
d = OD_1/sqrt(3);           % 离心偏距
uy_star = 1/r_1;            % tube预弯曲曲率
u_star = [0;uy_star;0];     % tube初始曲率矢量
I_1 = pi*(OD_1^4-ID_1^4)/64;% tube惯性矩
J_1 = pi*(OD_1^4-ID_1^4)/32;% tube极惯性矩
kb_1 = E*I_1;
kt_1 = G*J_1;
K_1 = [kb_1,0,0;0,kb_1,0;0,0,kt_1]; % sheah级tube刚度矩阵

% arm 1 其余物理参数
OD_12 = 1.80E-3;            % tube 12 外径(arm1,stage2)
ID_12 = 1.46E-3;            % tube 12 内径
length_12 = 25E-3;          % tube 12 弯曲段弧长（最大伸出长度）
r_12 = 80E-3;               % tube 12 预弯曲半径
u_12_star = [0;1/r_12;0];   % tube 12 初始曲率矢量
I_12 = pi*(OD_12^4-ID_12^4)/64;% tube 12 惯性矩
J_12 = pi*(OD_12^4-ID_12^4)/32;% tube 12 极惯性矩
kb_12 = E*I_12;
kt_12 = G*J_12;
K_12 = [kb_12,0,0;0,kb_12,0;0,0,kt_12]; % tube 12 刚度矩阵

OD_13 = 1.26E-3;            % tube 13 外径(arm1,stage3)
ID_13 = 1.10E-3;            % tube 13 内径
length_13 = 25E-3;          % tube 13 弯曲段弧长（最大伸出长度）
r_13 = 80E-3;               % tube 13 预弯曲半径
u_13_star = [0;1/r_13;0];   % tube 13 初始曲率矢量
I_13 = pi*(OD_13^4-ID_13^4)/64;% tube 13 惯性矩
J_13 = pi*(OD_13^4-ID_13^4)/32;% tube 13 极惯性矩
kb_13 = E*I_13;
kt_13 = G*J_13;
K_13 = [kb_13,0,0;0,kb_13,0;0,0,kt_13]; % tube 13 刚度矩阵

% arm 2 其余物理参数
OD_22 = 1.80E-3;            % tube 22 外径(arm2,stage2)
ID_22 = 1.46E-3;            % tube 22 内径
length_22 = 25E-3;          % tube 22 弯曲段弧长（最大伸出长度）
r_22 = 80E-3;               % tube 22 预弯曲半径
u_22_star = [0;1/r_22;0];   % tube 22 初始曲率矢量
I_22 = pi*(OD_22^4-ID_22^4)/64;% tube 22 惯性矩
J_22 = pi*(OD_22^4-ID_22^4)/32;% tube 22 极惯性矩
kb_22 = E*I_22;
kt_22 = G*J_22;
K_22 = [kb_22,0,0;0,kb_22,0;0,0,kt_22]; % tube 22 刚度矩阵

OD_23 = 1.26E-3;            % tube 23 外径(arm2,stage3)
ID_23 = 1.10E-3;            % tube 23 内径
length_23 = 25E-3;          % tube 23 弯曲段弧长（最大伸出长度）
r_23 = 80E-3;               % tube 23 预弯曲半径
u_23_star = [0;1/r_23;0];   % tube 23 初始曲率矢量
I_23 = pi*(OD_23^4-ID_23^4)/64;% tube 23 惯性矩
J_23 = pi*(OD_23^4-ID_23^4)/32;% tube 23 极惯性矩
kb_23 = E*I_23;
kt_23 = G*J_23;
K_23 = [kb_23,0,0;0,kb_23,0;0,0,kt_23]; % tube 23 刚度矩阵

% arm 3 其余物理参数
OD_32 = 1.80E-3;            % tube 32 外径(arm3,stage2)
ID_32 = 1.46E-3;            % tube 32 内径
length_32 = 25E-3;          % tube 32 弯曲段弧长（最大伸出长度）
r_32 = 80E-3;               % tube 32 预弯曲半径
u_32_star = [0;1/r_32;0];   % tube 32 初始曲率矢量
I_32 = pi*(OD_32^4-ID_32^4)/64;% tube 32 惯性矩
J_32 = pi*(OD_32^4-ID_32^4)/32;% tube 32 极惯性矩
kb_32 = E*I_32;
kt_32 = G*J_32;
K_32 = [kb_32,0,0;0,kb_32,0;0,0,kt_32]; % tube 32 刚度矩阵

OD_33 = 1.26E-3;            % tube 33 外径(arm3,stage3)
ID_33 = 1.10E-3;            % tube 33 内径
length_33 = 25E-3;          % tube 33 弯曲段弧长（最大伸出长度）
r_33 = 80E-3;               % tube 33 预弯曲半径
u_33_star = [0;1/r_33;0];   % tube 33 初始曲率矢量
I_33 = pi*(OD_33^4-ID_33^4)/64;% tube 33 惯性矩
J_33 = pi*(OD_33^4-ID_33^4)/32;% tube 33 极惯性矩
kb_33 = E*I_33;
kt_33 = G*J_33;
K_33 = [kb_33,0,0;0,kb_33,0;0,0,kt_33]; % tube 33 刚度矩阵
