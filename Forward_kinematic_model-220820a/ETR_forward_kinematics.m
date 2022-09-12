%% 3管ETR正运动学 - 基于扭转刚性模型
% 输入各管关节变量，输出各管末端位姿和整体形状
% by Zhang Chao
% Date：2022/7/12

% 首先运行ETR_params.m，加载物理参数

%% 运动学输入
arm_1_in = [0,0,25E-3,0,45E-3];         % arm 1 运动学输入 - 1级rotation；2级rotation和translation；3级rotation和translation
arm_2_in = [90,120,25E-3,120,45E-3];        % arm 2 运动学输入 - 1级rotation；2级rotation和translation；3级rotation和translation
arm_3_in = [110,-120,25E-3,-120,45E-3];       % arm 3 运动学输入 - 1级rotation；2级rotation和translation；3级rotation和translation

% arm 1的关节变量
theta_11 = deg2rad(arm_1_in(1));        % sheath级tube的rotation
theta_12 = deg2rad(arm_1_in(2));
transl_12 = arm_1_in(3);
theta_13 = deg2rad(arm_1_in(4));
transl_13 = arm_1_in(5);

% arm 2的关节变量
theta_21 = deg2rad(arm_2_in(1));
theta_22 = deg2rad(arm_2_in(2));
transl_22 = arm_2_in(3);
theta_23 = deg2rad(arm_2_in(4));
transl_23 = arm_2_in(5);

% arm 3的关节变量
theta_31 = deg2rad(arm_3_in(1));
theta_32 = deg2rad(arm_3_in(2));
transl_32 = arm_3_in(3);
theta_33 = deg2rad(arm_3_in(4));
transl_33 = arm_3_in(5);


%% sheath级运动学 - 初始输入
% 各tube初始输入示意 - 假设没有弹性相互作用时
p11_0 = [d*cos(phi_1);d*sin(phi_1);0];                   % 起始点
p21_0 = [d*cos(phi_2);d*sin(phi_2);0];
p31_0 = [d*cos(phi_3);d*sin(phi_3);0];

% 绘制形状
cal_cnt = 50;                   % 计算点数50（对应s_1 = 50mm）
p11_ini = zeros(3,cal_cnt);     % 数组预定义
p21_ini = zeros(3,cal_cnt);
p31_ini = zeros(3,cal_cnt);
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_1;
    % Dupont 2006年论文中矩阵T - 位置元素
    a14 = u_star(2)*(1-cos(s_temp*norm(u_star)))/(norm(u_star))^2;
    a24 = -u_star(1)*(1-cos(s_temp*norm(u_star)))/(norm(u_star))^2;
    a34 = sin(s_temp*norm(u_star))/norm(u_star);
    p11_ini(:,i) = p11_0+rotz(theta_11)*[a14;a24;a34];
    p21_ini(:,i) = p21_0+rotz(theta_21)*[a14;a24;a34];
    p31_ini(:,i) = p31_0+rotz(theta_31)*[a14;a24;a34];
end
% plot3(p11_ini(1,:),p11_ini(2,:),p11_ini(3,:),'--r','LineWidth',1);
% hold on
% plot3(p21_ini(1,:),p21_ini(2,:),p21_ini(3,:),'--g','LineWidth',1);
% hold on
% plot3(p31_ini(1,:),p31_ini(2,:),p31_ini(3,:),'--b','LineWidth',1);
% hold on


%% sheath级运动学 - 外鞘
syms u_x u_y
a_1 = sin(theta_11)+sin(theta_21)+sin(theta_31);
a_2 = cos(theta_11)+cos(theta_21)+cos(theta_31);
u_x = -a_1/a_2*u_y;

f1 = u_y*(1-d*cos(phi_2)*u_y+d*sin(phi_2)*u_x)*(1-d*cos(phi_3)*u_y+d*sin(phi_3)*u_x) + ...
     u_y*(1-d*cos(phi_1)*u_y+d*sin(phi_1)*u_x)*(1-d*cos(phi_3)*u_y+d*sin(phi_3)*u_x) + ...
     u_y*(1-d*cos(phi_1)*u_y+d*sin(phi_1)*u_x)*(1-d*cos(phi_2)*u_y+d*sin(phi_2)*u_x) - ...
     a_2*uy_star*(1-d*cos(phi_1)*u_y+d*sin(phi_1)*u_x)*(1-d*cos(phi_2)*u_y+d*sin(phi_2)*u_x)*(1-d*cos(phi_3)*u_y+d*sin(phi_3)*u_x);
f2 = vpa(simplify(f1));

% 求解上述方程（可化成一元三次方程，多解，需要选择符合实际情况的解）
Res = solve(f2,u_y)
u_y = Res(2);
u_x = eval(u_x);
u = [u_x;u_y;0]

% 绘制图形
cal_cnt = 50;                   % 计算点数50（对应s_1 = 50mm）
p01 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_1;
    a14 = u(2)*(1-cos(s_temp*norm(u)))/(norm(u))^2;
    a24 = -u(1)*(1-cos(s_temp*norm(u)))/(norm(u))^2;
    a34 = sin(s_temp*norm(u))/norm(u);
    p01(:,i) = [a14;a24;a34];
end
plot3(p01(1,:),p01(2,:),p01(3,:),'m','LineWidth',5);
hold on
T0_w1_g = T_u_s(u,s_1);          % sheath末端相对于global{0}系
hold on


%% sheath运动学 - 各管
T1_w0_g = [eye(3),p11_0;[0 0 0 1]];        % arm1 s=0起始点w0系相对于global{0}系
T2_w0_g = [eye(3),p21_0;[0 0 0 1]];        % arm2 s=0起始点w0系相对于global{0}系
T3_w0_g = [eye(3),p31_0;[0 0 0 1]];        % arm3 s=0起始点w0系相对于global{0}系

u11_w0 = u/(1-d*cos(phi_1)*u_y+d*sin(phi_1)*u_x);
p11 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_1;
    T_temp = T_u_s(u11_w0,s_temp);
    T11_w1_g = T1_w0_g*T_temp;
    p11(:,i) = transl(T11_w1_g);
end
plot3(p11(1,:),p11(2,:),p11(3,:),'-r','LineWidth',5);
hold on
T11_w1_w0 = T_u_s(u11_w0,s_1);          % sheath级tube 1末端w1相对于w0系
T11_w1_g = T1_w0_g*T11_w1_w0;           % sheath级tube 1末端w1相对于global{0}系

u21_w0 = u/(1-d*cos(phi_2)*u_y+d*sin(phi_2)*u_x);
p21 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_1;
    T_temp = T_u_s(u21_w0,s_temp);
    T21_w1_g = T2_w0_g*T_temp;
    p21(:,i) = transl(T21_w1_g);
end
plot3(p21(1,:),p21(2,:),p21(3,:),'-g','LineWidth',5);
hold on
T21_w1_w0 = T_u_s(u21_w0,s_1);           % sheath级tube 2末端相对于local {s=0}系
T21_w0_g = T2_w0_g*T21_w1_w0;           % sheath级tube 2末端相对于global{0}系

u31_w0 = u/(1-d*cos(phi_3)*u_y+d*sin(phi_3)*u_x);
p31 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_1;
    T_temp = T_u_s(u21_w0,s_temp);
    T31_w1_g = T3_w0_g*T_temp;
    p31(:,i) = transl(T31_w1_g);
end
plot3(p31(1,:),p31(2,:),p31(3,:),'-b','LineWidth',5);
hold on
T31_w1_w0 = T_u_s(u31_w0,s_1);           % sheath级tube 3末端相对于local {s=0}系
T31_w0_g = T3_w0_g*T31_w1_w0;           % sheath级tube 3末端相对于global{0}系

%% arm 1 运动学（基于扭转刚性模型）
% 由之前的计算，我们已经得到了sheath级整体的形状以及管间位姿
% T0_w1_g —— sheath末端w1系相对于global{0}系
% T11_w1_g —— sheath级tube 1末端w1系相对于global{0}系
% T21_w1_g —— sheath级tube 2末端w1系相对于global{0}系
% T31_w1_g —— sheath级tube 3末端w1系相对于global{0}系

% 按照tube重合情况分段(如果要更完善鲁棒，可增加判断条件)
s_14 = transl_13 - transl_12;
s_13 = length_13 - s_14;
s_12 = transl_12 - s_13;

% s_12段曲率
u_12_star_w1 = rotz(theta_12)*u_12_star;        % s12段初始曲率（相对于local initial frame）
u_12_w1 = (K_12+K_13)\(K_12*u_12_star_w1);      % s12段计算曲率（矩阵右除，即inv(K_12+K_13)*...）
% 绘制形状
cal_cnt = 25;                   % 计算点数50（对应s_1 = 50mm）
p12 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_12;
    T_temp = T_u_s(u_12_w1,s_temp);         % s12段各点相对于w1系
    T12_w2_g = T11_w1_g*T_temp;
    p12(:,i) = transl(T12_w2_g);
end
plot3(p12(1,:),p12(2,:),p12(3,:),'-','LineWidth',3.5);
hold on
T12_w2_w1 = T_u_s(u_12_w1,s_12);            % s12段末端相对于w1系
T12_w2_g = T11_w1_g*T12_w2_w1;              % s12段末端相对于global{0}系

% s_13段曲率
u_12_star_w2 = rotz(theta_12)*u_12_star;        % s12段初始曲率（相对于local initial frame）
u_13_star_w2 = rotz(theta_13)*u_13_star;        % s13段初始曲率（相对于local initial frame）
u_13_w2 = (K_12+K_13)\(K_12*u_12_star_w2+K_13*u_13_star_w2);      % s13段计算曲率（矩阵右除，即inv(K_12+K_13)*...）
% 绘制形状
cal_cnt = 25;                   % 计算点数50（对应s_1 = 50mm）
p13 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_13;
    T_temp = T_u_s(u_13_w2,s_temp);         % s13段各点相对于w2系
    T13_w3_g = T12_w2_g*T_temp;
    p13(:,i) = transl(T13_w3_g);
end
plot3(p13(1,:),p13(2,:),p13(3,:),'-','LineWidth',3.5);
hold on
T13_w3_w2 = T_u_s(u_13_w2,s_13);            % s13段末端w3相对于w2系
T13_w3_g = T12_w2_g*T13_w3_w2;              % s13段末端w3相对于global{0}系

% s_14段曲率
u_13_star_w3 = rotz(theta_13)*u_13_star;        % s14段初始曲率（相对于local initial frame）
u_14_w3 = K_13\(K_13*u_13_star_w3);             % s14段计算曲率（矩阵右除，即inv(K_12+K_13)*...）
% 绘制形状
cal_cnt = 25;                   % 计算点数50（对应s_1 = 50mm）
p14 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_14;
    T_temp = T_u_s(u_14_w3,s_temp);         % s13段各点相对于w2系
    T14_w4_g = T13_w3_g*T_temp;
    p14(:,i) = transl(T14_w4_g);
end
plot3(p14(1,:),p14(2,:),p14(3,:),'-','LineWidth',2);
hold on
T14_w4_w3 = T_u_s(u_14_w3,s_14);            % s14段末端相对于w3系
T14_w4_g = T13_w3_g*T14_w4_w3;              % s14段末端相对于global{0}系
% 最后一个tube是要考虑rotation的，需要修改添加


%% arm 2 运动学（基于扭转刚性模型）
% 由之前的计算，我们已经得到了sheath级整体的形状以及管间位姿
% T0_w1_g —— sheath末端w1系相对于global{0}系
% T11_w1_g —— sheath级tube 1末端w1系相对于global{0}系
% T21_w1_g —— sheath级tube 2末端w1系相对于global{0}系
% T31_w1_g —— sheath级tube 3末端w1系相对于global{0}系

% 按照tube重合情况分段(如果要更完善鲁棒，可增加判断条件)
s_24 = transl_23 - transl_22;
s_23 = length_23 - s_24;
s_22 = transl_22 - s_23;

% s_22段曲率
u_22_star_w1 = rotz(theta_22)*u_22_star;        % s22段初始曲率（相对于local initial frame）
u_22_w1 = (K_22+K_23)\(K_22*u_22_star_w1);      % s22段计算曲率（矩阵右除，即inv(K_22+K_23)*...）
% 绘制形状
cal_cnt = 25;                   % 计算点数50（对应s_1 = 50mm）
p22 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_22;
    T_temp = T_u_s(u_22_w1,s_temp);         % s22段各点相对于w1系
    T22_w2_g = T21_w1_g*T_temp;
    p22(:,i) = transl(T22_w2_g);
end
plot3(p22(1,:),p22(2,:),p22(3,:),'-','LineWidth',3.5);
hold on
T22_w2_w1 = T_u_s(u_22_w1,s_22);            % s22段末端相对于w1系
T22_w2_g = T21_w1_g*T22_w2_w1;              % s22段末端相对于global{0}系

% s_23段曲率
u_22_star_w2 = rotz(theta_22)*u_22_star;        % s22段初始曲率（相对于local initial frame）
u_23_star_w2 = rotz(theta_23)*u_23_star;        % s23段初始曲率（相对于local initial frame）
u_23_w2 = (K_22+K_23)\(K_22*u_22_star_w2+K_23*u_23_star_w2);      % s23段计算曲率（矩阵右除，即inv(K_12+K_13)*...）
% 绘制形状
cal_cnt = 25;                   % 计算点数50（对应s_1 = 50mm）
p23 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_23;
    T_temp = T_u_s(u_23_w2,s_temp);         % s23段各点相对于w2系
    T23_w3_g = T22_w2_g*T_temp;
    p23(:,i) = transl(T23_w3_g);
end
plot3(p23(1,:),p23(2,:),p23(3,:),'-','LineWidth',3.5);
hold on
T23_w3_w2 = T_u_s(u_23_w2,s_23);            % s23段末端w3相对于w2系
T23_w3_g = T22_w2_g*T23_w3_w2;              % s23段末端w3相对于global{0}系

% s_24段曲率
u_23_star_w3 = rotz(theta_23)*u_23_star;        % s24段初始曲率（相对于local initial frame）
u_24_w3 = K_23\(K_23*u_23_star_w3);             % s24段计算曲率（矩阵右除，即inv(K_12+K_13)*...）
% 绘制形状
cal_cnt = 25;                   % 计算点数50（对应s_1 = 50mm）
p24 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_24;
    T_temp = T_u_s(u_24_w3,s_temp);         % s24段各点相对于w3系
    T24_w4_g = T23_w3_g*T_temp;
    p24(:,i) = transl(T24_w4_g);
end
plot3(p24(1,:),p24(2,:),p24(3,:),'-','LineWidth',2);
hold on
T24_w4_w3 = T_u_s(u_24_w3,s_24);            % s14段末端相对于w3系
T24_w4_g = T23_w3_g*T24_w4_w3;              % s14段末端相对于global{0}系
% 最后一个tube是要考虑rotation的，需要修改添加


%% arm 3 运动学（基于扭转刚性模型）
% 由之前的计算，我们已经得到了sheath级整体的形状以及管间位姿
% T0_w1_g —— sheath末端w1系相对于global{0}系
% T11_w1_g —— sheath级tube 1末端w1系相对于global{0}系
% T21_w1_g —— sheath级tube 2末端w1系相对于global{0}系
% T31_w1_g —— sheath级tube 3末端w1系相对于global{0}系

% 按照tube重合情况分段(如果要更完善鲁棒，可增加判断条件)
s_34 = transl_33 - transl_32;
s_33 = length_33 - s_34;
s_32 = transl_32 - s_33;

% s_32段曲率
u_32_star_w1 = rotz(theta_32)*u_32_star;        % s32段初始曲率（相对于local initial frame）
u_32_w1 = (K_32+K_33)\(K_32*u_32_star_w1);      % s32段计算曲率（矩阵右除，即inv(K_12+K_13)*...）
% 绘制形状
cal_cnt = 25;                   % 计算点数50（对应s_1 = 50mm）
p32 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_32;
    T_temp = T_u_s(u_32_w1,s_temp);         % s32段各点相对于w1系
    T32_w2_g = T31_w1_g*T_temp;
    p32(:,i) = transl(T32_w2_g);
end
plot3(p32(1,:),p32(2,:),p32(3,:),'-','LineWidth',3.5);
hold on
T32_w2_w1 = T_u_s(u_32_w1,s_32);            % s32段末端相对于w1系
T32_w2_g = T31_w1_g*T32_w2_w1;              % s32段末端相对于global{0}系

% s_33段曲率
u_32_star_w2 = rotz(theta_32)*u_32_star;        % s32段初始曲率（相对于local initial frame）
u_33_star_w2 = rotz(theta_33)*u_33_star;        % s13段初始曲率（相对于local initial frame）
u_33_w2 = (K_32+K_33)\(K_32*u_32_star_w2+K_33*u_33_star_w2);      % s33段计算曲率（矩阵右除，即inv(K_12+K_13)*...）
% 绘制形状
cal_cnt = 25;                   % 计算点数50（对应s_1 = 50mm）
p33 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_33;
    T_temp = T_u_s(u_33_w2,s_temp);         % s33段各点相对于w2系
    T33_w3_g = T32_w2_g*T_temp;
    p33(:,i) = transl(T33_w3_g);
end
plot3(p33(1,:),p33(2,:),p33(3,:),'-','LineWidth',3.5);
hold on
T33_w3_w2 = T_u_s(u_33_w2,s_33);            % s33段末端w3相对于w2系
T33_w3_g = T32_w2_g*T33_w3_w2;              % s33段末端w3相对于global{0}系

% s_34段曲率
u_33_star_w3 = rotz(theta_33)*u_33_star;        % s34段初始曲率（相对于local initial frame）
u_34_w3 = K_33\(K_33*u_33_star_w3);             % s34段计算曲率（矩阵右除，即inv(K_12+K_13)*...）
% 绘制形状
cal_cnt = 25;                   % 计算点数50（对应s_1 = 50mm）
p34 = zeros(3,cal_cnt);         % 数组预定义
for i=1:1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_34;
    T_temp = T_u_s(u_34_w3,s_temp);         % s34段各点相对于w3系
    T34_w4_g = T33_w3_g*T_temp;
    p34(:,i) = transl(T34_w4_g);
end
plot3(p34(1,:),p34(2,:),p34(3,:),'-','LineWidth',2);
hold on
T34_w4_w3 = T_u_s(u_34_w3,s_34);            % s34段末端相对于w3系
T34_w4_g = T33_w3_g*T34_w4_w3;              % s34段末端相对于global{0}系
% 最后一个tube是要考虑rotation的，需要修改添加


%% 绘图标注
grid on;
axis equal;
% xlim([-0.020 0.020]);
% ylim([-0.020 0.020]);
% zlim([0 0.070]);
xlabel('x/m')
ylabel('y/m')
zlabel('z/m')
%legend(H([1,3]),'f1','f3');


