%% 3管ETR正运动学 - 基于扭转柔性模型
% 输入各管关节变量，输出各管末端位姿和整体形状
% by Zhang Chao
% Date：2022/8/22

% 首先运行ETR_params.m，加载物理参数

% 全局变量
global theta_11 theta_21 theta_31

%% 运动学输入
arm_1_in = [0,0,25E-3,0,45E-3];             % arm 1 运动学输入 - 1级rotation；2级rotation和translation；3级rotation和translation
arm_2_in = [90,120,25E-3,120,45E-3];        % arm 2 运动学输入 - 1级rotation；2级rotation和translation；3级rotation和translation
arm_3_in = [110,-120,25E-3,-120,45E-3];     % arm 3 运动学输入 - 1级rotation；2级rotation和translation；3级rotation和translation

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
for i=1:cal_cnt
    s_temp = (i-1)/(cal_cnt-1)*s_1;
    % Dupont 2006年论文中矩阵T - 位置元素
    a14 = u_star(2)*(1-cos(s_temp*norm(u_star)))/(norm(u_star))^2;
    a24 = -u_star(1)*(1-cos(s_temp*norm(u_star)))/(norm(u_star))^2;
    a34 = sin(s_temp*norm(u_star))/norm(u_star);
    p11_ini(:,i) = p11_0+rotz(theta_11)*[a14;a24;a34];
    p21_ini(:,i) = p21_0+rotz(theta_21)*[a14;a24;a34];
    p31_ini(:,i) = p31_0+rotz(theta_31)*[a14;a24;a34];
end
plot3(p11_ini(1,:),p11_ini(2,:),p11_ini(3,:),'--r','LineWidth',1);
hold on
plot3(p21_ini(1,:),p21_ini(2,:),p21_ini(3,:),'--g','LineWidth',1);
hold on
plot3(p31_ini(1,:),p31_ini(2,:),p31_ini(3,:),'--b','LineWidth',1);
hold on


%% sheath级运动学 - 外鞘

%-- 使用bvp5c求解微分方程 --%
cal_cnt = 11;                            % 计算201个点，这里即1/4mm (调整求解点数，程序运行时间会大幅度缩短)
s_mesh = linspace(0,s_1,cal_cnt);          % s各离散点弧长
solinit = bvpinit(s_mesh, @ETR_guess);
sol = bvp5c(@ETR_bvpfun, @ETR_bcfun, solinit);  % 这里可求解得到u_x(s)和u_y(s)

%plot(sol.x, sol.y)
%legend('theta_{11}','theta_{21}','theta_{31}','dot(theta_{11})','dot(theta_{21})','dot(theta_{31})');

%-- 计算曲率u(s)，并存储下来 --%
u_s = zeros(3,cal_cnt);
for i=1:cal_cnt
    %首先计算ψ和u_abs
    C = cos(sol.y(1,i))+cos(sol.y(2,i))+cos(sol.y(3,i));
    S = -sin(sol.y(1,i))-sin(sol.y(2,i))-sin(sol.y(3,i));
    sin_psi = C/sqrt(S^2+C^2);
    cos_psi = S/sqrt(S^2+C^2);
    
    a0 = 4*uy_star*sqrt(S^2+C^2);
    a1 = -12;
    a2 = -3/4*d^2*a0;
    a3 = 3*d^2+d^3*(-4*sin_psi^3+3*sin_psi)*a0/4;
    syms u_abs
    Res = vpasolve(a0+a1*u_abs+a2*u_abs^2+a3*u_abs^3,u_abs);
    u_abs = Res(2);
    
    u_s(1,i) = u_abs*cos_psi;   % u_x(s)
    u_s(2,i) = u_abs*sin_psi;   % u_y(s)
end

%-- 根据求解得到的曲率，计算绘制形状 --%
p_s = zeros(3,cal_cnt);
R_e3 = zeros(3,cal_cnt);
for i=1:cal_cnt
    s_temp = s_mesh(1:i);                   % s=0到当前弧长位置的所有点
    u_s_temp = u_s(:,1:i);                  % 各弧长s点对应曲率
    int_u_temp = trapz(s_temp,u_s_temp,2);  % 使用trapz做数值积分
    R_temp = expm(skew(int_u_temp));
    R_e3(:,i) = R_temp*[0;0;1];             % s=0到当前弧长位置点的中间值R*e3    
    p_s(:,i) = trapz(s_temp,R_e3(:,1:i),2); % 使用trapz做数值积分
end

plot3(p_s(1,:),p_s(2,:),p_s(3,:),'-m','LineWidth',5);
hold on

% sheath末端相对于global{0}系位姿
R0_w1_g = R_temp;
p_s_w1_g = p_s(:,end);
T0_w1_g = [R0_w1_g,p_s_w1_g;[0 0 0 1]];



%% sheath级运动学 - 各管
T1_w0_g = [eye(3),p11_0;[0 0 0 1]];        % arm1 s=0起始点w0系相对于global{0}系
T2_w0_g = [eye(3),p21_0;[0 0 0 1]];        % arm2 s=0起始点w0系相对于global{0}系
T3_w0_g = [eye(3),p31_0;[0 0 0 1]];        % arm3 s=0起始点w0系相对于global{0}系

s1_dot = 1+d*sin(phi_1)*u_s(1,:)-d*cos(phi_1)*u_s(2,:);
s2_dot = 1+d*sin(phi_2)*u_s(1,:)-d*cos(phi_2)*u_s(2,:);
s3_dot = 1+d*sin(phi_3)*u_s(1,:)-d*cos(phi_3)*u_s(2,:);

%-- 根据求解得到的曲率，计算绘制形状 --%
% 方法一：u_i是相对于tube截面贴体坐标系描述的，将u_i旋转回FS frame中描述，再进行计算
u1_F0 = (u_s+[0;0;1]*sol.y(4,:))./s1_dot;
u2_F0 = (u_s+[0;0;1]*sol.y(5,:))./s2_dot;
u3_F0 = (u_s+[0;0;1]*sol.y(6,:))./s3_dot;

% tube 11
p_11 = zeros(3,cal_cnt);
R1_e3 = zeros(3,cal_cnt);
for i=1:cal_cnt
    s_temp = s_mesh(1:i);                   % s=0到当前弧长位置的所有点
    u_temp = u1_F0(:,1:i);                  % 各弧长s点对应曲率
    int_u_temp = trapz(s_temp,u_temp,2);  % 使用trapz做数值积分
    R_temp = expm(skew(int_u_temp));
    R1_e3(:,i) = R_temp*[0;0;1];             % s=0到当前弧长位置点的中间值R*e3    
    p_11(:,i) = trapz(s_temp,R1_e3(:,1:i),2) + p11_0; % 使用trapz做数值积分
end
plot3(p_11(1,:),p_11(2,:),p_11(3,:),'-r','LineWidth',5);
hold on
% tube 11末端相对于global{0}系位姿
R0_w1_g = R_temp;
p_s_w1_g = p_s(:,end);
T0_w1_g = [R0_w1_g,p_s_w1_g;[0 0 0 1]];



% tube 21
p_21 = zeros(3,cal_cnt);
R2_e3 = zeros(3,cal_cnt);
for i=1:cal_cnt
    s_temp = s_mesh(1:i);                   % s=0到当前弧长位置的所有点
    u_temp = u2_F0(:,1:i);                  % 各弧长s点对应曲率
    int_u_temp = trapz(s_temp,u_temp,2);  % 使用trapz做数值积分
    R_temp = expm(skew(int_u_temp));
    R2_e3(:,i) = R_temp*[0;0;1];             % s=0到当前弧长位置点的中间值R*e3    
    p_21(:,i) = trapz(s_temp,R2_e3(:,1:i),2) + p21_0; % 使用trapz做数值积分
end
plot3(p_21(1,:),p_21(2,:),p_21(3,:),'-g','LineWidth',5);
hold on

% tube 31
p_31 = zeros(3,cal_cnt);
R3_e3 = zeros(3,cal_cnt);
for i=1:cal_cnt
    s_temp = s_mesh(1:i);                   % s=0到当前弧长位置的所有点
    u_temp = u3_F0(:,1:i);                  % 各弧长s点对应曲率
    int_u_temp = trapz(s_temp,u_temp,2);  % 使用trapz做数值积分
    R_temp = expm(skew(int_u_temp));
    R3_e3(:,i) = R_temp*[0;0;1];             % s=0到当前弧长位置点的中间值R*e3    
    p_31(:,i) = trapz(s_temp,R3_e3(:,1:i),2) + p31_0; % 使用trapz做数值积分
end
plot3(p_31(1,:),p_31(2,:),p_31(3,:),'-b','LineWidth',5);
hold on

% 方法二：根据公式Ri(si) =R(s)*Rz(theta_i(si))，求得Ri(si)；再根据pi(si)'=Ri*e3，积分求出pi(si)


%% arm 1 运动学（基于扭转柔性模型）
% 由之前的计算，我们已经得到了sheath级整体的形状以及管间位姿
% T0_w1_g —— sheath末端w1系相对于global{0}系
% T11_w1_g —— sheath级tube 1末端w1系相对于global{0}系
% T21_w1_g —— sheath级tube 2末端w1系相对于global{0}系
% T31_w1_g —— sheath级tube 3末端w1系相对于global{0}系

% 按照tube重合情况分段(如果要更完善鲁棒，可增加判断条件)
s_14 = transl_13 - transl_12;
s_13 = length_13 - s_14;
s_12 = transl_12 - s_13;








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


%% 函数定义 

% 微分方程
% y(1) = theta_1
% y(2) = theta_2
% y(3) = theta_3
% y(4) = theta_1'
% y(5) = theta_2'
% y(6) = theta_3'
function dydx = ETR_bvpfun(x,y)
    global poisson_rate uy_star d
    C = cos(y(1))+cos(y(2))+cos(y(3));
    S = -(sin(y(1))+sin(y(2))+sin(y(3)));
    C_dot = -(y(4)*sin(y(1))+y(5)*sin(y(2))+y(6)*sin(y(3)));
    S_dot = -(y(4)*cos(y(1))+y(5)*cos(y(2))+y(6)*cos(y(3)));
    % 用φ表示出中间参数ψ
    sin_psi = C/sqrt(S^2+C^2);
    cos_psi = S/sqrt(S^2+C^2);
    psi = atan(C/S);
    psi_dot = (S*C_dot-C*S_dot)/(S^2+C^2);

    % 求解一元三次方程，表示出u_abs
    a0 = 4*uy_star*sqrt(S^2+C^2);
    a1 = -12;
    a2 = -3/4*d^2*a0;
    %a3 = 3*d^2+d^3*sin(3*psi)*a0;
    a3 = 3*d^2+d^3*(-4*sin_psi^3+3*sin_psi)*a0/4;     %将sin(3ψ)写开
    syms u_abs
    Res = vpasolve(a0+a1*u_abs+a2*u_abs^2+a3*u_abs^3,u_abs);
    u_abs = Res(2);
    
    a0_dot = 4*uy_star*(S*S_dot+C*C_dot)/sqrt(S^2+C^2);
    %a1_dot = 0;
    a2_dot = -3/4*d^2*a0_dot;
    %a3_dot = d^3*sin(3*psi)*a0_dot/4 + 3/4*d^3*cos(3*psi)*a0_dot*psi_dot;
    a3_dot = d^3*(-4*sin_psi^3+3*sin_psi)*a0_dot/4 + 3/4*d^3*(4*cos_psi^3-3*cos_psi)*a0_dot*psi_dot;    %将sin(3ψ)和cos(3ψ)写开
    u_abs_dot = -(a0_dot+a2_dot*u_abs^2+a3_dot*u_abs^3)/(a1+2*a2*u_abs+3*a3*u_abs^2);
    
    dydx = [y(4);
            y(5);
            y(6);
            -(d*u_abs_dot*sin_psi+d*u_abs*cos_psi*psi_dot)/(1-d*u_abs*sin_psi)*y(4) + (1+poisson_rate)*uy_star*(1-d*u_abs*sin_psi)*u_abs*cos(psi-y(1));
            (d*u_abs_dot*(sqrt(3)/2*cos_psi+1/2*sin_psi)-d*u_abs*(-1/2*cos_psi+sqrt(3)/2*sin_psi)*psi_dot)/(1+d*u_abs*(sqrt(3)/2*cos_psi+1/2*sin_psi))*y(5) + (1+poisson_rate)*uy_star*(1+d*u_abs*(sqrt(3)/2*cos_psi+1/2*sin_psi))*u_abs*cos(psi-y(2));
            (d*u_abs_dot*(-sqrt(3)/2*cos_psi+1/2*sin_psi)+d*u_abs*(1/2*cos_psi+sqrt(3)/2*sin_psi)*psi_dot)/(1+d*u_abs*(-sqrt(3)/2*cos_psi+1/2*sin_psi))*y(6) + (1+poisson_rate)*uy_star*(1+d*u_abs*(-sqrt(3)/2*cos_psi+1/2*sin_psi))*u_abs*cos(psi-y(3))];
end


% 边界条件
% ya(1) = theta_1(s=0) = theta_11
% ya(2) = theta_2(s=0) = theta_21
% ya(3) = theta_2(s=0) = theta_31
% yb(4) = theta_1'(s=L) = 0
% yb(5) = theta_2'(s=L) = 0
% yb(6) = theta_3'(s=L) = 0
function res = ETR_bcfun(ya,yb)
    global theta_11 theta_21 theta_31
    res = [ya(1)-theta_11;
           ya(2)-theta_21;
           ya(3)-theta_31;
           yb(4);
           yb(5);
           yb(6)];
end

% 初始值估计
% 构建各x所在mesh点，对应y(i)的初始估计值，这对bvp5c求解结果影响较大
function y = ETR_guess(x)
    global theta_11 theta_21 theta_31
    y = [theta_11;
        theta_21;
        theta_31;
         0;
         0;
         0];
end