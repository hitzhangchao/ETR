%% 3管ETR sheath 扭转柔性模型 - 基于Lexuan Wang的求解方法
% 输入各管关节变量，输出sheath形状
% by Zhang Chao
% Date：2022/8/19

% 首先运行ETR_params.m，加载物理参数

% 全局变量
global theta_11 theta_21 theta_31

%% 运动学输入
sheath_in = [0,90,110];                        % sheath级各管运动学输入
theta_11 = deg2rad(sheath_in(1));               % tube 11的rotation
theta_21 = deg2rad(sheath_in(2));               % tube 21的rotation
theta_31 = deg2rad(sheath_in(3));               % tube 31的rotation


%% 微分方程求解

% 使用bvp5c求解
cal_point_cnt = 201;                            % 计算201个点，这里即1/4mm (调整求解点数，程序运行时间会大幅度缩短)
s_mesh = linspace(0,s_1,cal_point_cnt);          % s各离散点弧长
solinit = bvpinit(s_mesh, @ETR_guess);
sol = bvp5c(@ETR_bvpfun, @ETR_bcfun, solinit);  % 这里可求解得到u_x(s)和u_y(s)

%plot(sol.x, sol.y)
%legend('theta_{11}','theta_{21}','theta_{31}','dot(theta_{11})','dot(theta_{21})','dot(theta_{31})');

%% 计算曲率u(s)，并存储下来
u_s = zeros(3,cal_point_cnt);
for i=1:cal_point_cnt
    %首先计算ψ和u_abs
    C = cos(sol.y(1,i))+cos(sol.y(2,i))+cos(sol.y(3,i));
    S = -sin(sol.y(1,i))-sin(sol.y(2,i))-sin(sol.y(3,i));
    sin_psi = C/sqrt(S^2+C^2);
    cos_psi = S/sqrt(S^2+C^2);
    
    a0 = 4*uy_star*sqrt(S^2+C^2)
    a1 = -12;
    a2 = -3/4*d^2*a0;
    a3 = 3*d^2+d^3*(-4*sin_psi^3+3*sin_psi)*a0/4;
    syms u_abs
    Res = vpasolve(a0+a1*u_abs+a2*u_abs^2+a3*u_abs^3,u_abs);
    u_abs = Res(2);
    
    u_s(1,i) = u_abs*cos_psi;   % u_x(s)
    u_s(2,i) = u_abs*sin_psi;   % u_y(s)
end


%% 根据求解得到的曲率，计算绘制形状
p_s = zeros(3,cal_point_cnt);
for i=1:cal_point_cnt
    syms s_temp
    R_s_temp = s_temp*skew(u_s(:,i));       % vec2skew()不能用了，换成了skew()
    %test = expm(R_s_temp)*[0;0;1];
    p_s(:,i) = int(expm(R_s_temp)*[0;0;1],s_temp,0,s_mesh(i));
end

plot3(p_s(1,:),p_s(2,:),p_s(3,:),'-b','LineWidth',5);
hold on
% 经过测试，这个曲率计算和形状绘制的方法和程序应该是正确的，和扭转刚性模型相比，末端的位置偏差很小；


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
