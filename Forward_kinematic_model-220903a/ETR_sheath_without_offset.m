%% 假设三管无偏距情况下CTR扭转柔性模型
% 输入各管关节变量，输出sheath形状
% by Zhang Chao
% Date：2022/7/24

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
cal_point_cnt = 51;                            % 计算201个点，这里即1/4mm
s_mesh = linspace(0,s_1,cal_point_cnt);          % s各离散点弧长
solinit = bvpinit(s_mesh, @CTR_guess);
sol = bvp5c(@CTR_bvpfun, @CTR_bcfun, solinit);  % 这里可求解得到u_x(s)和u_y(s)

%plot(sol.x, sol.y)
%legend('theta_{11}','theta_{21}','theta_{31}','dot(theta_{11})','dot(theta_{21})','dot(theta_{31})');

% 计算曲率u(s)，并存储下来
u_s = zeros(3,cal_point_cnt);
for i=1:cal_point_cnt
    u_s(1,i) = -uy_star/3*(sin(sol.y(1,i))+sin(sol.y(2,i))+sin(sol.y(3,i)));    % u_x(s)
    u_s(2,i) = uy_star/3*(cos(sol.y(1,i))+cos(sol.y(2,i))+cos(sol.y(3,i)));     % u_y(s)
end


%% 根据求解得到的曲率，计算绘制形状
% p_s = zeros(3,cal_point_cnt);
% for i=1:cal_point_cnt
%     syms s_temp
%     R_s_temp = s_temp*skew(u_s(:,i));
%     %test = expm(R_s_temp)*[0;0;1];
%     p_s(:,i) = int(expm(R_s_temp)*[0;0;1],s_temp,0,s_mesh(i));
% end
% 
% plot3(p_s(1,:),p_s(2,:),p_s(3,:),'-g','LineWidth',5);
% hold on
% 经过测试，这个曲率计算和形状绘制的方法和程序应该是正确的，和扭转刚性模型相比，末端的位置偏差很小；

 
%% 根据求解得到的曲率，计算绘制形状(2)
%之前的编程中，对曲率积分存在问题
p_s = zeros(3,cal_point_cnt);
R_e3 = zeros(3,cal_point_cnt);
for i=1:cal_point_cnt
    s_temp = s_mesh(1:i);                   % s=0到当前弧长位置的所有点
    u_s_temp = u_s(:,1:i);                  % 各弧长s点对应曲率
    int_u_temp = trapz(s_temp,u_s_temp,2);  % 使用trapz做数值积分
    R_temp = expm(skew(int_u_temp));
    R_e3(:,i) = R_temp*[0;0;1];             % s=0到当前弧长位置点的中间值R*e3    
    p_s(:,i) = trapz(s_temp,R_e3(:,1:i),2);
end

plot3(p_s(1,:),p_s(2,:),p_s(3,:),'-g','LineWidth',5);
hold on



%% 函数定义 

% 微分方程
% y(1) = theta_1
% y(2) = theta_2
% y(3) = theta_3
% y(4) = theta_1'
% y(5) = theta_2'
% y(6) = theta_3'
function dydx = CTR_bvpfun(x,y)
    global poisson_rate uy_star 
    dydx = [y(4);
            y(5);
            y(6);
            (1+poisson_rate)/3*uy_star^2*(sin(y(1)-y(3))-sin(y(1)-y(3)));
            (1+poisson_rate)/3*uy_star^2*(sin(y(2)-y(1))-sin(y(2)-y(3)));
            (1+poisson_rate)/3*uy_star^2*(sin(y(3)-y(1))-sin(y(3)-y(2)))];
end

% 边界条件
% ya(1) = theta_1(s=0) = theta_11
% ya(2) = theta_2(s=0) = theta_21
% ya(3) = theta_2(s=0) = theta_31
% yb(4) = theta_1'(s=L) = 0
% yb(5) = theta_2'(s=L) = 0
% yb(6) = theta_3'(s=L) = 0
function res = CTR_bcfun(ya,yb)
    global theta_11 theta_21 theta_31
    res = [ya(1)-theta_11;
           ya(2)-theta_21;
           ya(3)-theta_31;
           yb(4);
           yb(5);
           yb(6)];
end

% 初始值估计
% 构建各x所在mesh点，对应y(i)的初始估计值，这\对bvp5c求解结果影响较大
function y = CTR_guess(x)
    global theta_11 theta_21 theta_31
    y = [theta_11;
        theta_21;
        theta_31;
         0;
         0;
         0];
end
