function T = T_u_s(u,s)
% func: 扭转刚性模型中的转换矩阵
% input: u - local系中的常曲率,输入格式为[u_x;u_y;0]; 
%        s - 弧长
% output: 曲线弧长s处local frame相对于伸出起始端的transformation matrix

% 先判断norm(u)是否为0
if norm(u)==0
    T = [eye(3),[0;0;s];[0 0 0 1]];
else
    % 旋转矩阵的元素
    a11 = (u(1)^2+u(2)^2*cos(s*norm(u)))/(norm(u))^2;
    a12 = u(1)*u(2)*(1-cos(s*norm(u)))/(norm(u))^2;
    a13 = u(2)*sin(s*norm(u))/norm(u);
    a21 = a12;
    a22 = (u(2)^2+u(1)^2*cos(s*norm(u)))/(norm(u))^2;
    a23 = -u(1)*sin(s*norm(u))/norm(u);
    a31 = -a13;
    a32 = -a23;
    a33 = cos(s*norm(u));
    R = [a11,a12,a13;a21,a22,a23;a31,a32,a33];

    % 位置矢量元素
    a14 = u(2)*(1-cos(s*norm(u)))/(norm(u))^2;
    a24 = -u(1)*(1-cos(s*norm(u)))/(norm(u))^2;
    a34 = sin(s*norm(u))/norm(u);
    p = [a14;a24;a34];

    T = [R,p;[0 0 0 1]];
end
end

