
clc;
trplot(eye(4),'frame','w0','color','black')
hold on
trplot(eval(T0_w1_g),'color','m');         % sheath末端相对于global{0}系
hold on
trplot(eval(T11_w1_g),'color','r');        % sheath级tube 1末端相对于global{0}系
hold on
trplot(eval(T21_w1_g),'color','g');        % sheath级tube 2末端相对于global{0}系
hold on
trplot(eval(T31_w1_g),'color','b');        % sheath级tube 2末端相对于global{0}系
hold on
