% syms a x
% f = sin(x + a);
% f1 = subs(f, a, 0);
% ezplot(f1, [0,4*pi]);
% figure;
% f2 = subs(f, a, pi/2);
% ezplot(f2, [0,4*pi]);
% 


clear;clc;

k=12;
s0=2;
s=s0:0.01:8;
s2=0:0.01:s0;
dmax=1-6/k;
d=s0*dmax./s;
plot(s,d);
grid;
% 移除坐标轴上的刻度标签
ax = gca;
set(ax, 'XTick', [],'YTick',[]);
% 
% % 可选：移除坐标轴上的刻度线
% set(ax, 'XColor', 'none', 'YColor', 'none');

% % 读取图像
% img = imread('image.jpg');
% imshow(img);
% hold on;
% 
% % 使用ginput函数选择点
% points = ginput();
% 
% % 在图像上绘制标记点
% plot(points(:, 1), points(:, 2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
% 
% % 关闭坐标轴
% axis off;
% 
% % 关闭图像绘制
% hold off;




% clear;clc;close all;
% lamda=[];
% lamda(1)=10000;
% N=[];
% N(1)=9000;
% a=0.2;
% c=0.00001;
% p=0.00001;
% k=12;
% s0=1;
% b=1;
% unep=normrnd(0,0.01,10);
% % unep=0*unep;
% y=1:31;
% for t=1:30
% s=lamda(t)/N(t);
% m=0.5;
% dmax=1-2/k;
% d=s0*dmax/(s^0.7);
% N(t+1)=k*(1-m)*(1-d)*N(t);
% if N(t+1)<0
%     N(t+1)=0;
% end
% lamda(t+1)=(1+(a+unep(t)))*lamda(t)-c*lamda(t)*N(t)-p*lamda(t)*lamda(t);
% end
% plot(y,N,'k',y,lamda,'r')
% grid;
% clear;clc;close all;
% lamda=[];
% lamda(1)=10000;
% N=[];
% N(1)=9000;
% a=0.2;
% c=0.00001;
% p=0.00001;
% k=12;
% s0=1;
% b=1;
% unep=normrnd(0,0.01,10);
% % unep=0*unep;
% y=1:31;
% for t=1:30
% s=lamda(t)/N(t);
% m=0.5;
% dmax=1-2/k;
% d=s0*dmax/(s^0.7);
% N(t+1)=k*(1-m)*(1-d)*N(t);
% if N(t+1)<0
%     N(t+1)=0;
% end
% lamda(t+1)=(1+(a+unep(t)))*lamda(t)-c*lamda(t)*N(t)-p*lamda(t)*lamda(t);
% end
% plot(y,N,'k',y,lamda,'r')
% grid;