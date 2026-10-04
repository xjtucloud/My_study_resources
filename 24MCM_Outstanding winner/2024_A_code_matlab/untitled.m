clear;clc;close all;
lamda=[];
lamda(1)=10000;
NN=[];
NN(1)=5000;
y=1:11;
a=0.2;
c=0.00001;
p=0.00001;
k=12;
s0=1;
b=1;
unep=normrnd(0,0.01,10);
% unep=0*unep;

M=[];
M(1)=5000;
cm=1*10^(-8);
rho=2.2;
alpha=0.04;
%%%%%%%%%%%%  M

for t=1:10
s=lamda(t)/NN(t);
m=0.5+1/(2+exp(b*(s-s0)));
dmax=1-6/k;
d=s0*dmax/(s^0.7);
NN(t+1)=k*(1-m)*(1-d)*NN(t);
lamda(t+1)=(1+(a+unep(t)))*lamda(t)-c*lamda(t)*NN(t)-p*lamda(t)*lamda(t);
M(t+1)=M(t)*exp(rho*(1-M(t))-NN(t)/(1+alpha*M(t)));
end

plot(y,NN,'k',y,lamda,'r')
grid;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画m-s图像

% clear;clc;
% 
% k=12;
% s0=2;
% s=s0:0.01:8;
% s2=0:0.01:s0;
% m=0.5+1./(2+exp(s-s0));
% m2=0.5+1./(2+exp(s2-s0));
% plot(s2,m2,'b--',s,m,'b-');
% grid;
% % 移除坐标轴上的刻度标签
% ax = gca;
% set(ax, 'XTick', []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%画d-s图像

% clear;clc;
% 
% k=12;
% s0=2;
% s=s0:0.01:8;
% s2=0:0.01:s0;
% dmax=1-6/k;
% d=s0*dmax./s;
% plot(s,d);
% grid;
% % 移除坐标轴上的刻度标签
% ax = gca;
% set(ax, 'XTick', [],'YTick',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











