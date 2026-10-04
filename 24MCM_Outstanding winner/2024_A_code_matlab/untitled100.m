clear;clc;close all;
X=[];Y=[];
xx=[];yy=[];
% for v=1:100
for j=1:100
turnmax=100;n=10;
lamda=zeros(1,turnmax);
lamda(1)=4500000;
N=zeros(n,turnmax);
N(:,1)=440000;

zr=0;

a=0.2;
c=normrnd(1*10^(-8),1*10^(-9),1,n);
p=3*10^(-8);
A=normrnd(0.0001,0.0001,1,n);
B=normrnd(0.0000001,0.0000001,1,n);
P=normrnd(0.000001,0.0000001,1,n);

% k=12;
% s0=1;
% b=1;

unep=normrnd(0,0.05,1,turnmax-1);
% unep=0*unep;

y=1:turnmax;
for t=1:turnmax-1
    for i=1:n
        N(i,t+1)=(1-A(i))*N(i,t)+B(i)*N(i,t)*lamda(t)-P(i)*N(i,t)*N(i,t);
        if N(i,t+1)<1 && N(i,t+1)~=0
            N(i,t+1)=0;
            zr=zr+1;
        end
    end
    lamda(t+1)=(1+a+unep(t))*lamda(t)-p*lamda(t)*lamda(t);
    for i=1:n
        lamda(t+1)=lamda(t+1)-c(i)*lamda(t)*N(i,t);
    end
end

% yyaxis right
% plot(y,lamda*10^3,'r','LineWidth',2);
% hold on;
% yyaxis left
% ylim([0,3*10^6]);
% for i=1:n
%     yyaxis left
%     plot(y,N(i,:),'LineWidth',1.5);
% end
% grid;

X(j)=n-zr;
Y(j)=lamda(n)*10^3;
end
Z=[X' Y'];
xx=mean(X);
yy=mean(Y);
% xx(v)=mean(X);
% yy(v)=mean(Y);
% end