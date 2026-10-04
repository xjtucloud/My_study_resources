clear;clc;close all;
% X=[];Y=[];
% xx1=[];yy1=[];
% for j=1:100
turnmax=100;n=7;m0=2;

lamda=zeros(1,turnmax);
lamda(1)=4500000;

N=zeros(n,turnmax);
N(:,1)=440000;

W=zeros(m0,turnmax);
W(:,1)=0.44*10^6;
cm0=normrnd(1*10^(-8),1*10^(-9),1,m0);
AM0=normrnd(0.0001,0.0001,1,m0);
BM0=normrnd(0.0000001,0.0000001,1,m0);
PM0=normrnd(0.000001,0.0000001,1,m0);

rho0=normrnd(2.2,0.1,1,m0);
alpha0=normrnd(0.04,0.005,1,m0);
phi0=normrnd(4,0.05,1,m0);


% M=zeros(m,turnmax);
% M(:,1)=0.44*10^6;
% cm=normrnd(1*10^(-8),1*10^(-9),1,m);
% AM=normrnd(0.0001,0.0001,1,m);
% BM=normrnd(0.0000001,0.0000001,1,m);
% PM=normrnd(0.000001,0.0000001,1,m);
% 
% rho=normrnd(2.2,0.1,1,m);
% alpha=normrnd(0.04,0.005,1,m);
% 
% %%%%%%%%%%%%  M

% NN=[];
% NN(1)=0.44*10^6;
% s0=0.5;
% b=1;
% cnn=1*10^(-8);
% phi=normrnd(1,0.05,1,m);
% %%%%%%%%%%%%  NN

VV=[];
VV(1)=0.44*10^6;
cvv=1*10^(-8);
%%%%%%%%%%%%  VV

zrn=0;
zr=0;
%%%%%   death number

a=0.2;
p=3*10^(-8);
% unepa0=normrnd(a,0.01,1,20);
% unepa1=normrnd(a,0.05,1,20);
% unepa2=normrnd(a,0.1,1,30);
% unepa3=normrnd(a,0.15,1,30);
% unepa=[unepa0 unepa1 unepa2 unepa3];

unepa=normrnd(a,0.1,1,turnmax-1);
%%%%%%%%%%%%  lamda

c=normrnd(1*10^(-8),1*10^(-9),1,n);

A=normrnd(0.0001,0.0001,1,n);
B=normrnd(0.0000001,0.0000001,1,n);
P=normrnd(0.000001,0.0000001,1,n);
%%%%%%%%%%%%  N
sunum=[];
s=[];
y=1:turnmax;

for t=1:turnmax-1
    m00=1/2;
    k00=3;
    VV(t+1)=sum(W(:,t)'.*phi0)*(1-exp(-VV(t)/(10^6)/(1+sum(alpha0.*W(:,t)')/(10^6))))*k00*(1-m00)^(5/4);
    %%%%%%%%%%%% VV


    for i=1:n
        N(i,t+1)=(1-A(i))*N(i,t)+B(i)*N(i,t)*lamda(t)-P(i)*N(i,t)*N(i,t);
        if N(i,t+1)<2 && N(i,t+1)~=0
            N(i,t+1)=0;
            zr=zr+1;
        end
    end
    %%%%%%%%  N


    for i=1:m0
        W(i,t+1)=0.5*W(i,t)*exp(rho0(i)*(1-W(i,t)/10^6)-VV(t)/(10^6+alpha0(i)*W(i,t)))...
        +0.5*((1-AM0(i))*W(i,t)+BM0(i)*W(i,t)*lamda(t)-PM0(i)*W(i,t)*W(i,t));
    end

    %%%%%%%%%%% W
    lamda(t+1)=(1+unepa(t))*lamda(t)-p*lamda(t)*lamda(t);
    for i=1:n
        lamda(t+1)=lamda(t+1)-c(i)*lamda(t)*N(i,t);
    end

    lamda(t+1)=lamda(t+1)-cvv*lamda(t)*VV(t);
    for i=1:m0
        lamda(t+1)=lamda(t+1)-cm0(i)*lamda(t)*W(i,t);
    end
    %%%%%%%%  lamda
end


yyaxis right
L1=plot(y,lamda*10^3,'r','LineWidth',0.5);hold on;
yyaxis left
L20=plot(y,VV,'k--','LineWidth',1.5);hold on;
L30=plot(y,W(1,:),'g-','LineWidth',1);hold on;
for i=2:m0
    yyaxis left
    plot(y,W(i,:),'g-','LineWidth',1);hold on;
end
L4=plot(y,N(1,:),'--','LineWidth',0.5);hold on;
for i=2:n
    yyaxis left
    plot(y,N(i,:),'--','LineWidth',0.5);hold on;
end
grid;
hold off;

legend([L1,L4,L20,L30],'resources',...
    'other species','the parasite','host of the parasite')


% X(j)=n-zr;
% for i=1:10-n
%     if zrm(i)==0
%         X(j)=X(j)+1;
%     end
% end
% Y(j)=lamda(n)*10^3;
% end
% X1=X';
% Y1=Y';
% Z1=[X1 Y1];