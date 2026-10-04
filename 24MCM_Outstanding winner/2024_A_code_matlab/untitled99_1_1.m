clear;clc;close all;
X=[];Y=[];
xx1=[];yy1=[];
% for v=1:100
% for j=1:100
turnmax=100;n=7;
lamda=zeros(1,turnmax);
lamda(1)=4500000;%4500000;
N=zeros(n,turnmax);
N(:,1)=440000;

M=zeros(10-n,turnmax);
M(:,1)=0.44*10^6;%40000+4000*v;
cm=normrnd(1*10^(-8),1*10^(-9),1,10-n);
rho=normrnd(2.2,0.1,1,10-n);
alpha=normrnd(0.04,0.005,1,10-n);
AM=normrnd(0.0001,0.0001,1,10-n);
BM=normrnd(0.0000001,0.0000001,1,10-n);
PM=normrnd(0.000001,0.0000001,1,10-n);
%%%%%%%%%%%%  M

NN=[];
NN(1)=0.44*10^6;%40000+4000*v;
k=12;
s0=0.5;
b=1;
cnn=1*10^(-8);
phi=normrnd(1,0.05,1,10-n);
%%%%%%%%%%%%  NN

zrn=0;
zr=0;
zrm=zeros(1,10-n);
%%%%%   death number

a=0.2;
p=3*10^(-8);
% unepa0=normrnd(a,0.01,1,20);
% unepa1=normrnd(a,0.05,1,20);
% unepa2=normrnd(a,0.1,1,30);
% unepa3=normrnd(a,0.15,1,30);
% unepa=[unepa0 unepa1 unepa2 unepa3];

unepa=normrnd(a,0.01,1,turnmax-1);
%%%%%%%%%%%%  lamda

c=normrnd(1*10^(-8),1*10^(-9),1,n);

A=normrnd(0.0001,0.0001,1,n);
B=normrnd(0.0000001,0.0000001,1,n);
P=normrnd(0.000001,0.0000001,1,n);
%%%%%%%%%%%%  N
sunum=[];
s=[];m=[];
y=1:turnmax;
for t=1:turnmax-1
    sunum(t)=NN(t)+sum(N(:,t))+sum(M(:,t));
    s(t)=lamda(t)/sunum(t);
    m(t)=0.5+1/(2+exp(b*(s(t)-s0)));
    %dmax=1-6/k;
    k(t)=12;%*(1-m(t))^(1/4);
    % d=2/(1+(exp(2*s(t))));  %s0*dmax/s;
    NN(t+1)=sum(M(:,t)'.*phi)*(1-exp(-NN(t)/(10^6)/(1+sum(alpha.*M(:,t)')/(10^6))))*k(t)*(1-m(t))^(5/4);
    %%%%%%%%  NN
    if NN(t+1)<2
            NN(t+1)=0;
            zrn=zrn+1;
    end
    for i=1:n
        N(i,t+1)=(1-A(i))*N(i,t)+B(i)*N(i,t)*lamda(t)-P(i)*N(i,t)*N(i,t);
        if N(i,t+1)<2 && N(i,t+1)~=0
            N(i,t+1)=0;
            zr=zr+1;
        end
    end
    if t==50
        NN(t+1)=NN(t+1)/2;
    end
    %%%%%%%%  N
    for i=1:10-n
        M(i,t+1)=0.5*M(i,t)*exp(rho(i)*(1-M(i,t)/10^6)-NN(t)/(10^6+alpha(i)*M(i,t)))...
            +0.5*((1-AM(i))*M(i,t)+BM(i)*M(i,t)*lamda(t)-PM(i)*M(i,t)*M(i,t));
            % +0.4*((1-1.2*AM(i))*M(i,t)+2*BM(i)*M(i,t)*lamda(t)-3*PM(i)*M(i,t)*M(i,t));
        if M(i,t+1)<2
                M(i,t+1)=0;
                zrm(i)=zrm(i)+1;
        end
    end
    lamda(t+1)=(1+unepa(t))*lamda(t)-p*lamda(t)*lamda(t)-cnn*lamda(t)*NN(t);
    for i=1:10-n
        lamda(t+1)=lamda(t+1)-cm(i)*lamda(t)*M(i,t);
    end
    for i=1:n
        lamda(t+1)=lamda(t+1)-c(i)*lamda(t)*N(i,t);
    end
    %%%%%%%%  lamda
end


yyaxis right
L1=plot(y,lamda*10^3,'r','LineWidth',1);
yyaxis right
ylabel('The total amount of resources/kJ','FontSize',24);
ylim([3*10^9,6*10^9]);
hold on;
yyaxis left
ylabel('The population size of species','FontSize',24);
ylim([0,3*10^6]);
L2=plot(y,NN,'k-','LineWidth',1.5);
L3=plot(y,M(1,:),'b-','LineWidth',1.5);
for i=2:10-n
    yyaxis left
    plot(y,M(i,:),'b-','LineWidth',1.5);
end
L4=plot(y,N(1,:),'--','LineWidth',0.5);
for i=2:n
    yyaxis left
    plot(y,N(i,:),'--','LineWidth',0.5);
end
grid;
hold off;
legend([L1,L2,L3,L4],'resources',...
    'the lamprey','host of the lamprey',...
    'other species')
xlabel('Number of reproductive rounds','FontSize',24);

% X(j)=n-zr;
% for i=1:10-n
%     if zrm(i)==0
%         X(j)=X(j)+1;
%     end
% end
% Y(j)=lamda(n)*10^3;
% end
% X1=X;
% Y1=Y;
% Z1=[X1' Y1'];
% xx1=mean(X1);
% yy1=mean(Y1);

% xx1(v)=mean(X1);
% yy1(v)=mean(Y1);
% end
% zz1=[xx1' yy1'];