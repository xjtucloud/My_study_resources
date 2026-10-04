% t=1:100;
% L1=plot(t,xx,'k-','LineWidth',1.5);
% hold on
% L2=plot(t,xx1,'r--','LineWidth',1.5);
% L3=plot(t,xx2,'b:','LineWidth',1.5);
% 
% L1=plot(t,yy,'k-','LineWidth',1.5);
% hold on
% L2=plot(t,yy1,'r--','LineWidth',1.5);
% L3=plot(t,yy2,'b:','LineWidth',1.5);
% yy1m=mean(yy1)
% yy2m=mean(yy2)
% legend([L1,L3,L2],'Ideal Reference Model','Constant Sex Ratio Model','Models of Changing Sex Ratios')
% grid;
% hold off

t=44000:4000:440000;
yyaxis left
L1=plot(t,xx1,'k-','LineWidth',1.5);hold on
ylabel('Total number of major species(Excluding lampreys)','FontSize',24);
yyaxis right
L2=plot(t,yy1,'r-','LineWidth',1.5);
ylabel('Total amount of resources','FontSize',24);
xlabel('Initial population size of species','FontSize',24);
hold off
grid