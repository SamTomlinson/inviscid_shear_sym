clear all; close all; clc;

%% delta=0.1

kd01=linspace(0.01,0.15,40); delta=0.1; store01=[]; erstore1=[]; sstore1=[];
for kk=1:40
    k=kd01(kk)
    point_norm_rsing
    store01=[store01,k*cn(end)]; erstore1=[erstore1,err]; sstore1=[sstore1,ssRn];
end
output = matfile('stability_output.mat','Writable',true);
output.d01=store01; output.s01=sstore1; output.e01=erstore1;

%% delta=0.3

kd03=linspace(0.01,0.2,40); delta=0.3; store03=[]; erstore3=[]; sstore3=[];
for kk=1:40
    k=kd03(kk)
    point_norm_rsing
    store03=[store03,k*cn(end)]; erstore3=[erstore3,err]; sstore3=[sstore3,ssRn];
end
output.d03=store03; output.s03=sstore3; output.e03=erstore3;

%% delta=0.5

kd05=linspace(0.01,0.25,40); delta=0.5; store05=[]; erstore5=[]; sstore5=[];
for kk=1:40
    k=kd05(kk)
    point_norm_rsing
    store05=[store05,k*cn(end)]; erstore5=[erstore5,err]; sstore5=[sstore5,ssRn];
end
output.d05=store05; output.s05=sstore5; output.e05=erstore5;

%% delta=0.7 

kd07=linspace(0.01,0.3,40); delta=0.7; store07=[]; erstore7=[]; sstore7=[];
for kk=1:40
    k=kd07(kk)
    point_norm_rsing
    store07=[store07,k*cn(end)]; erstore7=[erstore7,err]; sstore7=[sstore7,ssRn];
end
output.d07=store07; output.s07=sstore7; output.e07=erstore7;

%% delta=0.9 

kd09=linspace(0.01,0.35,40); delta=0.9; store09=[]; erstore9=[]; sstore9=[];
for kk=1:40
    k=kd09(kk)
    point_norm_rsing
    store09=[store09,k*cn(end)]; erstore9=[erstore9,err]; sstore9=[sstore9,ssRn];
end
output.d09=store09; output.s09=sstore9; output.e09=erstore9;

%% plotting 

% figure()
% plot(kd01(1:11),real(store01(1:11)),'LineWidth',2); hold on;
% plot(kd03(1:15),real(store03(1:15)),'LineWidth',2); hold on;
% plot(kd05(1:16),real(store05(1:16)),'LineWidth',2)
% plot(kd07(1:15),real(store07(1:15)),'LineWidth',2)
% plot(kd09(1:15),real(store09(1:15)),'LineWidth',2)
% set(gca,'fontsize',20)
% xlabel('$k$','FontSize',30,'Interpreter','latex'); 
% ylabel('$k \Re(c_m)$','FontSize',30,'Interpreter','latex'); 
% l1=legend('$\delta=0.1$','$\delta=0.3$','$\delta=0.5$','$\delta=0.7$','$\delta=0.9$');
% set(l1,'FontSize',30,'Interpreter','latex','location','northwest')
% set(gcf, 'Position', [0, 0, 600, 900])
% %ylim([0,0.45])
% 
% 
% figure()
% plot(kd01(1:11),imag(store01(1:11)),'LineWidth',2); hold on;
% plot(kd03(1:15),imag(store03(1:15)),'LineWidth',2); hold on;
% plot(kd05(1:16),imag(store05(1:16)),'LineWidth',2)
% plot(kd07(1:15),imag(store07(1:15)),'LineWidth',2)
% plot(kd09(1:15),imag(store09(1:15)),'LineWidth',2)
% set(gca,'fontsize',20)
% xlabel('$k$','FontSize',30,'Interpreter','latex'); 
% ylabel('$k \Im(c_m)$','FontSize',30,'Interpreter','latex'); 
% l1=legend('$\delta=0.1$','$\delta=0.3$','$\delta=0.5$','$\delta=0.7$','$\delta=0.9$');
% set(l1,'FontSize',30,'Interpreter','latex','location','northwest')
% set(gcf, 'Position', [0, 0, 600, 900])
% %ylim([0,0.0275])
% 
% figure()
% plot(kd01(1:11),abs(sstore1(1:11)),'LineWidth',2); hold on;
% plot(kd03(1:15),abs(sstore3(1:15)),'LineWidth',2); hold on;
% plot(kd05(1:16),abs(sstore5(1:16)),'LineWidth',2)
% plot(kd07(1:15),abs(sstore7(1:15)),'LineWidth',2)
% plot(kd09(1:15),abs(sstore9(1:15)),'LineWidth',2)
% set(gca,'fontsize',20)
% xlabel('$k$','FontSize',30,'Interpreter','latex'); 
% ylabel('$|P_0|$','FontSize',30,'Interpreter','latex'); 
% l1=legend('$\delta=0.1$','$\delta=0.3$','$\delta=0.5$','$\delta=0.7$','$\delta=0.9$');
% set(l1,'FontSize',30,'Interpreter','latex','location','northwest')
% set(gcf, 'Position', [0, 0, 600, 900])
%ylim([0,0.45]
% 
% %% interpolation and replotting 
% 
% kd01i=linspace(kd01(1),kd01(11),1000); kd03i=linspace(kd03(1),kd03(15),1000); kd05i=linspace(kd05(1),kd05(15),1000); kd07i=linspace(kd07(1),kd07(14),1000); kd09i=linspace(kd09(1),kd09(14),1000); 
% 
% kcr01=interp1(kd01(1:11),real(store01(1:11)),kd01i,'spline'); kcr03=interp1(kd03(1:15),real(store03(1:15)),kd03i,'spline'); kcr05=interp1(kd05(1:15),real(store05(1:15)),kd05i,'spline'); kcr07=interp1(kd07(1:14),real(store07(1:14)),kd07i,'spline'); kcr09=interp1(kd09(1:14),real(store09(1:14)),kd09i,'spline');
% 
% kci01=interp1(kd01(1:11),imag(store01(1:11)),kd01i,'spline'); kci03=interp1(kd03(1:15),imag(store03(1:15)),kd03i,'spline'); kci05=interp1(kd05(1:15),imag(store05(1:15)),kd05i,'spline'); kci07=interp1(kd07(1:14),imag(store07(1:14)),kd07i,'spline'); kci09=interp1(kd09(1:14),imag(store09(1:14)),kd09i,'spline');
% 
% ss01=interp1(kd01(1:11),abs(sstore1(1:11)),kd01i,'spline'); ss03=interp1(kd03(1:15),abs(sstore3(1:15)),kd03i,'spline'); ss05=interp1(kd05(1:15),abs(sstore5(1:15)),kd05i,'spline'); ss07=interp1(kd07(1:14),abs(sstore7(1:14)),kd07i,'spline'); ss09=interp1(kd09(1:14),abs(sstore9(1:14)),kd09i,'spline');
% 
% figure()
% plot(kd01i,kcr01,'LineWidth',2); hold on;
% plot(kd03i,kcr03,'LineWidth',2); hold on;
% plot(kd05i,kcr05,'LineWidth',2)
% plot(kd07i,kcr07,'LineWidth',2)
% plot(kd09i,kcr09,'LineWidth',2)
% set(gca,'fontsize',20)
% xlabel('$k$','FontSize',30,'Interpreter','latex'); 
% ylabel('$k \Re(c_m)$','FontSize',30,'Interpreter','latex'); 
% l1=legend('$\delta=0.1$','$\delta=0.3$','$\delta=0.5$','$\delta=0.7$','$\delta=0.9$');
% set(l1,'FontSize',30,'Interpreter','latex','location','northwest')
% set(gcf, 'Position', [0, 0, 600, 900])
% %ylim([0,0.45])
% 
% 
% figure()
% plot(kd01i,kci01,'LineWidth',2); hold on;
% plot(kd03i,kci03,'LineWidth',2); hold on;
% plot(kd05i,kci05,'LineWidth',2)
% plot(kd07i,kci07,'LineWidth',2)
% plot(kd09i,kci09,'LineWidth',2)
% set(gca,'fontsize',20)
% xlabel('$k$','FontSize',30,'Interpreter','latex'); 
% ylabel('$k \Im(c_m)$','FontSize',30,'Interpreter','latex'); 
% l1=legend('$\delta=0.1$','$\delta=0.3$','$\delta=0.5$','$\delta=0.7$','$\delta=0.9$');
% set(l1,'FontSize',30,'Interpreter','latex','location','northwest')
% set(gcf, 'Position', [0, 0, 600, 900])
% %ylim([0,0.0275])
% 
% figure()
% plot(kd01i,ss01,'LineWidth',2); hold on;
% plot(kd03i,ss03,'LineWidth',2); hold on;
% plot(kd05i,ss05,'LineWidth',2)
% plot(kd07i,ss07,'LineWidth',2)
% plot(kd09i,ss09,'LineWidth',2)
% set(gca,'fontsize',20)
% xlabel('$k$','FontSize',30,'Interpreter','latex'); 
% ylabel('$|P_0|$','FontSize',30,'Interpreter','latex'); 
% l1=legend('$\delta=0.1$','$\delta=0.3$','$\delta=0.5$','$\delta=0.7$','$\delta=0.9$');
% set(l1,'FontSize',30,'Interpreter','latex','location','northwest')
% set(gcf, 'Position', [0, 0, 600, 900])
