clear all
close all
clc

ifig =1

alpha0 = -4; %deg
n = 101;
b = 1;

Vinf= [0:1:60];; %dimensionless
theta = linspace(pi,0,n);
y = b/2*cos(theta); 

%% Q1 

%AR = [4 6 8 10];
AR = 6;
AOA = [0 5 10];
B=ones(n-1,1)*degtorad(AOA+4); %AOA - offset
c0 = b./AR;
S = b^2./AR;


for i = 1:(length(y)-1)
    yc(i) = y(i) + (y(i+1)-y(i))/2; %center of intervals
    dy(i) = y(i+1) - y(i);
end


% creation of the matrix A: A*gamma=wi
for i =1:n-1
    for j=2:n-2
        A(i,j) = 2*(1/(yc(i)-y(j))-1/(yc(i)-y(j+1)));
    end
    A(i,1) = 1/(yc(i)-y(1)) - 2/(yc(i)-y(2));
    A(i,n-1) = 2/(yc(i)-y(n-1)) - 1/(yc(i)-y(n));
end

A = A./(4*pi);
save('A.mat','A')
% resolution of Prandlt
%AOA1
% for j=1:length(AR) %loop over AR
%     WI.A1(:,j)=inv(eye(n-1)./(pi*Vinf*c0(j)).*inv(A)+1/Vinf)*B(:,1);
%     ALPH_i.A1= WI.A1/Vinf;
%     GAM.A1 = inv(A)*WI.A1;
%     CL(1,:) = 2./(Vinf.*S)*sum(GAM.A1(:,j).*dy');
%     CD(1,:) = 2./(Vinf.*S)*sum(ALPH_i.A1(:,j).*GAM.A1(:,j).*dy');
% end

for j=1:length(AR) %loop over AR
    GAM.A1(:,j) = inv(eye(n-1)./(pi*Vinf*c0(j))+ A./Vinf)*B(:,1);
    WI.A1(:,j) = A*GAM.A1(:,j);
    ALPH_i.A1(:,j) = WI.A1(:,j)./Vinf;
    CL(1,:) = 2./(Vinf.*S)*sum(GAM.A1(:,j).*dy');
    CD(1,:) = 2./(Vinf.*S)*sum(ALPH_i.A1(:,j).*GAM.A1(:,j).*dy');
end

%AOA2
% for j=1:length(AR) %loop over AR
%     WI.A2(:,j)=inv(eye(n-1)/(pi*Vinf*c0(j)).*inv(A)+1/Vinf)*B(:,2);
%     ALPH_i.A2=WI.A2/Vinf;
%     GAM.A2 = inv(A)*WI.A2;
%     CL(2,:) = 2./(Vinf.*S)*sum(GAM.A1(:,j).*dy');
%     CD(2,:) = 2./(Vinf.*S)*sum(ALPH_i.A1(:,j).*GAM.A1(:,j).*dy');
% end

for j=1:length(AR) %loop over AR
    GAM.A2(:,j) = inv(eye(n-1)./(pi*Vinf*c0(j))+ A./Vinf)*B(:,2);
    WI.A2(:,j) = A*GAM.A2(:,j);
    ALPH_i.A2(:,j) = WI.A2(:,j)./Vinf;
    CL(2,:) = 2./(Vinf.*S)*sum(GAM.A2(:,j).*dy');
    CD(2,:) = 2./(Vinf.*S)*sum(ALPH_i.A2(:,j).*GAM.A2(:,j).*dy');
end

%AOA3
% for j=1:length(AR) %loop over AR
%     WI.A3(:,j)=inv(eye(n-1)/(pi*Vinf*c0(j)).*inv(A)+1/Vinf)*B(:,3);
%     ALPH_i.A3= WI.A3/Vinf;
%     GAM.A3 = inv(A)*WI.A3;
%     CL(3,:) = 2./(Vinf.*S)*sum(GAM.A1(:,j).*dy');
%     CD(3,:) = 2./(Vinf.*S)*sum(ALPH_i.A1(:,j).*GAM.A1(:,j).*dy');
% end

for j=1:length(AR) %loop over AR
    GAM.A3(:,j) = inv(eye(n-1)./(pi*Vinf*c0(j))+ A./Vinf)*B(:,3);
    WI.A3(:,j) = A*GAM.A3(:,j);
    ALPH_i.A3(:,j) = WI.A3(:,j)./Vinf;
    CL(3,:) = 2./(Vinf.*S)*sum(GAM.A3(:,j).*dy');
    CD(3,:) = 2./(Vinf.*S)*sum(ALPH_i.A3(:,j).*GAM.A3(:,j).*dy');
end

fig1 = figure(1)
for i=1:length(AR)
    subplot(2,2,i)
    plot(yc,radtodeg(ALPH_i.A1(:,i)))
    hold on 
    plot(yc,radtodeg(ALPH_i.A2(:,i)))
    plot(yc,radtodeg(ALPH_i.A3(:,i)))
    legend('AOA = 0°','AOA = 5°','AOA = 10°')
    hold off
    title(['AR = ' num2str(AR(i)) ' (rectangular wing)'],'Fontsize',14)
    xlabel('y_c','Fontsize',14)
    ylabel('\alpha_i [°]','Fontsize',14)
end


fig2 = figure(2)
subplot(1,3,1) %AOA1
for i = 1:length(AR)
    hold on
    plot(yc,radtodeg(ALPH_i.A1(:,i)))
    cc(i) = cellstr(num2str(AR(i)));
end
hold off
lgd = legend(cc);
lgd.Title.String='Value of AR';
lcn='bestoutside';
grid minor
xlabel('y_c','Fontsize',14)
ylabel('\alpha_i [°]','Fontsize',14)
title('AOA = 0°','Fontsize',14)

subplot(1,3,2) %AOA2
for i = 1:length(AR)
    hold on
    plot(yc,radtodeg(ALPH_i.A2(:,i)))
    cc(i) = cellstr(num2str(AR(i)));
end
hold off
lgd = legend(cc);
lgd.Title.String='Value of AR';
lcn='bestoutside';
grid minor
xlabel('y_c','Fontsize',14)
ylabel('\alpha_i [°]','Fontsize',14)
title('AOA = 5°','Fontsize',14)

subplot(1,3,3) %AOA3
for i = 1:length(AR)
    hold on
    plot(yc,radtodeg(ALPH_i.A3(:,i)))
    cc(i) = cellstr(num2str(AR(i)));
end
hold off
lgd = legend(cc);
lgd.Title.String='Value of AR';
lcn='bestoutside';
grid minor
xlabel('y_c','Fontsize',14)
ylabel('\alpha_i [°]','Fontsize',14)
title('AOA = 10°','Fontsize',14)
saveas(fig2, 'alpha_Q1','epsc')


% for j=1:length(AR) %loop over AR
%     GAM.A4(:,j) = inv(eye(n-1)./(pi*Vinf*c0(j))+ A./Vinf)*B(:,4);
%     WI.A4(:,j) = A*GAM.A4(:,j);
%     ALPH_i.A4(:,j) = WI.A4(:,j)./Vinf;
%     CL(4,:) = 2./(Vinf.*S)*sum(GAM.A4(:,j).*dy');
%     CD(4,:) = 2./(Vinf.*S)*sum(ALPH_i.A4(:,j).*GAM.A4(:,j).*dy');
%     GAM_rect_tilde(:,j) = GAM.A4(:,j)/(Vinf.*c0(j));
% end



% %% Q2
% 
% c0=4*b./(pi.*AR); 
% 
% %AOA1
% for j=1:length(AR) %loop over AR
%     GAM_max(j) = 2*b*Vinf*(degtorad(4+AOA(1)))/(1+AR(j)/2);
%     %c.A1(:,j)=c0(j).*Vinf./GAM_max.*inv(A)*ones(n-1,1).*(1/(pi*Vinf*c0(j))*GAM_max-degtorad(4+AOA(1)));
%     c.A1(:,j)=((1-(yc./(b/2)).^2).^(1/2))';
%     %WI2.A1=A*GAM_max*c.A1  ;
%     WI2.A1(:,j) = -GAM_max(j)/(2*b)*ones(1,length(yc));
%     ALPH2_i.A1(:,j)= -WI2.A1(:,j)/Vinf;
%     CD2(1,:)= pi*GAM_max(j)^2./(4*S*Vinf^2);
%     CL2(1,:) = pi/2*b./S.*GAM_max(j)/Vinf;
%     GAM2.A1(:,j) = GAM_max(j)*sqrt(1-(2.*yc./b).^2);
% end
% 
% %AOA2
% for j=1:length(AR) %loop over AR
%     GAM_max=2*b*Vinf*(degtorad(4+AOA(2)))/(1+AR(j)/2);
%     %c.A2(:,j)=c0(j)*Vinf/GAM_max*inv(A)*ones(n-1,1)*(1/(pi*Vinf*c0(j))*GAM_max-degtorad(4+AOA(2)));
%     c.A2(:,j)=((1-(yc./(b/2)).^2).^(1/2))';
%     %WI2.A2=A*GAM_max*c.A2  ;
%     WI2.A2(:,j) = -GAM_max/(2*b)*ones(1,length(yc));
%     ALPH2_i.A2(:,j)= -WI2.A2(:,j)/Vinf;
%     CD2(2,:)= pi*GAM_max^2./(4*S*Vinf^2);
%     CL2(2,:) = pi/2*b./S.*GAM_max/Vinf;
%     GAM2.A2(:,j) = GAM_max*sqrt(1-(2.*yc./b).^2);
% end
% 
% %AOA3
% for j=1:length(AR) %loop over AR
%     GAM_max= 2*b*Vinf*(degtorad(4+AOA(3)))/(1+AR(j)/2);
%     %c.A3(:,j)=c0(j)*Vinf/GAM_max*inv(A)*ones(n-1,1)*(1/(pi*Vinf*c0(j))*GAM_max-degtorad(4+AOA(3)));
%     c.A3(:,j)=((1-(yc./(b/2)).^2).^(1/2))';
%     %WI2.A3=A*GAM_max*c.A3 ;
%     WI2.A3(:,j) = -GAM_max/(2*b)*ones(1,length(yc));
%     ALPH2_i.A3(:,j)= -WI2.A3(:,j)/Vinf;
%     CD2(3,:)=pi*GAM_max^2./(4*S*Vinf^2);
%     CL2(3,:) = pi/2*b./S.*GAM_max/Vinf;
%     GAM2.A3(:,j) = GAM_max*sqrt(1-(2.*yc./b).^2);
%     
% end
% 
% fig3 = figure(3)
% for i=1:length(AR)
%     plot(yc,radtodeg(ALPH2_i.A1(:,i)))
%     hold on
%     cc(i) = cellstr(num2str(AR(i)));
%     title('AOA = 0° (elliptic wing)','Fontsize',16)
%     xlabel('y_c','Fontsize',16)
%     ylabel('\alpha_i [°]','Fontsize',16)
% end
% hold off
% lgd = legend(cc);
% lgd.Title.String='Value of AR';
% lcn='bestoutside';
% grid minor
% saveas(fig3, 'alpha_cte','epsc')
% 
% fig4 = figure(4)
% 
% plot(AR,radtodeg(ALPH2_i.A1(1,:)))
% hold on
% plot(AR,radtodeg(ALPH2_i.A2(1,:)))
% plot(AR,radtodeg(ALPH2_i.A3(1,:)))
% legend('AOA = 0°', 'A0A = 5°','AOA = 10 °')
% xlabel('AR','Fontsize',16)
% ylabel('\alpha_i [°]','Fontsize',16)
% hold off
% grid minor
% saveas(fig4, 'alpha_AR','epsc')
% 
% %Plot gamma rectangular
% 
% fig5 = figure(5)
% subplot(1,3,1)
% hold on
% for i =1:length(AR)
%     plot(yc,GAM.A1(:,i))
%     cc(i)=cellstr(num2str(AR(i)));
% end
% lgd = legend(cc,'Location','Southeast');
% lgd.Title.String='Value of AR';
% grid minor
% title('Rectangular, AOA = 0°','Fontsize',14)
% xlabel('y_c','Fontsize',14)
% ylabel('\Gamma','Fontsize',14)
% hold off
% 
% subplot(1,3,2)
% hold on
% for i =1:length(AR)
%     plot(yc,GAM.A2(:,i))
% cc(i)=cellstr(num2str(AR(i)));
% end
% lgd = legend(cc,'Location','Southeast');
% lgd.Title.String='Value of AR';
% grid minor
% title('Rectangular, AOA = 5°','Fontsize',14)
% xlabel('y_c','Fontsize',14)
% ylabel('\Gamma','Fontsize',14)
% hold off
% 
% subplot(1,3,3)
% hold on
% for i =1:length(AR)
%     plot(yc,GAM.A3(:,i))
% cc(i)=cellstr(num2str(AR(i)));
% end
% lgd = legend(cc,'Location','Southeast');
% lgd.Title.String='Value of AR';
% lcn='bestoutside';
% grid minor
% title('Rectangular, AOA = 10°','Fontsize',14)
% xlabel('y_c','Fontsize',14)
% ylabel('\Gamma','Fontsize',14)
% hold off
% saveas(fig5, 'gamm_rect','epsc')
% 
% % Plot gamma eliptic
% fig6 = figure(6)
% subplot(1,3,1)
% hold on
% for i =1:length(AR)
%     plot(yc,GAM2.A1(:,i))
%     cc(i)=cellstr(num2str(AR(i)));
% end
% lgd = legend(cc,'Location','Southeast');
% lgd.Title.String='Value of AR';
% grid minor
% title('Elliptic, AOA = 0°','Fontsize',14)
% xlabel('y_c','Fontsize',14)
% ylabel('\Gamma','Fontsize',14)
% hold off
% 
% subplot(1,3,2)
% hold on
% for i =1:length(AR)
%     plot(yc,GAM2.A2(:,i))
% cc(i)=cellstr(num2str(AR(i)));
% end
% lgd = legend(cc,'Location','Southeast');
% lgd.Title.String='Value of AR';
% grid minor
% title('Elliptic, AOA = 5°','Fontsize',14)
% xlabel('y_c','Fontsize',14)
% ylabel('\Gamma','Fontsize',14)
% hold off
% 
% subplot(1,3,3)
% hold on
% for i =1:length(AR)
%     plot(yc,GAM2.A3(:,i))
% cc(i)=cellstr(num2str(AR(i)));
% end
% lgd = legend(cc,'Location','Southeast');
% lgd.Title.String='Value of AR';
% grid minor
% title('Elliptic, AOA = 10°','Fontsize',14)
% xlabel('y_c','Fontsize',14)
% ylabel('\Gamma','Fontsize',14)
% hold off
% saveas(fig6, 'gamm_ell','epsc')
% 
% %% Q3
% 
% AR = 6;
% AOA = 6;
% TR = [0.2 0.4 0.6 0.8 1];
% S=b^2/AR;
% b = 1;
% 
% c_r = 2*b./(1 + TR)*1/AR; %c_root
% c_t = c_r.*TR; %c_tip
% 
% for i=1:n-1
%     chord(:,i) = c_r - abs(yc(i))*(c_r - c_t)/b*2; %1line = 1 TR
% end
% 
% %plot tapered wing
% color = jet(length(TR));
% 
% fig7 = figure(7)
% hold on
% for i = 1:length(TR)  
%     p(i) = plot(yc,chord(i,:)/2,'color',color(i,:),'LineWidth',2);
%     plot(yc,-chord(i,:)/2,'color',color(i,:),'LineWidth',2)
%     cc(i) = cellstr(num2str(TR(i)));
% end
% hold off
% line([0 0], ylim);  %x-axis
% line(xlim, [0 0]);  %y-axis
% ylim([-0.2 0.2])
% lgd = legend(p,cc,'Location','Southeast');
% lgd.Title.String='Value of TR';
% grid minor
% xlabel('y_c','Fontsize',16)
% ylabel('chord','Fontsize',16)
% saveas(fig7, 'tap_wing','epsc')
% 
% for i =1:length(TR)
%     c_bar(i) = mean(chord(i,:)); %mean c for each TR value
% end
% 
% %Distribution gamma
% %alpha_i
% %Cl et Cdi 
% for i =1:length(TR)
%     C = diag(chord(i,:));
%     GAM3(:,i) = inv(inv(C)/(pi*Vinf)+A./Vinf)*(degtorad((4+AOA))*ones(n-1,1));
%     GAM_tilde(:,i) = GAM3(:,i)/(Vinf*c_bar(i));
%     Cl(:,i)=2*GAM3(:,i)./(Vinf*chord(i,:)');
%     Cdi(:,i)=Cl(:,i).*A*GAM3(:,i)./Vinf;
%     CDi(i)=1/S*sum(Cdi(:,i).*chord(i,:)'.*dy');
%     CL3(:,i) = 1/S*sum(Cl(:,i).*chord(i,:)'.*dy');
%     ALPHA_i3(:,i) = A*GAM3(:,i)./Vinf;
% end
% 
% fig8 = figure(8)
% hold on 
% for i =1:length(TR)
%     plot(yc,GAM_tilde(:,i))
%     cc(i) = cellstr(num2str(TR(i)));
% end
% %plot(yc,GAM_rect_tilde(:,2),'*g')
% hold off
% lgd=legend(cc,'Location','Southeast');
% lgd.Title.String='Value of TR';
% grid minor
% xlabel('y_c','Fontsize',16)
% ylabel( '$ \tilde \Gamma $','interpreter','latex','Fontsize',16)
% title('Distribution of the dimensionless circulation','Fontsize',16)
% saveas(fig8, 'gamm_tild','epsc')
% 
% fig9 = figure(9)
% subplot(1,2,1)
% hold on 
% for i =1:length(TR)
%     plot(yc,Cl(:,i))
%     cc(i) = cellstr(num2str(TR(i)));
% end
% hold off
% lgd=legend(cc,'Location','Southeast');
% lgd.Title.String='Value of TR';
% grid minor
% xlabel('y_c','Fontsize',16)
% ylabel('C_l','Fontsize',16)
% title('C_l','Fontsize',16)
% 
% 
% subplot(1,2,2)
% hold on
% for i =1:length(TR)
%     plot(yc,Cdi(:,i))
%     cc(i) = cellstr(num2str(TR(i)));
% end
% hold off
% lgd=legend(cc,'Location','Southeast');
% lgd.Title.String='Value of TR';
% grid minor
% xlabel('y_c','Fontsize',16)
% ylabel('C_d_i','Fontsize',16)
% title('C_d_i','Fontsize',16)
% saveas(fig9, 'Clcdi','epsc')
% 
% fig10 = figure(10)
% hold on
% for i =1:length(TR)
%     plot(yc,radtodeg(ALPHA_i3(:,i)))
%     cc(i) = cellstr(num2str(TR(i)));
% end
% hold off
% lgd=legend(cc);
% lgd.Title.String='Value of TR';
% grid minor
% title('Induced angle of attack','Fontsize',16)
% xlabel('y_c','Fontsize',16)
% ylabel('\alpha_i [°]','Fontsize',16)
% saveas(fig10, 'alpha_tild','epsc')
% 
% %% Q4
% 
% AOA=6;
% AR=6
% 
% GAM_max_elliptic=2*b*Vinf*(degtorad(4+AOA))/(1+AR/2);
% CDi_elliptic=pi*GAM_max_elliptic^2./(4*S*Vinf^2);
% 
% CDi_ell = S/(pi*b^2)*CL3.^2; %make sure same CL as tapered
% 
% Delta_CDi=(CDi-CDi_ell)./CDi_ell*100;
% 
% fig11 = figure(11)
% hold on
% for i =1:2
%     p(i) = plot(TR,Delta_CDi)
%     cc(i) = cellstr(num2str(AR));
% end
% lgd=legend(p(2),cc,'Location','Northeast');
% lgd.Title.String='Value of AR';
% hold off
% grid minor
% xlabel('TR','Fontsize',16)
% ylabel('\Delta C_D_i','Fontsize',16)
% title('\Delta C_D_i - AR = 6','Fontsize',16)
% saveas(fig11, 'delta','epsc')
% 
% %% Q5
% 
% AR = 6;
% AOA = 4;
% c0 = b./AR;
% S = b^2./AR;
% Theta=[-4 -2 0 2 4];
% B=ones(n-1,length(Theta))*degtorad(4+AOA);
% 
% for k=1:(n-1)
%     for j=1:length(Theta)
%         twist=[0 Theta(j)];
%         yc_twist=[0 max(yc)];
%         B(k,j)=B(k,j)+degtorad(interp1(yc_twist,twist,abs(yc(k))));
%     end
% end
% 
% for j=1:length(Theta)
%     GAM5(:,j) = inv(eye(n-1)./(pi*Vinf*c0)+ A./Vinf)*B(:,j);
%     WI5 (:,j)= A*GAM5(:,j);
%     ALPH_i5(:,j)= WI5(:,j)./Vinf;
%     Cl5(:,j)=2*GAM5(:,j)./(Vinf*c0);
%     Cdi5(:,j)=Cl(:,j).*A*GAM5(:,j)./Vinf;
%     CL5 (j)= 2./(Vinf.*S)*sum(GAM5(:,j).*dy');
%     CD5(j)= 2./(Vinf.*S)*sum(ALPH_i5(:,j).*GAM5(:,j).*dy');
% end
% 
% fig12 = figure(12)
% subplot(2,2,1)
% hold on 
% for i =1:length(Theta)
%     plot(yc,GAM5(:,i))
%     cc(i) = cellstr(num2str(Theta(i)));
% end
% hold off
% lgd=legend(cc,'Location','Southeast');
% lgd.Title.String='Value of \theta';
% grid minor
% xlabel('y_c')
% ylabel( '$ \Gamma $','interpreter','latex','Fontsize',14)
% title('Distribution of the circulation','Fontsize',14)
% 
% subplot(2,2,2)
% hold on 
% for i =1:length(Theta)
%     plot(yc,ALPH_i5(:,i))
%     cc(i) = cellstr(num2str(Theta(i)));
% end
% hold off
% lgd=legend(cc,'Location','Northeast');
% lgd.Title.String='Value of \theta';
% grid minor
% xlabel('y_c','Fontsize',14)
% ylabel('\alpha_i','Fontsize',14)
% title('\alpha_i','Fontsize',14)
% 
% 
% subplot(2,2,3)
% hold on 
% for i =1:length(Theta)
%     plot(yc,Cl5(:,i))
%     cc(i) = cellstr(num2str(Theta(i)));
% end
% hold off
% lgd=legend(cc,'Location','Southeast');
% lgd.Title.String='Value of \theta';
% grid minor
% xlabel('y_c','Fontsize',14)
% ylabel('Cl','Fontsize',14)
% title('Cl','Fontsize',14)
% 
% subplot(2,2,4)
% hold on 
% for i =1:length(Theta)
%     plot(yc,Cdi5(:,i))
%     cc(i) = cellstr(num2str(Theta(i)));
% end
% hold off
% lgd=legend(cc,'Location','Southeast');
% lgd.Title.String='Value of \theta';
% grid minor
% xlabel('y_c','Fontsize',14)
% ylabel('C_d_i','Fontsize',14)
% title('Cdi','Fontsize',14)
% saveas(fig12, 'twist','epsc')

