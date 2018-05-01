clear all
close all
clc

R=0.33;
A=pi*R^2;
rho=1.225; %kg/m^3

b=0.4; %Span (m)           
N=101;

alpha0 = -4; %deg
n = 101;


Vinf= [1:1:60];; %dimensionless
theta = linspace(pi,0,n);
y = b/2*cos(theta); 

%% Varying AOA

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


for j=1:length(Vinf) %loop over AR
    GAM.A1(:,j) = inv(eye(n-1)./(pi*Vinf(j)*c0)+ A./Vinf(j))*B(:,1);
    WI.A1(:,j) = A*GAM.A1(:,j);
    ALPH_i.A1(:,j) = WI.A1(:,j)./Vinf(j);
    CL(1,j) = 2./(Vinf(j).*S)*sum(GAM.A1(:,j).*dy');
    CD(1,j) = 2./(Vinf(j).*S)*sum(ALPH_i.A1(:,j).*GAM.A1(:,j).*dy');
    L(j,1) =  2*(rho*Vinf(j))*sum(GAM.A1(:,j).*dy');
    D(j,1) = 2*(rho*Vinf(j))*sum(ALPH_i.A1(:,j).*GAM.A1(:,j).*dy');
end


for j=1:length(Vinf) %loop over AR
    GAM.A2(:,j) = inv(eye(n-1)./(pi*Vinf(j)*c0)+ A./Vinf(j))*B(:,2);
    WI.A2(:,j) = A*GAM.A2(:,j);
    ALPH_i.A2(:,j) = WI.A2(:,j)./Vinf(j);
    CL(2,j) = 2./(Vinf(j).*S)*sum(GAM.A2(:,j).*dy');
    CD(2,j) = 2./(Vinf(j).*S)*sum(ALPH_i.A2(:,j).*GAM.A2(:,j).*dy');
    L(j,2) =  2*(rho*Vinf(j))*sum(GAM.A2(:,j).*dy');
    D(j,2) = 2*(rho*Vinf(j))*sum(ALPH_i.A2(:,j).*GAM.A2(:,j).*dy');
end


for j=1:length(Vinf) %loop over AR
    GAM.A3(:,j) = inv(eye(n-1)./(pi*Vinf(j)*c0)+ A./Vinf(j))*B(:,3);
    WI.A3(:,j) = A*GAM.A3(:,j);
    ALPH_i.A3(:,j) = WI.A3(:,j)./Vinf(j);
    CL(3,j) = 2./(Vinf(j).*S)*sum(GAM.A3(:,j).*dy');
    CD(3,j) = 2./(Vinf(j).*S)*sum(ALPH_i.A3(:,j).*GAM.A3(:,j).*dy');
    L(j,3) =  2*(rho*Vinf(j))*sum(GAM.A3(:,j).*dy');
    D(j,3) = 2*(rho*Vinf(j))*sum(ALPH_i.A3(:,j).*GAM.A3(:,j).*dy');
end
% 
% fig1 = figure(1)
% for i=1:length(AR)
%     subplot(2,2,i)
%     plot(yc,radtodeg(ALPH_i.A1(:,i)))
%     hold on 
%     plot(yc,radtodeg(ALPH_i.A2(:,i)))
%     plot(yc,radtodeg(ALPH_i.A3(:,i)))
%     legend('AOA = 0°','AOA = 5°','AOA = 10°')
%     hold off
%     title(['AR = ' num2str(AR(i)) ' (rectangular wing)'],'Fontsize',14)
%     xlabel('y_c','Fontsize',14)
%     ylabel('\alpha_i [°]','Fontsize',14)
% end

% 
% fig2 = figure(2)
% subplot(1,3,1) %AOA1
% for i = 1:length(AR)
%     hold on
%     plot(yc,radtodeg(ALPH_i.A1(:,i)))
%     cc(i) = cellstr(num2str(AR(i)));
% end
% hold off
% lgd = legend(cc);
% lgd.Title.String='Value of AR';
% lcn='bestoutside';
% grid minor
% xlabel('y_c','Fontsize',14)
% ylabel('\alpha_i [°]','Fontsize',14)
% title('AOA = 0°','Fontsize',14)

% subplot(1,3,2) %AOA2
% for i = 1:length(AR)
%     hold on
%     plot(yc,radtodeg(ALPH_i.A2(:,i)))
%     cc(i) = cellstr(num2str(AR(i)));
% end
% hold off
% lgd = legend(cc);
% lgd.Title.String='Value of AR';
% lcn='bestoutside';
% grid minor
% xlabel('y_c','Fontsize',14)
% ylabel('\alpha_i [°]','Fontsize',14)
% title('AOA = 5°','Fontsize',14)
% 
% subplot(1,3,3) %AOA3
% for i = 1:length(AR)
%     hold on
%     plot(yc,radtodeg(ALPH_i.A3(:,i)))
%     cc(i) = cellstr(num2str(AR(i)));
% end
% hold off
% lgd = legend(cc);
% lgd.Title.String='Value of AR';
% lcn='bestoutside';
% grid minor
% xlabel('y_c','Fontsize',14)
% ylabel('\alpha_i [°]','Fontsize',14)
% title('AOA = 10°','Fontsize',14)
% saveas(fig2, 'alpha_Q1','epsc')

% figure()
% subplot(1,3,1)
% plot(Vinf,L(:,1))
% 
% subplot(1,3,2)
% plot(Vinf,L(:,2))
% 
% subplot(1,3,3)
% plot(Vinf,L(:,3))
% 
% figure()
% subplot(1,3,1)
% plot(Vinf,D(:,1))
% 
% subplot(1,3,2)
% plot(Vinf,D(:,2))
% 
% subplot(1,3,3)
% plot(Vinf,D(:,3))

%% Varying AR

AR = [4 6 8 10];
AOA = 5;
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

for i =1:length(Vinf)
    for j=1:length(AR) %loop over AR
        GAM.A1(:,j) = inv(eye(n-1)./(pi*Vinf(i)*c0(j))+ A./Vinf(i))*B(:,1);
        WI.A1(:,j) = A*GAM.A1(:,j);
        ALPH_i.A1(:,j) = WI.A1(:,j)./Vinf(i);
        CL2(i,:) = 2./(Vinf(i).*S)*sum(GAM.A1(:,j).*dy');
        CD2(i,:) = 2./(Vinf(i).*S)*sum(ALPH_i.A1(:,j).*GAM.A1(:,j).*dy');
        L2(i,j) =  2*(rho*Vinf(i))*sum(GAM.A1(:,j).*dy');
        D2(i,j) = 2*(rho*Vinf(i))*sum(ALPH_i.A1(:,j).*GAM.A1(:,j).*dy');
    end
end



% fig1 = figure(1)
% for i=1:length(AR)
%     subplot(2,2,i)
%     plot(yc,radtodeg(ALPH_i.A1(:,i)))
%     hold on 
%     plot(yc,radtodeg(ALPH_i.A2(:,i)))
%     plot(yc,radtodeg(ALPH_i.A3(:,i)))
%     legend('AOA = 0°','AOA = 5°','AOA = 10°')
%     hold off
%     title(['AR = ' num2str(AR(i)) ' (rectangular wing)'],'Fontsize',14)
%     xlabel('y_c','Fontsize',14)
%     ylabel('\alpha_i [°]','Fontsize',14)
% end

% 
% fig2 = figure(2)%AOA1
% for i = 1:length(AR)
%     hold on
%     plot(yc,radtodeg(ALPH_i.A1(:,i)))
%     cc(i) = cellstr(num2str(AR(i)));
% end
% hold off
% lgd = legend(cc);
% lgd.Title.String='Value of AR';
% lcn='bestoutside';
% grid minor
% xlabel('y_c','Fontsize',14)
% ylabel('\alpha_i [°]','Fontsize',14)
% title('AOA = 0°','Fontsize',14)

%% Varying b

b = [0.25:0.25:1];
AR = 6;
AOA = 5;
B=ones(n-1,1)*degtorad(AOA+4); %AOA - offset
c0 = b./AR;
S = b.^2./AR;


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

for i =1:length(Vinf)
    for j=1:length(b) %loop over AR
        GAM.A1(:,j) = inv(eye(n-1)./(pi*Vinf(i)*c0(j))+ A./Vinf(i))*B(:,1);
        WI.A1(:,j) = A*GAM.A1(:,j);
        ALPH_i.A1(:,j) = WI.A1(:,j)./Vinf(i);
        CL3(i,j) = 2./(Vinf(i).*S(j))*sum(GAM.A1(:,j).*dy');
        CD3(i,j) = 2./(Vinf(i).*S(j))*sum(ALPH_i.A1(:,j).*GAM.A1(:,j).*dy');
        L3(i,j) =  2*(rho*Vinf(i))*sum(GAM.A1(:,j).*dy');
        D3(i,j) = 2*(rho*Vinf(i))*sum(ALPH_i.A1(:,j).*GAM.A1(:,j).*dy');
    end
end


% fig1 = figure(1)
% for i=1:length(AR)
%     plot(yc,radtodeg(ALPH_i.A1(:,i)))
%     hold on 
%     hold off
%     title(['AR = ' num2str(AR(i)) ' (rectangular wing)'],'Fontsize',14)
%     xlabel('y_c','Fontsize',14)
%     ylabel('\alpha_i [°]','Fontsize',14)
% end

% 
% fig2 = figure(2)
%  %AOA1
% for i = 1:length(AR)
%     hold on
%     plot(yc,radtodeg(ALPH_i.A1(:,i)))
%     cc(i) = cellstr(num2str(AR(i)));
% end
% hold off
% lgd = legend(cc);
% lgd.Title.String='Value of AR';
% lcn='bestoutside';
% grid minor
% xlabel('y_c','Fontsize',14)
% ylabel('\alpha_i [°]','Fontsize',14)
% title('AOA = 0°','Fontsize',14)

%% Changing twist 

AR = 6;
AOA = 5;
b = 0.4;
c0 = b./AR;
S = b^2./AR;
Theta=[-4 -2 0 2 4];
B=ones(n-1,length(Theta))*degtorad(4+AOA);

for k=1:(n-1)
    for j=1:length(Theta)
        twist=[0 Theta(j)];
        yc_twist=[0 max(yc)];
        B(k,j)=B(k,j)+degtorad(interp1(yc_twist,twist,abs(yc(k))));
    end
end

for i =1:length(Vinf)
for j=1:length(Theta)
    GAM5(:,j) = inv(eye(n-1)./(pi*Vinf(i)*c0)+ A./Vinf(i))*B(:,j);
    WI5 (:,j)= A*GAM5(:,j);
    ALPH_i5(:,j)= WI5(:,j)./Vinf(i);
%    Cl5(i,j)=2*GAM5(:,j)./(Vinf(i)*c0);
 %   Cdi5(:,j)=Cl(:,j).*A*GAM5(:,j)./Vinf(i);
    CL5 (i,j)= 2./(Vinf(i).*S)*sum(GAM5(:,j).*dy');
    CD5(i,j)= 2./(Vinf(i).*S)*sum(ALPH_i5(:,j).*GAM5(:,j).*dy');
    L5(i,j) =  2*(rho*Vinf(i))*sum(GAM5(:,j).*dy');
    D5(i,j) = 2*(rho*Vinf(i))*sum(ALPH_i5(:,j).*GAM5(:,j).*dy');
end
end

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
% 

save('data')


