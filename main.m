clear all
close all 
clc


%% data
m = 6; %kg
g= 9.81;
rho = 1.225; %rho_air`
E=650000; %J/kg

%% Q1 

V_tip = 170; %Mach = 0.4
N_b = [ 2 3 4 ];
%chord = [1:1:10]*10^(-2); %m
R = [10:5:90]*10^(-2); %m
CoR = [0.05:0.05:0.2];
omega = V_tip./R;
%omega = [200:50:450]; %rad/s
kappa = 1.2; %cf litterature
C_d0 = 0.011;
T= m*g;
mini_P = inf;

for i1=1:length(N_b)
    for i2=1:length(CoR)
        for i3=1:length(R)
            for i4 = 1:length(omega)
                chord = CoR(i2)*R(i3);
                %sigma = N_b(i1)*chord(i2)/(pi*R(i3));
                %C_T = m*g/(rho*pi*R(i3)^2*omega(i4)^2*R(i3)^2);
                %C_P{i1,i2}(i3,i4) = kappa*C_T^(3/2)/sqrt(2) + sigma*C_d0/8;
                P_ideal = T^(3/2)/sqrt(2*rho*pi*R(i3)^2);
                P_0 = 1/8*rho*N_b(i1)*omega(i4)^3*chord*C_d0*R(i3)^4;
                P{i1,i2}(i3,i4) = kappa*P_ideal + P_0;
                
                if P{i1,i2}(i3,i4) < mini_P
                    mini_P = P{i1,i2}(i3,i4);
                    result = [N_b(i1) CoR(i2) R(i3) omega(i4)];
                    index = [i1 i2 i3 i4];
                end
            end
        end
    end
end

N_b = result(1);
CoR = result(2);
%R = result(3);
omega = result(4);

   
figure()
plot(R,P{index(1),index(2)}(:,index(4)))
xlabel('Radius of the rotor')
ylabel('Power required')
title('Influence of the value of R on the power required to hover')

R = result(3);
chord = CoR*R;
time_hover=E*(m-4)/mini_P;
%% Q2

m_2 = [6:1:20]; %total mass with added batteries
clear P

for i=1:length(m_2)
    T= m_2(i)*g;
    
    P_ideal = T^(3/2)/sqrt(2*rho*pi*R^2);
    P_0 = 1/8*rho*N_b*omega^3*chord*C_d0*R^4;
    P(i) = kappa*P_ideal + P_0; %required power
    
    P_available(i) = 500*(m_2(i)-4); %mass of batteries
    time_hover_2(i)=E*(m_2(i)-4)/P(i);

end

figure()
plot(m_2-4,P_available-P)
title('Difference between the available and the required power')
xlabel('Number of batteries')
ylabel('Power difference [W]')

figure()
plot(m_2-4,time_hover_2)
title('Hovering time')
xlabel('Number of batteries')
ylabel('Hovering time [s]')
        

%% Q3 

N = 20;
dr = R/N;
r = [0:dr:R];

%alpha distrib = hyper stylé si on calcule avec assign 1 

%Cl et CD en fonction de alpha : http://airfoiltools.com/airfoil/details?airfoil=n0012-il

