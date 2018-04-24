clear all
close all 
clc


%% data
m = 6; %kg
g= 9.81;
rho = 1.225; %rho_air`
E=650000; %J/kg

%% Q1 

V_tip = 100; %Mach = 0.4
N_b = [ 2 3 4 ];
%chord = [1:1:10]*10^(-2); %m
R = [10:5:40]*10^(-2); %m
CoR = [0.05:0.05:0.2];
omega = V_tip./R;
%omega = [200:50:450]; %rad/s
kappa = 1.2; %cf litterature
C_d0 = 0.011;
T= m*g; %per rotor
mini_P = inf;

for i1=1:length(N_b)
    for i2=1:length(CoR)
        for i3=1:length(R)
            for i4 = 1:length(omega)
                chord = CoR(i2)*R(i3);
                %sigma = N_b(i1)*chord(i2)/(pi*R(i3));
                %C_T = m*g/(rho*pi*R(i3)^2*omega(i4)^2*R(i3)^2);
                %C_P{i1,i2}(i3,i4) = kappa*C_T^(3/2)/sqrt(2) + sigma*C_d0/8;
                P_ideal = T^(3/2)/sqrt(2*rho*pi*R(i3)^2); %/rotor
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

time_hover=E*(m-4)/(4*mini_P); %for all the rotors

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
dr = R/(N-1);
r = [0:dr:R]./R;
sigma = N_b*chord/(pi)./r;

%Data from source
load('data.txt');
alpha_ref = degtorad(data(:,1));
Cl_ref = data(:,2);
Cd_ref = data(:,3);

AOA_0 = degtorad(3); %10deg
AOA1 = AOA_0*ones(1,N); %constant
AOA2 = AOA_0.*linspace(1,0.2,N); %linear
%twist3 = ; %ideal ?

AOA3 = AOA2(end)./r; %not AOA_0 but theta_tip

AOA3(1) = AOA2(end) + degtorad(1);
AOA3_lin = (AOA3(end)-AOA3(1))*r + AOA3(1);
for i =1:length(AOA3)
   if AOA3(i) > AOA3_lin(i)
       AOA3(i) = AOA3_lin(i);
   end
end

figure()
plot(r,AOA3)
hold on 
plot(r,AOA2(end)./r);

% Cl and Cd values from http://airfoiltools.com/

Cl1 = interp1(alpha_ref,Cl_ref,AOA1);
Cl2 = interp1(alpha_ref,Cl_ref,AOA2);
Cl3 = interp1(alpha_ref,Cl_ref,AOA3);
Cd1 = interp1(alpha_ref,Cd_ref,AOA1);
Cd2 = interp1(alpha_ref,Cd_ref,AOA2);
Cd3 = interp1(alpha_ref,Cd_ref,AOA3);

C_Treq = T/(rho*pi*R^4*omega^2);
C_T = 0
F = 1;
max_it = 100; %maximum number of iterations

for i =1:max_it
    
    