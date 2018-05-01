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
R = [10:1:45]*10^(-2); %m
CoR = [0.05:0.05:0.2];
CoR = 0.1;
omega = V_tip./R;
%omega = [200:50:450]; %rad/s
kappa = 1.2; %cf litterature
C_d0 = 0.011;
T= m*g/4; %per rotor
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
xlabel('Radius of the rotor [m]')
ylabel('Power required [W]')
title('Influence of the value of R on the power required to hover')
enhance_plot('TIMES',14,2,8,0)

R = result(3);
chord = CoR*R;

time_hover=E*(m-4)/(4*mini_P); %for all the rotors

%% Q2

m_2 = [6:1:20]; %total mass with added batteries
clear P

for i=1:length(m_2)
    T= m_2(i)*g/4;
    
    P_ideal = T^(3/2)/sqrt(2*rho*pi*R^2);
    P_0 = 1/8*rho*N_b*omega^3*chord*C_d0*R^4;
    P(i) = kappa*P_ideal + P_0; %required power
    
    P_available(i) = 500*(m_2(i)-4); %mass of batteries
    time_hover_2(i)=E*(m_2(i)-4)/(4*P(i));

end

[time_max, i_max] = max(time_hover_2);
N_max = (m_2(i_max)-4)/2;


figure()
plot((m_2-4)/2,P_available-P)
title('Difference between the available and the required power')
xlabel('Number of batteries')
ylabel('Power difference [W]')
enhance_plot('TIMES',14,2,8,0)

figure()
plot((m_2-4)/2,time_hover_2)
title('Hovering time')
xlabel('Number of batteries')
ylabel('Hovering time [s]')
enhance_plot('TIMES',14,2,8,0)
        

%% Q3 

N = 20;
dr = R/(N-1);
r = [0:dr:R]./R;
sigma = N_b*chord/(pi*R);
T = m*g/4;
C_Treq = T/(4*rho*pi*R^4*omega^2);

theta_tip = 4*C_Treq/(sigma*2*pi) + sqrt(C_Treq/2);

theta_ref = degtorad(1); %4deg??

    
theta(1,:) = theta_ref*ones(1,N); %constant;
theta(2,:) = linspace(theta_ref,theta_tip,N); %linear
%twist3 = ; %ideal ?

theta(3,:) = theta_tip./r; %not AOA_0 but theta_tip

theta_tw = [0 ; (theta_tip - theta_ref);(theta_tip*(1 - 1/0.6)/0.4)]; %pente
 %theta 3/4R



%Simple version for ideal case

theta_simple = theta_tip./r;

%Complicated general version
max_it = 100; %maximum number of iterations

dC_T = zeros(max_it+1,N);
C_T = trapz(dC_T);
F = ones(max_it+1,N);

i1 = 3; %type of twist

theta0(1) = 6*C_Treq./(sigma*2*pi) - 3/4*theta_tw(i1) + 3/2*sqrt(C_Treq/2);
theta3_lin = (theta_tip*(1 - 1/0.6)/0.4)*r + theta0(1); %use that
theta(3,:) = theta3_lin; 

i2=1;
while i2<max_it && abs(C_T(i2)-C_Treq)>10^(-3)   
        theta0(i2+1) = theta0(i2) + (6*(C_Treq - C_T(i2))./(sigma*2*pi) ...
            + 3*sqrt(2)/4*(sqrt(C_Treq)-sqrt(C_T(i2))));
        theta(i1,:) = theta0(i2) + theta_tw(i1).*r;
        theta_75 = theta(:,15); %theta 3/4R
        lambda(i2,:) = sigma*2*pi./(16.*F(i2,:)).*(sqrt(1 +  ... 
            32.*F(i2,:)./(sigma*2*pi).*theta(i1,:).*r) - 1);
        phi = lambda(i2,:)./r;
        f(i2,:) = N_b/2.*(1-r)./(r.*phi);
        F(i2+1,:) = 2/pi .* acos(exp(-f(i2)));
        %dC_T(i2+1,:) = 1/2*sigma*2*pi.*r.^2;
        dC_T(i2+1,:) = 1/2*sigma.*2*pi.*(theta_75(i1)./3 - lambda(i2,:)/2);
        C_T(i2+1) = trapz(r,dC_T(i2+1,:));
        T_2(i2) = (4*rho*pi*R^4*omega^2)*C_T(i2+1);
        
        C_T_sol = C_T(i2+1);
        theta_sol = theta(i1,:);
        T_sol = C_T_sol *pi* 4* rho * omega^2 * R^4;
i2=i2+1;
end

figure()
%plot(r,radtodeg(theta(3,:)))
hold on 
plot(r,radtodeg(theta(2,end)./r));
plot(r,radtodeg(theta(i1,:)));
xlabel('r')
ylabel('Twist [°]')
enhance_plot('TIMES',14,2,8,0)

%theta_ref = theta_ref + degtorad(1);


%Data from source
load('data.txt');
alpha_ref = degtorad(data(:,1));
Cl_ref = data(:,2);
Cd_ref = data(:,3);

% Cl and Cd values from http://airfoiltools.com/
for i =1:3
    Cl(i,:) = interp1(alpha_ref,Cl_ref,theta(i,:));
    Cd(i,:) = interp1(alpha_ref,Cd_ref,theta(i,:));
end

%% Q4

A_front = 0.1^2;
V_inf = [1:1:60];
W = m*g;
D_0 = 1/2*rho*A_front*V_inf.^2*1.28;
alpha = atan(D_0/W);

v_h = sqrt(T_sol*4/(2*rho*pi*R^2)); %hovering for one rotor

for j = 1:length(V_inf)
f= @(v_i) v_i - v_h^2/(sqrt((cos(alpha(j))*V_inf(j))^2 + (V_inf(j)*sin(alpha(j)) + v_i)^2)) ;
V_i(j) = fzero(f,0.001);
end

P = 4*T*(V_inf.*sin(alpha) + V_i);



%close all

figure()
plot(V_inf, P)
hold on 

xlabel('Wind  speed [m/s]')
ylabel('Power required [W]')
enhance_plot('TIMES',14,2,8,0)
plot(V_inf, 1000*ones(1,length(V_inf)),'--')

%% Q5

load('data.mat');
L_rotor = T*cos(alpha);
L_wing = L5;

b=[0.25:0.25:1];
AOA = [0 5 10];
AR = [ 4 6 8 10 ];
twist = [-4:2:4];

D_tot = D_0' + D5;

alpha_tot = atan(D_tot./W);
T_wing = L_wing./cos(alpha_tot);


v_h = sqrt((T_sol + T_wing).*4/(2*rho*pi*R^2)); %hovering for one rotor

for i =1:size(L_wing,2)
    for j = 1:length(V_inf)
    f= @(v_i) v_i - v_h(j,i)^2/(sqrt((cos(alpha_tot(j,i))*V_inf(j))^2 + (V_inf(j)*sin(alpha_tot(j,i)) + v_i)^2)) ;
    V_i_tot(j,i) = fzero(f,0.1);
    P_tot(j,i) = 4*T*(V_inf(j).*sin(alpha_tot(j,i)) + V_i_tot(j,i));
    end
end

figure()
hold on 
for i =1:size(L_wing,2)
    p(i) = plot(V_inf, P_tot(:,i))
    cc(i) = cellstr(num2str(twist(i)));
end
lgd=legend(p,cc,'Location','Northwest');
lgd.Title.String='Value of twist';
xlabel('Wind  speed [m/s]')
ylabel('Power required [W]')
enhance_plot('TIMES',14,2,8,0)
plot(V_inf, 1000*ones(1,length(V_inf)),'k--')


