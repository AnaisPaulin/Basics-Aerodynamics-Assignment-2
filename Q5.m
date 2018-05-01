clear all
close all
clc

%% Varying AR
R=0.13;
A_rotor=pi*R^2;
rho=1.225; %kg/m^3

b=0.4; %Span (m)           
N=101;
AR=[0 2 4 6];              
AOA=6;                      
S=b^2./AR;
c=b./AR;
Vinf=linspace(0.01,55,70); 
alpha0=ones(N-1,1)*degtorad(AOA+4); %AOA - offset     
Theta=linspace(pi,0,N);
y=b/2*cos(Theta); %y axis vector

%Creation of the intervals center vector
for i = 1:(N-1)
    yc(i) = y(i) + (y(i+1)-y(i))/2; 
end

% Matrix A: A*circulation=wi
for z=1:length(Vinf)
for k=1:length(AR)
for i =1:N-1
    for j=2:N-2
        A(i,j,k,z) =1/(4*pi*Vinf(z))*(2*(1/(yc(i)-y(j))-1/(yc(i)-y(j+1))));
        if i==j
            A(i,j,k,z)= A(i,j,k,z)+1/(pi*Vinf(z)*c(k));
        end
    end
    A(i,1,k,z) =1/(4*pi*Vinf(z))*(1/(yc(i)-y(1)) - 2/(yc(i)-y(2)));
    A(i,N-1,k,z) =1/(4*pi*Vinf(z))*(2/(yc(i)-y(N-1)) - 1/(yc(i)-y(N)));
end
end
end

% Gamma, cl,cd,induced alpha:

for i=1:length(Vinf)
for j=1:length(AR) 
    B=inv(A(:,:,j,i))*alpha0;
    wi=(A(:,:,j,i)-eye(N-1)./(pi*Vinf(i)*c(j)))*B;
    induced_AOA=wi/Vinf(i);
    Cl(i,j,:)=2/(Vinf(i)*S(j))*trapz(y(1:end-1), B);
    L(i,j)=2*(rho*Vinf(i)*trapz(y(1:end-1), B));
    Cd(i,j,:)=2/(Vinf(i)*S(j))*trapz(y(1:end-1), B.*induced_AOA);
    D(i,j)=2*(rho*Vinf(i)*trapz(y(1:end-1), B.*induced_AOA));
end
end

A_box=0.1^2;
Nr=4;
g=9.81;
W=6*g;

for j=1:length(AR)
for i=1:length(Vinf)
    Drag(i,j)=(0.5*rho*A_box*Vinf(i)^2+D(i,j))/4;
    alpha(i,j)=abs(atan(Drag(i,j)/(W/4-L(i,j)/4)));
    T5(i,j)=Drag(i,j)/sin(alpha(i,j));
    B(i,j)=2*Vinf(i)*sin(alpha(i,j));
    X(i)=Vinf(i)^2;
    K(i,j)=T5(i,j)^2/(4*rho^2*A_rotor^2);
    p=[1/K(i,j) B(i,j)/K(i,j) X(i)/K(i,j) 0 -1];
    r=roots(p);
    vi_5=max(real(r));
    P5(i,j)=4*(T5(i,j)*Vinf(i)*sin(alpha(i,j))+T5(i,j)*vi_5);
end
end

figure ()
subplot(1,2,1)
plot(Vinf,D(:,1)) %drag two airfoils
hold on
plot(Vinf,D(:,2))
plot(Vinf,D(:,3))
plot(Vinf,D(:,4))
ylabel('Drag');
xlabel('wind speed');
hold off

subplot(1,2,2)
plot(Vinf,L(:,1))
hold on
plot(Vinf,L(:,2))
plot(Vinf,L(:,3))
plot(Vinf,L(:,4))
ylabel('Lift');
xlabel('Wind speed');
hold off

figure()
plot(Vinf,P5(:,1))
hold on
plot(Vinf,P5(:,2))
plot(Vinf,P5(:,3))
plot(Vinf,P5(:,4))
ylabel('Total power (W)');
xlabel('Wind speed (m/s)');
hline = refline([0 1000]);
hline.Color = 'k';
legend('AR=0','AR=2','AR=4','AR=8')
enhance_plot('TIMES',14,1.5,8,0)
hold off


%% Varying AOA
R=0.13;
A_rotor=pi*R^2;
rho=1.225; %kg/m^3

b=0.4; %Span (m)           
N=101;
AR=4;              
AOA=[0 5 10];                      
S=b^2./AR;
c=b./AR;
Vinf=linspace(0.01,55,70); 
alpha0=ones(N-1,1)*degtorad(AOA+4); %AOA - offset     
Theta=linspace(pi,0,N);
y=b/2*cos(Theta); %y axis vector

%Creation of the intervals center vector
for i = 1:(N-1)
    yc(i) = y(i) + (y(i+1)-y(i))/2; 
end

% Matrix A: A*circulation=wi
for z=1:length(Vinf)
for k=1:length(AR)
for i =1:N-1
    for j=2:N-2
        A(i,j,k,z) =1/(4*pi*Vinf(z))*(2*(1/(yc(i)-y(j))-1/(yc(i)-y(j+1))));
        if i==j
            A(i,j,k,z)= A(i,j,k,z)+1/(pi*Vinf(z)*c(k));
        end
    end
    A(i,1,k,z) =1/(4*pi*Vinf(z))*(1/(yc(i)-y(1)) - 2/(yc(i)-y(2)));
    A(i,N-1,k,z) =1/(4*pi*Vinf(z))*(2/(yc(i)-y(N-1)) - 1/(yc(i)-y(N)));
end
end
end

% Gamma, cl,cd,induced alpha:

for i=1:length(Vinf)
for j=1:length(AOA) 
    B(:,j)=inv(A(:,:,j,i))*alpha0(:,j);
    wi=(A(:,:,j,i)-eye(N-1)./(pi*Vinf(i)*c))*B(:,j);
    induced_AOA=wi/Vinf(i);
    Cl(i,j,:)=2/(Vinf(i)*S)*trapz(y(1:end-1), B(:,j));
    L(i,j)=2*(rho*Vinf(i)*trapz(y(1:end-1), B(:,j)));
    Cd(i,j,:)=2/(Vinf(i)*S)*trapz(y(1:end-1), B(:,j).*induced_AOA);
    D(i,j)=2*(rho*Vinf(i)*trapz(y(1:end-1), B(:,j).*induced_AOA));
end
end

A_box=0.1^2;
Nr=4;
g=9.81;
W=6*g;

for j=1:length(AOA)
for i=1:length(Vinf)
    Drag(i,j)=(0.5*rho*A_box*Vinf(i)^2+D(i,j))/4;
    alpha(i,j)=abs(atan(Drag(i,j)/(W/4-L(i,j)/4)));
    T5(i,j)=Drag(i,j)/sin(alpha(i,j));
    B(i,j)=2*Vinf(i)*sin(alpha(i,j));
    X(i)=Vinf(i)^2;
    K(i,j)=T5(i,j)^2/(4*rho^2*A_rotor^2);
    p=[1/K(i,j) B(i,j)/K(i,j) X(i)/K(i,j) 0 -1];
    r=roots(p);
    vi_5=max(real(r));
    P4(i,j)=4*(T5(i,j)*Vinf(i)*sin(alpha(i,j))+T5(i,j)*vi_5);
end
end


figure()
plot(Vinf,P4(:,1))
hold on
plot(Vinf,P4(:,2))
plot(Vinf,P4(:,3))
ylabel('Total power (W)');
xlabel('Wind speed (m/s)');
refline([0 1000]);hline = refline([0 1000]);
hline.Color = 'k';
legend('AOA=0','AOA=5','AOA=10')
enhance_plot('TIMES',14,1.5,8,0)

hold off

%% Varying b
R=0.33;
A_rotor=pi*R^2;
rho=1.225; %kg/m^3

b=[0.25:0.25:1]; %Span (m)           
N=101;
AR=4;              
AOA=5;                      
S=b.^2./AR;
c=b./AR;
Vinf=linspace(0.01,55,70); 
alpha0=ones(N-1,1)*degtorad(AOA+4); %AOA - offset     
Theta=linspace(pi,0,N);

y=b'./2*cos(Theta) ;%y axis vector

%Creation of the intervals center vector
for j=1:length(b)
for i = 1:(N-1)
    yc(j,i) = y(j,i) + (y(j,i+1)-y(j,i))/2; 
end
end

% Matrix A: A*circulation=wi

for z=1:length(Vinf)
for k=1:length(b)
for i =1:N-1
    for j=2:N-2
        A(i,j,k,z) =1/(4*pi*Vinf(z))*(2*(1/(yc(k,i)-y(k,j))-1/(yc(k,i)-y(k,j+1))));
        if i==j
            A(i,j,k,z)= A(i,j,k,z)+1/(pi*Vinf(z)*c(k));
        end
    end
    A(i,1,k,z) =1/(4*pi*Vinf(z))*(1/(yc(k,i)-y(1)) - 2/(yc(k,i)-y(k,2)));
    A(i,N-1,k,z) =1/(4*pi*Vinf(z))*(2/(yc(k,i)-y(k,N-1)) - 1/(yc(k,i)-y(k,N)));
end
end
end

% Gamma, cl,cd,induced alpha:

for i=1:length(Vinf)
for j=1:length(b) 
    B=inv(A(:,:,j,i))*alpha0;
    wi=(A(:,:,j,i)-eye(N-1)./(pi*Vinf(i)*c(j)))*B;
    induced_AOA=wi/Vinf(i);
    Cl(i,j,:)=2/(Vinf(i)*S(j))*trapz(y(j,1:end-1), B);
    L(i,j)=2*(rho*Vinf(i)*trapz(y(j,1:end-1), B));
    Cd(i,j,:)=2/(Vinf(i)*S(j))*trapz(y(j,1:end-1), B.*induced_AOA);
    D(i,j)=2*(rho*Vinf(i)*trapz(y(j,1:end-1), B.*induced_AOA));
end
end

A_box=0.1^2;
Nr=4;
g=9.81;
W=6*g + 0.5*b;

for j=1:length(b)
for i=1:length(Vinf)
    Drag(i,j)=(0.5*rho*A_box*Vinf(i)^2+D(i,j))/4;
    alpha(i,j)=abs(atan(Drag(i,j)/(W(j)/4-L(i,j)/4)));
    T5(i,j)=Drag(i,j)/sin(alpha(i,j));
    B(i,j)=2*Vinf(i)*sin(alpha(i,j));
    X(i)=Vinf(i)^2;
    K(i,j)=T5(i,j)^2/(4*rho^2*A_rotor^2);
    p=[1/K(i,j) B(i,j)/K(i,j) X(i)/K(i,j) 0 -1];
    r=roots(p);
    vi_5=max(real(r));
    P3(i,j)=4*(T5(i,j)*Vinf(i)*sin(alpha(i,j))+T5(i,j)*vi_5);
end
end


figure()
plot(Vinf,P3(:,1))
hold on
plot(Vinf,P3(:,2))
plot(Vinf,P3(:,3))
plot(Vinf,P3(:,4))
ylabel('Total power (W)');
xlabel('Wind speed (m/s)');
hline = refline([0 1000]);
hline.Color = 'k';
legend('b=0.25','b=0.5','b=0.75', 'b=1')
enhance_plot('TIMES',14,1.5,8,0)

hold off

%% Varying twist
R=0.33;
A_rotor=pi*R^2;
rho=1.225; %kg/m^3

b=0.4; %Span (m)           
N=101;
AR=4;              
AOA=5;                      
S=b.^2./AR;
c=b./AR;
Theta=[-4 -2 0 2 4];
Vinf=linspace(0.01,55,70); 
alpha0=ones(N-1,1)*degtorad(AOA+4); %AOA - offset     
Theta=linspace(pi,0,N);
B=ones(N-1,length(Theta))*degtorad(4+AOA);

y=b'./2*cos(Theta) ;%y axis vector

%Creation of the intervals center vector
for j=1:length(b)
for i = 1:(N-1)
    yc(j,i) = y(j,i) + (y(j,i+1)-y(j,i))/2; 
end
end

% Matrix A: A*circulation=wi

for z=1:length(Vinf)
for k=1:length(b)
for i =1:N-1
    for j=2:N-2
        A(i,j,k,z) =1/(4*pi*Vinf(z))*(2*(1/(yc(k,i)-y(k,j))-1/(yc(k,i)-y(k,j+1))));
        if i==j
            A(i,j,k,z)= A(i,j,k,z)+1/(pi*Vinf(z)*c(k));
        end
    end
    A(i,1,k,z) =1/(4*pi*Vinf(z))*(1/(yc(k,i)-y(1)) - 2/(yc(k,i)-y(k,2)));
    A(i,N-1,k,z) =1/(4*pi*Vinf(z))*(2/(yc(k,i)-y(k,N-1)) - 1/(yc(k,i)-y(k,N)));
end
end
end

for k=1:(N-1)
    for j=1:length(Theta)
        twist=[0 Theta(j)]
        yc_twist=[0 max(yc)]
        yc(k)
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
    Cl (i,j)= 2./(Vinf(i).*S)*sum(GAM5(:,j).*dy');
    Cd(i,j)= 2./(Vinf(i).*S)*sum(ALPH_i5(:,j).*GAM5(:,j).*dy');
    L(i,j) =  2*(rho*Vinf(i))*sum(GAM5(:,j).*dy');
    D(i,j) = 2*(rho*Vinf(i))*sum(ALPH_i5(:,j).*GAM5(:,j).*dy');
end
end



A_box=0.1^2;
Nr=4;
g=9.81;
W=6*g;

for j=1:length(twist)
for i=1:length(Vinf)
    Drag(i,j)=(0.5*rho*A_box*Vinf(i)^2+D(i,j))/4;
    alpha(i,j)=abs(atan(Drag(i,j)/(W/4-L(i,j)/4)));
    T5(i,j)=Drag(i,j)/sin(alpha(i,j));
    B(i,j)=2*Vinf(i)*sin(alpha(i,j));
    X(i)=Vinf(i)^2;
    K(i,j)=T5(i,j)^2/(4*rho^2*A_rotor^2);
    p=[1/K(i,j) B(i,j)/K(i,j) X(i)/K(i,j) 0 -1];
    r=roots(p);
    vi_5=max(real(r));
    P3(i,j)=4*(T5(i,j)*Vinf(i)*sin(alpha(i,j))+T5(i,j)*vi_5);
end
end


figure()
plot(Vinf,P3(:,1))
hold on
plot(Vinf,P3(:,2))
plot(Vinf,P3(:,3))
plot(Vinf,P3(:,4))
ylabel('Total power (W)');
xlabel('Wind speed (m/s');
hline = refline([0 1000]);
hline.Color = 'k';
legend('b=0.25','b=0.5','b=0.75', 'b=1')
enhance_plot('TIMES',14,1.5,8,0)

hold off