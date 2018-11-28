clear;
clc;
close all;


%타임도메인
SimulTime=200;
N=10000/2;
dt=SimulTime/N;

%제원
rho=1025;
disp=24801; %(ton)
m=disp*1000; %(kg)
draft=9.5;
U=10.29; %20knot
L=175;
k=0.48*L; Iz=m*(k^2); 
%k값 산정방법에 따라 그래프 바뀜
%통상값:0.25, 그래프와 비슷한값:0.42

%미계수
U_=U/U;
m_=m/((1/2)*rho*(L^2)*draft);
mx_=0.0049;
my_=0.1451;
Yv_=-0.2478;
Yr_=0.0605;
Yd_=-0.0531;
Nv_=-0.0791;
Nr_=-0.05;
Nd_=0.0259;
Iz_=Iz/((1/2)*rho*(L^4)*draft);
Jz_=0.0086;

%초기값
u(1)=0; udot(1)=0;
v(1)=0; vdot(1)=0; v_(1)=0;
r(1)=0; rdot(1)=0; r_(1)=0;
psi(1)=0;
X(1)=0; Xdot(1)=0; Xdotdot(1)=0; 
Y(1)=0; Ydot(1)=0; Ydotdot(1)=0;
t(1)=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ZigZag 실험조건
d0_dgr=10; %(degree)
d0=10*pi/180; %(radian)
t_RudChange=1;
d(1)=d0;%(추가)
d00=d0;
d_Targeting(1)=0;

%ZigZag_Turn Rate:10deg/1sec
for i=2:N
    
    t(i)=dt*(i-1);
    v_(i-1)=v(i-1)/U; 
    r_(i-1)=r(i-1)/(U/L);
    
    if psi(i-1)>d00;
        d(i)=-d0;
        d00=-d0;
    else
        d(i)=d0;
        d00=d0;
    end
    
    if d_Targeting(i-1)<d(i);
        d_Targeting(i)=d_Targeting(i-1)+(d0/t_RudChange)*dt;
    else
        d_Targeting(i)=d_Targeting(i-1)-(d0/t_RudChange)*dt;
    end
    
    vdot(i)=((-(m_+mx_)*U_*r_(i-1)+Yv_*v_(i-1)+Yr_*r_(i-1)+Yd_*d_Targeting(i))/(m_+my_))*((U^2)/L);
    rdot(i)=((Nv_*v_(i-1)+Nr_*r_(i-1)+Nd_*d_Targeting(i))/(Iz_+Jz_))*((U^2)/(L^2));
    
    v(i)=v(i-1)+vdot(i-1)*dt;
    r(i)=r(i-1)+rdot(i-1)*dt;
    psi(i)=psi(i-1)+r(i-1)*dt+(1/2)*rdot(i-1)*(dt^2);
        
    u(i)=sqrt((U^2)-(v(i)^2));
    udot(i)=(u(i)-u(i-1))/dt;
    
    Xdot(i)=u(i)*cos(psi(i))-v(i)*sin(psi(i));
    Ydot(i)=u(i)*sin(psi(i))+v(i)*cos(psi(i));
        
    Xdotdot(i)=udot(i)*cos(psi(i))-vdot(i)*sin(psi(i));
    Ydotdot(i)=udot(i)*sin(psi(i))+vdot(i)*cos(psi(i));
    
    X(i)=X(i-1)+Xdot(i-1)*dt+(1/2)*Xdotdot(i)*(dt^2);
    Y(i)=Y(i-1)+Ydot(i-1)*dt+(1/2)*Ydotdot(i)*(dt^2);
       
end

figure(1); 
plot(t,psi*180/pi,'r:','LineWidth',1.8); grid on; 
axis([0, 200, -35, 35]);
set(gca,'XTick',[0:20:200]);
set(gca,'YTick',[-30:10:30]);
hold on; 
plot(t,d_Targeting*180/pi,'r:','LineWidth',1.8); xlabel('t[s]'); ylabel('Ψ,δ[deg]')
% legend('Ψ','δ');

figure(2);
plot(t,r*180/pi,'r:','LineWidth',1.8); grid on; 
axis([0, 200, -1.4, 1.4]);
set(gca,'XTick',[0:20:200]);
set(gca,'YTick',[-1.4:0.2:1.4]);
hold on;
xlabel('t[s]'); ylabel('r[deg/s]');

%ZigZag_Turn Rate:10deg/3sec
d0_dgr=10; %(degree)
d0=10*pi/180; %(radian)
t_RudChange=3;
d(1)=d0;%(추가)
d00=d0;
d_Targeting(1)=0;
for i=2:N
    
    t(i)=dt*(i-1);
    v_(i-1)=v(i-1)/U; 
    r_(i-1)=r(i-1)/(U/L);
    
    if psi(i-1)>d00;
        d(i)=-d0;
        d00=-d0;
    else
        d(i)=d0;
        d00=d0;
    end
    
    if d_Targeting(i-1)<d(i);
        d_Targeting(i)=d_Targeting(i-1)+(d0/t_RudChange)*dt;
    else
        d_Targeting(i)=d_Targeting(i-1)-(d0/t_RudChange)*dt;
    end
    
    vdot(i)=((-(m_+mx_)*U_*r_(i-1)+Yv_*v_(i-1)+Yr_*r_(i-1)+Yd_*d_Targeting(i))/(m_+my_))*((U^2)/L);
    rdot(i)=((Nv_*v_(i-1)+Nr_*r_(i-1)+Nd_*d_Targeting(i))/(Iz_+Jz_))*((U^2)/(L^2));
    
    v(i)=v(i-1)+vdot(i-1)*dt;
    r(i)=r(i-1)+rdot(i-1)*dt;
    psi(i)=psi(i-1)+r(i-1)*dt+(1/2)*rdot(i-1)*(dt^2);
        
    u(i)=sqrt((U^2)-(v(i)^2));
    udot(i)=(u(i)-u(i-1))/dt;
    
    Xdot(i)=u(i)*cos(psi(i))-v(i)*sin(psi(i));
    Ydot(i)=u(i)*sin(psi(i))+v(i)*cos(psi(i));
        
    Xdotdot(i)=udot(i)*cos(psi(i))-vdot(i)*sin(psi(i));
    Ydotdot(i)=udot(i)*sin(psi(i))+vdot(i)*cos(psi(i));
    
    X(i)=X(i-1)+Xdot(i-1)*dt+(1/2)*Xdotdot(i)*(dt^2);
    Y(i)=Y(i-1)+Ydot(i-1)*dt+(1/2)*Ydotdot(i)*(dt^2);
       
end

figure(1); 
plot(t,psi*180/pi,'r','LineWidth',1.8); grid on; 
axis([0, 200, -35, 35]);
set(gca,'XTick',[0:20:200]);
set(gca,'YTick',[-30:10:30]);
hold on; 
plot(t,d_Targeting*180/pi,'r','LineWidth',1.8); xlabel('t[s]'); ylabel('Ψ,δ[deg]')
% legend('Ψ','δ');

figure(2);
plot(t,r*180/pi,'r','LineWidth',1.8); grid on; 
axis([0, 200, -1.4, 1.4]);
set(gca,'XTick',[0:20:200]);
set(gca,'YTick',[-1.4:0.2:1.4]);
hold on;
xlabel('t[s]'); ylabel('r[deg/s]');

%ZigZag_Turn Rate:10deg/5sec
d0_dgr=10; %(degree)
d0=10*pi/180; %(radian)
t_RudChange=5;
d(1)=d0;%(추가)
d00=d0;
d_Targeting(1)=0;
for i=2:N
    
    t(i)=dt*(i-1);
    v_(i-1)=v(i-1)/U; 
    r_(i-1)=r(i-1)/(U/L);
    
    if psi(i-1)>d00;
        d(i)=-d0;
        d00=-d0;
    else
        d(i)=d0;
        d00=d0;
    end
    
    if d_Targeting(i-1)<d(i);
        d_Targeting(i)=d_Targeting(i-1)+(d0/t_RudChange)*dt;
    else
        d_Targeting(i)=d_Targeting(i-1)-(d0/t_RudChange)*dt;
    end
    
    vdot(i)=((-(m_+mx_)*U_*r_(i-1)+Yv_*v_(i-1)+Yr_*r_(i-1)+Yd_*d_Targeting(i))/(m_+my_))*((U^2)/L);
    rdot(i)=((Nv_*v_(i-1)+Nr_*r_(i-1)+Nd_*d_Targeting(i))/(Iz_+Jz_))*((U^2)/(L^2));
    
    v(i)=v(i-1)+vdot(i-1)*dt;
    r(i)=r(i-1)+rdot(i-1)*dt;
    psi(i)=psi(i-1)+r(i-1)*dt+(1/2)*rdot(i-1)*(dt^2);
        
    u(i)=sqrt((U^2)-(v(i)^2));
    udot(i)=(u(i)-u(i-1))/dt;
    
    Xdot(i)=u(i)*cos(psi(i))-v(i)*sin(psi(i));
    Ydot(i)=u(i)*sin(psi(i))+v(i)*cos(psi(i));
        
    Xdotdot(i)=udot(i)*cos(psi(i))-vdot(i)*sin(psi(i));
    Ydotdot(i)=udot(i)*sin(psi(i))+vdot(i)*cos(psi(i));
    
    X(i)=X(i-1)+Xdot(i-1)*dt+(1/2)*Xdotdot(i)*(dt^2);
    Y(i)=Y(i-1)+Ydot(i-1)*dt+(1/2)*Ydotdot(i)*(dt^2);
       
end

figure(1); 
plot(t,psi*180/pi,'r--','LineWidth',1.8); grid on; 
axis([0, 200, -35, 35]);
set(gca,'XTick',[0:20:200]);
set(gca,'YTick',[-30:10:30]);
hold on; 
plot(t,d_Targeting*180/pi,'r--','LineWidth',1.8); xlabel('t[s]'); ylabel('Ψ,δ[deg]')
% legend(psi,'hide','hide','g');

figure(2);
plot(t,r*180/pi,'r--','LineWidth',1.8); grid on; 
axis([0, 200, -1.4, 1.4]);
set(gca,'XTick',[0:20:200]);
set(gca,'YTick',[-1.4:0.2:1.4]);
hold on;
xlabel('t[s]'); ylabel('r[deg/s]');
