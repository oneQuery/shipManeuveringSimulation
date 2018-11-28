clc;
clear all;
close all;

%Turning 실험조건
d_right_deg=5; d_left_deg=-d_right_deg; %타각설정
d_r_rad=d_right_deg*pi/180; d_l_rad=d_left_deg*pi/180;
cycle=10;
StopAngle=360*cycle; %Stop Heading Angle(degree)

%타임도메인
SimulTime=1000;
N=10000;
dt=SimulTime/N;

%제원
rho=1025;
disp=24801; %(ton)
m=disp*1000; %(kg)
draft=9.5;
U=10.29; %20knot
L=175;
k=0.48*L; Iz=m*(k^2); %k값 산정방법에 따라 결과 바뀜

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

%Turning 변위
for i=2:N
    
    t(i)=dt*(i-1);
    v_(i-1)=v(i-1)/U; 
    r_(i-1)=r(i-1)/(U/L);
            
    vdot(i)=((-(m_+mx_)*U_*r_(i-1)+Yv_*v_(i-1)+Yr_*r_(i-1)+Yd_*d_r_rad)/(m_+my_))*((U^2)/L);
    rdot(i)=((Nv_*v_(i-1)+Nr_*r_(i-1)+Nd_*d_r_rad)/(Iz_+Jz_))*((U^2)/(L^2));
    
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
    
    X0(i)=X(i)/L;
    Y0(i)=Y(i)/L;

    if psi(i)>StopAngle*pi/180;
        break;
    end
    
end

%결과
figure(1); plot(Y,X,'r','LineWidth',3); axis([-100, 2000, -100, 2000]);
set(gca,'XTick',[0:500:2000]);
set(gca,'YTick',[0:500:2000]);
grid on; xlabel('Y0'); ylabel('X0'); hold on;
fprintf('좌선회\nAD/L=%4.2f\n', max(X0))
fprintf('DT/L=%4.2f\n\n', max(Y0))

figure(2); plot(t,r,'r','LineWidth',3); axis([0, 1000, 0, 0.015]);
set(gca,'XTick',[0:200:1000]);
set(gca,'YTick',[0:0.005:0.015]);
grid on; xlabel('t[s]'); ylabel('r[rad/s]');


% for i=2:N
%     
%     t(i)=dt*i;
%     v_(i-1)=v(i-1)/U; 
%     r_(i-1)=r(i-1)/(U/L);
%             
%     vdot(i)=((-(m_+mx_)*U_*r_(i-1)+Yv_*v_(i-1)+Yr_*r_(i-1)+Yd_*d_l_rad)/(m_+my_))*((U^2)/L);
%     rdot(i)=((Nv_*v_(i-1)+Nr_*r_(i-1)+Nd_*d_l_rad)/(Iz_+Jz_))*((U^2)/(L^2));
%     
%     v(i)=v(i-1)+vdot(i-1)*dt;
%     r(i)=r(i-1)+rdot(i-1)*dt;
%     psi(i)=psi(i-1)+r(i-1)*dt+(1/2)*rdot(i-1)*(dt^2);
%     
%     u(i)=sqrt((U^2)-(v(i)^2));
%     udot(i)=(u(i)-u(i-1))/dt;
%     
%     Xdot(i)=u(i)*cos(psi(i))-v(i)*sin(psi(i));
%     Ydot(i)=u(i)*sin(psi(i))+v(i)*cos(psi(i));
%         
%     Xdotdot(i)=-vdot(i)*sin(psi(i));
%     Ydotdot(i)=vdot(i)*cos(psi(i));
%     
%     Xdotdot(i)=udot(i)*cos(psi(i))-vdot(i)*sin(psi(i));
%     Ydotdot(i)=udot(i)*sin(psi(i))+vdot(i)*cos(psi(i));
%     
%     X(i)=X(i-1)+Xdot(i-1)*dt+(1/2)*Xdotdot(i)*(dt^2);
%     Y(i)=Y(i-1)+Ydot(i-1)*dt+(1/2)*Ydotdot(i)*(dt^2);
%     
%     X0(i)=X(i)/L;
%     Y0(i)=Y(i)/L;
% 
%     if psi(i)<-StopAngle*pi/180;
%         break;
%     end
% end
% 
% %결과
% figure(1); plot(Y0,X0); grid on;
% fprintf('우선회\nAD/L=%4.2f\n', max(X0))
% fprintf('DT/L=%4.2f\n', abs(min(Y0)))






