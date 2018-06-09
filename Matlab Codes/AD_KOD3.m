clc; clear all; close all
format long
% Constants %
i=1;
J2      = 0; %1.0826359e-3;
Re      = 6378;
% Initial Values %
% Quaternions
q = [1 0 0 0];
%Euler Angles in rad %
Euler = quatern2euler(q);
Roll = Euler(1);
Pitch = Euler(2);
Yaw = Euler(3);
% Angular Velocities  in rad/sec%
wx(1)   = 0.5*pi/180;
wy(1)   = 0.5*pi/180;
wz(1)   = 0.5*pi/180;
W(:,i)  = [wx(i),wy(i),wz(i)];                                             % Angular Velocities Vector %
% Positions in km %
rx(1)   = -5634.985593;
ry(1)   = 3944.062506;
rz(1)   = 8.769715;
R(:,i)  = [rx(i),ry(i),rz(i)];                                             % Position Vector %
% Velocities in km/sec %
vx(1)   = 0.574218;
vy(1)   = 0.803617;
vz(1)   = 7.548263;
V(:,i)  = [vx(i),vy(i),vz(i)];                                             % Velocity Vector %
% RK Conditions %
f       = @RKINPUTm;
y0      = [[Roll; Pitch; Yaw]; W; R; V];
% 4th Order Runge Kutta Calculations %
dt_sim  = 1;                                                               % Simulation Time  %
step      = 0.1;                                                             %
tf      = dt_sim; 
t0      = 0;
m       = 3000/dt_sim;
ySave   = zeros(12,1);
Coe     = zeros(7,m);
t_gost     = tf:tf:tf*m;
Mp   = 0.5;                                                               % Maximum Overshoot %
Ts   = 100;                                                                 % Settling Time %
[Kp,Kd]    = Control(Ts,Mp);                                               % Control function %
for i          = 1:m
    t          = (i-1)*dt_sim;
    ySave(:,i) = y0;
    dv         = [0;0;0];
    % Torque Control %
    if t < 700
         Ref_Angle  = [40*sin((t_gost(i)*pi/150));0;1*sin((t_gost(i)*pi/150))]*pi/180;
    elseif t> 700 && t <= 1500
         Ref_Angle  = [40;0;1]*pi/180;
    elseif t > 1500 
        Ref_Angle  = [40*sin((t_gost(i)*pi/150));0;1*sin((t_gost(i)*pi/150))]*pi/180;
    end
    Ref_AnVel  = [0;0;0];                                                              % Referance Angular Velocity %
    trqx        = (Ref_Angle(1)-y0(1))*Kp(1,1)+((Ref_AnVel(1)-y0(4)).*Kd(1,1));        % Toque Function %
    trqy        = (Ref_Angle(2)-y0(2))*Kp(2,2)+((Ref_AnVel(2)-y0(5)).*Kd(2,2));
    trqz        = (Ref_Angle(3)-y0(3))*Kp(3,3)+((Ref_AnVel(3)-y0(6)).*Kd(3,3));
    trq         = [trqx;trqy;trqz];
    trqMax = 0.2;
    for k = 1:3
        if abs(trq(k)) > trqMax
            trq(k) = sign(trq(k))*trqMax;
        end
     end
     tDist = [1e-1*sin(t/3000*pi); 1e-1*cos(t/3000*pi); -1e-1*sin(t/3000*pi)];
%      tDist = [0;0;0];
     trq   = trq+tDist;
     trqSave(:,i)=trq;
     if t>700 && t<1500
        v_uydu = y0(10:12);
        v_uydu = v_uydu/norm(v_uydu);
        dv     = v_uydu*0.000008;
    end
    [T,Y]      = RK4m(f,t0,tf,y0,step,dv,trq);
    y0         = Y(:,end);  
    Coe(:,i)   = Kepler(y0(7:9),y0(10:12)); 
end
% Subplot Graphs %
% --------------- %0
figure(1)
% Plot3D%
earth_sphere;
grid on
hold on
plot3(ySave(7,:),ySave(8,:),ySave(9,:),'Color','r','linewidth',3);
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');
title('Orbit Figure');
% Graphs Euler Angler
figure(2)
subplot(3,1,1),plot(t_gost,ySave(1,:)*180/pi)
title('Roll')
legend('Roll')
subplot(3,1,2),plot(t_gost,ySave(2,:)*180/pi)
title('Pitch')
legend('Pitch')
subplot(3,1,3),plot(t_gost,ySave(3,:)*180/pi)
title('Yaw')
legend('Yaw')
% Graphs Angular Velocity %
figure(3)
subplot(3,1,1),plot(t_gost,ySave(4,:)*180/pi)
title('Wx')
legend('Wx')
subplot(3,1,2),plot(t_gost,ySave(5,:)*180/pi)
title('Wy')
legend('Wy')
subplot(3,1,3),plot(t_gost,ySave(6,:)*180/pi)
title('Wz')
legend('Wz')
% Graphs Kepler Elements %
figure(4)
subplot(6,1,1),plot(t_gost,Coe(7,:))
title('Semimajor Axis')
legend('a')
subplot(6,1,2),plot(t_gost,Coe(2,:))
title('Eccentricity')
legend('e')
subplot(6,1,3),plot(t_gost,Coe(4,:))
title('Inclination')
legend('incl')
subplot(6,1,4),plot(t_gost,Coe(3,:))
title('RAAN')
legend('RA')
subplot(6,1,5),plot(t_gost,Coe(5,:))
title('Argument of Perigee')
legend('Wa')
subplot(6,1,6),plot(t_gost,Coe(6,:))
title('True Anomaly')
legend('TA')
figure(5)
subplot(3,1,1),plot(t_gost,trqSave(1,:))
title('Torque x')
legend('Torque x')
subplot(3,1,2),plot(t_gost,trqSave(2,:))
title('Torque y')
legend('Torque y')
subplot(3,1,3),plot(t_gost,trqSave(3,:))
title('Torque z')
legend('Torque z')