clc; clear all; close all
% Mustafa Çağatay ŞAHİN %
% 15/Summer Intern in TAI %
% Created : 02.09.2015 in TAI-TUSAŞ Turkish Aerospace Industries. INC %
% Introduction %
% This code can be used determination of Attitude, velocity, and other specific
% information about satellite
% Information about Constant;
% Mu : Earth's gravitional constant.
% J2 : Coefficent.
% Re : Radius of Earth.
% Roll,Pitch,Yaw : Euler Angles of Satellite in rad.
% Wx,Wy,Wz : Angular velocities in satellite in rad/sec.
% Rx,Ry,Rz : Position of satellite on ECI in km.
% Vx,Vy,Vz : Velocities of satellite in km/sec.
% y0 : Initial State Vector.
% f : Function in Used RK4 solver
% t0 : Initial time in sec
% tf : Final time in sec
% h : Step size used in RK4 solver.
format long
% Constants %
i=1;
J2 = 0; %1.0826359e-3;
Re = 6378;
% Initial Values %
%Euler Angles in rad %
Roll(1) = 0;
Pitch(1)= 0;
Yaw(1) = 0;
% Angular Velocities in rad/sec%
wx(1) = 0;
wy(1) = 0;
wz(1) = 0;
W(:,i) = [wx(i),wy(i),wz(i)]; %
% Angular Velocities Vector %
% Positions in km %
rx(1) = -5634.985593;
ry(1) = 3944.062506;
rz(1) = 8.769715;
R(:,i) = [rx(i),ry(i),rz(i)]; %
% Position Vector %
% Velocities in km/sec %
vx(1) = 0.574218;
vy(1) = 0.803617;2
vz(1) = 7.548263;
V(:,i) = [vx(i),vy(i),vz(i)]; %
% Velocity Vector %
% RK Conditions %
f = @RKINPUTm;
y0 = [[Roll; Pitch; Yaw]; W; R; V];
% 4th Order Runge Kutta Calculations %
dt_sim = 1; % Simulation Time %
step = 0.1;
tf = dt_sim;
t0 = 0;
m = 4000/dt_sim;
ySave = zeros(12,1);
Coe = zeros(7,m);
trq = [0;0;0];
t_gost = tf:tf:tf*m;
for i = 1:m
t = (i-1)*dt_sim;
ySave(:,i) = y0;
dv = [0;0;0];
    if t>700 && t<800
        v_uydu = y0(10:12);
        v_uydu = v_uydu/norm(v_uydu);
        dv = v_uydu*0.000008;
    end
    tDist = [1e-4*sin(t/5800*pi); 1e-3*cos(t/5800*pi); -1e-4*sin(t/5800*pi)];
    trq = trq+tDist;
    [T,Y] = RK4m(f,t0,tf,y0,step,dv,trq);
    y0 = Y(:,end);
    Coe(:,i) = Kepler(y0(7:9),y0(10:12));
end
% Plot3D%
earth_sphere;
grid on
hold on
plot3(ySave(7,:),ySave(8,:),ySave(9,:),'Color','r','linewidth',3);
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');
title('Orbit Figure');
% Subplot Graphs %
% --------------- %
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