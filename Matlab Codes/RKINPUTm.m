function Y = RKINPUTm(t,y0,dv,trq);
% This code is 4th order Runge Kutta solver function.
% Used in Runge Kutta solver.
% Information about Constant;
%
% Mu              : Earth's gravitional constant.
% J2              : Coefficent.
% Re              : Radius of Earth.
% I               : Spacecraft moment of Inertia
% phi,theta,yaw   : Euler Angles of Satellite in rad.
% W               : Angular velocity vector
% R               : Position vector
% V               : Velocity vector
% y0              : Initial State Vector.
%
% Information about Results;
% Edot            : Attitude Kinematics
% Cbg             : Rotation Matrix
% Rb              : Spacecraft orbital position in spacecraft body frame
% Rbcross         : Cross Matrix on Rb
% Wcross          : Cross Matrix on W
% Tg              : Disturbance Torque
% Wdot            : Angular accelarations
% Rdot            : Velocities
% Vdot            : Accelarations

% Constants %
%-----------%
Mu    = 398600.4418;
J2    = 1.0826359e-3; %0; 
Re    = 6378;
zg    = [0 0 1]';
I     = [300 0 0; 0 320 0; 0 0 250];
%State Defination%
phi   = y0(1);
theta = y0(2);
yaw   = y0(3);
W     = y0(4:6);
R     = y0(7:9);
V     = y0(10:12);
% Functions and Calculations %
%----------------------------%
% Euler Functions and Calculations %
Edot1   = [1, sin(phi)*tan(theta), cos(phi)*tan(theta)                                                                       % Roll %
           0, cos(phi)           , -sin(phi)                                                                                 % Pitch %
           0, sin(phi)*sec(theta), cos(phi)*sec(theta)];                                                                     % Yaw %
Edot    = Edot1*W;                                                                                                           % Calculation Matrix %    
%Angular Velocities%
Cbg     = [                           cos(theta)*cos(yaw), cos(theta)*sin(yaw)                           , -sin(theta)       % Rotation Matrix %         
           sin(phi)*sin(theta)*cos(yaw)-cos(phi)*sin(yaw), sin(phi)*sin(theta)*sin(yaw)+cos(phi)*cos(yaw), sin(phi)*cos(theta)
           cos(phi)*sin(theta)*cos(yaw)+sin(phi)*sin(yaw), cos(phi)*sin(theta)*sin(yaw)-sin(phi)*cos(yaw), cos(phi)*cos(theta)];       
Rb      = Cbg*R;                                                                                                             % Body Frame Position %
r       = norm(Rb);                                                                                                          % Magnititude %
Wcross  = [0 -W(3) W(2)                                                                                                      % Cross Vectors %
           W(3) 0 -W(1)
           -W(2) W(1) 0];       
Rbcross = [0 -Rb(3) Rb(2)                                                                                                    % Cross Vectors %
           Rb(3) 0 -Rb(1)
           -Rb(2) Rb(1) 0];   
Tg      = ((3*Mu)/(r^5))*(Rbcross*I*Rb);                                                                                     % Disturbance Torque %
% Tg      = 0;
Wdot    = (inv(I))*(-Wcross*I*W+Tg+trq);
% Position Vector %
Rdot    = V;
% Velocities %
Vdot    = (((-Mu)/(r^3))*R)+((3*Mu*J2*Re^2)/2/r^5)*(((5*((R'*zg)^2)/r^2)-1)*R-2*(R'*zg)*zg);
% Conculision %
Y       = [Edot; Wdot; Rdot; Vdot+dv];