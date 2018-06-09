function [Kp,Kd] = Control(Ts,Mp)
% This function solves which Kp and Kd constant on PD Controller with given
% Maximum overshoot and settling time condition
% Mp    : Maximum overshoot
% Ts    : Settling Time
% Kp    : Proportianl Constant
% Kd    : Derivative Constant
% I     : Inertial Matrix
% Wn    : Undamped Natural Frequency
% Damp  : Damping Ratio
% C     : Constant on Calculation in Maximum Overshoot
% -------------------

% Constnant
I    = [300 0 0; 0 300 0; 0 0 300];

% Calculations
C    = -log(Mp)/pi;
Damp = abs(sqrt(C/(1+C)));
A    = -log(0.02*sqrt(1-Damp^2));
Wn   = (A/Ts/Damp);
Kp = Wn^2*I;
Kd = 2*Damp*Wn*I;
