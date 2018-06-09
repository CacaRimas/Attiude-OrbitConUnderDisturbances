function [T,S]= RK4m(f,t0,tf,y0,step,dv,trq)
% RK4 This function solves RK4 routine for any given F function.
% [T,Y]= rk40(f,t0,tf,y0,h)
%
% INPUTS:
% f     : Function Name
% t0    : Initial "t" Value
% tf    : Final "t" Value
% y0    : Inital "y" Value
% h     : Step size

T          = t0:step:tf;
S          = zeros(12,length(T)-1);
Y          = y0;
for i      = 1:(length(T)-1)
    k1     = f(T(i),Y,dv,trq);
    k2     = f(T(i)+0.5*step,Y+0.5*step*k1,dv,trq);
    k3     = f((T(i)+0.5*step),(Y+0.5*step*k2),dv,trq);
    k4     = f((T(i)+step),(Y+k3*step),dv,trq);
    Y      = Y+(1/6)*(k1+2*k2+2*k3+k4)*step;
    S(:,i) = Y;
end
