function Coe = Kepler(R,V)
% This function computes the classical orbital elements which is called COE
% from the state vector (R,V) using Algorithm
% Constants %
% Mu    : Gravitiational Parameter in km^3/s^2
% R     : Position vector in geocentric equatorial frame in km
% V     : Velocity vector in geocentric equatorial frame in km/s
% r & v : The magnitudes of R and V
% vr    : Radial velocity component in km/s
% H     : The angular momentum vector in km^2/s
% h     : The magnitude of  H in km^2/s
% incl  : Inclination of the orbit in rad
% N     : The node line vector km^2/s
% n     : The magnitude of N
% cp    : Cross product of N and R
% RA    : Right ascension of the ascending node in rad (RAAN)
% E     : Eccentricity vector
% e     : Eccentricity 
% eps   : a small number below which the eccentricity is considired be zero
% wa    : argument of perigee in rad
% TA    : True anomaly in rad
% a     : Semi major axis in km
% pi    : 3.1415926
% Coe   : Vector of orbital elements which includes [h e Ra incl wa TA a]

% Kepler Elements %
%----------------%

% Constants
Mu   = 398600.4418;
eps  = 1.e-10;
r    = norm(R);
v    = norm(V);
vr   = dot(R,V)/r;
H    = cross(R,V);
h    = norm(H);
incl = acos(H(3)/h);
N    = cross([0 0 1],H);
n    = norm(N);

% Calculatinos %
if n ~= 0
    RA = acos(N(1)/n);
    if N(2) < 0
        RA = 2*pi - RA;
    end
else
    RA = 0;
end
E = 1/Mu*((v^2 - Mu/r)*R - r*vr*V);
e = norm(E);
if n ~= 0
    if e > eps
        wa = acos(dot(N,E)/n/e);
        if E(3) < 0
            wa = 2*pi - wa;
        end
    else
        wa = 0;
    end
else
    wa = 0;
end
if e > eps
    TA = acos(dot(E,R)/e/r);
    if vr < 0
        TA = 2*pi - TA;
    end
else
    cp = cross(N,R);
    if cp(3) >= 0
        TA = 2*pi - acos(dot(N,R)/n/r);
    end
end
a = h^2/Mu/(1 - e^2);
Coe = [h e RA incl wa TA a];
Coe(3:6)=Coe(3:6)*180/pi;

% End of Kepler Calculations %
%----------------------------%

