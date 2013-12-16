% orbitalElements(position, velocity)
% Takes 3D position, 3D velocity, and mass of parent body and returns a 
% vector of orbital elements as follows:
% [orbital energy, eccentricity, period]

% TO DO: 
%   -extend orbital elements to include:
%       inclination, longitude of ascending node, argument of periapsis

function elts = orbitalElements(pos,vel, priMass)
    G=6.67384e-11;
    
    r=norm(pos);
    v=norm(vel);

    h_vec=cross(pos,vel);
    h=norm(h_vec);
    mu=G*priMass;

    a=1/(2/r-(v*v)/(G*priMass));
    en_orb=(v*v/2)-((G*priMass)/r);
    ecc=sqrt(1+(2*en_orb*h*h)/(mu*mu));
    period=(2*pi*(a^(3/2)))/sqrt(mu);
    
    elts=[en_orb,ecc,period];