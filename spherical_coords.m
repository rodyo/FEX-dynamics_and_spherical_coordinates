% Please report bugs and inquiries to:
%
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com    (personal)
%              oldenhuis@luxspace.lu  (professional)
% Affiliation: LuxSpace s�rl
% Licence    : BSD


% If you find this work useful, please consider a donation:
% https://www.paypal.me/RodyO/3.5




% System I ("Terran"): 
% ------------------------------------------------------------
% XYZ right-handed system
% +Z North pole
% -Z South pole
% azimuth from +X

% System I-1a: elevation, theta := azimuth
% System I-1b: elevation, phi := azimuth

% System I-2a: Northern colatitude, theta := azimuth
% System I-2b: Northern colatitude, phi := azimuth

% System I-3a: Southern colatitude, theta := azimuth
% System I-3b: Southern colatitude, phi := azimuth



% System II ("Uranian"): 
% ------------------------------------------------------------
% XYZ right-handed system
% +X North pole
% -X South pole
% azimuth from +Y

% System II-1a: Uranian +X elevation: phi := azimuth
% System II-1b: Uranian +X elevation: theta := azimuth

% System II-2a: Uranian +X colatitude: positive X-colatitude ("northern"), phi := azimuth
% System II-2b: Uranian +X colatitude: positive X-colatitude ("northern"), theta := azimuth

% System II-3a: Uranian -X colatitude: negative X-colatitude ("southern"), phi := azimuth
% System II-3b: Uranian -X colatitude: negative X-colatitude ("southern"), theta := azimuth




% System III ("Uranian"): 
% ------------------------------------------------------------
% XYZ right-handed system
% +Y North pole
% -Y South pole
% azimuth from +Z
%
% System III-1a: Uranian +Y elevation: phi := azimuth
% System III-1b: Uranian +Y elevation: theta := azimuth

% System III-2a: Uranian +Y: positive Y-colatitude ("northern"), phi := azimuth
% System III-2b: Uranian +Y: positive Y-colatitude ("northern"), theta := azimuth

% System III-3a: Uranian -Y: negative Y-colatitude ("southern"), phi := azimuth
% System III-3b: Uranian -Y: negative Y-colatitude ("southern"), theta := azimuth






%% Simple test case

clc

% the Target vector
rT = [2 3 4];

% its velocity
vT = [-3 8 6];

% and acceleration
aT = [0.4 -0.3 -1];




% Rename these to easier-to-recognize quantities
x = rT(1);   x_dot = vT(1);   x_dot_dot = aT(1);
y = rT(2);   y_dot = vT(2);   y_dot_dot = aT(2);
z = rT(3);   z_dot = vT(3);   z_dot_dot = aT(3); 





%% System I-1a: elevation, theta := azimuth


% theta: azimuth (angle between positive X-axis and projection of target vector onto the XY plane, positive counterclockwise)
% phi  : elevation (angle between XY plane and target vector, positive counterclockwise)

% Compute angles
R     = norm(rT);
phi   = asin(z/R); % ASSERT: must equal atan2(z, hypot(x,y))
theta = atan2(y,x);

% Radial and angular velocities
R_dot     = vT*rT.'/R;
phi_dot   = (R*z_dot - z*R_dot) / (R^2*sqrt(1-(z/R)^2));
theta_dot = (x*y_dot - y*x_dot) / (x^2 + y^2);

% Radial and angular accelerations
R_dot_dot = ( R*(vT*vT.' + rT*aT.') - R_dot*(vT*rT.') ) / R^2;
phi_dot_dot = ( ...
    R_dot_dot*R*z^3 + R^2*z*(2*R_dot^2 - z*z_dot_dot + z_dot^2) - ...
    z^3*R_dot^2 - R^3*(R_dot_dot*z + 2*R_dot*z_dot) + R^4*z_dot_dot ) / ...
    (R^3*(R^2-z^2)*sqrt(1-(z/R)^2));
theta_dot_dot = (...
    x*y*(2*x_dot^2 + y*y_dot_dot - 2*y_dot^2) - x^2*(y*x_dot_dot + 2*x_dot*y_dot) + ...
    y^2*(2*x_dot*y_dot - y*x_dot_dot) + x^3*y_dot_dot) / (x^2 + y^2)^2;


% pre-compute sines, cosines
cT = cos(theta);    cP = cos(phi);
sT = sin(theta);    sP = sin(phi);

% Unit vectors
eR     = [+cP*cT; +cP*sT; +sP];
ePhi   = [-sP*cT; -sP*sT; +cP];
eTheta = [-sT; +cT; +0];
% ASSERT: 
%
% eR = rT/R
% |eTheta| = |ePhi| = |eR| = 1
% cross(eR,     eTheta) = ePhi
% cross(eTheta,   ePhi) = eR       (orthonormal basis)
% cross(ePhi,       eR) = eTheta


% Time derivatives of unit vectors
%  eR_dot     =       0           +  phi_dot *    ePhi +  theta_dot*cP*eTheta;
%  ePhi_dot   = -phi_dot *eR      +         0          -  theta_dot*sP*eTheta;
%  eTheta_dot = -theta_dot*cP*eR  +  theta_dot*sP*ePhi +         0           ;
% %
% % ASSERT:
% % 
% deRdph = [-sP*cT; -sP*sT; cP];
% deRdth = [-cP*sT; +cP*cT; 0];
% deRdth*theta_dot + deRdph*phi_dot; % Should equal eR_dot
% 
% dePhidph = [-cP*cT; -cP*sT; -sP];
% dePhidth = [+sP*sT; -sP*cT; +0]; 
% dePhidth*theta_dot + dePhidph*phi_dot; % Should equal dePhi_dot
% 
% deThetadth = [-cT; -sT; 0];
% deThetadth*theta_dot % Should equal eTheta_dot



% Radius vector
r = R*eR;
% ASSERT: should equal rT

% Velocity
v = R_dot*eR + R*phi_dot*ePhi + R*theta_dot*cP*eTheta;
% ASSERT: v must equal:
% v = R_dot*eR + R*eR_dot
% which equals vT.

% Acceleration
a = ( R_dot_dot - R*(phi_dot^2 + theta_dot^2*cP^2) ) * eR + ...
    ( 2*R_dot*phi_dot + R*(phi_dot_dot + theta_dot^2*cP*sP) ) * ePhi + ...
    ( 2*R_dot*theta_dot*cP + R*(theta_dot_dot*cP - 2*theta_dot*phi_dot*sP) ) * eTheta;
% ASSERT: a must equal:
% a = R_dot_dot*eR + 2*R_dot*eR_dot + R*eR_dot_dot
% which, in terms of the unit vector derivatives, is
% a = R_dot_dot*eR + 2*R_dot*eR_dot + R*(...
%     phi_dot_dot*ePhi + phi_dot*ePhi_dot + ...
%     theta_dot_dot*cP*eTheta - theta_dot*phi_dot*sP*eTheta + theta_dot*cP*eTheta_dot)
% which equals aT.
 


%% System I-1b: elevation, phi := azimuth

% phi: azimuth (angle between positive X-axis and projection of target vector onto the XY plane, positive counterclockwise)
% theta: elevation (angle between XY plane and target vector, positive counterclockwise)

% Compute angles
R     = norm(rT);
theta = asin(z/R); % ASSERT: must equal atan2(z, hypot(x,y))
phi   = atan2(y,x);

% Radial and angular velocities
R_dot     = vT*rT.'/R;
theta_dot = (R*z_dot - z*R_dot) / (R^2*sqrt(1-(z/R)^2));
phi_dot   = (x*y_dot - y*x_dot) / (x^2 + y^2);

% Radial and angular accelerations
R_dot_dot = ( R*(vT*vT.' + rT*aT.') - R_dot*(vT*rT.') ) / R^2;
theta_dot_dot = ( ...
    R_dot_dot*R*z^3 + R^2*z*(2*R_dot^2 - z*z_dot_dot + z_dot^2) - ...
    z^3*R_dot^2 - R^3*(R_dot_dot*z + 2*R_dot*z_dot) + R^4*z_dot_dot ) / ...
    (R^3*(R^2-z^2)*sqrt(1-(z/R)^2));
phi_dot_dot  = (...
    x*y*(2*x_dot^2 + y*y_dot_dot - 2*y_dot^2) - x^2*(y*x_dot_dot + 2*x_dot*y_dot) + ...
    y^2*(2*x_dot*y_dot - y*x_dot_dot) + x^3*y_dot_dot) / (x^2 + y^2)^2;


% pre-compute sines, cosines
cT = cos(theta);    cP = cos(phi);
sT = sin(theta);    sP = sin(phi);


% Unit vectors
eR     = [+cT*cP; +cT*sP; +sT];
eTheta = [-sT*cP; -sT*sP; +cT];
ePhi   = [-sP; +cP; +0];
% ASSERT: 
%
% eR = rT/R
% |eTheta| = |ePhi| = |eR| = 1
% cross(eR,     ePhi)   = eTheta
% cross(ePhi,   eTheta) = eR       (orthonormal basis)
% cross(eTheta, eR)     = ePhi


% Time derivatives of unit vectors
%  eR_dot     =       0         +  theta_dot *eTheta  +  phi_dot*cT*ePhi;
%  eTheta_dot = -theta_dot *eR  +         0           -  phi_dot*sT*ePhi;
%  ePhi_dot   = -phi_dot*cT*eR  +  phi_dot*sT*eTheta  +         0       ;
%
% ASSERT:
% 
% deRdth = [-sT*cP; -sT*sP; cT];
% deRdph = [-cT*sP; +cT*cP; 0];
% deRdth*theta_dot + deRdph*phi_dot % Should equal eR_dot
% 
% deThetadth = [-cT*cP; -cT*sP; -sT];
% deThetadph = [+sT*sP; -sT*cP; +0]; 
% deThetadth*theta_dot + deThetadph*phi_dot % Should equal eTheta_dot
% 
% dePhidph = [-cP; -sP; 0];
% dePhidph*phi_dot % Should equal ePhi_dot




% Radius vector
r = R*eR;
% ASSERT: should equal rT

% Velocity
v = R_dot*eR + R*theta_dot*eTheta + R*phi_dot*cT*ePhi;
% ASSERT: v must equal:
% v = R_dot*eR + R*eR_dot
% which equals vT.

% Acceleration
a = ( R_dot_dot - R*(theta_dot^2 + phi_dot^2*cT^2) ) * eR + ...
    ( 2*R_dot*theta_dot + R*(theta_dot_dot + phi_dot^2*cT*sT) ) * eTheta + ...
    ( 2*R_dot*phi_dot*cT + R*(phi_dot_dot*cT - 2*phi_dot*theta_dot*sT) ) * ePhi;
% ASSERT: a must equal:
% a = R_dot_dot*eR + 2*R_dot*eR_dot + R*eR_dot_dot
% which, in terms of the unit vector derivatives, is
% a = R_dot_dot*eR + 2*R_dot*eR_dot + R*(...
%     theta_dot_dot*eTheta  + theta_dot*eTheta_dot + ...
%      phi_dot_dot*cT*ePhi - phi_dot*theta_dot*sT*ePhi + phi_dot*cT*ePhi_dot)
% which equals aT.
 


%% System I-2a: Northern Colatitude, theta := azimuth

% theta: azimuth (angle between positive X-axis and projection of target vector onto the XY plane, positive counterclockwise)
% phi: northern colatitude (angle between positive Z-axis and target vector, positive clockwise)

% Compute angles
R     = norm(rT);
%phi   = 
theta = atan2(y, x);

% Radial and angular velocities
R_dot     = vT*rT.'/R;
% phi_dot   = 
theta_dot = (x*y_dot - y*x_dot) / (x^2 + y^2);

% Radial and angular accelerations
R_dot_dot = ( R*(vT*vT.' + rT*aT.') - R_dot*(vT*rT.') ) / R^2;
% phi_dot_dot   = 
% theta_dot_dot = 


% pre-compute sines, cosines
cT = cos(theta);    cP = cos(phi);
sT = sin(theta);    sP = sin(phi);



% eR     = r/norm(r);
ePhi   = [];
eTheta = [];



%% System I-2b: Northern Colatitude, phi := azimuth

% phi: azimuth (angle between positive X-axis and projection of target vector onto the XY plane, positive counterclockwise)
% theta: northern colatitude (angle between positive Z-axis and target vector, positive clockwise)

% Compute angles
R     = norm(rT);
% theta =
phi   = atan2(y, x);

% Radial and angular velocities
% TODO

% Radial and angular accelerations
% TODO

% pre-compute sines, cosines
cT = cos(theta);    cP = cos(phi);
sT = sin(theta);    sP = sin(phi);


% eR     = r/norm(r);
ePhi   = [];
eTheta = [];



%% System I-3a: Southern Colatitude, theta := azimuth

% theta: azimuth (angle between positive X-axis and projection of target vector onto the XY plane, positive counterclockwise)
% phi: southern colatitude (angle between negative Z-axis and target vector, positive counterclockwise)



%% System I-3b: Southern Colatitude, phi := azimuth

% phi  : azimuth (angle between positive X-axis and projection of target vector onto the XY plane, positive counterclockwise)
% theta: southern colatitude (angle between negative Z-axis and target vector, positive counterclockwise)

% Compute angles
R     = norm(rT);
theta = acos(-z/R); % ASSERT: should equal atan2(hypot(rT(1),rT(2)), -rT(3));
phi   = atan2(y, x);

% Radial and angular velocities
R_dot     = vT*rT.'/R;
theta_dot = (R*z_dot - R_dot*z) / (R^2*sqrt(1-(z/R)^2));
phi_dot   = (x*y_dot - y*x_dot) / (x^2 + y^2);

% Radial and angular accelerations
R_dot_dot = ( R*(vT*vT.' + rT*aT.') - R_dot*(vT*rT.') ) / R^2;
theta_dot_dot = ( ...
    R_dot_dot*R*z^3 + R^2*z*(2*R_dot^2 - z*z_dot_dot + z_dot^2) - ...
    z^3*R_dot^2 - R^3*(R_dot_dot*z + 2*R_dot*z_dot) + R^4*z_dot_dot ) / ...
    (R^3*(R^2-z^2)*sqrt(1-(z/R)^2));
phi_dot_dot  = (...
    x*y*(2*x_dot^2 + y*y_dot_dot - 2*y_dot^2) - x^2*(y*x_dot_dot + 2*x_dot*y_dot) + ...
    y^2*(2*x_dot*y_dot - y*x_dot_dot) + x^3*y_dot_dot) / (x^2 + y^2)^2;

% pre-compute sines, cosines
cT = cos(theta);    cP = cos(phi);
sT = sin(theta);    sP = sin(phi);

% The unit vectors
eR     = [sT*cP; sT*sP; -cT];
ePhi   = [-sP;  cP;  0];
eTheta = [cT*cP;  cT*sP;  sT];
% ASSERT:
%
% eR = r/R
% |eTheta| = |ePhi| = |eR| = 1
% cross(eR,     ePhi)   = eTheta
% cross(ePhi,   eTheta) = eR       (orthonormal basis)
% cross(eTheta, eR)     = ePhi


% Time derivatives of unit vectors
% eR_dot     =       0         +  theta_dot *eTheta  +  phi_dot*sT*ePhi
% eTheta_dot = -theta_dot *eR  +         0           +  phi_dot*cT*ePhi
% ePhi_dot   = -phi_dot*sT*eR  -  phi_dot*cT*eTheta  +         0       
%
% ASSERT:
%
% deRdth = [cT*cP; cT*sP; sT];
% deRdph = [-sT*sP; sT*cP; 0];
% deRdth*theta_dot + deRdph*phi_dot % Should equal eR_dot
% 
% deThetadth = [-sT*cP;  -sT*sP;  cT];
% deThetadph = [-cT*sP;  cT*cP;  0];
% deThetadth*theta_dot + deThetadph*phi_dot % Should equal eTheta_dot
% 
% dePhidph = [-cP; -sP; 0];
% dePhidph*phi_dot % Should equal ePhi_dot


% Radius vector
r = R*eR;
% ASSERT: should equal rT

% Velocity
v = R_dot*eR + R*theta_dot*eTheta + R*phi_dot*sT*ePhi;
% ASSERT: v must equal:
%  v = R_dot*eR + R*eR_dot

% Acceleration
a = ( R_dot_dot - R*(theta_dot^2 + phi_dot^2*sT^2) ) * eR + ...
    ( 2*R_dot*theta_dot + R*(theta_dot_dot - phi_dot^2*sT*cT) ) * eTheta + ...
    ( 2*R_dot*phi_dot*sT + R*(2*phi_dot*theta_dot*cT + phi_dot_dot*sT) ) * ePhi;
% ASSERT: a must equal:
% a = R_dot_dot*eR + 2*R_dot*eR_dot + R*eR_dot_dot
% which, in terms of the unit vector derivatives, is
% a = R_dot_dot*eR + 2*R_dot*eR_dot + R*(...
%     theta_dot_dot*eTheta  + theta_dot*eTheta_dot + ...
%     phi_dot_dot*sT*ePhi + phi_dot*theta_dot*cT*ePhi + phi_dot*sT*ePhi_dot)







%% System II-1a: Uranian +X elevation: phi := azimuth

% theta: elevation from YZ plane (angle between projection of target vector onto the YZ plane and the target vector, positive towards positive X-axis)
% phi  : azimuth (angle between positive Y-axis and projection of target vector onto the YZ plane, positive counterclockwise)


%% System II-1b: Uranian +X elevation: theta := azimuth

% phi  : elevation from YZ plane (angle between projection of target vector onto the YZ plane and the target vector, positive towards positive X-axis)
% theta: azimuth (angle between positive Y-axis and projection of target vector onto the YZ plane, positive counterclockwise)



%% System II-2a: Uranian +X colatitude: positive X-colatitude, phi := azimuth

% theta: +X colatitude (angle between positive X-axis and target vector, positive clockwise)
% phi  : azimuth (angle between positive Y-axis and projection of target vector onto the YZ plane, positive counterclockwise)



%% System II-2b: Uranian +X colatitude: positive X-colatitude, theta := azimuth

% phi  : +X colatitude (angle between positive X-axis and target vector, positive clockwise)
% theta: azimuth (angle between positive Y-axis and projection of target vector onto the YZ plane, positive counterclockwise)



%% System II-3a: Uranian -X colatitude: negative X-colatitude, phi := azimuth

% theta: -X colatitude (angle between negative X-axis and target vector, positive counterclockwise)
% phi  : azimuth (angle between positive Y-axis and projection of target vector onto the YZ plane, positive counterclockwise)

% Compute radius, angles
R     = norm(rT);
theta = acos(-x/R); % ASSERT: should equal atan2(hypot(y,z), -x);
phi   = atan2(z, y);

% Radial and angular velocities
R_dot     = vT*rT.'/R;
theta_dot = (R*x_dot - R_dot*x) / (R^2*sqrt(1-(x/R)^2));
phi_dot   = (y*z_dot - z*y_dot) / (y^2 + z^2);

% Radial and angular accelerations
R_dot_dot = ( R*(vT*vT.' + rT*aT.') - R_dot*(vT*rT.') ) / R^2;
theta_dot_dot = ( ...
    R_dot_dot*R*x^3 + R^2*x*(2*R_dot^2 - x*x_dot_dot + x_dot^2) - ...
    x^3*R_dot^2 - R^3*(R_dot_dot*x + 2*R_dot*x_dot) + R^4*x_dot_dot) /...
    (R^3*(R^2-x^2)*sqrt(1-(x/R)^2));
phi_dot_dot = (...
    y*z*(2*y_dot^2 + z*z_dot_dot - 2*z_dot^2) - y^2*(z*y_dot_dot + 2*y_dot*z_dot) + ...
    z^2*(2*y_dot*z_dot - z*y_dot_dot) + y^3*z_dot_dot) / (y^2 + z^2)^2;


% pre-compute sines, cosines
cT = cos(theta);    cP = cos(phi);
sT = sin(theta);    sP = sin(phi);

% The unit vectors
eR     = [-cT; sT*cP; sT*sP];
ePhi   = [0; -sP; cP];
eTheta = [sT; cT*cP;  cT*sP];

% ASSERT:
%
% eR = r/R
% |eTheta| = |ePhi| = |eR| = 1
% cross(eR,     ePhi)   = eTheta
% cross(ePhi,   eTheta) = eR       (orthonormal basis)
% cross(eTheta, eR)     = ePhi
%
% FIXME: this does not prove that it's a valid basis for these coordinates

% Time derivatives of unit vectors
% eR_dot     =       0         +  theta_dot *eTheta  +  phi_dot*sT*ePhi;
% eTheta_dot = -theta_dot *eR  +         0           +  phi_dot*cT*ePhi;
% ePhi_dot   = -phi_dot*sT*eR  -  phi_dot*cT*eTheta  +         0       ;
%
% ASSERT:
%
% deRdth = [sT; cT*cP; cT*sP];
% deRdph = [0; -sT*sP; sT*cP];
% deRdth*theta_dot + deRdph*phi_dot % Should equal eR_dot
% 
% deThetadth = [cT; -sT*cP;  -sT*sP];
% deThetadph = [0; -cT*sP;  cT*cP];
% deThetadth*theta_dot + deThetadph*phi_dot % Should equal eTheta_dot
% 
% dePhidph = [0; -cP; -sP];
% dePhidph*phi_dot % Should equal ePhi_dot



% Radius vector
r = R*eR;
% ASSERT: should equal rT

% Velocity
v = R_dot*eR + R*theta_dot*eTheta + R*phi_dot*sT*ePhi;
% ASSERT: v must equal:
% v = R_dot*eR + R*eR_dot

% Acceleration
a = ( R_dot_dot - R*(theta_dot^2 + phi_dot^2*sT^2) ) * eR + ...
    ( 2*R_dot*theta_dot + R*(theta_dot_dot - phi_dot^2*sT*cT) ) * eTheta + ...
    ( 2*R_dot*phi_dot*sT + R*(2*phi_dot*theta_dot*cT + phi_dot_dot*sT) ) * ePhi;
% ASSERT: a must equal:
% a = R_dot_dot*eR + 2*R_dot*eR_dot + R*eR_dot_dot
% which, in terms of the unit vector derivatives, is
% a = R_dot_dot*eR + 2*R_dot*eR_dot + R*(...
%     theta_dot_dot*eTheta  + theta_dot*eTheta_dot + ...
%      phi_dot_dot*sT*ePhi + phi_dot*theta_dot*cT*ePhi + phi_dot*sT*ePhi_dot)
 


%% System II-3b: Uranian -X colatitude: negative X-colatitude, theta := azimuth

% phi  : -X colatitude (angle between negative X-axis and target vector, positive counterclockwise)
% theta: azimuth (angle between positive Y-axis and projection of target vector onto the YZ plane, positive counterclockwise)











%% System III-1a: Uranian +Y elevation: phi := azimuth

% theta: elevation from XZ plane (angle between projection of target vector onto the XZ plane and the target vector, positive towards positive Y-axis)
% phi  : azimuth (angle between positive Z-axis and projection of target vector onto the XZ plane, positive counterclockwise)



%% System III-1b: Uranian +Y elevation: theta := azimuth

% phi  : elevation from XZ plane (angle between projection of target vector onto the XZ plane and the target vector, positive towards positive Y-axis)
% theta: azimuth (angle between positive Z-axis and projection of target vector onto the XZ plane, positive counterclockwise)



%% System III-2a: Uranian +Y: positive Y-colatitude, phi := azimuth

% theta: +Y colatitude (angle between positive Y-axis and target vector, positive clockwise)
% phi  : azimuth (angle between positive Z-axis and projection of target vector onto the XZ plane, positive counterclockwise)



%% System III-2b: Uranian +Y: positive Y-colatitude, theta := azimuth

% phi  : +Y colatitude (angle between positive Y-axis and target vector, positive clockwise)
% theta: azimuth (angle between positive Z-axis and projection of target vector onto the XZ plane, positive counterclockwise)



%% System III-3a: Uranian -Y: negative Y-colatitude, phi := azimuth

% theta: -Y colatitude (angle between negative Y-axis and target vector, positive counterclockwise)
% phi  : azimuth (angle between positive Z-axis and projection of target vector onto the XZ plane, positive counterclockwise)



%% System III-3b: Uranian -Y: negative Y-colatitude, theta := azimuth

% phi  : -Y colatitude (angle between negative Y-axis and target vector, positive counterclockwise)
% theta: azimuth (angle between positive Z-axis and projection of target vector onto the XZ plane, positive counterclockwise)







