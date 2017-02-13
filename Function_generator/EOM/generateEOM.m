%% This code generates the equations of motion for a two-knee hopper

%% Some Housekeeping 
clear; clc; 
%% State Variable declarations
syms x y phi alphaR betaR alphaL betaL vx vy vphi valphaR vbetaR valphaL vbetaL
q = [x; y; phi; alphaR; betaR; alphaL; betaL];
vq = [vx; vy; vphi; valphaR; vbetaR; valphaL; vbetaL];
%% Parameters
% System
syms m1 J1 lH m2 J2 l2 lL2 m3 J3 l3 lL3
% Environmental 
syms g
param = [m1 J1 lH m2 J2 l2 lL2 m3 J3 l3 lL3 g];
%% Calculations
% Mass & Joint position calculations
    % Body
CoGB = [x;y];
Hip = CoGB + [lH*sin(phi);
    - lH*cos(phi)];
    % Right Foot 
CoGTR = Hip + [l2*sin(phi+alphaR);
    -l2*cos(phi+alphaR)];
KneeR = Hip + [lL2*sin(phi+alphaR);
    -lL2*cos(phi+alphaR)];
CoGSR = KneeR + [l3*sin(phi+alphaR+betaR);
    - l3*cos(phi+alphaR+betaR)];
FootR = KneeR + [lL3*sin(phi+alphaR+betaR);
    - lL3*cos(phi+alphaR+betaR)];
    % Left Foot
CoGTL = Hip + [l2*sin(phi+alphaL);
    -l2*cos(phi+alphaL)];
KneeL = Hip + [lL2*sin(phi+alphaL);
    -lL2*cos(phi+alphaL)];
CoGSL = KneeL + [l3*sin(phi+alphaL+betaL);
    - l3*cos(phi+alphaL+betaL)];
FootL = KneeL + [lL3*sin(phi+alphaL+betaL);
    - lL3*cos(phi+alphaL+betaL)];

% Mass velocity calculations
vCoGB = jacobian(CoGB,q)*vq;
vCoGTR = jacobian(CoGTR,q)*vq;
vCoGSR = jacobian(CoGSR,q)*vq;
vCoGTL = jacobian(CoGTL,q)*vq;
vCoGSL = jacobian(CoGSL,q)*vq;

% Lagrangian Calculations
T = 0.5*( m1*sum(vCoGB.^2) + J1*vphi^2 ...
        + m2*sum(vCoGTR.^2) + J2*(vphi+valphaR)^2 ...
        + m3*sum(vCoGSR.^2) + J3*(vphi+valphaR+vbetaR)^2 ...
        + m2*sum(vCoGTL.^2) + J2*(vphi+valphaL)^2 ...
        + m3*sum(vCoGSL.^2) + J3*(vphi+valphaL+vbetaL)^2);
U = m1*CoGB(2)*g + m2*CoGTR(2)*g + m3*CoGSR(2)*g + m2*CoGTL(2)*g + m3*CoGSL(2)*g;
L = simplify(T-U);

dL_dq = jacobian(L,q).';
dL_dqdt = jacobian(L,vq).';
ddL_dqdt_dt = jacobian(dL_dqdt,q)*vq;

% Getting EOM terms
M = jacobian(dL_dqdt,vq).';
M = simplify(M);
fCG = simplify(dL_dq - ddL_dqdt_dt);

%% Calculating contact point jacobian
% Right foot
JFootR = jacobian(FootR,q);
JFootR = simplify(JFootR);
dJFootR = sym(zeros(size(JFootR)));
for i = 1:size(JFootR,2)
    dJFootR(:,i) = jacobian(JFootR(:,i),q)*vq;
end
dJFootR = simplify(dJFootR);
% Left foot
JFootL = jacobian(FootL,q);
JFootL = simplify(JFootL);
dJFootL = sym(zeros(size(JFootL)));
for i = 1:size(JFootL,2)
    dJFootL(:,i) = jacobian(JFootL(:,i),q)*vq;
end
dJFootL = simplify(dJFootL);

%% Energy
E = simplify(T+U);

%% Center of mass of the whole system
CoG_tot = (m1*CoGB+m2*CoGTR+m3*CoGSR+m2*CoGTL+m3*CoGSL)/(m1+m2+m3+m2+m3);
CoG_tot = simplify(CoG_tot);

%% Theta and dTheta (angle between vertical line and spring)
% The spring here is placed between CoG and foot.
% GF stands for CoG-Foot
% Right Foot: 
thetaR_GF = atan2((FootR(1)-CoG_tot(1)),(CoG_tot(2)-FootR(2)));
dthetaR_GF = jacobian(thetaR_GF,q)*vq;
% Left Foot: 
thetaL_GF = atan2((FootL(1)-CoG_tot(1)),(CoG_tot(2)-FootL(2)));
dthetaL_GF = jacobian(thetaL_GF,q)*vq;
%% Theta and dTheta (angle between vertical line and thigh)
% HK stands for Hip-Knee
% Right Foot: 
thetaR_HK = atan2((KneeR(1)-Hip(1)),(Hip(2)-KneeR(2)));
dthetaR_HK = jacobian(thetaR_HK,q)*vq;
% Left Foot: 
thetaL_HK = atan2((KneeL(1)-Hip(1)),(Hip(2)-KneeL(2)));
dthetaL_HK = jacobian(thetaL_HK,q)*vq;

%% L and dL (Virtual spring length and its changing speed)
% Right Foot: 
LengthR = sum((CoGB-FootR).^2)^0.5;
dLengthR = jacobian(LengthR,q)*vq;
% Left Foot: 
LengthL = sum((CoGB-FootL).^2)^0.5;
dLengthL = jacobian(LengthL,q)*vq;

%% Create functions
if ~exist('Functions','dir')
    mkdir('Functions');
end

matlabFunction(M,'file','Functions\MassMatrix','vars',{q,param});
matlabFunction(fCG,'file','Functions\FCorGrav','vars',{[q;vq],param});
matlabFunction(JFootR,'file','Functions\JcontPointR','vars',{q,param});
matlabFunction(dJFootR,'file','Functions\dJcontPointR','vars',{[q;vq],param});
matlabFunction(JFootL,'file','Functions\JcontPointL','vars',{q,param});
matlabFunction(dJFootL,'file','Functions\dJcontPointL','vars',{[q;vq],param});
matlabFunction(CoGB,'file','Functions\posBody','vars',{q,param});
matlabFunction(Hip,'file','Functions\posHip','vars',{q,param});
matlabFunction(CoGTR,'file','Functions\CoGThighR','vars',{q,param});
matlabFunction(KneeR,'file','Functions\posKneeR','vars',{q,param});
matlabFunction(CoGSR,'file','Functions\CoGShankR','vars',{q,param});
matlabFunction(FootR,'file','Functions\posFootR','vars',{q,param});
matlabFunction(CoGTL,'file','Functions\CoGThighL','vars',{q,param});
matlabFunction(KneeL,'file','Functions\posKneeL','vars',{q,param});
matlabFunction(CoGSL,'file','Functions\CoGShankL','vars',{q,param});
matlabFunction(FootL,'file','Functions\posFootL','vars',{q,param});
matlabFunction(E,'file','Functions\energy','vars',{[q;vq],param});
matlabFunction(CoG_tot,'file','Functions\CoG_tot','vars',{q,param});
matlabFunction(thetaR_GF,'file','Functions\ThetaR','vars',{q,param});
matlabFunction(dthetaR_GF,'file','Functions\dThetaR','vars',{[q;vq],param});
matlabFunction(LengthR,'file','Functions\SpringLengthR','vars',{q,param});
matlabFunction(dLengthR,'file','Functions\dSpringLengthR','vars',{[q;vq],param});
matlabFunction(thetaL_GF,'file','Functions\ThetaL','vars',{q,param});
matlabFunction(dthetaL_GF,'file','Functions\dThetaL','vars',{[q;vq],param});
matlabFunction(LengthL,'file','Functions\SpringLengthL','vars',{q,param});
matlabFunction(dLengthL,'file','Functions\dSpringLengthL','vars',{[q;vq],param});

matlabFunction(thetaR_HK,'file','Functions\ThetaR_HK','vars',{q,param});
matlabFunction(dthetaR_HK,'file','Functions\dThetaR_HK','vars',{[q;vq],param});
matlabFunction(thetaL_HK,'file','Functions\ThetaL_HK','vars',{q,param});
matlabFunction(dthetaL_HK,'file','Functions\dThetaL_HK','vars',{[q;vq],param});

