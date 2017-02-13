function fCG = FCorGrav(in1,in2)
%FCORGRAV
%    FCG = FCORGRAV(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    12-Feb-2017 19:57:20

alphaL = in1(6,:);
alphaR = in1(4,:);
betaL = in1(7,:);
betaR = in1(5,:);
g = in2(:,12);
l2 = in2(:,6);
l3 = in2(:,10);
lH = in2(:,3);
lL2 = in2(:,7);
m1 = in2(:,1);
m2 = in2(:,4);
m3 = in2(:,8);
phi = in1(3,:);
valphaL = in1(13,:);
valphaR = in1(11,:);
vbetaL = in1(14,:);
vbetaR = in1(12,:);
vphi = in1(10,:);
t2 = alphaL+betaL+phi;
t3 = sin(t2);
t4 = alphaR+betaR+phi;
t5 = sin(t4);
t6 = vphi.^2;
t7 = valphaL.^2;
t8 = valphaR.^2;
t9 = alphaL+phi;
t10 = sin(t9);
t11 = alphaR+phi;
t12 = sin(t11);
t13 = sin(phi);
t14 = vbetaL.^2;
t15 = cos(t2);
t16 = vbetaR.^2;
t17 = cos(t4);
t18 = cos(t9);
t19 = cos(t11);
t20 = cos(phi);
t21 = alphaL+betaL;
t22 = sin(t21);
t23 = alphaR+betaR;
t24 = sin(t23);
t25 = sin(alphaL);
t26 = sin(alphaR);
t27 = sin(betaL);
t28 = sin(betaR);
t29 = l3.*lL2.*m3.*t16.*t28;
t30 = l3.*lL2.*m3.*t28.*valphaR.*vbetaR.*2.0;
t31 = l3.*lL2.*m3.*t28.*vbetaR.*vphi.*2.0;
t32 = l3.*lL2.*m3.*t14.*t27;
t33 = l3.*lL2.*m3.*t27.*valphaL.*vbetaL.*2.0;
t34 = l3.*lL2.*m3.*t27.*vbetaL.*vphi.*2.0;
fCG = [l3.*m3.*t3.*t6+l3.*m3.*t3.*t7+l3.*m3.*t5.*t6+l3.*m3.*t5.*t8+l2.*m2.*t6.*t10+l2.*m2.*t7.*t10+l2.*m2.*t6.*t12+l3.*m3.*t3.*t14+l2.*m2.*t8.*t12+l3.*m3.*t5.*t16+lH.*m2.*t6.*t13.*2.0+lH.*m3.*t6.*t13.*2.0+lL2.*m3.*t6.*t10+lL2.*m3.*t7.*t10+lL2.*m3.*t6.*t12+lL2.*m3.*t8.*t12+l3.*m3.*t3.*valphaL.*vbetaL.*2.0+l3.*m3.*t5.*valphaR.*vbetaR.*2.0+l3.*m3.*t3.*valphaL.*vphi.*2.0+l2.*m2.*t10.*valphaL.*vphi.*2.0+l3.*m3.*t5.*valphaR.*vphi.*2.0+l2.*m2.*t12.*valphaR.*vphi.*2.0+l3.*m3.*t3.*vbetaL.*vphi.*2.0+l3.*m3.*t5.*vbetaR.*vphi.*2.0+lL2.*m3.*t10.*valphaL.*vphi.*2.0+lL2.*m3.*t12.*valphaR.*vphi.*2.0;-g.*m1-g.*m2.*2.0-g.*m3.*2.0-l3.*m3.*t6.*t15-l2.*m2.*t6.*t18-l3.*m3.*t7.*t15-l2.*m2.*t6.*t19-l2.*m2.*t7.*t18-l3.*m3.*t6.*t17-l2.*m2.*t8.*t19-l3.*m3.*t8.*t17-l3.*m3.*t14.*t15-l3.*m3.*t16.*t17-lH.*m2.*t6.*t20.*2.0-lH.*m3.*t6.*t20.*2.0-lL2.*m3.*t6.*t18-lL2.*m3.*t6.*t19-lL2.*m3.*t7.*t18-lL2.*m3.*t8.*t19-l3.*m3.*t15.*valphaL.*vbetaL.*2.0-l3.*m3.*t17.*valphaR.*vbetaR.*2.0-l3.*m3.*t15.*valphaL.*vphi.*2.0-l2.*m2.*t18.*valphaL.*vphi.*2.0-l2.*m2.*t19.*valphaR.*vphi.*2.0-l3.*m3.*t17.*valphaR.*vphi.*2.0-l3.*m3.*t15.*vbetaL.*vphi.*2.0-l3.*m3.*t17.*vbetaR.*vphi.*2.0-lL2.*m3.*t18.*valphaL.*vphi.*2.0-lL2.*m3.*t19.*valphaR.*vphi.*2.0;t29+t30+t31+t32+t33+t34-g.*l3.*m3.*t3-g.*l3.*m3.*t5-g.*l2.*m2.*t10-g.*l2.*m2.*t12-g.*lH.*m2.*t13.*2.0-g.*lH.*m3.*t13.*2.0-g.*lL2.*m3.*t10-g.*lL2.*m3.*t12+l3.*lH.*m3.*t7.*t22+l2.*lH.*m2.*t7.*t25+l2.*lH.*m2.*t8.*t26+l3.*lH.*m3.*t8.*t24+l3.*lH.*m3.*t14.*t22+l3.*lH.*m3.*t16.*t24+lH.*lL2.*m3.*t7.*t25+lH.*lL2.*m3.*t8.*t26+l3.*lH.*m3.*t22.*valphaL.*vbetaL.*2.0+l3.*lH.*m3.*t24.*valphaR.*vbetaR.*2.0+l3.*lH.*m3.*t22.*valphaL.*vphi.*2.0+l2.*lH.*m2.*t25.*valphaL.*vphi.*2.0+l2.*lH.*m2.*t26.*valphaR.*vphi.*2.0+l3.*lH.*m3.*t24.*valphaR.*vphi.*2.0+l3.*lH.*m3.*t22.*vbetaL.*vphi.*2.0+l3.*lH.*m3.*t24.*vbetaR.*vphi.*2.0+lH.*lL2.*m3.*t25.*valphaL.*vphi.*2.0+lH.*lL2.*m3.*t26.*valphaR.*vphi.*2.0;t29+t30+t31-g.*l3.*m3.*t5-g.*l2.*m2.*t12-g.*lL2.*m3.*t12-l2.*lH.*m2.*t6.*t26-l3.*lH.*m3.*t6.*t24-lH.*lL2.*m3.*t6.*t26;-l3.*m3.*(g.*t5+lH.*t6.*t24+lL2.*t6.*t28+lL2.*t8.*t28+lL2.*t28.*valphaR.*vphi.*2.0);t32+t33+t34-g.*l3.*m3.*t3-g.*l2.*m2.*t10-g.*lL2.*m3.*t10-l3.*lH.*m3.*t6.*t22-l2.*lH.*m2.*t6.*t25-lH.*lL2.*m3.*t6.*t25;-l3.*m3.*(g.*t3+lH.*t6.*t22+lL2.*t6.*t27+lL2.*t7.*t27+lL2.*t27.*valphaL.*vphi.*2.0)];
