function dJFootR = dJcontPointR(in1,in2)
%DJCONTPOINTR
%    DJFOOTR = DJCONTPOINTR(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    12-Feb-2017 19:57:21

alphaR = in1(4,:);
betaR = in1(5,:);
lH = in2(:,3);
lL2 = in2(:,7);
lL3 = in2(:,11);
phi = in1(3,:);
valphaR = in1(11,:);
vbetaR = in1(12,:);
vphi = in1(10,:);
t2 = alphaR+phi;
t3 = sin(t2);
t4 = lL2.*t3;
t5 = alphaR+betaR+phi;
t6 = sin(t5);
t7 = lL3.*t6;
t8 = t4+t7;
t9 = cos(t2);
t10 = lL2.*t9;
t11 = cos(t5);
t12 = lL3.*t11;
t13 = t10+t12;
t14 = t13.*valphaR;
t15 = lL3.*t11.*vbetaR;
t16 = valphaR+vbetaR+vphi;
dJFootR = reshape([0.0,0.0,0.0,0.0,-t8.*valphaR-vphi.*(t4+t7+lH.*sin(phi))-lL3.*t6.*vbetaR,t14+t15+vphi.*(t10+t12+lH.*cos(phi)),-t8.*valphaR-t8.*vphi-lL3.*t6.*vbetaR,t14+t15+t13.*vphi,-lL3.*t6.*t16,lL3.*t11.*t16,0.0,0.0,0.0,0.0],[2,7]);
