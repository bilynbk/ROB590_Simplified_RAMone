function JFootR = JcontPointR(in1,in2)
%JCONTPOINTR
%    JFOOTR = JCONTPOINTR(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    12-Feb-2017 17:50:27

alphaR = in1(4,:);
betaR = in1(5,:);
lH = in2(:,3);
lL2 = in2(:,7);
lL3 = in2(:,11);
phi = in1(3,:);
t2 = alphaR+phi;
t3 = cos(t2);
t4 = lL2.*t3;
t5 = alphaR+betaR+phi;
t6 = cos(t5);
t7 = lL3.*t6;
t8 = sin(t2);
t9 = lL2.*t8;
t10 = sin(t5);
t11 = lL3.*t10;
JFootR = reshape([1.0,0.0,0.0,1.0,t4+t7+lH.*cos(phi),t9+t11+lH.*sin(phi),t4+t7,t9+t11,t7,t11,0.0,0.0,0.0,0.0],[2,7]);
