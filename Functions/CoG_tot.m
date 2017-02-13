function CoG_tot = CoG_tot(in1,in2)
%COG_TOT
%    COG_TOT = COG_TOT(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    12-Feb-2017 23:35:15

alphaL = in1(6,:);
alphaR = in1(4,:);
betaL = in1(7,:);
betaR = in1(5,:);
l2 = in2(:,6);
l3 = in2(:,10);
lH = in2(:,3);
lL2 = in2(:,7);
m1 = in2(:,1);
m2 = in2(:,4);
m3 = in2(:,8);
phi = in1(3,:);
x = in1(1,:);
y = in1(2,:);
t2 = sin(phi);
t3 = alphaL+phi;
t4 = sin(t3);
t5 = alphaR+phi;
t6 = sin(t5);
t7 = m2.*2.0;
t8 = m3.*2.0;
t9 = m1+t7+t8;
t10 = 1.0./t9;
t11 = cos(phi);
t12 = alphaL+betaL+phi;
t13 = alphaR+betaR+phi;
t14 = cos(t3);
t15 = cos(t5);
CoG_tot = [t10.*(m1.*x+m2.*x.*2.0+m3.*x.*2.0+l3.*m3.*sin(t12)+l3.*m3.*sin(t13)+l2.*m2.*t4+l2.*m2.*t6+lH.*m2.*t2.*2.0+lH.*m3.*t2.*2.0+lL2.*m3.*t4+lL2.*m3.*t6);-t10.*(-m1.*y-m2.*y.*2.0-m3.*y.*2.0+l3.*m3.*cos(t12)+l3.*m3.*cos(t13)+l2.*m2.*t14+l2.*m2.*t15+lH.*m2.*t11.*2.0+lH.*m3.*t11.*2.0+lL2.*m3.*t14+lL2.*m3.*t15)];
