function E = energy(in1,in2)
%ENERGY
%    E = ENERGY(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    12-Feb-2017 19:57:23

J1 = in2(:,2);
J2 = in2(:,5);
J3 = in2(:,9);
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
vx = in1(8,:);
vy = in1(9,:);
y = in1(2,:);
t2 = alphaL+phi;
t3 = sin(t2);
t4 = lL2.*t3;
t5 = alphaL+betaL+phi;
t6 = sin(t5);
t7 = l3.*t6;
t14 = sin(phi);
t15 = lH.*t14;
t8 = vy+valphaL.*(t4+t7)+vphi.*(t4+t7+t15)+l3.*t6.*vbetaL;
t9 = cos(t2);
t10 = lL2.*t9;
t11 = cos(t5);
t12 = l3.*t11;
t25 = cos(phi);
t26 = lH.*t25;
t13 = vx+valphaL.*(t10+t12)+vphi.*(t10+t12+t26)+l3.*t11.*vbetaL;
t16 = alphaR+phi;
t17 = sin(t16);
t18 = lL2.*t17;
t19 = alphaR+betaR+phi;
t20 = sin(t19);
t21 = l3.*t20;
t22 = vy+valphaR.*(t18+t21)+vphi.*(t15+t18+t21)+l3.*t20.*vbetaR;
t23 = cos(t16);
t24 = lL2.*t23;
t27 = cos(t19);
t28 = l3.*t27;
t29 = vx+valphaR.*(t24+t28)+vphi.*(t24+t26+t28)+l3.*t27.*vbetaR;
t30 = valphaL+vphi;
t31 = valphaR+vphi;
t32 = valphaL+vbetaL+vphi;
t33 = valphaR+vbetaR+vphi;
t38 = l2.*t9;
t34 = vx+vphi.*(t26+t38)+l2.*t9.*valphaL;
t35 = vy+vphi.*(t15+l2.*t3)+l2.*t3.*valphaL;
t39 = l2.*t23;
t36 = vx+vphi.*(t26+t39)+l2.*t23.*valphaR;
t37 = vy+vphi.*(t15+l2.*t17)+l2.*t17.*valphaR;
E = J2.*t30.^2.*(1.0./2.0)+J2.*t31.^2.*(1.0./2.0)+J3.*t32.^2.*(1.0./2.0)+J3.*t33.^2.*(1.0./2.0)+J1.*vphi.^2.*(1.0./2.0)+m3.*(t8.^2+t13.^2).*(1.0./2.0)+m3.*(t22.^2+t29.^2).*(1.0./2.0)+m2.*(t34.^2+t35.^2).*(1.0./2.0)+m2.*(t36.^2+t37.^2).*(1.0./2.0)+m1.*(vx.^2+vy.^2).*(1.0./2.0)-g.*m2.*(t26+t38-y)-g.*m2.*(t26+t39-y)+g.*m1.*y-g.*m3.*(t10+t12+t26-y)-g.*m3.*(t24+t26+t28-y);
