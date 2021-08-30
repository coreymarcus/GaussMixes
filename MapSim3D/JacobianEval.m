function J = JacobianEval(bhat_1,bhat_2,bhat_3,d,m0_1,m0_2,m0_3,nhat_1,nhat_2,nhat_3)
%JACOBIANEVAL
%    J = JACOBIANEVAL(BHAT_1,BHAT_2,BHAT_3,D,M0_1,M0_2,M0_3,NHAT_1,NHAT_2,NHAT_3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    11-Aug-2021 19:28:24

t2 = nhat_3+1.0;
t3 = 1.0./t2;
t4 = nhat_1.^2;
t5 = t3.*t4;
t6 = t5-1.0;
t7 = bhat_1.*nhat_1;
t8 = bhat_2.*nhat_2;
t9 = bhat_3.*nhat_3;
t10 = t7+t8+t9;
t11 = 1.0./t10;
t12 = nhat_2.^2;
t13 = t3.*t12;
t14 = t13-1.0;
t15 = 1.0./t10.^2;
t16 = m0_1.*nhat_1;
t17 = m0_2.*nhat_2;
t18 = m0_3.*nhat_3;
t19 = -d+t16+t17+t18;
t20 = m0_3.*nhat_1;
t21 = m0_1.*t6;
t22 = m0_2.*nhat_1.*nhat_2.*t3;
t23 = t20+t21+t22;
t24 = bhat_3.*nhat_1;
t25 = bhat_1.*t6;
t26 = bhat_2.*nhat_1.*nhat_2.*t3;
t27 = t24+t25+t26;
t28 = m0_3.*nhat_2;
t29 = m0_2.*t14;
t30 = m0_1.*nhat_1.*nhat_2.*t3;
t31 = t28+t29+t30;
t32 = bhat_3.*nhat_2;
t33 = bhat_2.*t14;
t34 = bhat_1.*nhat_1.*nhat_2.*t3;
t35 = t32+t33+t34;
J = reshape([bhat_1.*t11.*t23-bhat_1.*t15.*t19.*t27,bhat_2.*t11.*t23-bhat_2.*t15.*t19.*t27,bhat_3.*t11.*t23-bhat_3.*t15.*t19.*t27,bhat_1.*t11.*t31-bhat_1.*t15.*t19.*t35,bhat_2.*t11.*t31-bhat_2.*t15.*t19.*t35,bhat_3.*t11.*t31-bhat_3.*t15.*t19.*t35,bhat_1.*t11,bhat_2.*t11,bhat_3.*t11],[3,3]);