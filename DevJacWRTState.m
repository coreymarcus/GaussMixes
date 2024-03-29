function J_wrt_pt = DevJacWRTState(dphihat_1,dphihat_2,dphihat_3,nhat_1,nhat_2,nhat_3)
%DevJacWRTState
%    J_wrt_pt = DevJacWRTState(DPHIHAT_1,DPHIHAT_2,DPHIHAT_3,NHAT_1,NHAT_2,NHAT_3)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    17-Nov-2022 12:24:58

t2 = nhat_3+1.0;
t3 = 1.0./t2;
J_wrt_pt = [dphihat_3.*nhat_1-dphihat_1.*(nhat_1.^2.*t3-1.0)-dphihat_2.*nhat_1.*nhat_2.*t3,dphihat_3.*nhat_2-dphihat_2.*(nhat_2.^2.*t3-1.0)-dphihat_1.*nhat_1.*nhat_2.*t3,-dphihat_1.*nhat_1-dphihat_2.*nhat_2+dphihat_3.*nhat_3];
