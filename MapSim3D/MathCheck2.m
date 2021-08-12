clear
close all
clc

% Define symbolics
m0 = sym('m0_',[3,1],'real'); % sensor position
bhat = sym('bhat_',[3,1],'real'); % measurement bearing
syms d dd 'real' % plane(4) and peturbation
nhat = sym('nhat_',[3,1],'real'); %plane normal
dphi = sym('dphi_',[2 1],'real'); %plane normal error
dphihat = sym('dphihat_',[3 1],'real'); % dphi as a unit vector

%state
state = [dphihat; dd];

%range between sensor and plane
t = (d + dd - (UnitVectorToPlus(nhat)*dphihat)'*m0) / ((UnitVectorToPlus(nhat)*dphihat)'*bhat);

%predicted mearurement
mPred = m0 + t*bhat;

% Jacobian
J1 = jacobian(mPred,state);
J1 = subs(J1,dd,0); % Expected dd = 0
J1 = subs(J1,dphihat,[0 0 1]'); % Expected dphihat

% This is the jacobian of the augmented state with respect to the true
% state
J2 = [1 0 0
    0 1 0
    0 0 0
    0 0 1];

% Final J
J = J1*J2;

% Convert to a matlab function
Jfunc = matlabFunction(J,'File','JacobianEval');

% Use the jacobian to predict changes in mPred
nhat_eval = [.1 .1 1]';
nhat_eval = nhat_eval/norm(nhat_eval);
d_eval = 1.5;
m0_eval = [2 2 10]';
bhat_eval = [-.1 .1 -1]';
bhat_eval = bhat_eval/norm(bhat_eval);

% Evaluate jacobian
J_eval = JacobianEval(bhat_eval(1),bhat_eval(2),bhat_eval(3),...
    d_eval,...
    m0_eval(1),m0_eval(2),m0_eval(3),...
    nhat_eval(1),nhat_eval(2),nhat_eval(3));

% find first predicted point
t1 = (d_eval - nhat_eval'*m0_eval)/(nhat_eval'*bhat_eval);
mPred1 = m0_eval + t1*bhat_eval;

% create a small change in phi and d
tol = 1E-4;

% peturb d
% dd_eval = tol;
% dphi_eval = [0 0]';

% peturb phi(1)
% dd_eval = 0;
% dphi_eval = [tol 0]';

% peturb phi(2)
dd_eval = 0;
dphi_eval = [0 tol]';

% concatenate into state
state_eval = [dphi_eval; dd_eval];

% find the second predicted point
t2 = (d_eval + dd_eval - UnitVectorAdd(nhat_eval,dphi_eval)'*m0_eval) / (UnitVectorAdd(nhat_eval,dphi_eval)'*bhat_eval);
mPred2 = m0_eval + t2*bhat_eval

% predict second point with jacobian
mPred2Jac = mPred1 + J_eval*state_eval

% Find the error
errJac = mPred2 - mPred2Jac






