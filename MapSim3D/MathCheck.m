clear
close all
clc



% checking some random jacobians
syms a b c real

n = [a b c]';
kappanorm = sqrt(a^2 + b^2);

phi = atan2(kappanorm,c)/kappanorm*[a; b];

J = jacobian(phi,n);

Jeval = subs(subs(subs(J,b,-0.00000001),a,-0.000000001),c,1);

vpa(Jeval,3)

clear

%linearizing mPerp
syms dd dhat real
nhat = sym('nhat',[3,1],'real');
dphi = sym('dphi',[2,1],'real');
e_nhat = sym('e_nhat',[3,1],'real');
m = sym('m',[3,1],'real');
state = [e_nhat; dd];

d = dhat + dd;
n = simplify(UnitVectorToPlus(nhat)*e_nhat);

mPerp = simplify((n'*m - d)*n);

J = simplify(jacobian(mPerp,state));

J2 = [1 0 0
    0 1 0
    0 0 0
    0 0 1];

Jtotal = J*J2;

Jtotaleval = subs(Jtotal,e_nhat,[0 0 1]');
Jtotaleval = simplify(subs(Jtotaleval, dd, 0));

nhateval = [0 0 1]';
dhateval = 2;
Jtotaltest = subs(Jtotaleval,nhat,nhateval);
Jtotaltest = subs(Jtotaltest,dhat,dhateval);

% test jacobian
mEval = [2 2 3]';
Jeval = double(subs(Jtotaltest,m,mEval));

% evaluate the nomial mPerp
mPerp1 = (nhateval'*mEval - dhateval)*nhateval;

% create a small change in phi and d
tol = 1E-4;

% peturb d
% ddeval = tol;
% dphieval = [0 0]';

% peturb phi(1)
ddeval = 0;
dphieval = [tol 0]';

% peturb phi(2)
% ddeval = 0;
% dphieval = [0 tol]';

stateeval = [dphieval; ddeval];

% predict mPerp2
mPerp2predict = mPerp1 + Jeval*stateeval

% find the real mPerp2
d2 = dhateval + ddeval;
n2 = UnitVectorAdd(nhateval,dphieval);
mPerp2 = (n2'*mEval - d2)*n2