function [xi, theta, alpha, pixis, csi, fii] = M_H_solve_int_r( sig,dsig,dHapp,dsig11, dsig22, Happ, dHapp1, dHapp2, xi0, theta0, C, eps_rmc,S1,S2,S3,tol,beta)
% This function solves for all internal variables assuming rotation.
Y = C(1);
% A is just a function that restricts all values in an array or matrix to
% be between 'a' and 'b'
A = @(X, a, b) arrayfun(@(x) max(min(x,b),a), X);
ulxi=1;% upper & lower limit of xi
llxi=1-ulxi;
%upper & lower limit of theta
ultheta2=pi/2;
lltheta2=0;
if beta<0 || beta>0
  ultheta1=pi/2-beta;
else
  ultheta1=pi/2;
end
lltheta1=0;
%Set fsolve options
options = optimset('Display', 'off');

% Alpha is a function that finds alpha based on Happ
alpha = [alpha_i_r(Happ,1);alpha_i_r(Happ,2)];

% solves for initial rotation 
theta00=@(thet)[theta_i(Happ,xi0,thet,1,alpha,beta);theta_i(Happ,xi0,thet,2,alpha,beta)];
theta=fsolve(@(thet) theta00(thet),theta0,options);
theta(1)=A(theta(1),lltheta1,ultheta1);
theta(2)=A(theta(2),lltheta2,ultheta2);
theta=[theta(1);theta(2)];

% Find pixi0 based on above quantities
pixi0=[pixi_i_r( sig, Happ, 1, xi0, theta, alpha,C,eps_rmc,S1,S2,S3,beta);...
       pixi_i_r( sig, Happ, 2, xi0, theta, alpha,C,eps_rmc,S1,S2,S3,beta)];

[pixi_p, p] = max(pixi0);  % p = max
[pixi_q, q] = min(pixi0);  % q = min
% % if reorientation not initiated, xi=xi00=xi0
xi=xi0;
%elaticity check    [fi_f_n=pi1-pi2-Y;   fi_r_n=-(pi1-pi2)-Y]
fii=(pixi_i_r( sig, Happ, p, xi, theta, alpha,C,eps_rmc,S1,S2,S3,beta)-...
       pixi_i_r( sig, Happ, q, xi, theta, alpha,C,eps_rmc,S1,S2,S3,beta)-Y);

%load check      sii_n=pi1-pi2-y    without xi term
diff_load_1=(pixi_i_r(sig, Happ, 1, [0;0], [0;0], [0;0],C,eps_rmc,S1,S2,S3,beta)-...
       pixi_i_r(sig, Happ, 2, [0;0], [0;0], [0;0],C,eps_rmc,S1,S2,S3,beta)-Y);
   % sii_(n-1)=pi1-pi2-Y    without xi term
diff_load_0=(pixi_i_r((sig-dsig), (Happ-dHapp), 1, [0;0], [0;0], [0;0],C,eps_rmc,S1,S2,S3,beta)-...
       pixi_i_r((sig-dsig), (Happ-dHapp), 2, [0;0], [0;0], [0;0],C,eps_rmc,S1,S2,S3,beta)-Y);
 
csi=diff_load_1-diff_load_0; %loading index %sii_n-sii_(n-1) 

pixis = arrayfun(@(n) pixi_i_r(sig,Happ,n,xi,theta,alpha,C,eps_rmc,S1,S2,S3,beta), [1;2]);
end