function [ P ] = pixi_i_r(sig, Happ, i, xi, theta, alpha,C,eps_rmc,S1,S2,S3,beta)
% Calculate variant driving force including rotation

[mu0, Msat, rhok1, D, Halpha] = Mat_consts();

% j = 2*(i~=2)+1*(i==2);
% 
% % hardening function
% f_xi=C2/6*xi1^6+C3/2*xi1^2;

% weighted compliance tensor
%S=xi(1)*S1+xi(2)*S2;
S_xi1_1111 = S1(1,1);
S_xi2_1111 = S2(1,1);
S_xi1=S_xi1_1111;
S_xi2=S_xi2_1111;
D11=D(1,1);
D22=D(2,2);
C2=C(2);
C3=C(3);
theta1=theta(1);
theta2=theta(2);
alpha1=alpha(1);
alpha2=alpha(2);
xi1=xi(1);
xi2=xi(2);
sig11=sig(1);
sig22=sig(2);
Happ1=Happ(1);
Happ2=Happ(2);
k1=1;
rho=rhok1;
bet1=beta;
% %% Resume code
if i==1  %pixi_1
 P=eps_rmc*sig11 - C3*xi1 - rho*(k1*sin(theta1)^2 - k1*sin(theta2)^2) + mu0*(Happ1*Msat*(cos(bet1 - theta1)*(alpha1 - 1) - sin(theta2) + alpha1*cos(bet1 + theta1)) + Happ2*Msat*(sin(bet1 - theta1)*(alpha1 - 1) - cos(theta2)*(2*alpha2 - 1) + alpha1*sin(bet1 + theta1)) - 2*D11*Msat^2*(cos(bet1 - theta1)*(alpha1 - 1) - sin(theta2) + alpha1*cos(bet1 + theta1))*(xi1*cos(bet1 - theta1)*(alpha1 - 1) - sin(theta2)*(xi1 - 1) + alpha1*xi1*cos(bet1 + theta1)) - 2*D22*Msat^2*(xi1*sin(bet1 - theta1)*(alpha1 - 1) - cos(theta2)*(2*alpha2 - 1)*(xi1 - 1) + alpha1*xi1*sin(bet1 + theta1))*(sin(bet1 - theta1)*(alpha1 - 1) - cos(theta2)*(2*alpha2 - 1) + alpha1*sin(bet1 + theta1))) - C2*xi1^5 + (sig11^2*(S_xi1 - S_xi2))/2;


else  %pixi_2
P=eps_rmc*sig22 - mu0*(Happ1*Msat*(cos(bet1 - theta1)*(alpha1 - 1) - sin(theta2) + alpha1*cos(bet1 + theta1)) + Happ2*Msat*(sin(bet1 - theta1)*(alpha1 - 1) - cos(theta2)*(2*alpha2 - 1) + alpha1*sin(bet1 + theta1)) + 2*D22*Msat^2*(sin(bet1 - theta1)*(alpha1 - 1) - cos(theta2)*(2*alpha2 - 1) + alpha1*sin(bet1 + theta1))*(sin(bet1 - theta1)*(alpha1 - 1)*(xi2 - 1) + alpha1*sin(bet1 + theta1)*(xi2 - 1) - xi2*cos(theta2)*(2*alpha2 - 1)) + 2*D11*Msat^2*(cos(bet1 - theta1)*(alpha1 - 1) - sin(theta2) + alpha1*cos(bet1 + theta1))*(alpha1*cos(bet1 + theta1)*(xi2 - 1) - xi2*sin(theta2) + cos(bet1 - theta1)*(alpha1 - 1)*(xi2 - 1))) + rho*(k1*sin(theta1)^2 - k1*sin(theta2)^2) - (sig11^2*(S_xi1 - S_xi2))/2 - (C3*(2*xi2 - 2))/2 - C2*(xi2 - 1)^5;

 
end

end

