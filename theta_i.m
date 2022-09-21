function [U]=theta_i(Happ,xi,theta,i,alpha,beta)
% theta is much more complicated than alpha. Solving the pi^theta equations
% yields a transcendental equation. The follwoing code finds the zeros
% based on pi^theta.

[mu0, Msat, rhok1, D, Halpha] = Mat_consts();

D11=D(1,1);
D22=D(2,2);
theta1=theta(1);
theta2=theta(2);
alpha1=alpha(1);
alpha2=alpha(2);
xi1=xi(1);
xi2=xi(2);
Happ1=Happ(1);
Happ2=Happ(2);
k1=1;
rho=rhok1;
bet1=beta;

if i==1  %pi_theta1
    U=- mu0*(Happ2*Msat*(xi1*cos(bet1 - theta1)*(alpha1 - 1) - alpha1*xi1*cos(bet1 + theta1)) ...
        - Happ1*Msat*(xi1*sin(bet1 - theta1)*(alpha1 - 1) - alpha1*xi1*sin(bet1 + theta1)) ...
        + 2*D11*Msat^2*(xi1*sin(bet1 - theta1)*(alpha1 - 1) - alpha1*xi1*sin(bet1 + theta1))*(xi2*sin(theta2) ...
        + xi1*cos(bet1 - theta1)*(alpha1 - 1) + alpha1*xi1*cos(bet1 + theta1)) - 2*D22*Msat^2*(xi1*cos(bet1 - theta1)*(alpha1 - 1) ...
        - alpha1*xi1*cos(bet1 + theta1))*(xi1*sin(bet1 - theta1)*(alpha1 - 1) + xi2*cos(theta2)*(2*alpha2 - 1) ...
        + alpha1*xi1*sin(bet1 + theta1))) - 2*k1*rho*xi1*cos(theta1)*sin(theta1);
 
else %pi_theta2
    U=mu0*(Happ1*Msat*xi2*cos(theta2) - 2*D11*Msat^2*xi2*cos(theta2)*(xi2*sin(theta2) + xi1*cos(bet1 - theta1)*(alpha1 - 1) ...
        + alpha1*xi1*cos(bet1 + theta1)) - Happ2*Msat*xi2*sin(theta2)*(2*alpha2 - 1) ...
        + 2*D22*Msat^2*xi2*sin(theta2)*(2*alpha2 - 1)*(xi1*sin(bet1 - theta1)*(alpha1 - 1) + xi2*cos(theta2)*(2*alpha2 - 1) ...
        + alpha1*xi1*sin(bet1 + theta1))) - 2*k1*rho*xi2*cos(theta2)*sin(theta2);
 
end
   
end