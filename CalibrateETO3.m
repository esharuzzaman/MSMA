% This is the main script for calibrating a 2-var MSMA specimen for Model 8
%
% Inputs needed:
%   calibration file
%   - The key points to pick:
%           - Location of a point in initial elastic region
%           - Location of two points in the final elastic region
%           - Location of 3 points during reorientation

clear;  
close all;
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mu0, Msat, rhok1, D] = Mat_consts();
D11 = D(1,1);
D22 = D(2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read in data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stress2 = xlsread('102719 calibration.xlsx', 'Processed data', 'B3:B3046');
stress2 = 1e6*stress2;
strain2 = xlsread('102719 calibration.xlsx', 'Processed data', 'A3:A3046');
strain2 = strain2/100;

%%%%%%%%%%%%%%%%%%%%%%%%% Selecting Data Points %%%%%%%%%%%%%%%%%%%%%%%%%
% Total size of calibration data
N2 = length(stress2);

% pick two points during initial and final elastic regions
% for determining elastic moduli
% ~0-60, ~1600-3137
E2_1=264;
E2_2=274;
E1_1=2470; 
E1_2=2483;
    
% pick three points during reorientation
% ~100-700
dot1 = 393;
dot2 = 455;
dot3 = 825;

% points in calibration for moduli
in_pt = [stress2(E2_1),strain2(E2_1)];
E2_pt = [stress2(E2_2),strain2(E2_2)];
E1_pt1 = [stress2(E1_1),strain2(E1_1)];
E1_pt2 = [stress2(E1_2),strain2(E1_2)];

% points in calibration for determination of constants
point1 = [stress2(dot1),strain2(dot1)]; %stress strain in cal. data 2 to 1
point2 = [stress2(dot2),strain2(dot2)]; %stress strain in cal. data 2 to 1
point3 = [stress2(dot3),strain2(dot3)]; %stress strain in cal. data 2 to 1
% point4 = [stress2(dot4),strain2(dot4)]; %stress strain in cal. data 2 to 1

% plot above points against calibration data
figure;
plot(strain2, stress2)
hold on
scatter([point1(2) point2(2) point3(2)], [point1(1) point2(1) point3(1)],'o')
hold on
scatter([in_pt(2) E2_pt(2) E1_pt1(2) E1_pt2(2)], [in_pt(1) E2_pt(1) E1_pt1(1) E1_pt2(1)], 'k', 'o')
hold off

%%%%%%%%%%%%%%%%%%%%%% Calculate Elastic Modulii %%%%%%%%%%%%%%%%%%%%%%%%%

% From above points
E1_1 = -(E1_pt2(1)-E1_pt1(1))/(E1_pt2(2)-E1_pt1(2));  %Pa
E1_2 = ((E2_pt(1)-in_pt(1))/(E2_pt(2)-in_pt(2)));
E1_3 = E1_2; % placeholder, not actually tested here

% E1_1=abs(E1_1);
% E1_2=abs(E1_2);
% E1_3=abs(E1_3);

nu1_12 = .3; %not measured yet, just a guess...
nu1_13 = nu1_12;
nu1_23 = .3; %Again, this is just a guess here
nu1_21 = E1_2/E1_1*nu1_12; %Doug's code
nu1_31 = E1_3/E1_1*nu1_13; %Doug's code
nu1_32 = E1_3/E1_2*nu1_23; %Doug's code
S1 = [1/E1_1, -nu1_21/E1_2;
    -nu1_12/E1_1,       1/E1_2];
S2 = [S1(2,2) S1(2,1);
    S1(1,2) S1(1,1)];
S3=0;

S_xi1_1111 = 1/E1_1; % [m^2/N] %3var
S_xi2_1111 = 1/E1_2; % [m^2/N] %3var

S_xi1 = S1;
S_xi2 = S2;

sig_pt1 = [point1(1) 0; 0 0];
sig_pt2 = [point2(1) 0; 0 0];
sig_pt3 = [point3(1) 0; 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%% Calculate eps_rmc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_rmax = min(strain2)  - min(stress2)/E1_1

%%%%%%%%%%%%%%%%%%%%%%%% Calculate C's & Y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi1_1 = (point1(2)-S_xi2_1111*point1(1))/((S_xi1_1111-S_xi2_1111)*point1(1)+eps_rmax)
xi1_2 = (point2(2)-S_xi2_1111*point2(1))/((S_xi1_1111-S_xi2_1111)*point2(1)+eps_rmax)
xi1_3 = (point3(2)-S_xi2_1111*point3(1))/((S_xi1_1111-S_xi2_1111)*point3(1)+eps_rmax)

prompt = 'Are these xi1s OK?';
s = input(prompt)

% known values
b =    [point1(1)*eps_rmax + point1(1)^2*(S_xi1_1111-S_xi2_1111);
        point2(1)*eps_rmax + point2(1)^2*(S_xi1_1111-S_xi2_1111);
        point3(1)*eps_rmax + point3(1)^2*(S_xi1_1111-S_xi2_1111)];    

% coefficients of the constants (Y, C1, C2)
A = [ 1, 2*xi1_1^5, 2*xi1_1;
      1, 2*xi1_2^5, 2*xi1_2;
      1, 2*xi1_3^5, 2*xi1_3];

% solves for [Y; C1; C2]
C = A\b;

disp('The calbration results in: ')
C
prompt = 'Is this OK?';
    s = input(prompt)

%%%%%%%%%%%%%%%%%%%%%%%%% Check Calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = C(1)

if Y < 0
    Y
    disp('Error Y negative');
    prompt = 'Better stop now - OK?';
    s = input(prompt)
end

%%%%%%%%%%%%%%%%%% Write calibration constants to file %%%%%%%%%%%%%%%%%%%

csvwrite('output_calibrate_s2.csv',[C', eps_rmax, E1_1, E1_2, E1_3]);


%%%%%%%%%%%%%%%%%%%%% Simulate the Calibration Data %%%%%%%%%%%%%%%%%%%%%


% Resolution & Tolerance
res = 1e4
tol =1 ;
A=10.5625e-6;    %cross section area
nt=1000; %number of turns

calibrate_results = csvread('output_calibrate_s2.csv');

C = calibrate_results(1:3)
Y = calibrate_results(1)
eps_rmc = calibrate_results(4)
E1_1 = calibrate_results(5)
E1_2 = calibrate_results(6)
E1_3 = calibrate_results(7)


% loadtypes =    18;
% lat_stresses = 0;
% fields =       0;
% axials =       4;
%fi=21;
        %fields=0.5;
        f2=0
        f1=0
beta=(pi/180)*(8)
        rf=f2/f1;
loadtypes=18
lat_stresses=0
fields=f2;
axials=4 

for counter = 1:length(loadtypes) 
    loadtype = loadtypes(counter)
    lat_stress = lat_stresses(counter)
    field = fields(counter)/(4*pi*10^-7);
    axial =  axials(counter)
   
    %Start timer
    tic
    
    % Get loading profile
    [ sig, Happ, time, xi00,dt ] = loading(loadtype, res, lat_stress, field, axial,A,rf);
    
    % Create empty matrices for internal variables
    N = size(Happ,2);
    % N is the number of steps in the loading.  In loading.m sig & Happ to be
    % the applied stress & field in each step.
    
    % Initialize values
    pixis = zeros(2,N);
    alphas = zeros(2,N);
    xis = zeros(2,N);
    thetas = zeros(2,N);
    M = zeros(2,N);
    eps_tot = zeros(2,N);
    eps_tot(:,1)=[eps_rmc;0];
    
    dM1=zeros(1,N);
    emf=zeros(1,N);
    
    dsig11=zeros(1,N);
    dsig22=zeros(1,N);
    dsig=zeros(2,N);
    
    dHapp1=zeros(1,N);
    dHapp2=zeros(1,N);
    dHapp=zeros(2,N);
    
    dxi=zeros(2,N);
    dff=zeros(1,N);
    
    pixis_pp=zeros(1,N);
    pixis_qq=zeros(1,N);
    pq=zeros(1,N);
    df=zeros(1,N);
    dpixi=zeros(1,N);
    save_f=zeros(1,N);

    csi=zeros(1,N);
    fii=zeros(1,N);

    % Set initial xi and U values
    xi0 = xi00;
    theta0 = [0;0];
    xis(:,1)=xi0; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Find initial pixi and alpha values
    n=1;
     dsig11(:,n)=0;
     dsig22(:,n)=0;
     dsig(:,n)=[dsig11(:,n);dsig22(:,n)];
        
     dHapp1(:,n)=0;
     dHapp2(:,n)=0;
     dHapp(:,n)=[dHapp1(:,n);dHapp2(:,n)];
    
    [xi0, theta, alphas(:,n), pixis(:,n),csi(:,n),fii(:,n)] = ...
        solve_int_r(sig(:,n),dsig(:,n),dHapp(:,n),dsig11(:,n),dsig22(:,n),Happ(:,n),dHapp1(:,n),dHapp2(:,n),xi0,theta0,C,eps_rmc,S1,S2,S3,tol,beta);
    xis(:,n) = xi0;
    thetas(:,n) = theta;
    theta0=theta;
    save_f(:,n) = (max(pixis(:,n))-min(pixis(:,n)))/Y - 1; %ploting max min pixis
    
    eps_tot(:,n) = xi0(1)*S1*sig(:,n) + xi0(2)*S2*sig(:,n) + eps_rmc * xi0;
        
%     M(:,n)=[xi0(1)*cos(theta(1))*(2*alphas(1,n)-1)+xi0(2)*sin(theta(2));...
%                 xi0(1)*sin(theta(1))+xi0(2)*cos(theta(2))*(2*alphas(2,n)-1)].*Msat;
M(:,n)=[xi0(1)*alphas(1,n)*cos(theta(1)+beta)-xi0(1)*(1-alphas(1,n))*cos(theta(1)-beta)+xi0(2)*sin(theta(2));...
    xi0(1)*alphas(1,n)*sin(theta(1)+beta)+xi0(1)*(1-alphas(1,n))*sin(theta(1)-beta)+xi0(2)*(2*alphas(2,n)-1)*cos(theta(2))].*Msat;
 %emf
     dM1(:,n)=M(1,n)-0;
     emf(:,n)=(mu0*(1-D(1,1))*A*nt/dt)*(M(1,n)-0);
    % This is the main loop that solves for internal variables at
    % timestep n
    for n = 2:N
        % If we are 1%, 2%, etc. done with the model, say so
        if(mod(n,round(N/100))==0)
            fprintf('\n%2.0f %% %G min\n', n/N*100, round((toc/60),2));
        end
        % We outsource all the conditionals to the solve_int_r (solve
        % internal variables) function
        dsig11(:,n)=sig(1,n)-sig(1,n-1);
        dsig22(:,n)=sig(2,n)-sig(2,n-1);
        dsig(:,n)=[dsig11(:,n);dsig22(:,n)];
        
        dHapp1(:,n)=Happ(1,n)-Happ(1,n-1);
        dHapp2(:,n)=Happ(2,n)-Happ(2,n-1); 
        dHapp(:,n)=[dHapp1(:,n);dHapp2(:,n)];

[xi0, theta, alphas(:,n), pixis(:,n),csi(:,n),fii(:,n)] = ...
            solve_int_r(sig(:,n),dsig(:,n),dHapp(:,n),dsig11(:,n),dsig22(:,n),Happ(:,n),dHapp1(:,n),dHapp2(:,n),xi0,theta0,C,eps_rmc,S1,S2,S3,tol,beta);
        
        xis(:,n) = xi0;
        thetas(:,n) = theta;
        theta0=theta;
              
save_f(:,n) = (max(pixis(:,n))-min(pixis(:,n)))/Y - 1; %ploting max min pixis
      
        eps_tot(:,n) = xi0(1)*S1*sig(:,n) + xi0(2)*S2*sig(:,n) + eps_rmc * xi0;
        
%         M(:,n)=[xi0(1)*cos(theta(1))*(2*alphas(1,n)-1)+xi0(2)*sin(theta(2));...
%                 xi0(1)*sin(theta(1))+xi0(2)*cos(theta(2))*(2*alphas(2,n)-1)].*Msat;
M(:,n)=[xi0(1)*alphas(1,n)*cos(theta(1)+beta)-xi0(1)*(1-alphas(1,n))*cos(theta(1)-beta)+xi0(2)*sin(theta(2));...
    xi0(1)*alphas(1,n)*sin(theta(1)+beta)+xi0(1)*(1-alphas(1,n))*sin(theta(1)-beta)+xi0(2)*(2*alphas(2,n)-1)*cos(theta(2))].*Msat;

%     %emf
     dM1(:,n)=M(1,n)-M(1,n-1);
     emf(:,n)=(mu0*(1-D(1,1))*A*nt/dt)*(M(1,n)-M(1,n-1));
    end
    % Write values to a CSV file for plotting later
%     csvwrite(sprintf('output %i - %i - %.1f - %.1f - %.1f - %i time_emf_Sig1_sig2_Happ1_Happ2_xi1_xi2_theta1_theta2_pixi1_pixi2.csv',res, loadtype, lat_stress, axial, field*mu0, j),[time;emf;sig;Happ;xis;thetas;pixis]');
expe=xlsread('exp176pkpk'); %read experimental emf and time
time_exp=(expe(:,1))';
emf_exp=(expe(:,2))';
  

fprintf('\n Simulation ran in %G minutes\n',round((toc/60),2))
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    sp3 = subplot(4,1,3);
    p3 = plot(time, (Happ(1,:)*mu0)', 'g');
    hold on
    p3 = plot(time, (Happ(2,:)*mu0)', 'b');
    ylabel('Happ');
    xlabel('time');
    
    sp2 = subplot(4,1,2);
    hold on;
    p2 = plot(time, (sig(1,:))', 'g');
    hold on;
    p2 = plot(time, (sig(2,:))', 'b');
    ylabel('sigs');
    xlabel('time');
    
    
    sp1 = subplot(4,1,1);
    hold on;
    p1 = plot(time, (xis(1,:))', 'g');
    hold on;
    p1 = plot(time, (xis(2,:))', 'b');
    hold on;
    ylabel('xis');
    xlabel('time');
    l2=legend('green is 1-dir','blue is 2-dir');
    title(sprintf('case %d', loadtype), 'fontsize',16);
    
save_f = [];
    for i=1:N
        save_f = [save_f, (max(pixis(:,i))-min(pixis(:,i)))/Y - 1];
    end
    sp4 = subplot(4,1,4);
    p4 = plot(time, save_f, 'k');
    ylabel('f = (pi_p - pi_q)/Y - 1');
    xlabel('time');    
% simulate calibration test
    figure;
    plot(-100*eps_tot(1,:),-1e-6*sig(1,:),'r','LineWidth',1.5)
    hold on
    plot(-100.*strain2,-1e-6.*stress2,'b:','LineWidth',1.5)
    hold on
    plot(-100*point1(2),-1e-6*point1(1),'ro')
    hold on
    plot(-100*point2(2),-1e-6*point2(1),'ro')
    hold on
    plot(-100*point3(2),-1e-6*point3(1),'ro')

    title('Calibration Direction-2','FontSize',16)
    xlabel('-\epsilon_1_1 (%)','fontWeight','bold','FontSize',14)
    ylabel('-\sigma_1_1 (MPa)','fontWeight','bold','FontSize',14)
    h=legend('Model','Experimental Data','Calibration points','location','northwest');
    set(h,'fontsize',10)
    axis([0,7,0,5]) 
    
end


