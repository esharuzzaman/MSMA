% This is the main script for simulating the behavior of MSMAs assuming
% full evolution and normal demagnetization (Hd=f(M))
% Inputs needed:
%   - calibration_results.csv
%   - Resolution & Tolerance
%   - Loadtypes, lat_stresses, axials & fields

clear; clc
tic
%%%%%%%%%%%%%%%%%%% Define Resolution & Tolerance & sample %%%%%%%%%%%%%%%%
sample_no = 3;
res = 1e3;
tol =1 ;
A=8.26e-6;    %cross section area
nt=2010; %number of turns
emfll=-0.2; %upper & lower limit of emf
emful=0.2;

texpll=0;%experimental limit
texpul=.4;

%load types
% 50-experimental data 10HZ 5.85 deg H2=0.537 H1=0.055
% 51-experimental data 02HZ 5.85 deg H2=0.537 H1=0.055
% 52-experimental data 02HZ 0 deg H2=0.54 H1=0
% 54-experimental data 02HZ 0 deg H2=0.5 H1=0
%200-experimental data 10HZ 2.87 deg H2=0.578 H1=0.029
%201-experimental data 10HZ 5.61 deg H2=0.575 H1=0.057
%202-experimental data 10HZ 5.61 deg H2=0.575 H1=0.057

%inputs
H=0.578;
H1=0.029;
beta_d=-8;
% H2=.3
H2=sqrt(H^2-H1^2);        
beta=(pi/180)*(beta_d);
rf=H2/H1;
angle=atand(H1/H2);

loadtypes=202;
lat_stresses=0;
fields=H2;
axials=5;

% current const. field tests use 15 axial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load constants %%%%%%%%%%%%%%%%%%%%%%%%%%
[mu0, Msat, rhok1, D] = Mat_consts();
D11 = D(1,1);
D22 = D(2,2);

calibrate_results = csvread('output_calibrate_s2.csv');

    
C = calibrate_results(1:3);
Y = calibrate_results(1);
eps_rmc = calibrate_results(4);
E1_1 = calibrate_results(5);
E1_2 = calibrate_results(6);
E1_3 = calibrate_results(7);


nu1_12 = .3; %not measured yet, just a guess...
nu1_13 = nu1_12;
nu1_23 = .3; %Again, this is just a guess here
nu1_21 = E1_2/E1_1*nu1_12; %Doug's code
nu1_32 = E1_3/E1_2*nu1_23; %Doug's code
S1 = [1/E1_1, -nu1_21/E1_2;
    -nu1_12/E1_1, 1/E1_2];
S2 = [S1(2,2) S1(2,1);
    S1(1,2) S1(1,1)];

S3 = 0; % not used

for counter = 1:length(loadtypes) 
    loadtype = loadtypes(counter);
    lat_stress = lat_stresses(counter);
    field = fields(counter)/(4*pi*10^-7);
    axial =  axials(counter);
   
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
     emf(:,n)=-(mu0*(1-D(1,1))*A*nt/dt)*(M(1,n)-0);
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
     emf(:,n)=-(mu0*(1-D(1,1))*A*nt/dt)*(M(1,n)-M(1,n-1));
    end
    % Write values to a CSV file for plotting later
%     csvwrite(sprintf('output %i - %i - %.1f - %.1f - %.1f - %i time_emf_Sig1_sig2_Happ1_Happ2_xi1_xi2_theta1_theta2_pixi1_pixi2.csv',res, loadtype, lat_stress, axial, field*mu0, j),[time;emf;sig;Happ;xis;thetas;pixis]');
% expe=xlsread('exp176pkpk'); %read experimental emf and time
% time_exp=(expe(:,1))';
% emf_exp=(expe(:,2))';
%  emf_exp = (xlsread('dtr_10hz_0.055H1_0.537H2.xlsx', 'Sheet1', 'J2:J26'))'; 
%  time_exp = (xlsread('dtr_10hz_0.055H1_0.537H2.xlsx', 'Sheet1', 'K2:K26'))'; 
%  stress11 = (xlsread('dtr_10hz_0.055H1_0.537H2.xlsx', 'Sheet1', 'E1:E82'))'; 
%  time_exp1 = (xlsread('dtr_10hz_0.055H1_0.537H2.xlsx', 'Sheet1', 'D1:D82'))'; 

fprintf('\n Simulation ran in %G minutes\n',round((toc/60),2))
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	% output some plots to ensure code is operating correctly
    figure
    plot(eps_tot(1,:),sig(1,:))
    title(sprintf('Beta= %G deg; H_1= %G T; H_2= %G T; Phi= %G deg; rhok1= %e', beta_d, H1,H2,angle,rhok1), 'fontsize',10)
    xlabel('$\varepsilon_{11}$','interpreter','latex')
    ylabel('$\sigma_{11}$ (Pa)','interpreter','latex')
    grid minor
% %     
% %     figure
% %     plot(time,xis)
% %     legend('xi1','xi2')
% %     title('xi vs time')
% %     xlabel('time(sec)')
% %     ylabel('xi')
% %     grid minor
% %     
% %     figure
% %     plot(time,pixis)
% %     legend('pixi1','pixi2')
% %     title('pixi vs time')
% %     xlabel('time(sec)')
% %     ylabel('pixi')
% %     grid minor
% %     
    figure
    plot(eps_tot(1,:),Happ(2,:).*4*pi*10^-7)
    title(sprintf('Beta= %G deg; H_1= %G T; H_2= %G T; Phi= %G deg; rhok1= %e', beta_d, H1,H2,angle,rhok1), 'fontsize',10)
    xlabel('$\varepsilon_{11}$','interpreter','latex')
    ylabel('H_2 (T)')
    grid minor
    
%     figure
%     plot(time,thetas)
%     legend('theta1','theta2')
%     title('Theta vs Time')
%     xlabel('time(sec)')
%     ylabel('theta')
%     grid minor
%     
% %     figure
% %     plot(time,M)
% %     legend('M1','M2')
% %     title('M vs time')
% %     xlabel('time(sec)')
% %     ylabel('M (A/m)')
% %     grid minor
%     
%     figure
%     plot(time,(mu0*Happ))
%     legend('Happ1','Happ2')
%     title('Happ vs time')
%     xlabel('time(sec)')
%     ylabel('Happ(Tesla)')
%     grid minor
%     
    figure
    plot(time,sig)
    legend('sig1','sig2')
    title('Sigma vs Time')
    xlabel('time(sec)')
    ylabel('sig(Pascals)')
    grid minor
%     
figure
    plot(time,emf)
    ylim([emfll emful])
    title('EMF vs Time')
    xlabel('time(sec)')
    ylabel('emf(volt)')
    grid minor
%     
%     figure
%     plot(time,dM1)
%     title('dM1 vs time')
%     xlabel('time(sec)')
%     ylabel('dM1(A/m)')
%     grid minor
%     
% %     figure
% %     plot(time,emf,time,sig(1,:),time,Happ)
% %     legend('emf','sig1','Happ1','Happ2')
% %     title('emf,sig1,Happ vs time')
% %     xlabel('time(sec)')
% %     ylabel('emf(volt),sig1(pascals),Happ(A/m)')
% %     grid minor
%     
    figure
    yyaxis left
    plot(time,sig)
    %legend('sig1','sig2')
    title('Sigma & Happ vs Time')
    xlabel('time(sec)')
    ylabel('sig(Pascals)')
    grid minor
    
    yyaxis right
    plot(time,(mu0*Happ))
    %legend('Happ1','Happ2')
    title('Sig & Happ vs time')
    xlabel('time(sec)')
    ylabel('Happ(Tesla)')
    grid minor
    legend('sig1','sig2','Happ1','Happ2')
    
    figure
    yyaxis left
    plot(time,sig)
    legend('sig1','sig2')
    title('Sigma & EMF vs Time')
    xlabel('time (sec)')
    ylabel('sig (Pascals)')
    grid minor
    
    yyaxis right
    plot(time,emf)
    %title('EMF vs Time')
    xlabel('time (sec)')
    ylabel('EMF (volt)')
    grid minor
    legend('sig1','sig2','EMF')
    
    figure
    subplot(5,1,1)
    plot(time,(mu0*Happ))
    legend('Happ1','Happ2')
    xlabel('time(sec)')
    ylabel('Happ(Tesla)')
    title(sprintf('case %d ; run time %G min, res %G', loadtype, round((toc/60),2),res), 'fontsize',12);
    grid minor
    
    subplot(5,1,2)
    plot(time,sig)
    legend('sig1','sig2')
    xlabel('time(sec)')
    ylabel('sig(Pascals)')
    grid minor
    
    subplot(5,1,4)
    plot(time,alphas)
    legend('alpha1','alpha2')
    title('Alpha vs Time')
    xlabel('time(sec)')
    ylabel('alpha')
    grid minor
    
    subplot(5,1,5)
    plot(time,thetas)
    legend('theta1','theta2')
    xlabel('time(sec)')
    ylabel('theta')
    grid minor
    
    subplot(5,1,3)
    plot(time,xis)
    legend('xi1','xi2')
    xlabel('time(sec)')
    ylabel('xi')
    grid minor
% %     
% %     figure
% %     subplot(4,1,1)
% %     plot(time,csi)
% %     xlabel('time(sec)')
% %     ylabel('csi')
% %     title(sprintf('case %d ; run time %G min, res %G', loadtype, round((toc/60),2),res), 'fontsize',12);
% %     grid minor
% %     
% %     subplot(4,1,2)
% %     plot(time,fii)
% %     xlabel('time(sec)')
% %     ylabel('fii')
% %     grid minor
% %     
% %     subplot(4,1,3)
% %     plot(time,xis)
% %     legend('xi1','xi2')
% %     xlabel('time(sec)')
% %     ylabel('xi')
% %     grid minor
% %     
% %     subplot(4,1,4)
% %     plot(time,thetas)
% %     legend('theta1','theta2')
% %     xlabel('time(sec)')
% %     ylabel('theat(rad)')
% %     grid minor
% %  
% %     figure;
% %     sp2 = subplot(5,1,1);
% %     yyaxis left
% %     hold on;
% %     p2 = plot(time, (sig(1,:))', 'g');
% %     hold on;
% %     p2 = plot(time, (sig(2,:))', 'b');
% %     ylabel('sigs');
% %     xlabel('time');
% %     title(sprintf('case %d ; run time %G min, res %G', loadtype, round((toc/60),2),res), 'fontsize',12);
% %     grid on;
% %     grid minor;
% %     
% %     yyaxis right
% %    plot(time,Happ.*4*pi*10^-7)
% %     ylabel('Happ');
% %     xlabel('time');
% %     legend('sig1','sig2','Happ1','Happ2')
% %     grid on;
% %     grid minor;
% %     
% %     sp1 = subplot(5,1,3);
% %     hold on;
% %     p1 = plot(time, (xis(1,:))', 'g');
% %     hold on;
% %     p1 = plot(time, (xis(2,:))', 'b');
% %     hold on;
% %     ylabel('xis');
% %     xlabel('time');
% %     l2=legend('xi1','xi2');
% %     grid on;
% %     grid minor;
% % 
% %     sp4 = subplot(5,1,2);
% %     p3 = plot(time, save_f, 'r');
% %     ylabel('(pi^p- pi^q)/Y - 1');
% %     xlabel('time');
% %     grid on;
% %     grid minor;
% %     
% %     subplot(5,1,4)
% %     plot(time,alphas)
% %     legend('alpha1','alpha2')
% %     xlabel('time(sec)')
% %     ylabel('alpha')
% %     grid minor
% %     
% %     subplot(5,1,5)
% %     plot(time,thetas)
% %     legend('theta1','theta2')
% %     xlabel('time(sec)')
% %     ylabel('theta(rad)')
% %     grid minor
% %     
  figure;
    sp2 = subplot(5,1,1);
    p2 = plot(time, (sig(1,:))', 'g');
    hold on;
    p2 = plot(time, (sig(2,:))', 'b');
    ylabel('sigs');
    xlabel('time');
    title(sprintf('case %d ; run time %G min, res %G', loadtype, round((toc/60),2),res), 'fontsize',12)
    grid on;
    grid minor;
    
    subplot (5,1,2)
    plot(time,emf)
    ylim([emfll emful])
    title('EMF vs Time')
    xlabel('time(sec)')
    ylabel('emf(volt)')
    grid minor
    
    sp1 = subplot(5,1,4);
    hold on;
    p1 = plot(time, (xis(1,:))', 'g');
    hold on;
    p1 = plot(time, (xis(2,:))', 'b');
    hold on;
    ylabel('xis');
    xlabel('time');
    l2=legend('xi1','xi2');
    grid on;
    grid minor;

    sp4 = subplot(5,1,3);
    p3 = plot(time, save_f, 'r');
    ylabel('(pi^p- pi^q)/Y - 1');
    xlabel('time');
    grid on;
    grid minor;
    
    subplot(5,1,5)
    plot(time,thetas)
    legend('theta1','theta2')
    xlabel('time(sec)')
    ylabel('theta(rad)')
    grid minor  
    
    
    figure;
    subplot (3,1,1)
    plot(time,M)
    title(sprintf('beta= %G deg; H1= %G T; H2= %G T; angle= %G deg; rhok1= %e', beta_d, H1,H2,angle,rhok1), 'fontsize',8)
%     title('M vs Time')
    legend('M1','M2')
    xlabel('time(sec)')
    ylabel('M')
    grid minor
    
%     subplot (5,1,1)
%     plot(time,emf)
%     ylim([emfll emful])
%     title(sprintf('beta= %G rad; H1= %G T; H2= %G T; rhok1= %e', beta, f1,f2,rhok1), 'fontsize',11)
% %     title('EMF vs Time')
%     xlabel('time(sec)')
%     ylabel('emf(volt)')
%     grid minor
    
    sp1 = subplot(3,1,2);
    hold on;
    plot(time, (xis(1,:))', 'g');
    hold on;
    plot(time, (xis(2,:))', 'b');
    hold on;
    ylabel('xis');
    xlabel('time');
    l2=legend('xi1','xi2');
    grid on;
    grid minor;
    
    subplot(3,1,3)
    plot(time,thetas)
    legend('theta1','theta2')
    xlabel('time(sec)')
    ylabel('theta(rad)')
    grid minor  
%     
%     subplot(5,1,5)
%     plot(time,alphas)
%     legend('alpha1','alpha2')
%     xlabel('time(sec)')
%     ylabel('alpha')
%     grid minor  
    
%     figure;
%     subplot (3,1,1)
%     plot(time,sig(1,:))
%     title(sprintf('beta= %G rad; H1= %G T; H2= %G T', beta, f1,f2), 'fontsize',12)
%     title('Stress vs Time')
%     hold on
%     plot(time_exp1,stress11/A)
%     legend('Sig11-model','Sig11-exp')
%     xlabel('time(sec)')
%     ylabel('Sig11 (Pa)')
%     grid minor
%     
%     subplot (3,1,2)
%     plot(time,emf)
%     ylim([emfll emful])
%     
%     xlabel('time(sec)')
%     ylabel('Model-emf(volt)')
%     grid minor
%     
%     subplot (3,1,3)
%     plot(emf_exp,time_exp)
% %     xlim([texpll texpul])
% %     legend('model','exp')
% title('Experimental emf vs Time')
% xlabel('time(sec)')
%     ylabel('exp-emf(volt)')
%     grid minor
%     
figure;
    subplot (2,1,1)
    plot(time,emf)
    ylim([emfll emful])
    title('EMF vs Time')
    xlabel('time(sec)')
    ylabel('emf(volt)')
    grid minor
    
    sp1 = subplot(2,1,2);
    hold on;
    plot(time, (xis(1,:))', 'g');
    hold on;
    plot(time, (xis(2,:))', 'b');
    hold on;
    ylabel('xis');
    xlabel('time');
    l2=legend('xi1','xi2');
    grid on;
    grid minor;
    
    figure;
    subplot (3,1,1)
    plot(time,sig)
    legend('$\sigma_{11}$','$\sigma_{22}$','interpreter','latex')
    title(sprintf('Beta= %G deg; H_1= %G T; H_2= %G T; Phi= %G deg; rhok1= %e', beta_d, H1,H2,angle,rhok1), 'fontsize',10)
    xlabel('time (sec)')
    ylabel('$\sigma_{ii}$ (Pa)','interpreter','latex')
    grid minor
    
    subplot (3,1,2)
    plot(time,emf)
    ylim([emfll emful])
%     title('EMF vs Time')
    xlabel('time (sec)')
    ylabel('Emf (V)')
    grid minor
    
    sp1 = subplot(3,1,3);
    hold on;
    plot(time, (xis(1,:))', 'g');
    hold on;
    plot(time, (xis(2,:))', 'b');
    hold on;
    ylabel('$\xi_i$','interpreter','latex');
    xlabel('time (sec)');
    l2=legend('$\xi_1$','$\xi_2$','interpreter','latex');
    grid on;
    grid minor;
end





