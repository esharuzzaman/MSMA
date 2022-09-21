function [ sig, Happ, time, xi00,dt] = loading( loadtype, res, lat_stress, field, axial,A,rf)
%case 1: 2D constant field loading
%case 2: 2D constant stress loading
%case 3: No field, stress in 3-direction
%case 4: Testing to see at what stress all variants could be put into xi3;
%        Alternating sig_33 while everything else is 0
%case 5: Constant stress in 1-direction & rotating field about 1-direction
%case 6: Constant stress in 1-direction & rotating stress about 1-direction
%case 7: 2D comparison
%case 8: Constant field in the 3-direction; Tri-axial actuator
%case 9: See fig 3.5 of (Doug's???) thesis
%        Predicted response of an MSMA to constant stress in the 3-direction, 
%        constant stress in the 1-direction, and a linearly varying field in the 2-direction.
%case 10: See fig 3.6 of (Doug's???) thesis
%         Predicted response of an MSMA to a constant applied stress in the 3- direction, 
%         constant applied magnetic field in the 2-direction, and a linearly varying stress in the 1-direction.
%case 11: See fig 3.7 of the thesis
%         Predicted response of an MSMA beginning at ?1 = ?2 = ?3 = 1/3, and continuing to evolve under 2D loading with a constant stress of 2 [MPa] in the 1-direction.
%case 12: See fig 3.8 of the thesis
%         Predicted response of an MSMA beginning at ?1 = ?2 = ?3 = 1/3, and continuing to evolve under 2D loading with a constant field in the 2-direction of .7 [T].
%case 13: H vs M - Hard axis
%case 14: H vs M - Easy axis
%case 15: .8T constant field - For comparison with 2D loading #5
%case 16: .6T constant field in the 3-direction, 
%         lateral load in the 2-direction, and 
%         linearly varying stress in the 1-direction
%case 17: .6T constant field in the 2-direction, 
%         lateral load in the 3-direction, and 
%         linearly varying stress in the 1-direction
%case 18: Calibration dir 2
%case 19: Calibration dir 3
%case 20: the axial and lateral stresses are constant and the field varies 
%         lateral stress is applied in 2 direction and the field is applied in 3 direction
%case 21: the axial and lateral stresses are constant and the field varies 
%         lateral stress is applied in 3 direction and the field is applied in 2 direction
dt = 1e-4;

%global mu0 D Msat

[mu0, Msat, rhok1, D] = Mat_consts();

switch loadtype
    case 202
        % experimental data 10HZ 2.877 deg H2=0.577
                % INPUTS:
       Happ_c=field;
        sig_max = -axial*10^6;
        %sig_min = 0;
        sig_lat_max = -lat_stress*10^6;
        
        res_pre = res/1000;
       %start of 1st cycle
 %1st piece
t1=[0:dt:0.047];
sigx1 =- 833.6*t1 - 0.72;
%2nd piece
t2=[.047:dt:.101];
sigx2 =730.2*t2 - 74.22;
%3rd piece
t3=[0.101:dt:0.147];
sigx3 =87.66 - 872.6*t3;
%end of 1st cycle

%start of 2nd cycle
 %1st piece
t4=[.147:dt:.2];
sigx4 =761.1*t4 - 152.5;
% % %2nd piece
% t5=[0.143:dt:0.189];
% sigx5 =676.9*t5 - 128.2;
% %3rd piece
% t6=[0.111:dt:0.147];
% sigx6 =121.2 - 1094.0*t6;
% %end of 2nd cycle
% 
% %start of 3rd cycle
%  %1st piece
% t7=[0.147:dt:0.161];
% sigx7 =1707.0*t7 - 290.5;
% %2nd piece
% t8=[0.161:dt:0.191];
% sigx8 =516.3*t8 - 98.74;

sig11 =[sigx1,sigx2,sigx3,sigx4]/A; %,sigx32,sigx33,sigx34...
   
 sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)));
%         sig_lat = [[sig_lat_pre] , [sig_lat_load]];
        sig_lat = [sig_lat_load];        
        sig = [sig11;sig_lat];
       %Next pre-load is apply field
        H_pre = [0:Happ_c/(res_pre):Happ_c];
        H_load = Happ_c*ones(1, (size(sig11,2)-size(H_pre,2)));

        H2 = [[H_pre] , [H_load]];
%         H2 =[H_load];
        %field at axial direction
H_pre1 = [[zeros(1, res_pre+1)], [0:Happ_c/(rf*res_pre):Happ_c/rf]];
        H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)-size(H_pre1,2)));

        H1 = [[H_pre1] , [H_load1]];

Happ = [ H1;H2];
time = [t1,t2,t3,t4]; %,t32,t33,t34,t35,t36,t37];
xi00 = [1;0]; 
    
    case 201
        % experimental data 10HZ 0 deg H2=0.575
                % INPUTS:
       Happ_c=field;
        sig_max = -axial*10^6;
        %sig_min = 0;
        sig_lat_max = -lat_stress*10^6;
        
        res_pre = res/1000;
       %start of 1st cycle
 %1st piece
t1=[0:dt:0.033];
sigx1 =- 941.9*t1 - 0.501;
%2nd piece
t2=[.033:dt:.089];
sigx2 =554.9*t2 - 49.89;
%3rd piece
t3=[0.089:dt:0.1];
sigx3 =21.91*t3 - 2.46;
%end of 1st cycle

%start of 2nd cycle
 %1st piece
t4=[.1:dt:.143];
sigx4 =72.2 - 724.7*t4;
% %2nd piece
t5=[0.143:dt:0.189];
sigx5 =676.9*t5 - 128.2;
% %3rd piece
% t6=[0.111:dt:0.147];
% sigx6 =121.2 - 1094.0*t6;
% %end of 2nd cycle
% 
% %start of 3rd cycle
%  %1st piece
% t7=[0.147:dt:0.161];
% sigx7 =1707.0*t7 - 290.5;
% %2nd piece
% t8=[0.161:dt:0.191];
% sigx8 =516.3*t8 - 98.74;

sig11 =[sigx1,sigx2,sigx3,sigx4,sigx5]/A; %,sigx32,sigx33,sigx34...
   
 sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)));
%         sig_lat = [[sig_lat_pre] , [sig_lat_load]];
        sig_lat = [sig_lat_load];        
        sig = [sig11;sig_lat];
       %Next pre-load is apply field
        H_pre = [0:Happ_c/(res_pre):Happ_c];
        H_load = Happ_c*ones(1, (size(sig11,2)-size(H_pre,2)));

        H2 = [[H_pre] , [H_load]];
%         H2 =[H_load];
        %field at axial direction
H_pre1 = [[zeros(1, res_pre+1)], [0:Happ_c/(rf*res_pre):Happ_c/rf]];
        H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)-size(H_pre1,2)));

        H1 = [[H_pre1] , [H_load1]];

Happ = [ H1;H2];
time = [t1,t2,t3,t4,t5]; %,t32,t33,t34,t35,t36,t37];
xi00 = [1;0]; 
    
        case 200
        % experimental data 10HZ 2.877 deg H2=0.577
                % INPUTS: 6.7 degree
       Happ_c=field;
        sig_max = -axial*10^6;
        %sig_min = 0;
        sig_lat_max = -lat_stress*10^6;
        
        res_pre = res/1000;
       %start of 1st cycle
 %1st piece
t1=[0:dt:0.014];
sigx1 =- 4.171*t1 - 0.0636;
%2nd piece
t2=[.014:dt:.047];
sigx2 =16.65 - 1198.0*t2;
%3rd piece
t3=[0.047:dt:0.062];
sigx3 =1666.0*t3 - 117.9;
%end of 1st cycle

%start of 2nd cycle
 %1st piece
t4=[.062:dt:.093];
sigx4 =473.3*t4 - 44.03;
% %2nd piece
t5=[0.093:dt:0.111];
sigx5 =0.7008 - 7.611*t5;
%3rd piece
t6=[0.111:dt:0.147];
sigx6 =121.2 - 1094.0*t6;
%end of 2nd cycle

%start of 3rd cycle
 %1st piece
t7=[0.147:dt:0.161];
sigx7 =1707.0*t7 - 290.5;
%2nd piece
t8=[0.161:dt:0.191];
sigx8 =516.3*t8 - 98.74;

sig11 =[sigx1,sigx2,sigx3,sigx4,sigx5,sigx6,sigx7,sigx8]/A; %,sigx32,sigx33,sigx34...
   
 sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)));
%         sig_lat = [[sig_lat_pre] , [sig_lat_load]];
        sig_lat = [sig_lat_load];        
        sig = [sig11;sig_lat];
       %Next pre-load is apply field
        H_pre = [0:Happ_c/(res_pre):Happ_c];
        H_load = Happ_c*ones(1, (size(sig11,2)-size(H_pre,2)));

        H2 = [[H_pre] , [H_load]];
%         H2 =[H_load];
        %field at axial direction
H_pre1 = [[zeros(1, res_pre+1)], [0:Happ_c/(rf*res_pre):Happ_c/rf]];
        H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)-size(H_pre1,2)));

        H1 = [[H_pre1] , [H_load1]];

Happ = [ H1;H2];
time = [t1,t2,t3,t4,t5,t6,t7,t8]; %,t32,t33,t34,t35,t36,t37];
xi00 = [1;0]; 
    
    case 54
        % experimental data 2HZ 0 deg H2=0.5 H1=0.0, dl
                % INPUTS:
       Happ_c=field;
        sig_max = -axial*10^6;
        %sig_min = 0;
        sig_lat_max = -lat_stress*10^6;
        
        res_pre = res/1000;
       %start of 1st cycle
 %1st piece
t1=[0:dt:0.061];
sigx1 =423.1*t1 - 44.1;
%2nd piece
t2=[0.061:dt:0.196];
sigx2 =132.7*t2 - 26.39;
%3rd piece
t3=[0.196:dt:0.347];
sigx3 =0.1808 - 2.821*t3; 
 %4th piece
t4=[0.347:dt:0.433];
sigx4 =84.74 - 246.5*t4;
%5th piece
t5=[0.433:dt:0.5];
sigx5=126.4 - 342.6*t5;
%6th piece
t6=[0.5:dt:.552];
sigx6 =457.1*t6 - 273.5;
 %7th piece
t7=[.552:dt:.704];
sigx7 =137.8*t7 - 97.26;
%8th piece
t8=[.707:dt:.841];
sigx8 =1.624 - 2.642*t8;
%9th piece
t9=[.841:dt:.938];
sigx9 =192.7 - 229.8*t9;
 %10th piece
t10=[0.938:dt:1];
sigx10 =319.0 - 364.5*t10;

sig11 =[sigx1,sigx2,sigx3,sigx4,sigx5,sigx6,sigx7,sigx8,sigx9,sigx10]/A;
 
 
        %first pre-load is apply lateral stress
        %sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
        % once lateral stress is applied it stays constant
%         sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
 sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)));
%         sig_lat = [[sig_lat_pre] , [sig_lat_load]];
        sig_lat = [sig_lat_load];        
        sig = [sig11;sig_lat];
        %Next pre-load is apply field
        H_pre = [0:Happ_c/(res_pre):Happ_c];
        H_load = Happ_c*ones(1, (size(sig11,2)-size(H_pre,2)));
% H_load = Happ_c*ones(1, (size(sig11,2)));
        H2 = [[H_pre] , [H_load]];
%         H2 =[H_load];
        %field at axial direction
H_pre1 = [[zeros(1, res_pre+1)], [0:Happ_c/(rf*res_pre):Happ_c/rf]];
        H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)-size(H_pre1,2)));
% H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)));
        H1 = [[H_pre1] , [H_load1]];
%         H1 = [H_load1];
%         Happ = [ zeros(1,size(sig11,2));
%                  H];
Happ = [ H1;H2];
time = [t1,t2,t3,t4,t5,t6,t7,t8,t9,t10];
xi00 = [1;0]; 
    
    case 50
        % experimental data 10HZ 5.85 deg H2=0.537 H1=0.055
                % INPUTS:
       Happ_c=field;
        sig_max = -axial*10^6;
        %sig_min = 0;
        sig_lat_max = -lat_stress*10^6;
        
        res_pre = res/1000;
       %start of 1st cycle
 %1st piece
t1=[0:dt:.011];
sigx1 =- 1301.0*t1 - 30.06;
%2nd piece
t2=[0.011:dt:0.02];
sigx2 =457.7*t2 - 49.41;
%3rd piece
t3=[0.02:dt:0.033];
sigx3 =2571.0*t3 - 91.67; 
 %4th piece
t4=[0.033:dt:0.059];
sigx4 =205.3*t4 - 13.61;
%5th piece
t5=[0.059:dt:0.068];
sigx5=5.773 - 123.2*t5;
%6th piece
t6=[0.068:dt:0.112];
sigx6 =60.7 - 930.9*t6;
 %7th piece
t7=[0.112:dt:0.132];
sigx7 =1820.0*t7 - 247.4;
%8th piece
t8=[0.132:dt:0.157];
sigx8 =223.1*t8 - 36.62;
%9th piece
t9=[0.157:dt:0.168];
sigx9 =14.34 - 101.5*t9;
 %10th piece
t10=[0.168:dt:0.211];
sigx10 =155.7 - 943.0*t10;
%11th piece
t11=[0.211:dt:0.218];
sigx11 =343.3*t11 - 115.7;
%12th piece
t12=[0.218:dt:0.234];
sigx12 =2167.0*t12 - 513.3;
 %13th piece
t13=[0.234:dt:0.259];
sigx13 =188.1*t13 - 50.21;
%14th piece
t14=[0.259:dt:0.268];
sigx14 =31.78 - 128.4*t14;
%15th piece
t15=[0.268:dt:0.312];
sigx15 =244.9 - 923.5*t15;
%16th piece
t16=[0.312:dt:0.316];
sigx16 =344.5*t16 - 150.8;
%17th piece
t17=[0.316:dt:0.333];
sigx17 =2076.0*t17 - 697.8;
%18th piece
t18=[0.333:dt:0.36];
sigx18 =190.5*t18 - 70.04;
%19th piece
t19=[0.36:dt:0.368];
sigx19 =37.5 - 108.2*t19;
%20th piece
t20=[0.368:dt:0.393];
sigx20 =309.7 - 848.0*t20;

sig11 =[sigx1,sigx2,sigx3,sigx4,sigx5,sigx6,sigx7,sigx8,sigx9,sigx10,sigx11,sigx12...
   ,sigx13,sigx14,sigx15,sigx16,sigx17,sigx18,sigx19,sigx20]/A;
 
 
        %first pre-load is apply lateral stress
        %sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
        % once lateral stress is applied it stays constant
%         sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
 sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)));
%         sig_lat = [[sig_lat_pre] , [sig_lat_load]];
        sig_lat = [sig_lat_load];        
        sig = [sig11;sig_lat];
       %Next pre-load is apply field
        H_pre = [0:Happ_c/(res_pre):Happ_c];
        H_load = Happ_c*ones(1, (size(sig11,2)-size(H_pre,2)));
% H_load = Happ_c*ones(1, (size(sig11,2)));
        H2 = [[H_pre] , [H_load]];
%         H2 =[H_load];
        %field at axial direction
H_pre1 = [[zeros(1, res_pre+1)], [0:Happ_c/(rf*res_pre):Happ_c/rf]];
        H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)-size(H_pre1,2)));
% H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)));
        H1 = [[H_pre1] , [H_load1]];
%         H1 = [H_load1];
%         Happ = [ zeros(1,size(sig11,2));
%                  H];
Happ = [ H1;H2];
time = [t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20];
xi00 = [1;0]; 
    
    case 52
        % experimental data 2HZ 5.85 deg H2=0.537 H1=0.055
                % INPUTS:
       Happ_c=field;
        sig_max = -axial*10^6;
        %sig_min = 0;
        sig_lat_max = -lat_stress*10^6;
        
        res_pre = res/1000;
       %start of 1st cycle
 %1st piece
t1=[0:dt:.127];
sigx1 =314.0*t1 - 45.62;
%2nd piece
t2=[0.127:dt:0.307];
sigx2 =27.42*t2 - 9.232;
%3rd piece
t3=[0.307:dt:0.5];
sigx3 =70.0 - 230.7*t3; 
 %4th piece
t4=[0.5:dt:0.639];
sigx4 =292.4*t4 - 191.5;
%5th piece
t5=[0.639:dt:0.797];
sigx5=26.16*t5 - 21.41;
%6th piece
t6=[0.797:dt:1.002];
sigx6 =170.8 - 215.1*t6;
 %7th piece
t7=[1.002:dt:1.134];
sigx7 =298.2*t7 - 343.5;
%8th piece
t8=[1.134:dt:1.301];
sigx8 =27.78*t8 - 36.78;
%9th piece
t9=[1.301:dt:1.504];
sigx9 =282.8 - 217.9*t9;

sig11 =[sigx1,sigx2,sigx3,sigx4,sigx5,sigx6,sigx7,sigx8,sigx9]/A;
 
 
        %first pre-load is apply lateral stress
        %sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
        % once lateral stress is applied it stays constant
%         sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
 sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)));
%         sig_lat = [[sig_lat_pre] , [sig_lat_load]];
        sig_lat = [sig_lat_load];        
        sig = [sig11;sig_lat];
      %Next pre-load is apply field
        H_pre = [0:Happ_c/(res_pre):Happ_c];
        H_load = Happ_c*ones(1, (size(sig11,2)-size(H_pre,2)));
% H_load = Happ_c*ones(1, (size(sig11,2)));
        H2 = [[H_pre] , [H_load]];
%         H2 =[H_load];
        %field at axial direction
H_pre1 = [[zeros(1, res_pre+1)], [0:Happ_c/(rf*res_pre):Happ_c/rf]];
        H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)-size(H_pre1,2)));
% H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)));
        H1 = [[H_pre1] , [H_load1]];
%         H1 = [H_load1];
%         Happ = [ zeros(1,size(sig11,2));
%                  H];
Happ = [ H1;H2];
time = [t1,t2,t3,t4,t5,t6,t7,t8,t9];
xi00 = [1;0]; 
   case 51
        % experimental data 2HZ 0 deg H2=0.54
                % INPUTS:
       Happ_c=field;
        sig_max = -axial*10^6;
        %sig_min = 0;
        sig_lat_max = -lat_stress*10^6;
        
        res_pre = res/1000;
       %start of 1st cycle
 %1st piece
t1=[0:dt:.192];
sigx1 =225.7*t1 - 44.22;
%2nd piece
t2=[.192:dt:.335];
sigx2 =- 2.724*t2 - 0.3632;
%3rd piece
t3=[0.335:dt:.5];
sigx3 =84.94 - 257.4*t3;
%end of 1st cycle

%start of 2nd cycle
 %1st piece
t4=[.5:dt:.692];
sigx4 =223.3*t4 - 155.4;
% %2nd piece
t5=[0.692:dt:0.835];
sigx5 =1.051 - 2.786*t5;
%3rd piece
t6=[0.835:dt:1];
sigx6 =215.9 - 260.1*t6;
%end of 2nd cycle

%start of 3rd cycle
 %1st piece
t29=[1:dt:1.192];
sigx29 =229.7*t29 - 273.9;
%2nd piece
t30=[1.192:dt:1.343];
sigx30 =9.413 - 7.972*t30;
%3rd piece
t31=[1.343:dt:1.5];
sigx31 =367.1 - 274.3*t31;
%end of 3rd cycle

% %start of 4th cycle
%  %1st piece
% t32=[1.5:dt:1.692];
% sigx32 =226.4*t32 - 384.0;
% %2nd piece
% t33=[1.692:dt:1.835];
% sigx33 =3.562 - 2.63*t33;
% %3rd piece
% t34=[1.835:dt:2];
% sigx34 =477.1 - 260.7*t34;
% %end of 4th cycle
% 
% %start of 5th cycle
%  %1st piece
% t35=[2:dt:2.185];
% sigx35 =234.1*t35 - 512.5;
% %2nd piece
% t36=[2.185:dt:2.335];
% sigx36 =3.507 - 2.045*t36;
% %3rd piece
% t37=[2.335:dt:2.5];
% sigx37 =607.3 - 260.6*t37;
% %end of 5th cycle

sig11 =[sigx1,sigx2,sigx3,sigx4,sigx5,sigx6,sigx29,sigx30,sigx31]/A; %,sigx32,sigx33,sigx34...
   %,sigx35,sigx36,sigx37]/A;
 
 
        %first pre-load is apply lateral stress
        %sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
        % once lateral stress is applied it stays constant
%         sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
 sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)));
%         sig_lat = [[sig_lat_pre] , [sig_lat_load]];
        sig_lat = [sig_lat_load];        
        sig = [sig11;sig_lat];
       %Next pre-load is apply field
        H_pre = [0:Happ_c/(res_pre):Happ_c];
        H_load = Happ_c*ones(1, (size(sig11,2)-size(H_pre,2)));
% H_load = Happ_c*ones(1, (size(sig11,2)));
        H2 = [[H_pre] , [H_load]];
%         H2 =[H_load];
        %field at axial direction
H_pre1 = [[zeros(1, res_pre+1)], [0:Happ_c/(rf*res_pre):Happ_c/rf]];
        H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)-size(H_pre1,2)));
% H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)));
        H1 = [[H_pre1] , [H_load1]];
%         H1 = [H_load1];
%         Happ = [ zeros(1,size(sig11,2));
%                  H];
Happ = [ H1;H2];
time = [t1,t2,t3,t4,t5,t6,t29,t30,t31]; %,t32,t33,t34,t35,t36,t37];
xi00 = [1;0]; 
    
    
case 95
        % experimental data
                % INPUTS:
       Happ_c=field;
        sig_max = -axial*10^6;
        %sig_min = 0;
        sig_lat_max = -lat_stress*10^6;
        
        res_pre = res/1000;
       %1st piece
t1=[0:dt:0.01];
sigx01 =(69037*t1)/10 - 63;
%2nd piece
t2=[.01:dt:.04];
sigx02 =1767/250 - (1031*t2)/10;
%3rd piece
t3=[0.04:dt:0.056];
sigx03 =7315147616616973/35184372088832 - (40993*t3)/8;
%4th piece
t4=[0.056:dt:0.065];
sigx04 =(784*t4)/9 - 369085040555517/4398046511104;
%5th piece
t5=[0.065:dt:0.077];
sigx05 =(42115*t5)/6 - 2350772719363381/4398046511104;
%6th piece
t6=[0.077:dt:0.108];
sigx06 =1284738052168039/70368744177664 - (4946*t6)/31;
%7th piece
t7=[0.108:dt:0.125];
sigx07 =2293875665948361/4398046511104 - (81937*t7)/17;
%8th piece
t8=[0.125:dt:0.131];
sigx08=(2165*t8)/6 - 8867529025672365/70368744177664;
%9th piece
t9=[0.131:dt:0.143];
sigx09=(42269*t9)/6 - 4405167681413235/4398046511104;
%10th piece
t10=[0.143:dt:0.175];
sigx10=4757855094223929/140737488355328 - (6269*t10)/32;
%11th piece
t11=[0.175:dt:0.187];
sigx11=2503963826170711/2199023255552 - (78113*t11)/12;
%12th piece
t12=[0.187:dt:0.192];
sigx12=- 382*t12 - 503558733335377/70368744177664;
%13th piece
t13=[0.192:dt:0.198];
sigx13=364*t13 - 37597/250;
%14th piece
t14=[0.198:dt:0.2];
sigx14=7658*t14 - 7973/5;        
sig11 =[sigx01,sigx02,sigx03,sigx04,sigx05,sigx06,sigx07,sigx08,sigx09,sigx10,sigx11,sigx12,sigx13,sigx14]/A;
 
        %first pre-load is apply lateral stress
        %sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
        % once lateral stress is applied it stays constant
%         sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
 sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)));
%         sig_lat = [[sig_lat_pre] , [sig_lat_load]];
        sig_lat = [sig_lat_load];        
        sig = [sig11;sig_lat];
       %Next pre-load is apply field
        H_pre = [0:Happ_c/(res_pre):Happ_c];
        H_load = Happ_c*ones(1, (size(sig11,2)-size(H_pre,2)));
% H_load = Happ_c*ones(1, (size(sig11,2)));
        H2 = [[H_pre] , [H_load]];
%         H2 =[H_load];
        %field at axial direction
H_pre1 = [[zeros(1, res_pre+1)], [0:Happ_c/(rf*res_pre):Happ_c/rf]];
        H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)-size(H_pre1,2)));
% H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)));
        H1 = [[H_pre1] , [H_load1]];
%         H1 = [H_load1];
%         Happ = [ zeros(1,size(sig11,2));
%                  H];
Happ = [ H1;H2];
time = [t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14];
xi00 = [1;0];

    case 89
        %zero stress , varrying field
        % 2D constant stress loading
        % INPUTS:
        Happ_c = field;
        sig_max = -axial*10^6;
        sig_min = 0;
        sig_lat_max = -lat_stress*10^6;

        res_pre = res/10;
        
        % applied field
        H_pre = [zeros(1, res_pre+1), zeros(1, res_pre+1)];
        H_load = [[0:(Happ_c)/(res):Happ_c,Happ_c:-Happ_c/(res):0],[0:(Happ_c)/(res):Happ_c,Happ_c:-Happ_c/(res):0]];
        Happ = [H_pre , H_load];
        Happ=[zeros(1,length(Happ));Happ];
        
%         %first pre-load is apply lateral stress
%         sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
%         %next pre-load is apply axial stress
%         sig11 = [zeros(1, res_pre+1),0:(sig_max)/(res_pre):sig_max, sig_max.*ones(1,length(H_load))];
%          
%         % once lateral stress is applied it stays constant
%         sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
%         sig_lat = [[sig_lat_pre] , [sig_lat_load]];
                
%         sig = [sig11;sig_lat];
        sig=[zeros(1,length(Happ));zeros(1,length(Happ))];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [1;0];
    
    case 99
        % 2D constant field loading
                % INPUTS:
        Happ_c = field;
        sig_max = -axial*10^6;
        sig_min = 0;
        sig_lat_max = -lat_stress*10^6;
        
        res_pre = res/10;
        %res_pre2= res/100;
        
        sig11 = [[zeros(1, res_pre+1)],[zeros(1, res_pre+1)],[0:(sig_max)/(res):sig_max,sig_max:-(sig_max-sig_min)/(res):sig_min],[sig_min:(sig_max-sig_min)/(res):sig_max,sig_max:-(sig_max-sig_min)/(res):sig_min]];
        
        %first pre-load is apply lateral stress
        sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
        % once lateral stress is applied it stays constant
        sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
        sig_lat = [[sig_lat_pre] , [sig_lat_load]];
              
        sig = [sig11;sig_lat];
%         disp('sig: ');
%         disp(sig);
        %Next pre-load is apply field
        H_pre = [[zeros(1, res_pre+1)], [0:Happ_c/(res_pre):Happ_c]];
        H_load = Happ_c*ones(1, (size(sig11,2)-size(H_pre,2)));
        H2 = [[H_pre] , [H_load]];
%         

%field at axial direction
H_pre1 = [[zeros(1, res_pre+1)], [0:Happ_c/(rf*res_pre):Happ_c/rf]];
        H_load1 = (Happ_c)/rf*ones(1, (size(sig11,2)-size(H_pre1,2)));
        H1 = [[H_pre1] , [H_load1]];

%         Happ = [ zeros(1,size(sig11,2));
%                  H];
Happ = [ H1;H2];

%Happ = [H;H];


        time = (1:size(Happ,2))*dt;
        xi00 = [1;0];
%         disp('time: ');
%         disp(time);

case 90
        % 2D variable stress & field
                % INPUTS:
        Happ_c = field;
        sig_max = -axial*10^6;
        sig_min = 0;
        sig_lat_max = -lat_stress*10^6;
        
        
        res_pre = res/10;
        
        sig11 = [[zeros(1, res_pre+1)],[zeros(1, res_pre+1)],[0:(sig_max)/(res):sig_max,sig_max:-(sig_max-sig_min)/(res):sig_min],[sig_min:(sig_max-sig_min)/(res):sig_max,sig_max:-(sig_max-sig_min)/(res):sig_min]];
        
        %first pre-load is apply lateral stress
        sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
        % once lateral stress is applied it stays constant
        sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
        sig_lat = [[sig_lat_pre] , [sig_lat_load]];
            
        sig = [sig11;sig_lat];
        
        % applied field
        H_pre = [zeros(1, res_pre+1), zeros(1, res_pre+1)];
        H_load = [[0:(Happ_c)/(res):Happ_c,Happ_c:-Happ_c/(res):0],[0:(Happ_c)/(res):Happ_c,Happ_c:-Happ_c/(res):0]];
        Happ = [H_pre , H_load];
        Happ=[zeros(1,length(Happ));Happ];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [1;0];
        
    case 1
        % 2D constant field loading
                % INPUTS:
        Happ_c = field;
        sig_max = -axial*10^6;
        sig_min = 0;
        sig_lat_max = -lat_stress*10^6;
        
        
        res_pre = res/10;
        
        sig11 = [[zeros(1, res_pre+1)],[zeros(1, res_pre+1)],[0:(sig_max)/(res):sig_max,sig_max:-(sig_max-sig_min)/(res):sig_min],[sig_min:(sig_max-sig_min)/(res):sig_max,sig_max:-(sig_max-sig_min)/(res):sig_min]];
        
        %first pre-load is apply lateral stress
        sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
        % once lateral stress is applied it stays constant
        sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
        sig_lat = [[sig_lat_pre] , [sig_lat_load]];
        
                
        sig = [sig11;sig_lat];
        
        %Next pre-load is apply field
        H_pre = [[zeros(1, res_pre+1)], [0:Happ_c/(res_pre):Happ_c]];
        H_load = Happ_c*ones(1, (size(sig11,2)-size(H_pre,2)));
        H = [[H_pre] , [H_load]];
        
        Happ = [ zeros(1,size(sig11,2));
                 H];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [1;0];
        
    case 2
        % 2D constant stress loading
        % INPUTS:
        Happ_c = field;
        sig_max = -axial*10^6;
        sig_min = 0;
        sig_lat_max = -lat_stress*10^6;

        res_pre = res/10;
        
        % applied field
        H_pre = [zeros(1, res_pre+1), zeros(1, res_pre+1)];
        H_load = [[0:(Happ_c)/(res):Happ_c,Happ_c:-Happ_c/(res):0],[0:(Happ_c)/(res):Happ_c,Happ_c:-Happ_c/(res):0]];
        Happ = [H_pre , H_load];
        Happ=[zeros(1,length(Happ));Happ];
        
        
        %first pre-load is apply lateral stress
        sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
        %next pre-load is apply axial stress
        sig11 = [zeros(1, res_pre+1),0:(sig_max)/(res_pre):sig_max, sig_max.*ones(1,length(H_load))];
         
        % once lateral stress is applied it stays constant
        sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
        sig_lat = [[sig_lat_pre] , [sig_lat_load]];
                
        sig = [sig11;sig_lat];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [0;1];
        
         case 222
        % 2D constant stress loading H1 varries
        % INPUTS:
        Happ_c = field;
        sig_max = -axial*10^6;
        sig_min = 0;
        sig_lat_max = -lat_stress*10^6;

        res_pre = res/10;
        
        % applied field
        H_pre = [zeros(1, res_pre+1), zeros(1, res_pre+1)];
        H_load = [[0:(Happ_c)/(res):Happ_c,Happ_c:-Happ_c/(res):0],[0:(Happ_c)/(res):Happ_c,Happ_c:-Happ_c/(res):0]];
        Happ = [H_pre , H_load];
        Happ=[Happ;zeros(1,length(Happ))];
        
        
        %first pre-load is apply lateral stress
        sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
        %next pre-load is apply axial stress
        sig11 = [zeros(1, res_pre+1),0:(sig_max)/(res_pre):sig_max, sig_max.*ones(1,length(H_load))];
         
        % once lateral stress is applied it stays constant
        sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
        sig_lat = [[sig_lat_pre] , [sig_lat_load]];
                
        sig = [sig11;sig_lat];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [0;1];
     case 3
        % No field, stress in 3-direction
        sig33_max = -8e6;
        Happ1_c = 0/mu0;
        Happ2_c = 0/mu0;
        
        sig33 = [0:sig33_max/res:sig33_max,sig33_max:-sig33_max/res:0];
        
        sig = [0*sig33; 0*sig33; sig33];
        Happ = [ones(1,2*res+2)*Happ1_c;...
                ones(1,2*res+2)*Happ2_c;...
                zeros(1,2*res+2)];
        time = dt*(1:2*res+2);
        
        xi00 = [.75;.25;0];
    case 4
        % Testing to see at what stress all variants could be put into xi3
        % Alternating sig_33 while everything else is 0
        sig33_max = -50e6;
        Happ1_c = 0/mu0;
        Happ2_c = 0/mu0;
        
        sig33 = [0:sig33_max/res:sig33_max,sig33_max:-sig33_max/res:0];
        
        sig = [sig33; sig33; sig33];
        Happ = [ones(1,2*res+2)*Happ1_c;...
                ones(1,2*res+2)*Happ2_c;...
                zeros(1,2*res+2)];
        time = dt*(1:2*res+2);
        
        xi00 = [.7;.3;0];
    case 5
        % Constant stress in 1-direction
        % rotating field about 1-direction
        sig11_c = -2e6;
        
        Happ_c = 1/mu0;
        gamma = 0:2*pi/res:4*pi;
        
        sig11 = sig11_c * ones(1,2*res+1);
        sig22 = zeros(1,2*res+1);
        sig33 = zeros(1,2*res+1);
        
        Happ = [zeros(1,2*res+1);cos(gamma);sin(gamma)]*Happ_c;
        
        sig = [sig11; sig22; sig33];
        
        time = dt*(1:size(sig,2));

        xi00 = [0;1;0];
    case 6
        % Constant stress in 1-direction
        % rotating stress about 1-direction
        sig_c = -4.5e6;
        Happ_c = .4/mu0;
        gamma = 0:2*pi/res:4*pi;
        
        Happ = [ones(1,2*res+1);
                zeros(2,2*res+1)]*Happ_c;
        
        sig = [zeros(1,2*res+1);.5+.5*cos(gamma);.5-.5*cos(gamma)]*sig_c;
        
        time = 1:size(sig,2);
        xi00 = [0;1;0];
    case 7
        % 2D comparison
        sig_c=-1e6;
        %sig_c = -2.0418606314948e6;
        Happ_max = 0.8/mu0;
        %Happ_max = 1/mu0;
        phi = 0*pi/180;
        
        type=2;
        
        preload_sig = [0:sig_c*1/res:sig_c,sig_c];
        unload_sig  = [sig_c,sig_c:-sig_c*1/res:0];
        Happ2 = [0:(Happ_max)/(res):Happ_max,Happ_max:-(Happ_max)/(res):0] * cos(phi);
        Happ1 = [0:(Happ_max)/(res):Happ_max,Happ_max:-(Happ_max)/(res):0] * sin(phi);
        sig   = [  preload_sig,sig_c*ones(1,size(Happ2,2)*2),unload_sig;
                  [zeros(2,size(preload_sig,2)),zeros(2,size(Happ2,2)*2),zeros(2,size(unload_sig,2))]];
        Happ = [zeros(3,size(preload_sig,2)),[Happ1,Happ1;Happ2,Happ2;zeros(1,size(Happ1,2)*2)],zeros(3,size(unload_sig,2))];
        time = (1:size(Happ,2))*dt;
        
        xi00 = [0;1;0];
    case 8
        % Constant field in the 3-direction
        % Tri-axial actuator
        sig_c = -.5e6;
        Happ_max = 1/mu0;
        
        Happ_tent = [0:Happ_max/res:Happ_max,Happ_max:-Happ_max/res:0];
        
        Happ = [Happ_tent, zeros(1,size(Happ_tent,2)), zeros(1,size(Happ_tent,2));
                zeros(1,size(Happ_tent,2)), Happ_tent, zeros(1,size(Happ_tent,2));
                zeros(1,size(Happ_tent,2)), zeros(1,size(Happ_tent,2)), Happ_tent;];
        
        sig = [zeros(2,size(Happ,2));
                ones(1,size(Happ,2))*sig_c;];
        
        time = (1:size(Happ,2))*dt;
        
        xi00 = [1/3;1/3;1/3];
    case 9
        % See fig 3.5 of the thesis
        sig3_max = -2e6;
        sig3_pre = [0:sig3_max/res:sig3_max];
        
        sig1_max = -2e6;
        sig1_pre = [0:sig1_max/res:sig1_max];
        
        Happ2_max = 1/mu0;
        Happ2 = [0:Happ2_max/res:Happ2_max, Happ2_max:-Happ2_max/res:0];
        
        sig = [ zeros(1,size(sig3_pre,2)), sig1_pre, ones(1,size(Happ2,2)*2)*sig1_max, fliplr(sig1_pre), zeros(1,size(sig3_pre,2))
                                            zeros(1, size(sig3_pre,2)*2+size(sig1_pre,2)*2+size(Happ2,2)*2);
                                 sig3_pre,          ones(1,size(Happ2,2)*2+size(sig1_pre,2)*2)*sig3_max,          fliplr(sig3_pre)];
        Happ = [zeros(1,size(sig,2));
                zeros(1,size(sig1_pre,2)+size(sig3_pre,2)),Happ2,Happ2,zeros(1,size(sig1_pre,2)+size(sig3_pre,2));
                zeros(1,size(sig,2));];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [0;1;0];
    case 10
        % See fig 3.6 of the thesis
        %    Predicted response of an MSMA to a constant applied stress in the 3- direction, 
        %    constant applied magnetic field in the 2-direction, and a linearly varying stress in the 1-direction.

        % INPUTS:
        Happ2_max = field;
        sig1_max = -axial*10^6;
        sig3_max = -lat_stress*10^6;
        
        
        %sig3_max = -2.0e6;
        sig3_pre = [0:sig3_max/res:sig3_max];
        
        %sig1_max = -4.752e6;
        sig1_loading = [0:sig1_max/res:sig1_max];
        
        
        Happ2_pre = 0:Happ2_max/res:Happ2_max;
        
        sig = [   zeros(1,size(sig3_pre,2)+size(Happ2_pre,2)),     sig1_loading, fliplr(sig1_loading),sig1_loading, fliplr(sig1_loading), zeros(1,size(sig3_pre,2)+size(Happ2_pre,2)); 
                                    zeros(1,(size(sig3_pre,2)+size(Happ2_pre,2)+size(sig1_loading,2)*2)*2)
                        sig3_pre,          ones(1,size(Happ2_pre,2)*2+size(sig1_loading,2)*4)*sig3_max,          fliplr(sig3_pre)];
                    
                    
                    
        Happ = [zeros(1,size(sig,2));
                zeros(1,size(sig3_pre,2)),Happ2_pre,ones(1,size(sig1_loading,2)*4)*Happ2_max,fliplr(Happ2_pre),zeros(1,size(sig3_pre,2));
                zeros(1,size(sig,2));];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [0;1;0];
    case 11
        % See fig 3.7 of the thesis
        sig3_max = -1e6;
        sig3_pre = fliplr([0:sig3_max/res:sig3_max]);
        
        sig1_max = -2e6;
        sig1_pre = [0:sig1_max/res:sig1_max];
        
        Happ2_max = 1/mu0;
        Happ2 = [0:Happ2_max/res:Happ2_max, Happ2_max:-Happ2_max/res:0];
        
        sig = [ zeros(1,size(sig3_pre,2)), sig1_pre, ones(1,size(Happ2,2)*2)*sig1_max, fliplr(sig1_pre), zeros(1,size(sig3_pre,2))
                                            zeros(1, size(sig3_pre,2)*2+size(sig1_pre,2)*2+size(Happ2,2)*2);
                                 sig3_pre,                  zeros(1,size(Happ2,2)*2+size(sig1_pre,2)*2),          fliplr(sig3_pre)];
        Happ = [zeros(1,size(sig,2));
                zeros(1,size(sig1_pre,2)+size(sig3_pre,2)),Happ2,Happ2,zeros(1,size(sig1_pre,2)+size(sig3_pre,2));
                zeros(1,size(sig,2));];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [1/3;1/3;1/3];
    case 12
        % See fig 3.8 of the thesis
        sig3_max = -1e6;
        sig3_pre = fliplr([0:sig3_max/res:sig3_max]);
        
        sig1_max = -8e6;
        sig1_loading = [0:sig1_max/res:sig1_max];
        
        Happ2_max = .7/mu0;
        Happ2_pre = 0:Happ2_max/res:Happ2_max;
        
        sig = [   zeros(1,size(sig3_pre,2)+size(Happ2_pre,2)),     sig1_loading, fliplr(sig1_loading),sig1_loading, fliplr(sig1_loading), zeros(1,size(sig3_pre,2)+size(Happ2_pre,2)); 
                                    zeros(1,(size(sig3_pre,2)+size(Happ2_pre,2)+size(sig1_loading,2)*2)*2)
                        sig3_pre,          zeros(1,size(Happ2_pre,2)*2+size(sig1_loading,2)*4),          fliplr(sig3_pre)];
        Happ = [zeros(1,size(sig,2));
                zeros(1,size(sig3_pre,2)),Happ2_pre,ones(1,size(sig1_loading,2)*4)*Happ2_max,fliplr(Happ2_pre),zeros(1,size(sig3_pre,2));
                zeros(1,size(sig,2));];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [1/3;1/3;1/3];
    case 13
        % Used to generate H vs M curves
        % Hard axis
        
        % INPUTS:
        Happ_c = field;
        sig_max = -axial*10^6;
        
        res_pre = res/10;
        
        sig11 = ones(1,res).*sig_max;
        sig11(1:res_pre)=linspace(0,sig_max,res_pre);
        sig = [sig11;zeros(1,length(sig11))];
        %Next pre-load is apply field
        Happ = [zeros(1,res);
                zeros(1,res_pre), linspace(0,Happ_c,res-res_pre)];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [1;0];
    case 14
        % Used to generate H vs M curves
        % Easy axis
        
        % INPUTS:
        % INPUTS:
        Happ_c = field;
        sig_max = -axial*10^6;
        
        res_pre = res/10;
        
        sig11 = ones(1,res).*sig_max;
        sig11(1:res_pre)=linspace(0,sig_max,res_pre);
        sig = [sig11;zeros(1,length(sig11))];
        %Next pre-load is apply field
        Happ = [zeros(1,res);
                zeros(1,res_pre), linspace(0,Happ_c,res-res_pre)];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [0;1];
    case 15
        % .8T constant field
        % For comparison with 2D loading #5
        Happ_c = 0.787534978165939/mu0;
        sig_max = -8.6e6;
        phi = 0*pi/180;
        %mytitle = sprintf('Constant B^{app} = %.2f T',Happ_c*mu0);
        
        sig11 = [[0:sig_max/(res):sig_max,sig_max:-sig_max/(res):0],[0:sig_max/(res):sig_max,sig_max:-sig_max/(res):0]];
        sig = [sig11;zeros(2,size(sig11,2))];
        Happ = [ Happ_c*ones(1,size(sig11,2))*sin(phi);
                 Happ_c*ones(1,size(sig11,2))*cos(phi);
                 zeros(1,size(sig11,2))];
        time = (1:size(Happ,2))*dt;
        xi00 = [0;1;0];
    case 16
        % constant field in 2-dir, lateral stress in 3-dir, linearly
        % varying stres (up do a defined max) in the 1-dir
        
        % INPUTS:
        Happ_c = field;
        sig_max = -5.33*10^6+axial*10^6;
        sig_min = -axial*10^6;
        sig_lat_max = -lat_stress*10^6;
        
        
        res_pre = res/10;
        
        sig11 = [[zeros(1, res_pre+1)],[zeros(1, res_pre+1)],[0:(sig_max)/(res):sig_max,sig_max:-(sig_max-sig_min)/(res):sig_min],[sig_min:(sig_max-sig_min)/(res):sig_max,sig_max:-(sig_max-sig_min)/(res):sig_min]];
        
        %first pre-load is apply lateral stress
        sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
        % once lateral stress is applied it stays constant
        sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
        sig_lat = [[sig_lat_pre] , [sig_lat_load]];
        
                
        sig = [sig11;zeros(1,size(sig11,2));sig_lat];
        
        %Next pre-load is apply field
        H_pre = [[zeros(1, res_pre+1)], [0:Happ_c/(res_pre):Happ_c]];
        H_load = Happ_c*ones(1, (size(sig11,2)-size(H_pre,2)));
        H = [[H_pre] , [H_load]];
        
        Happ = [ Happ_c*ones(1,size(sig11,2))*0;
                 H;
                 Happ_c*ones(1,size(sig11,2))*0];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [1;0;0];
    
    case 17
        % constant field in 3-dir, lateral stress in 2-dir, linearly
        % varying stres (up do a defined max) in the 1-dir
        
                % INPUTS:
        Happ_c = field;
        sig_max = -5.33*10^6+axial*10^6;
        sig_min = -axial*10^6;
        sig_lat_max = -lat_stress*10^6;
        
        
        res_pre = res/10;
        
        sig11 = [[zeros(1, res_pre+1)],[zeros(1, res_pre+1)],[0:(sig_max)/(res):sig_max,sig_max:-(sig_max-sig_min)/(res):sig_min],[sig_min:(sig_max-sig_min)/(res):sig_max,sig_max:-(sig_max-sig_min)/(res):sig_min]];
        
        %first pre-load is apply lateral stress
        sig_lat_pre = 0:sig_lat_max/(res_pre):sig_lat_max;
        % once lateral stress is applied it stays constant
        sig_lat_load = sig_lat_max*ones(1, (size(sig11,2)-size(sig_lat_pre,2)));
        sig_lat = [[sig_lat_pre] , [sig_lat_load]];
        
                
        sig = [sig11;sig_lat; zeros(1,size(sig11,2))];
        
        %Next pre-load is apply field
        H_pre = [[zeros(1, res_pre+1)], [0:Happ_c/(res_pre):Happ_c]];
        H_load = Happ_c*ones(1, (size(sig11,2)-size(H_pre,2)));
        H = [[H_pre] , [H_load]];
        
        Happ = [ Happ_c*ones(1,size(sig11,2))*0;
                 Happ_c*ones(1,size(sig11,2))*0;
                 H];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [1;0;0];
        
        
        
    case 18
        %Calibration dir 2
        Happ_c = field;
        sig_max = -axial*10^6;
        phi = 0*pi/180;
        %mytitle = sprintf('Constant B^{app} = %.2f T',Happ_c*mu0);
        
        %sig11 = [[0:sig_max/(res):sig_max,sig_max:-sig_max/(res):0],[0:sig_max/(res):sig_max,sig_max:-sig_max/(res):0]];
        %sig11 =(1e6*expdata_dir2(:,1))';
        sig11 = [0:sig_max/(res/2):sig_max sig_max:-sig_max/(res/2):0];
        siglat = zeros(1,length(sig11));
        sig = [sig11;siglat];
        Happ = [ Happ_c*ones(1,size(sig11,2))*0;
                 Happ_c*ones(1,size(sig11,2))*0];
        time = (1:size(Happ,2))*dt;
        xi00 = [0;1];
    case 19
        %Calibration dir 3
        Happ_c = field;
        sig_max = -axial*10^6;
        phi = 0*pi/180;
        %mytitle = sprintf('Constant B^{app} = %.2f T',Happ_c*mu0);
        
        sig11 = 0:sig_max/(res):sig_max;
        %sig11 =(1e6*expdata_dir3(:,1))';
        siglat = zeros(1,length(sig11));
        sig = [sig11;siglat;siglat];
        Happ = [ Happ_c*ones(1,size(sig11,2))*0;
                 Happ_c*ones(1,size(sig11,2))*0;
                 Happ_c*ones(1,size(sig11,2))*0];
        time = (1:size(Happ,2))*dt;
        xi00 = [0;0;1];
        
        
    case 20
        % the axial and lateral stresses are constant and the field varies 
        % lateral stress is applied in 2 direction and the field is applied in 3 direction
        % Note this is "90-deg"
        
        
        % INPUTS:
        Happ_max = field;
        sig1_c = -axial*10^6;
        sig2_c = -lat_stress*10^6;
        
        
        
%           res_pre = res/10;
%           res_load = res-res_pre;
%              
        
        
%         sig11 = [[zeros(1, res_pre)],[0:sig1_max/(res_load):sig1_max,sig1_max:-sig1_max/(res_load):0],[0:sig1_max/(res_load):sig1_max,sig1_max:-sig1_max/(res_load):0]];
%          
%         sig33_pre = 0:sig3_max/(res_pre):sig3_max;
%         sig33_load = sig3_max*ones(1, (size(sig11,2)-size(sig33_pre,2)));
%         sig33 = [[sig33_pre] , [sig33_load]];
%         

        
        Happ33 = [[zeros(1, res+1)], [0:Happ_max/(res):Happ_max,Happ_max:-Happ_max/(res):0],[0:Happ_max/(res):Happ_max,Happ_max:-Happ_max/(res):0],[0:Happ_max/(res):Happ_max,Happ_max:-Happ_max/(res):0]];
        Happ = [zeros(1, size(Happ33,2));
                zeros(1, size(Happ33,2));
                Happ33];
        
        
        % first preload is the lateral stress
        sig22_pre = 0:sig2_c/(res):sig2_c;
        sig22 = [[sig22_pre], [sig2_c*ones(1, size(Happ33,2)-size(sig22_pre,2))]];
        
        % the axial stress is just always there
        sig11 = [sig1_c*ones(1, size(Happ33,2))];
        
        sig = [sig11;sig22;zeros(1,size(sig11,2))];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [1;0;0];
        
        
                
    case 21
        % the axial and lateral stresses are constant and the field varies 
        % lateral stress is applied in 3 direction and the field is applied in 2 direction
        % Note this is "0-deg"
        
        
        % INPUTS:
        Happ_max = field;
        sig1_c = -axial*10^6;
        sig3_c = -lat_stress*10^6;
        
%         res_pre = res/10;
%         res_load = res-res_pre;
        
        
        Happ22 = [[zeros(1, res+1)], [0:Happ_max/(res):Happ_max,Happ_max:-Happ_max/(res):0],[0:Happ_max/(res):Happ_max,Happ_max:-Happ_max/(res):0],[0:Happ_max/(res):Happ_max,Happ_max:-Happ_max/(res):0]];
        Happ = [zeros(1, size(Happ22,2));
                Happ22;
                zeros(1, size(Happ22,2))];
        
        
        % first preload is the lateral stress
        sig33_pre = 0:sig3_c/(res):sig3_c;
        sig33 = [[sig33_pre], [sig3_c*ones(1, size(Happ22,2)-size(sig33_pre,2))]];
        
        % the axial stress is just always there
        sig11 = [sig1_c*ones(1, size(Happ22,2))];
        
        sig = [sig11;zeros(1,size(sig11,2));sig33];
        
        time = (1:size(Happ,2))*dt;
        xi00 = [1;0;0];
        
end


end

