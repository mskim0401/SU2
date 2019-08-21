clear all
close all

BEST_DESIGN = 9

history_project = csvread('history_project.dat',1,0)
history_project_header = {"EVALUATION"    , "LIFT"          , ...
    "DRAG"          , "SIDEFORCE"     , "MOMENT_X"      , "MOMENT_Y"...
    , "MOMENT_Z"      , "FORCE_X"       , "FORCE_Y"       , "FORCE_Z"  ...
    , "EFFICIENCY"    , "TOTAL_HEATFLUX", "MAXIMUM_HEATFLUX", ...
    "INVERSE_DESIGN_LIFT_FIML", "INVERSE_DESIGN_PRESSURE_FIML",...
    "INVERSE_DESIGN_LIFT", "INVERSE_DESIGN_DRAG", ...
    "INVERSE_DESIGN_DRAG_FIML", "ITERATION"    ...
    , "Res_Flow[0]"   , "Res_Flow[1]"   , "Res_Flow[2]" ...
    , "Res_Flow[3]"   , "Res_Flow[4]"   , "Res_Turb[0]"   , "Loss"  ...
    , "Linear_Solver_Iterations", "CFL_Number"    , "TIME"    };


%NOTE: Lift matches S809Delft experiment and Drag matches NREL/Ohio State Experiment pg 51 ? This is weird
%Table B2 S809 Clean Re = 1 Million
AoA2 = [-6.2 -4.1 -2.1 0 2.1 4.1 6.1 8.2 10.1 11.2 12.2 13.3 14.2 15.2 16.2 17.2 18.1 19.2 20 22.1 24];
CL2 = [-0.61 -0.40 -0.16 0.07 0.3 0.55 0.79 0.9 0.94 0.93 0.97 1 1.02 1.03 1.01 0.95 0.9 0.78 0.67 0.7 0.77];
CD2 = [0.0183 0.0004 0.0009 0.0022 0.0037 0.0050 0.0063 0.0096 0.0231 0.0236 0.0368 0.0551 0.0618 0.0705 0.0880 0.1043 0.1325 0.3474 0.3211 0.3699 0.4348];
%https://wind.nrel.gov/airfoils/Coefficients/S809Delft.TXT
AoA = [-3.09 -2.06 -1.01 -0.01 1.02 2.05 3.08 4.1 5.13 6.16 7.17 8.2 9.22 10.21 11.21 12.22 13.24 14.24 15.24 16.24 17.23];
CL = [-0.222 -0.1 0.019 0.139 0.258 0.378 0.497 0.617 0.736 0.851 0.913 0.952 0.973 0.952 0.947 1.007 1.031 1.055 1.062 1.043 0.969];
CD = [0.0071 0.007 0.0095 0.0094 0.0096 0.0099 0.01 0.01 0.0097 0.0095 0.0127 0.0169 0.0247 0.0375 0.0725 0.0636 0.0703 0.0828 0.1081 0.1425 0.1853];

%CD = interp1(AoA2,CD2,AoA,'spline')

%Plot Objective Functions vs Evaluation

i_CL_FIML = 14;
i_CL = 16;

figure
plot(history_project(:,1),history_project(:,i_CL),'-o','LineWidth',2)
hold on
grid on
plot(history_project(:,1),history_project(:,i_CL_FIML),'-ro','LineWidth',2)
set(gca,'FontSize',18)
xlabel('Evaluation')
ylabel('C_L Objective Function')
legend('Without Regularization','With Regularization')

%Plot Initial 
Baseline0 = csvread('Baseline/DIRECT_0/history_direct.dat',3,0);
Baseline1 = csvread('Baseline/DIRECT_1/history_direct.dat',3,0);
Baseline2 = csvread('Baseline/DIRECT_2/history_direct.dat',3,0);
Baseline3 = csvread('Baseline/A5p13_Baseline/history_direct.dat',3,0);
Baseline4 = csvread('Baseline/A11p21_Baseline/history_direct.dat',3,0);
Baseline5 = csvread('Baseline/A12p22_Baseline/history_direct.dat',3,0);
Baseline6 = csvread('Baseline/A15p24_Baseline/history_direct.dat',3,0);

AoA_Base = [1.02 8.2 14.24 5.13 11.21 12.22 15.24];
CL_Base = [Baseline0(length(Baseline0),2) Baseline1(length(Baseline1),2) Baseline2(length(Baseline2),2) Baseline3(length(Baseline3),2) Baseline4(length(Baseline4),2) Baseline5(length(Baseline5),2) Baseline6(length(Baseline6),2)];
CD_Base = [Baseline0(length(Baseline0),3) Baseline1(length(Baseline1),3) Baseline2(length(Baseline2),3) Baseline3(length(Baseline3),3) Baseline4(length(Baseline4),3) Baseline5(length(Baseline5),3) Baseline6(length(Baseline6),3)]; 

if BEST_DESIGN < 10
    Corrected0 = csvread(['DESIGNS/DSN_00' num2str(BEST_DESIGN) '/DIRECT_0/history_direct.dat'],3,0);
    Corrected1 = csvread(['DESIGNS/DSN_00' num2str(BEST_DESIGN) '/DIRECT_1/history_direct.dat'],3,0);
    Corrected2 = csvread(['DESIGNS/DSN_00' num2str(BEST_DESIGN) '/DIRECT_2/history_direct.dat'],3,0);
	Corrected3 = csvread(['DESIGNS/DSN_00' num2str(BEST_DESIGN) '/DIRECT_3/history_direct_direct.dat'],3,0);
	Corrected4 = csvread(['DESIGNS/DSN_00' num2str(BEST_DESIGN) '/DIRECT_4/history_direct_direct.dat'],3,0);
	Corrected5 = csvread(['DESIGNS/DSN_00' num2str(BEST_DESIGN) '/DIRECT_5/history_direct_direct.dat'],3,0);
	Corrected6 = csvread(['DESIGNS/DSN_00' num2str(BEST_DESIGN) '/DIRECT_6/history_direct_direct.dat'],3,0);
else
    if BEST_DESIGN < 100
        Corrected0 = csvread(['DESIGNS/DSN_0' num2str(BEST_DESIGN) '/DIRECT_0/history_direct.dat'],3,0);
        Corrected1 = csvread(['DESIGNS/DSN_0' num2str(BEST_DESIGN) '/DIRECT_1/history_direct.dat'],3,0);
        Corrected2 = csvread(['DESIGNS/DSN_0' num2str(BEST_DESIGN) '/DIRECT_2/history_direct.dat'],3,0);
Corrected3 = csvread(['DESIGNS/DSN_0' num2str(BEST_DESIGN) '/DIRECT_3/history_direct.dat'],3,0);
Corrected4 = csvread(['DESIGNS/DSN_0' num2str(BEST_DESIGN) '/DIRECT_4/history_direct.dat'],3,0);
Corrected5 = csvread(['DESIGNS/DSN_0' num2str(BEST_DESIGN) '/DIRECT_5/history_direct.dat'],3,0);
Corrected6 = csvread(['DESIGNS/DSN_0' num2str(BEST_DESIGN) '/DIRECT_6/history_direct.dat'],3,0);
    else
        Corrected0 = csvread(['DESIGNS/DSN_' num2str(BEST_DESIGN) '/DIRECT_0/history_direct.dat'],3,0);
        Corrected1 = csvread(['DESIGNS/DSN_' num2str(BEST_DESIGN) '/DIRECT_1/history_direct.dat'],3,0);
        Corrected2 = csvread(['DESIGNS/DSN_' num2str(BEST_DESIGN) '/DIRECT_2/history_direct.dat'],3,0);
Corrected3 = csvread(['DESIGNS/DSN_' num2str(BEST_DESIGN) '/DIRECT_3/history_direct.dat'],3,0);  
Corrected4 = csvread(['DESIGNS/DSN_' num2str(BEST_DESIGN) '/DIRECT_4/history_direct.dat'],3,0);  
Corrected5 = csvread(['DESIGNS/DSN_' num2str(BEST_DESIGN) '/DIRECT_5/history_direct.dat'],3,0);  
Corrected6 = csvread(['DESIGNS/DSN_' num2str(BEST_DESIGN) '/DIRECT_6/history_direct.dat'],3,0);      
    end
end

AoA_Corrected = AoA_Base;
CL_Corrected = [Corrected0(length(Corrected0),2) Corrected1(length(Corrected1),2) Corrected2(length(Corrected2),2) Corrected3(length(Corrected3),2) Corrected4(length(Corrected4),2) Corrected5(length(Corrected5),2) Corrected6(length(Corrected6),2)];
CD_Corrected = [Corrected0(length(Corrected0),3) Corrected1(length(Corrected1),3) Corrected2(length(Corrected2),3) Corrected3(length(Corrected3),3) Corrected4(length(Corrected4),3) Corrected5(length(Corrected5),3) Corrected6(length(Corrected6),3)];

%Holdout0 = csvread(['Holdout/A5p13_Baseline/history_direct.dat'],3,0);
%Holdout1 = csvread(['Holdout/A11p21_Baseline/history_direct.dat'],3,0);
%Holdout2 = csvread(['Holdout/A12p22_Baseline/history_direct.dat'],3,0);
%Holdout3 = csvread(['Holdout/A15p24_Baseline/history_direct.dat'],3,0);

%AoA_Holdout = [5.13 11.21 12.22 15.24];

%CL_Holdout = [Holdout0(length(Holdout0),2) Holdout1(length(Holdout1),2) Holdout2(length(Holdout2),2) Holdout3(length(Holdout3),2)];
%CD_Holdout = [Holdout0(length(Holdout0),3) Holdout1(length(Holdout1),3) Holdout2(length(Holdout2),3) Holdout3(length(Holdout3),3)];

%Augmented0 = csvread(['Holdout/A5p13_Augmented/history_direct.dat'],3,0);
%Augmented1 = csvread(['Holdout/A11p21_Augmented/history_direct.dat'],3,0);
%Augmented2 = csvread(['Holdout/A12p22_Augmented/history_direct.dat'],3,0);
%Augmented3 = csvread(['Holdout/A15p24_Augmented/history_direct.dat'],3,0);

%CL_Augmented = [Augmented0(length(Augmented0),2) Augmented1(length(Augmented1),2) Augmented2(length(Augmented2),2) Augmented3(length(Augmented3),2)];
%CD_Augmented = [Augmented0(length(Augmented0),3) Augmented1(length(Augmented1),3) Augmented2(length(Augmented2),3) Augmented3(length(Augmented3),3)];

%AoA_Augmented = AoA_Holdout;

figure
plot(AoA, CL,'-+','LineWidth',2)
hold on
grid on
%plot(AoA2,CL2,'-x','LineWidth',2)
plot(AoA_Base, CL_Base,'bo','LineWidth',2,'MarkerSize',12)
set(gca,'FontSize',18)
xlabel('Angle of Attack (deg)')
ylabel('Lift Coefficient')
legend('Experiment','SA-BC')


figure
plot(AoA, CL,'-+','LineWidth',2)
hold on
grid on
%plot(AoA2,CL2,'-x','LineWidth',2)
plot(AoA_Base, CL_Base,'bo','LineWidth',2,'MarkerSize',12)
plot(AoA_Corrected, CL_Corrected,'b*','LineWidth',2,'MarkerSize',12)
%plot(AoA_Holdout, CL_Holdout,'ko','LineWidth',2,'MarkerSize',12)
%plot(AoA_Augmented, CL_Augmented,'k*','LineWidth',2,'MarkerSize',12)
set(gca,'FontSize',18)
xlabel('Angle of Attack (deg)')
ylabel('Lift Coefficient')
legend('Experiment','No Correction','Augmented')

figure
plot(AoA, CD,'-+','LineWidth',2)
hold on
grid on
%plot(AoA2, CD2,'-x','LineWidth',2)
plot(AoA_Base, CD_Base,'bo','LineWidth',2,'MarkerSize',12)
plot(AoA_Corrected, CD_Corrected,'b*','LineWidth',2,'MarkerSize',12)
%plot(AoA_Holdout, CD_Holdout,'ko','LineWidth',2,'MarkerSize',12)
%plot(AoA_Augmented, CD_Augmented,'k*','LineWidth',2,'MarkerSize',12)
set(gca,'FontSize',18)
xlabel('Angle of Attack (deg)')
ylabel('Drag Coefficient')
legend('Experiment','No Correction','Augmented')

figure
plot(CL, CD,'-+','LineWidth',2)
hold on
grid on
%plot(CL2, CD2,'-+','LineWidth',2)
plot(CL_Base, CD_Base,'bo','LineWidth',2,'MarkerSize',12)
plot(CL_Corrected, CD_Corrected,'b*','LineWidth',2,'MarkerSize',12)
% plot(CL_Holdout, CD_Holdout,'ko','LineWidth',2,'MarkerSize',12)
%plot(CL_Augmented, CD_Augmented,'k*','LineWidth',2,'MarkerSize',12)
set(gca,'FontSize',18)
xlabel('Lift Coefficient')
ylabel('Drag Coefficient')
legend('Experiment','No Correction','Augmented')
