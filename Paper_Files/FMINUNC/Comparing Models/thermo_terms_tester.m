clear;
clc;

%%

r = 8.617333262E-5; % eV/K;

load('crossfit_3term_model.mat');
load('Corners_3term_model.mat');
load('No_DFT_3term_model.mat');

d0_dirs = ["../model3_d0s.mat","../DFT_end_mem_only/model_dft_on_endmem_use_endd0s_ballpark_d0s.mat",...
    "../No_DFT/model_no_dft_d0s.mat"];

fit_col = 3;

%final_Fulleqn = three_term_noDFT_model;
final_Fulleqn = Corners_3term_model;
%final_Fulleqn = crossfit_3term_model;

load(d0_dirs(1));

%% experimental data range

exp_data = readtable('Transformed_exp_data.xlsx');

Press_exp = table2array(exp_data(:,2));  % partial pressure O2
Temp_exp = table2array(exp_data(:,3))+273.15;  % K
x_exp = table2array(exp_data(:,4));  % mol fraciton Ba
dd_exp = table2array(exp_data(:,5));  % delta of delta

min_x = min(x_exp);
max_x = max(x_exp);
min_dd = min(dd_exp);
max_dd = max(dd_exp);

%% Equations from Gibbs

syms T x y
G = simplify(final_Fulleqn);
% S = -dG/dT
S = simplify(-(diff(G,T)));
% H where H = G + TS  
H = simplify(G + T*S);

%% Check self consistant d0

x_vals = [0,0.05,0.1,0.15,0.2];
for k = 1:length(x_vals)
    d = d0_progress(k,end-(fit_col)+1);
    muh_check = double(subs(-diff(G,y),[T x y],[573.15 x_vals(k) d]))
end


%% Check the heat capcity

cp1 = double(subs(diff(H*96.487,T),[x y],[0 0]))*1000/8.314;
fprintf("Here is the heat capacity of SrFeO3 (R): %f \n",cp1)

cp2 = double(subs(diff(H*96.487,T),[x y],[1 0]))*1000/8.314;
fprintf("Here is the heat capacity of BaFeO3 (R): %f \n",cp2)

cp3 = double(subs(diff(H*96.487,T),[x y],[0 0.5]))*1000/8.314;
fprintf("Here is the heat capacity of SrFeO2.5 (R): %f \n",cp3)

cp4 = double(subs(diff(H*96.487,T),[x y],[1 0.5]))*1000/8.314;
fprintf("Here is the heat capacity of BaFeO2.5 (R): %f \n",cp4)

%% H and S Plots

temp1 = 800; % K
temp2 = 1000;
temp3 = 1400;
temp4 = 1800;

%% Entropy comapres between G_end and G_excess
% 
% S_excess = simplify(-(diff(final_G_ex,T)));
% H_excess = simplify(final_G_ex + T*S_excess);
% S_end = simplify(-(diff(final_G_end,T)));
% H_end = simplify(final_G_end + T*S_end);
% 
% 
% figure
% fsurf(subs(S_excess/r,[T],[temp1]),[0 1 0 0.5])
% %title (['Gibbs Free Energy @ ' num2str(temp1) ' K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' S Excess Correcton (R)')
% 
% figure
% fsurf(subs(S_end/r,[T],[temp1]),[0 1 0 0.5])
% %title (['Gibbs Free Energy @ ' num2str(temp1) ' K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' S End Members (R)')
% 
% S_total = simplify(S_end + S_excess);
% % isequal(S,S_total)
% 
% figure
% fsurf(subs(S_total/r,[T],[temp1]),[0 1 0 0.5])
% %title (['Gibbs Free Energy @ ' num2str(temp1) ' K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' S_end + S_excess (R)')


%% -dG/ds = muO

figure
fsurf(-subs(diff(G,y),[T],[temp1]),[0 1 0 0.5])
title (['\partialG/\partial\delta = \mu^O @ ' num2str(temp1) ' K']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel(' \mu^O')

%% S

figure
fsurf(subs(S/r,[T],[temp1]),[0 1 0 0.5])
title (['Entropy @ ' num2str(temp1) ' K']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel(' S (R)')

%% dS/ds

figure
fsurf(diff(subs(S/r,[T],[temp1]),y),[0 1 0.01 0.48])
title (['\partialS/\partial\delta @ ' num2str(temp1) ' K']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel(' \partialS/\partial\delta (R)')


%% dH/ds 

H_red_1 = diff(subs(H,[T],[temp1]),y);
H_red_2 = diff(subs(H,[T],[temp2]),y);
H_red_3 = diff(subs(H,[T],[temp3]),y);
H_red_4 = diff(subs(H,[T],[temp4]),y);


figure
hold on
fs1 = fsurf(H_red_1*96.487,[0 1 0 0.5]);
fs2 = fsurf(H_red_2*96.487,[0 1 0 0.5]);
fs3 = fsurf(H_red_3*96.487,[0 1 0 0.5]);
fs4 = fsurf(H_red_4*96.487,[0 1 0 0.5]);
fs1.FaceColor = '#FE0000';
fs2.FaceColor = '#00CB00';
fs3.FaceColor = '#009898';
fs4.FaceColor = '#FE7300';
fs1.FaceAlpha = 0.75;
fs2.FaceAlpha = 0.75;
fs3.FaceAlpha = 0.75;
fs4.FaceAlpha = 0.75;
fs1.EdgeColor = 'none';
fs2.EdgeColor = 'none';
fs3.EdgeColor = 'none';
fs4.EdgeColor = 'none';

title (['\partialH/\partial\delta']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel('\partialH/\partial\delta (kJ/mol)')
legend(strcat(num2str(temp1), ' K'), strcat(num2str(temp2), ' K'),...
    strcat(num2str(temp3), ' K'),strcat(num2str(temp4), ' K'))
grid on
box on
hold off

%% dH/ds Delta H reduction??

[Ho_1 So_1] = get_O2_thermo(temp1);
[Ho_2 So_2] = get_O2_thermo(temp2);
[Ho_3 So_3] = get_O2_thermo(temp3);
[Ho_4 So_4] = get_O2_thermo(temp4);

H_red_1 = simplify(diff(subs(H,[T],[temp1]),y) + 0.5*Ho_1/98.4875);
H_red_2 = simplify(diff(subs(H,[T],[temp2]),y) + 0.5*Ho_2/98.4875);
H_red_3 = simplify(diff(subs(H,[T],[temp3]),y) + 0.5*Ho_3/98.4875);
H_red_4 = simplify(diff(subs(H,[T],[temp4]),y) + 0.5*Ho_4/98.4875);


figure
hold on
fs1 = fsurf(H_red_1*96.487,[0 1 0 0.5]);
fs2 = fsurf(H_red_2*96.487,[0 1 0 0.5]);
fs3 = fsurf(H_red_3*96.487,[0 1 0 0.5]);
fs4 = fsurf(H_red_4*96.487,[0 1 0 0.5]);
fs1.FaceColor = '#FE0000';
fs2.FaceColor = '#00CB00';
fs3.FaceColor = '#009898';
fs4.FaceColor = '#FE7300';
fs1.FaceAlpha = 0.75;
fs2.FaceAlpha = 0.75;
fs3.FaceAlpha = 0.75;
fs4.FaceAlpha = 0.75;
fs1.EdgeColor = 'none';
fs2.EdgeColor = 'none';
fs3.EdgeColor = 'none';
fs4.EdgeColor = 'none';

title (['\partialH/\partial\delta + 0.5*H_{O_{2}}']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel(' \DeltaH_{red} (kJ/mol)')
legend(strcat(num2str(temp1), ' K'), strcat(num2str(temp2), ' K'),...
    strcat(num2str(temp3), ' K'),strcat(num2str(temp4), ' K'))
grid on
box on
hold off

%% look at T dependence fo dH/ds for x = 0.15

x_frac = 0.15;

H_red_5 = simplify(diff(subs(H,[x],[x_frac]),y));

figure
hold on
fs1 = fsurf(H_red_5*96.487,[500 1500 0 0.5]);


title (['\partialH/\partial\delta @ x = ' num2str(x_frac)]);
xlabel('Temperature (K)')
ylabel('Extent of Reduction (\delta)')
zlabel(' \partialH/\partial\delta (kJ/mol)')
%legend('400 K', '600 K', '800 K', '1000 K')
grid on
box on
hold off


%%  dH/ds for x = 0 + O enthalpy

x_frac = 0;

[Ho_1 So_1] = get_O2_thermo(400+273.15);
[Ho_2 So_2] = get_O2_thermo(575+273.15);
[Ho_3 So_3] = get_O2_thermo(750+273.15);
[Ho_4 So_4] = get_O2_thermo(925+273.15);
[Ho_5 So_5] = get_O2_thermo(1100+273.15);


H_red_400 = simplify(diff(subs(H,[x T],[x_frac (400+273.15)]),y)+ 0.5*Ho_1/98.4875);
H_red_575 = simplify(diff(subs(H,[x T],[x_frac (575+273.15)]),y)+ 0.5*Ho_2/98.4875);
H_red_750 = simplify(diff(subs(H,[x T],[x_frac (750+273.15)]),y)+ 0.5*Ho_3/98.4875);
H_red_925 = simplify(diff(subs(H,[x T],[x_frac (925+273.15)]),y)+ 0.5*Ho_4/98.4875);
H_red_1100 = simplify(diff(subs(H,[x T],[x_frac (1100+273.15)]),y)+ 0.5*Ho_5/98.4875);


figure
hold on
fs1 = fplot(H_red_400*96.487*2,[0 0.5],'linewidth',2.0);
fs2 = fplot(H_red_575*96.487*2,[0 0.5],'linewidth',2.0);
fs3 = fplot(H_red_750*96.487*2,[0 0.5],'linewidth',2.0);
fs4 = fplot(H_red_925*96.487*2,[0 0.5],'linewidth',2.0);
fs5 = fplot(H_red_1100*96.487*2,[0 0.5],'linewidth',2.0);

ylim([0 250])
title (['Full Model SrFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
xlabel('Extent of Reduction (\delta)')
ylabel(' \partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1})')
legend('400 C', '575 C', '750 C', '925 C','1100 C')
box on
hold off

%% look at T dependence fo dS/ds for x = 0

x_frac = 0.00000000000000000000000000001;

S_red_400 = simplify(diff(subs(S,[x T],[x_frac (400+273.15)]),y) + 0.5*So_1/98.4875);
S_red_575 = simplify(diff(subs(S,[x T],[x_frac (575+273.15)]),y) + 0.5*So_2/98.4875);
S_red_750 = simplify(diff(subs(S,[x T],[x_frac (750+273.15)]),y) + 0.5*So_3/98.4875);
S_red_925 = simplify(diff(subs(S,[x T],[x_frac (925+273.15)]),y) + 0.5*So_4/98.4875);
S_red_1100 = simplify(diff(subs(S,[x T],[x_frac (1100+273.15)]),y) + 0.5*So_5/98.4875);


figure
hold on
fs1 = fplot(S_red_400*96.487*2*1000,[0 0.5],'linewidth',2.0); %+205.152
fs2 = fplot(S_red_575*96.487*2*1000,[0 0.5],'linewidth',2.0);
fs3 = fplot(S_red_750*96.487*2*1000,[0 0.5],'linewidth',2.0);
fs4 = fplot(S_red_925*96.487*2*1000,[0 0.5],'linewidth',2.0);
fs5 = fplot(S_red_1100*96.487*2*1000,[0 0.5],'linewidth',2.0);

ylim([0 250])
title (['Full Model SrFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
xlabel('Extent of Reduction (\delta)')
ylabel(' \partialS/\partial\delta (J \bullet (mol O_{2} \bullet K)^{-1})')
legend('400 C', '575 C', '750 C', '925 C','1100 C')
box on
hold off







