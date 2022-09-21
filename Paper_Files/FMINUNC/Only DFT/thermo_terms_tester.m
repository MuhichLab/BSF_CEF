clear;
clc;

%%

expansion = 2;
CEF_build = "phase_shift";
% Constants
r = 8.617333262E-5; % eV/K;

if CEF_build == "phase_shift"
    [G_soln,G_excess,G_end_mem] =CEF_for_plots(expansion);
elseif CEF_build == "no_phase_shift"
    [G_soln,G_excess,G_end_mem] =CEF_for_plots(expansion);
elseif CEF_build == "no_dft"
    [G_soln,G_excess,G_end_mem] =CEF_for_plots(expansion);
end

syms x y T real
syms L [1 24] real
syms gDiffA [1 2] real
syms gDiffB [1 2] real
syms gDiffC [1 2] real
syms gDiffD [1 2] real
syms A [1 24] real
syms B [1 24] real
syms C [1 24] real
syms D [1 24] real


%% We found that L9 = 6*L21 and L11 = 6*L23 thus we combine them into one term
% Note the equiotns no longer have L so the A B C and D terms depeding on
% expansion are set equal this is only for the eqautions that contain the excess terms. 

 G_soln = subs(G_soln,[L21 L23],[L9 L11]);
 G_excess = subs(G_excess,[L21 L23],[L9 L11]);


%% Expand L terms in e1100+273.xcess to be G(T) terms
syms A [1 24] real
syms B [1 24] real
syms C [1 24] real
syms D [1 24] real


sub_me = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13...
    L14 L15 L16 L17 L18 L19 L20 L22 L24];

sub_value = ...
    [A1 A2 ...
    A3 (A4) ...
    (A5) (A6) ...
    (A7) (A8) ...
    (A9) (A10) ...
    (A11) (A12) ...
    (A13) (A14) ...
    (A15) (A16) ...
    (A17) (A18) ...
    (A19) (A20) ...
    (A22) ...
    (A24)];

sub_value_dft = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13...
        A14 A15 A16 A17 A18 A19 A20 A22 A24];

G_excess = subs(G_excess,sub_me,sub_value);
G_soln = subs(G_soln,sub_me,sub_value);


% How many terms we want

load('model_guess0001_no_intial_fits_allfits.mat') % This is all_fits
load('model_guess0001_no_intial_fits_d0s.mat') % this is d0_progress  % only needed for check
%self consisten d0s

fit_col = 3; % number of terms 

d0s = d0_progress(:,end-(fit_col)+1);

% % For monte saves
% all_fits = compiled_coeff(:,:,4);

a1          = all_fits(1,end-(fit_col)+1);
a2          = all_fits(2,end-(fit_col)+1);
a3          = all_fits(3,end-(fit_col)+1);
a4          = all_fits(4,end-(fit_col)+1);
a5          = all_fits(5,end-(fit_col)+1);
a6          = all_fits(6,end-(fit_col)+1);
a7          = all_fits(7,end-(fit_col)+1);
a8          = all_fits(8,end-(fit_col)+1);
a9          = all_fits(9,end-(fit_col)+1);
a10         = all_fits(10,end-(fit_col)+1);
a11         = all_fits(11,end-(fit_col)+1);
a12         = all_fits(12,end-(fit_col)+1);
a13         = all_fits(13,end-(fit_col)+1);
a14         = all_fits(14,end-(fit_col)+1);
a15         = all_fits(15,end-(fit_col)+1);
a16         = all_fits(16,end-(fit_col)+1);
a17         = all_fits(17,end-(fit_col)+1);
a18         = all_fits(18,end-(fit_col)+1);
a19         = all_fits(19,end-(fit_col)+1);
a20         = all_fits(20,end-(fit_col)+1);
a22         = all_fits(21,end-(fit_col)+1);
a24         = all_fits(22,end-(fit_col)+1);
gDiffa1     = all_fits(23,end-(fit_col)+1);
gDiffa2     = all_fits(24,end-(fit_col)+1);    


sub_me = [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                A15,A16,A17,A18,A19,A20,A22,A24,gDiffA1,gDiffA2];
            
sub_value = [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,...
                a15,a16,a17,a18,a19,a20,a22,a24,gDiffa1,gDiffa2];


final_Fulleqn = simplify(expand(subs(G_soln,[sub_me],[sub_value])));
final_G_ex = simplify(expand(subs(G_excess,[sub_me(1:22)],[sub_value(1:22)])));
final_G_end = simplify(expand(subs(G_end_mem,[sub_me(23:24)],[sub_value(23:24)])));

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

syms T
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

S_excess = simplify(-(diff(final_G_ex,T)));
H_excess = simplify(final_G_ex + T*S_excess);
S_end = simplify(-(diff(final_G_end,T)));
H_end = simplify(final_G_end + T*S_end);


figure
fsurf(subs(S_excess/r,[T],[temp1]),[0 1 0 0.5])
%title (['Gibbs Free Energy @ ' num2str(temp1) ' K']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel(' S Excess Correcton (R)')

figure
fsurf(subs(S_end/r,[T],[temp1]),[0 1 0 0.5])
%title (['Gibbs Free Energy @ ' num2str(temp1) ' K']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel(' S End Members (R)')

S_total = simplify(S_end + S_excess);
% isequal(S,S_total)

figure
fsurf(subs(S_total/r,[T],[temp1]),[0 1 0 0.5])
%title (['Gibbs Free Energy @ ' num2str(temp1) ' K']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel(' S_{end} + S_{excess} (R)')



%% G

figure
fsurf(subs(G,[T],[temp1]),[0 1 0 0.5])
title (['Gibbs Free Energy @ ' num2str(temp1) ' K']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel(' G (eV/BSF)')

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

%% S

figure
fsurf(subs(S/r,[x],[1e-6]),[100 1500 0 0.5])
title (['Entropy @  x = 0']);
xlabel('Temp (K)')
ylabel('Extent of Reduction (\delta)')
zlabel(' S (R)')

%% dS/ds

figure
fsurf(diff(subs(S/r,[T],[temp1]),y),[0 1 0.01 0.48])
title (['\partialS/\partial\delta @ ' num2str(temp1) ' K']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel(' \partialS/\partial\delta (R)')

%% dS/ds @ x = 0

figure
fsurf(diff(subs(S/r,[x],[1e-6]),y),[100 1500 0.02 0.48])
title (['\partialS/\partial\delta @ x = 0']);
xlabel('Temp (K)')
ylabel('Extent of Reduction (\delta)')
zlabel(' \partialS/\partial\delta (R)')
%% H where H = G + TS  where S = -dG/dT

figure
fsurf(subs(H,[T],[temp1]),[0 1 0 0.5])
title (['Enthalpy @ ' num2str(temp1) ' K']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel(' H (eV/BSF)')

%% Compare H to DFT data

% DFT data
dft_data = table2array(readtable('BSF_data_phase_shift.xlsx'));
dft_data_copy = dft_data;
rows = any(isnan(dft_data),2);
dft_data(rows,:) = [];

X_dft = dft_data(:,2); % mol fract Ba
Y_dft = dft_data(:,1); % delta
Z_dft = dft_data(:,3)/8; % E eV per ABO3

figure
hold on
fsurf(subs(H,[T],[0]),[0 1 0 0.5])
scatter3(X_dft,Y_dft,Z_dft,100,'or')
title (['Enthalpy @ 0 K']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel(' ~G @ 0 K (eV/ mol BSF)')

%%

% 
% %x = 0
% line0 = dft_data(dft_data(:,2)==0,:);
% %x = 0.125
% line0125 = dft_data(dft_data(:,2)==0.125,:);
% %x = 0.5
% line05 = dft_data(dft_data(:,2)==0.5,:);
% % x = 1.0
% line1 = dft_data(dft_data(:,2)==1,:);
% 
% 
% 
% H_0 = simplify(subs(H,[x T],[0 0]));
% H_0125 = simplify(subs(H,[x T],[0.125 0]));
% H_05 = simplify(subs(H,[x T],[0.5 0]));
% H_1 = simplify(subs(H,[x T],[1.0 0]));
% 
% 
% 
% figure
% hold on
% fplot(H_0,[0 0.5],'linewidth',2.0);
% fplot(H_0125,[0 0.5],'linewidth',2.0);
% fplot(H_05,[0 0.5],'linewidth',2.0);
% fplot(H_1,[0 0.5],'linewidth',2.0);
% scatter(line0(:,1),line0(:,3)/8,100,'xr')
% scatter(line0125(:,1),line0125(:,3)/8,100,'xr')
% scatter(line05(:,1),line05(:,3)/8,100,'xr')
% scatter(line1(:,1),line1(:,3)/8,100,'xr')
% %ylim([0 10])
% title (['Ba_{x}Sr_{1-x}FeO_{3-\delta} @ 0 K']);
% xlabel('Extent of Reduction (\delta)')
% ylabel(' H (eV/mol BSF)')
% legend('X = 0','X = 0.125', 'X = 0.5', 'X = 1.0', 'DFT Data')
% box on
% hold off

%%
% figure
% hold on
% fplot(H_0,[0 0.5],'linewidth',2.0);
% scatter(line0(:,1),line0(:,3)/8,100,'xr')
% %ylim([0 10])
% title (['Ba_{x}Sr_{1-x}FeO_{3-\delta} @ 0 K']);
% xlabel('Extent of Reduction (\delta)')
% ylabel(' H (eV/mol BSF)')
% legend('X = 0','DFT Data')
% box on
% hold off

%%
% %x = 0
% line0 = dft_data(dft_data(:,2)==0,:);
% %x = 0.125
% line0125 = dft_data(dft_data(:,2)==0.125,:);
% %x = 0.5
% line05 = dft_data(dft_data(:,2)==0.5,:);
% % x = 1.0
% line1 = dft_data(dft_data(:,2)==1,:);
% 
% 
% 
% H_0 = simplify(diff(subs(H,[x T],[0 0]),y));
% H_0125 = simplify(diff(subs(H,[x T],[0.125 0]),y));
% H_05 = simplify(diff(subs(H,[x T],[0.5 0]),y));
% H_1 = simplify(diff(subs(H,[x T],[1.0 0]),y));
% 
% figure
% hold on
% fplot(H_0,[0 0.5],'linewidth',2.0);
% fplot(H_0125,[0 0.5],'linewidth',2.0);
% fplot(H_05,[0 0.5],'linewidth',2.0);
% fplot(H_1,[0 0.5],'linewidth',2.0);
% % scatter(line0(:,1),line0(:,3)/8,100,'xr')
% % scatter(line0125(:,1),line0125(:,3)/8,100,'xr')
% % scatter(line05(:,1),line05(:,3)/8,100,'xr')
% % scatter(line1(:,1),line1(:,3)/8,100,'xr')
% ylim([0 1])
% title (['Ba_{x}Sr_{1-x}FeO_{3-\delta} @ 0 K']);
% xlabel('Extent of Reduction (\delta)')
% ylabel(' \partialH/\partial\delta (eV/mol BSF)')
% legend('X = 0','X = 0.125', 'X = 0.5', 'X = 1.0', 'DFT Data')
% box on
% hold off


%% dH/ds 

% H_red_1 = diff(subs(H,[T],[temp1]),y);
% H_red_2 = diff(subs(H,[T],[temp2]),y);
% H_red_3 = diff(subs(H,[T],[temp3]),y);
% H_red_4 = diff(subs(H,[T],[temp4]),y);
% 
% 
% figure
% hold on
% fs1 = fsurf(H_red_1*96.487,[0 1 0 0.5]);
% fs2 = fsurf(H_red_2*96.487,[0 1 0 0.5]);
% fs3 = fsurf(H_red_3*96.487,[0 1 0 0.5]);
% fs4 = fsurf(H_red_4*96.487,[0 1 0 0.5]);
% fs1.FaceColor = '#FE0000';
% fs2.FaceColor = '#00CB00';
% fs3.FaceColor = '#009898';
% fs4.FaceColor = '#FE7300';
% fs1.FaceAlpha = 0.75;
% fs2.FaceAlpha = 0.75;
% fs3.FaceAlpha = 0.75;
% fs4.FaceAlpha = 0.75;
% fs1.EdgeColor = 'none';
% fs2.EdgeColor = 'none';
% fs3.EdgeColor = 'none';
% fs4.EdgeColor = 'none';
% 
% title (['\partialH/\partial\delta']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel('\partialH/\partial\delta (kJ/mol)')
% legend(strcat(num2str(temp1), ' K'), strcat(num2str(temp2), ' K'),...
%     strcat(num2str(temp3), ' K'),strcat(num2str(temp4), ' K'))
% grid on
% box on
% hold off

%% dH/ds Delta H reduction??
% 
% [Ho_1 So_1] = get_O2_thermo(temp1);
% [Ho_2 So_2] = get_O2_thermo(temp2);
% [Ho_3 So_3] = get_O2_thermo(temp3);
% [Ho_4 So_4] = get_O2_thermo(temp4);
% 
% H_red_1 = simplify(diff(subs(H,[T],[temp1]),y) + 0.5*Ho_1/98.4875);
% H_red_2 = simplify(diff(subs(H,[T],[temp2]),y) + 0.5*Ho_2/98.4875);
% H_red_3 = simplify(diff(subs(H,[T],[temp3]),y) + 0.5*Ho_3/98.4875);
% H_red_4 = simplify(diff(subs(H,[T],[temp4]),y) + 0.5*Ho_4/98.4875);
% 
% 
% figure
% hold on
% fs1 = fsurf(H_red_1*96.487,[0 1 0 0.5]);
% fs2 = fsurf(H_red_2*96.487,[0 1 0 0.5]);
% fs3 = fsurf(H_red_3*96.487,[0 1 0 0.5]);
% fs4 = fsurf(H_red_4*96.487,[0 1 0 0.5]);
% fs1.FaceColor = '#FE0000';
% fs2.FaceColor = '#00CB00';
% fs3.FaceColor = '#009898';
% fs4.FaceColor = '#FE7300';
% fs1.FaceAlpha = 0.75;
% fs2.FaceAlpha = 0.75;
% fs3.FaceAlpha = 0.75;
% fs4.FaceAlpha = 0.75;
% fs1.EdgeColor = 'none';
% fs2.EdgeColor = 'none';
% fs3.EdgeColor = 'none';
% fs4.EdgeColor = 'none';
% 
% title (['\partialH/\partial\delta + 0.5*H_{O_{2}}']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' \DeltaH_{red} (kJ/mol)')
% legend(strcat(num2str(temp1), ' K'), strcat(num2str(temp2), ' K'),...
%     strcat(num2str(temp3), ' K'),strcat(num2str(temp4), ' K'))
% grid on
% box on
% hold off

%% look at T dependence fo dH/ds for x = 0.15

% x_frac = 0.15;
% 
% H_red_5 = simplify(diff(subs(H,[x],[x_frac]),y));
% 
% figure
% hold on
% fs1 = fsurf(H_red_5*96.487,[500 1500 0 0.5]);
% 
% 
% title (['\partialH/\partial\delta @ x = ' num2str(x_frac)]);
% xlabel('Temperature (K)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' \partialH/\partial\delta (kJ/mol)')
% %legend('400 K', '600 K', '800 K', '1000 K')
% grid on
% box on
% hold off


%%  dH/ds for x = 0 + O enthalpy

x_frac = 1e-6;

[Ho_1 So_1] = get_O2_thermo(400+273.15);
[Ho_2 So_2] = get_O2_thermo(575+273.15);
[Ho_3 So_3] = get_O2_thermo(750+273.15);
[Ho_4 So_4] = get_O2_thermo(925+273.15);
[Ho_5 So_5] = get_O2_thermo(1100+273.15);
[Ho_6 So_6] = get_O2_thermo(1200+273.15);


H_red_400 = simplify(diff(subs(H,[x T],[x_frac (400+273.15)]),y)+ 0.5*Ho_1/98.4875);
H_red_575 = simplify(diff(subs(H,[x T],[x_frac (575+273.15)]),y)+ 0.5*Ho_2/98.4875);
H_red_750 = simplify(diff(subs(H,[x T],[x_frac (750+273.15)]),y)+ 0.5*Ho_3/98.4875);
H_red_925 = simplify(diff(subs(H,[x T],[x_frac (925+273.15)]),y)+ 0.5*Ho_4/98.4875);
H_red_1100 = simplify(diff(subs(H,[x T],[x_frac (1100+273.15)]),y)+ 0.5*Ho_5/98.4875);
%H_red_1200 = simplify(diff(subs(H,[x T],[x_frac (1200+273.15)]),y)+ 0.5*Ho_6/98.4875);


figure
hold on
%fs6 = fplot(H_red_1200*96.487*2,[0 0.5],'linewidth',2.0);
fs5 = fplot(H_red_1100*96.487*2,[0 0.5],'linewidth',2.0);
fs4 = fplot(H_red_925*96.487*2,[0 0.5],'linewidth',2.0);
fs3 = fplot(H_red_750*96.487*2,[0 0.5],'linewidth',2.0);
fs2 = fplot(H_red_575*96.487*2,[0 0.5],'linewidth',2.0);
fs1 = fplot(H_red_400*96.487*2,[0 0.5],'linewidth',2.0);

% plot([0,0.5],[211,211],'--r');
% plot([0,0.5],[191,191],'--r');
% plot([0,0.5],[168,168],'--b');

ylim([80 220])
%title (['Full Model Ba_xSr_1-xFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
xlabel('Extent of Reduction (\delta)')
ylabel(' \partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1})')
legend('1100 ^\circC','925   ^\circC','750   ^\circC','575   ^\circC','400   ^\circC')
box on
hold off


%%  dH/ds heat maps over x and delta at a specific temperture

T_map = 800 + 273.15;

[Ho_map So_map] = get_O2_thermo(T_map);


H_red_map = simplify(diff(subs(H,[T],[(T_map)]),y)+ 0.5*Ho_map/98.4875);


figure
hold on

fs5 = fsurf(H_red_map*96.487*2,[0 1.0 0 0.5],'edgecolor','none');

title (['T = ' num2str(T_map-273.15) ' ^{\circ}C']);
xlabel('Mol Fraction Ba (x)');
ylabel('Extent of Reduction (\delta)');
zlabel(' \partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1})');
colormap jet
cb = colorbar
cb.Label.String = ['\partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1}'];
cb.Limits = [80, 200];
cb.Location = 'northoutside'
%legend('1100 ^\circC','925   ^\circC','750   ^\circC','575   ^\circC','400   ^\circC')
box on
hold off


%% 3D temperture x s heatmaps

syms T
a=30.03235;b=8.772972;c=-3.988133;d=0.788313;e=-0.741599;f=-11.32468;g=236.1663;h=0;
Ho_T = a*T/1000 + b*(T/1000)^2/2 + c*(T/1000)^3/3 + d*(T/1000)^4/4 - e/(T/1000) + f - h; %kJ/mol
H_red_Tmap = matlabFunction((simplify(diff(H,y) + 0.5*Ho_T/98.4875))*96.487*2);

N = 50;
i = 1;
MAP = [];
for mol = linspace(0,1,N)
    for del = linspace(0,0.5,N)
        for TEMP = linspace(500,1200,N)
            val = H_red_Tmap(TEMP+273.15,mol,del);
            MAP(i,:) = [mol,del,TEMP,val];
            i = i + 1;
        end
    end
end


figure
box on
ax = gca;
ax.BoxStyle = 'full';
scatter3(MAP(:,1),MAP(:,2),MAP(:,3),[],MAP(:,4),'filled');
cb = colorbar
cb.Label.String = ['\partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1}'];
cb.Location = 'northoutside'
colormap jet
cb.Label.String = ['\partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1}'];
xlabel('Mol Fraction Ba (x)');
ylabel('Extent of Reduction (\delta)');
zlabel('Temperture (C^\circ)');
xlim([-0.01 1])
ylim([-0.01 0.5])



%% look at T dependence fo dS/ds for x = 0

x_frac = 1e-6;

S_red_400 = simplify(diff(subs(S,[x T],[x_frac (400+273.15)]),y) + 0.5*So_1/98.4875);
S_red_575 = simplify(diff(subs(S,[x T],[x_frac (575+273.15)]),y) + 0.5*So_2/98.4875);
S_red_750 = simplify(diff(subs(S,[x T],[x_frac (750+273.15)]),y) + 0.5*So_3/98.4875);
S_red_925 = simplify(diff(subs(S,[x T],[x_frac (925+273.15)]),y) + 0.5*So_4/98.4875);
S_red_1100 = simplify(diff(subs(S,[x T],[x_frac (1100+273.15)]),y) + 0.5*So_5/98.4875);

figure
hold on

fs5 = fplot(S_red_1100*96.487*2*1000,[0 0.5],'linewidth',2.0);
fs4 = fplot(S_red_925*96.487*2*1000,[0 0.5],'linewidth',2.0);
fs3 = fplot(S_red_750*96.487*2*1000,[0 0.5],'linewidth',2.0);
fs2 = fplot(S_red_575*96.487*2*1000,[0 0.5],'linewidth',2.0);
fs1 = fplot(S_red_400*96.487*2*1000,[0 0.5],'linewidth',2.0);

%ylim([0 500])
%title (['Full Model SrFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
xlabel('Extent of Reduction (\delta)')
ylabel(' \partialS/\partial\delta (J \bullet (mol O_{2} \bullet K)^{-1})')
legend('1100 ^\circC','925   ^\circC','750   ^\circC','575   ^\circC','400   ^\circC')
box on
hold off

%%  dH/ds heat maps over x and delta at a specific temperture

T_map = 800 + 273.15;

[Ho_map So_map] = get_O2_thermo(T_map);


S_red_map = simplify(diff(subs(S,[T],[(T_map)]),y)+ 0.5*So_map/98.4875);


figure
hold on
fs5 = fsurf(S_red_map*96.487*1000*2,[0 1.0 0.01 0.49],'edgecolor','none');


title (['T = ' num2str(T_map-273.15) ' ^{\circ}C']);
xlabel('Mol Fraction Ba (x)');
ylabel('Extent of Reduction (\delta)');
zlabel(' \partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1})');
colormap jet
cb = colorbar
cb.Label.String = ['\partialS/\partial\delta (J \bullet (mol O_{2})^{-1}'];
cb.Limits = [0, 300];
cb.Location = "northoutside";
%legend('1100 ^\circC','925   ^\circC','750   ^\circC','575   ^\circC','400   ^\circC')
box on
hold off


%% 3D temperture x s heatmaps

syms T
a=30.03235;b=8.772972;c=-3.988133;d=0.788313;e=-0.741599;f=-11.32468;g=236.1663;h=0;
So_T = (a*log((T/1000)) + b*(T/1000) + c*(T/1000)^2/2 + d*(T/1000)^3/3 - e/(2*(T/1000)^2) + g)/1000; %kJ/mol
S_red_Tmap = matlabFunction((simplify(diff(S,y) + 0.5*So_T/98.4875))*96.487*1000*2);

N = 50;
i = 1;
MAP = [];
for del = linspace(0.015,0.485,N)
    for TEMP = linspace(500,1200,N)
        val = S_red_Tmap(TEMP+273.15,del);
        MAP(i,:) = [del,TEMP,val];
        i = i + 1;
    end
end
    


figure()

scatter(MAP(:,1),MAP(:,2),[],MAP(:,3),'filled');
cb = colorbar
colormap jet
cb.Label.String = ['\partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1}'];
ylabel('Temp (K)');
xlabel('Extent of Reduction (\delta)');
xlim([0.05 0.45])
ylim([500 1200])
%cb.Limits = [0, 300];
cb.Label.String = ['\partialS/\partial\delta (J \bullet (mol O_{2})^{-1}'];
cb.Location = 'northoutside'

box on
hold off

%%  dS/ds heat maps over x and delta and tmeperature with colorbar for reduction energy

%T_map = 800 + 273.15;

[Ho_map So_map] = get_O2_thermo(T_map);


S_red_map = simplify(diff(subs(S,[T],[(T_map)]),y)+ 0.5*So_map/98.4875);


figure
hold on
fs5 = fsurf(S_red_map*96.487*1000*2,[0 1.0 0.01 0.49],'edgecolor','none');


title (['T = ' num2str(T_map-273.15) ' ^{\circ}C']);
xlabel('Mol Fraction Ba (x)');
ylabel('Extent of Reduction (\delta)');
zlabel(' \partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1})');
colormap jet
cb = colorbar
cb.Label.String = ['\partialS/\partial\delta (J \bullet (mol O_{2})^{-1}'];
cb.Limits = [0, 300];
cb.Location = "northoutside";
%legend('1100 ^\circC','925   ^\circC','750   ^\circC','575   ^\circC','400   ^\circC')
box on
hold off
%% reduction energy compares

% DFT data
dft_data = table2array(readtable('BSF_data_phase_shift.xlsx'));
dft_data_copy = dft_data;
rows = any(isnan(dft_data),2);
dft_data(rows,:) = [];

X_dft = dft_data(:,2); % mol fract Ba
Y_dft = dft_data(:,1); % delta
Z_dft = dft_data(:,3)/8; % E eV per ABO3


test_T = 400+273.15;
test_T2 = 800+273.15;


%x = 0
line0 = dft_data(dft_data(:,2)==0,:);
(line0(end,3) - line0(1,3))/8
(line0(end-1,3) - line0(1,3))/8
(line0(end-2,3) - line0(1,3))/8
(line0(end-3,3) - line0(1,3))/8
(line0(end-4,3) - line0(1,3))/8
(line0(end-5,3) - line0(1,3))/8
%x = 0.125
line0125 = dft_data(dft_data(:,2)==0.125,:);
%x = 0.5
line05 = dft_data(dft_data(:,2)==0.5,:);
% x = 1.0
line1 = dft_data(dft_data(:,2)==1,:);
(line1(end,3) - line1(1,3))/8;

H_0 = subs(H,[x T],[0 test_T]);
%H_0_T = subs(H,[x T],[0 test_T2]);
dH_0 = diff(subs(H,[x T],[0 test_T]),y);
%dH_0_T = diff(subs(H,[x T],[0 test_T2]),y);
dH_0125 = diff(subs(H,[x T],[0.125 test_T]),y);
dH_05 = diff(subs(H,[x T],[0.5 test_T]),y);
dH_1 = diff(subs(H,[x T],[0.999 test_T]),y);


double(int(dH_0,0,0.5))
dH_ds0 = matlabFunction(dH_0);
matH_0 = matlabFunction(H_0);
double(int(dH_05,0,0.25));
dH_ds05 = matlabFunction(dH_05);

s = 0:0.01:0.5;

fig = figure;
left_color = [0 0 0];
right_color = [140/255 29/255 64/255];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
ax = gca;
hold on
fplot(dH_0,[0 0.5],'-k','linewidth',2.0);
%fplot(dH_0_T,[0 0.5],'-','Color','#528D1D','linewidth',2.0);
%area(s,dH_ds0(s),'FaceColor','#8C1D40','Linewidth',4,'FaceAlpha',0.75)
ylabel('\partialH/\partial\delta (eV \bullet (mol BSF)^{-1})')
ylim([0,1.2])
yyaxis right
fplot(H_0-subs(H_0,y,0),[0 0.5],'Color','#FFC627','Linewidth',2)
%fplot(H_0_T-subs(H_0_T,y,0),[0 0.5],'Color','#528D1D','Linewidth',2)

scatter(0.0630,(line0(end-4,3) - line0(1,3))/8,100,'*','MarkerFaceColor','#8C1D40')
scatter(0.1250,(line0(end-3,3) - line0(1,3))/8,100,'*','MarkerFaceColor','#8C1D40')
scatter(0.1880,(line0(end-2,3) - line0(1,3))/8,100,'*','MarkerFaceColor','#8C1D40')
scatter(0.2500,(line0(end-1,3) - line0(1,3))/8,100,'*','MarkerFaceColor','#8C1D40')
scatter(0.5000,(line0(end  ,3) - line0(1,3))/8,100,'*','MarkerFaceColor','#8C1D40')

scatter(0.0630,matH_0(0.0630)-subs(H_0,y,0),75,'d','MarkerFaceColor','#FFC627')
scatter(0.1250,matH_0(0.1250)-subs(H_0,y,0),75,'d','MarkerFaceColor','#FFC627')
scatter(0.1880,matH_0(0.1880)-subs(H_0,y,0),75,'d','MarkerFaceColor','#FFC627')
scatter(0.2500,matH_0(0.2500)-subs(H_0,y,0),75,'d','MarkerFaceColor','#FFC627')
scatter(0.5000,matH_0(0.5000)-subs(H_0,y,0),75,'d','MarkerFaceColor','#FFC627')

plot([0.0630 0.5],[matH_0(0.0630)-subs(H_0,y,0) matH_0(0.0630)-subs(H_0,y,0)],'--','Color','#FFC627','Linewidth',2)
plot([0.1250 0.5],[matH_0(0.1250)-subs(H_0,y,0) matH_0(0.1250)-subs(H_0,y,0)],'--','Color','#FFC627','Linewidth',2)
plot([0.1880 0.5],[matH_0(0.1880)-subs(H_0,y,0) matH_0(0.1880)-subs(H_0,y,0)],'--','Color','#FFC627','Linewidth',2)
plot([0.2500 0.5],[matH_0(0.2500)-subs(H_0,y,0) matH_0(0.2500)-subs(H_0,y,0)],'--','Color','#FFC627','Linewidth',2)

plot([0.0630 0.5],[(line0(end-4,3) - line0(1,3))/8 (line0(end-4,3) - line0(1,3))/8],'--','Color','#8C1D40','Linewidth',2)
plot([0.1250 0.5],[(line0(end-3,3) - line0(1,3))/8 (line0(end-3,3) - line0(1,3))/8],'--','Color','#8C1D40','Linewidth',2)
plot([0.1880 0.5],[(line0(end-2,3) - line0(1,3))/8 (line0(end-2,3) - line0(1,3))/8],'--','Color','#8C1D40','Linewidth',2)
plot([0.2500 0.5],[(line0(end-1,3) - line0(1,3))/8 (line0(end-1,3) - line0(1,3))/8],'--','Color','#8C1D40','Linewidth',2)

ylim([0,0.5])
title (['SrFeO_{3-\delta} @ 800 ^\circC']);
xlabel('Extent of Reduction (\delta)')
ylabel('H - H(\delta=0) (eV \bullet (mol BSF)^{-1})')
legend('Model \partialH/\partial\delta', 'Model H - H(\delta=0)','DFT')
box on
hold off

%% reduction energy compares for end-member only

% %x = 0
% line0 = dft_data(dft_data(:,2)==0,:);
% (line0(end,3) - line0(1,3))/8
% (line0(end-1,3) - line0(1,3))/8
% (line0(end-2,3) - line0(1,3))/8
% (line0(end-3,3) - line0(1,3))/8
% (line0(end-4,3) - line0(1,3))/8
% (line0(end-5,3) - line0(1,3))/8
% %x = 0.125
% line0125 = dft_data(dft_data(:,2)==0.125,:);
% %x = 0.5
% line05 = dft_data(dft_data(:,2)==0.5,:);
% % x = 1.0
% line1 = dft_data(dft_data(:,2)==1,:);
% (line1(end,3) - line1(1,3))/8;
% 
% H_0 = subs(H_end,[x T],[0 0]);
% dH_0 = diff(subs(H_end,[x T],[0 0]),y);
% dH_0125 = diff(subs(H_end,[x T],[0.125 0]),y);
% dH_05 = diff(subs(H_end,[x T],[0.5 0]),y);
% dH_1 = diff(subs(H_end,[x T],[0.999 0]),y);
% 
% 
% double(int(dH_0,0,0.5));
% dH_ds0 = matlabFunction(dH_0);
% matH_0 = matlabFunction(H_0);
% double(int(dH_05,0,0.25));
% dH_ds05 = matlabFunction(dH_05);
% 
% s = 0:0.01:0.5;
% 
% fig = figure;
% left_color = [0 0 0];
% right_color = [140/255 29/255 64/255];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% ax = gca;
% hold on
% fplot(dH_0,[0 0.5],'-k','linewidth',2.0);
% %area(s,dH_ds0(s),'FaceColor','#8C1D40','Linewidth',4,'FaceAlpha',0.75)
% ylabel('\partialH/\partial\delta (eV \bullet (mol BSF)^{-1})')
% ylim([0,1.2])
% yyaxis right
% fplot(H_0-subs(H_0,y,0),[0 0.5],'Color','#FFC627','Linewidth',2)
% 
% scatter(0.0630,(line0(end-4,3) - line0(1,3))/8,100,'*','MarkerFaceColor','#8C1D40')
% scatter(0.1250,(line0(end-3,3) - line0(1,3))/8,100,'*','MarkerFaceColor','#8C1D40')
% scatter(0.1880,(line0(end-2,3) - line0(1,3))/8,100,'*','MarkerFaceColor','#8C1D40')
% scatter(0.2500,(line0(end-1,3) - line0(1,3))/8,100,'*','MarkerFaceColor','#8C1D40')
% scatter(0.5000,(line0(end  ,3) - line0(1,3))/8,100,'*','MarkerFaceColor','#8C1D40')
% 
% scatter(0.0630,matH_0(0.0630)-subs(H_0,y,0),75,'d','MarkerFaceColor','#FFC627')
% scatter(0.1250,matH_0(0.1250)-subs(H_0,y,0),75,'d','MarkerFaceColor','#FFC627')
% scatter(0.1880,matH_0(0.1880)-subs(H_0,y,0),75,'d','MarkerFaceColor','#FFC627')
% scatter(0.2500,matH_0(0.2500)-subs(H_0,y,0),75,'d','MarkerFaceColor','#FFC627')
% scatter(0.5000,matH_0(0.5000)-subs(H_0,y,0),75,'d','MarkerFaceColor','#FFC627')
% 
% plot([0.0630 0.5],[matH_0(0.0630)-subs(H_0,y,0) matH_0(0.0630)-subs(H_0,y,0)],'--','Color','#FFC627','Linewidth',2)
% plot([0.1250 0.5],[matH_0(0.1250)-subs(H_0,y,0) matH_0(0.1250)-subs(H_0,y,0)],'--','Color','#FFC627','Linewidth',2)
% plot([0.1880 0.5],[matH_0(0.1880)-subs(H_0,y,0) matH_0(0.1880)-subs(H_0,y,0)],'--','Color','#FFC627','Linewidth',2)
% plot([0.2500 0.5],[matH_0(0.2500)-subs(H_0,y,0) matH_0(0.2500)-subs(H_0,y,0)],'--','Color','#FFC627','Linewidth',2)
% 
% plot([0.0630 0.5],[(line0(end-4,3) - line0(1,3))/8 (line0(end-4,3) - line0(1,3))/8],'--','Color','#8C1D40','Linewidth',2)
% plot([0.1250 0.5],[(line0(end-3,3) - line0(1,3))/8 (line0(end-3,3) - line0(1,3))/8],'--','Color','#8C1D40','Linewidth',2)
% plot([0.1880 0.5],[(line0(end-2,3) - line0(1,3))/8 (line0(end-2,3) - line0(1,3))/8],'--','Color','#8C1D40','Linewidth',2)
% plot([0.2500 0.5],[(line0(end-1,3) - line0(1,3))/8 (line0(end-1,3) - line0(1,3))/8],'--','Color','#8C1D40','Linewidth',2)
% 
% ylim([0,0.5])
% title (['SrFeO_{3-\delta} @ 0 K']);
% xlabel('Extent of Reduction (\delta)')
% ylabel('H - H(\delta=0) (eV \bullet (mol BSF)^{-1})')
% legend('Model \partialH/\partial\delta', 'Model H - H(\delta=0)','DFT')
% box on
% hold off

%% dH/dT = Cp 

% Cp = simplify(diff(H,T));
% 
% figure
% hold on
% fs5 = fsurf(Cp*96.487*1000,[0 1.0 0 0.5]);
% 
% %fs5.Color = 'r';
% 
% 
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' \partialH/\partialT (J/mol*K)')
% title ('C^O_{P}');
% grid on
% box on
% hold off

%% T*dS/dT = Cp @ x = 0

% Cp = T*simplify(diff(S,T));
% 
% figure
% hold on
% fs5 = fsurf(Cp*96.487*1000,[0 1 0 0.5]);
% 
% %fs5.Color = 'r';
% 
% 
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' T*\partialS/\partialT (J/mol*K)')
% title ('C^O_{P}');
% grid on
% box on
% hold off

%% dH/dT = Cp + 1/2 H_o2  <700K

% new_H = H + 0.5*(31.32234*T/1000 - 20.23531*(T/1000)^2/2 + ...
%     57.86644*(T/1000)^3/3 - 36.50624*(T/1000)^4/4 + 0.007374/(T/1000) - 8.903471)/98.4875; %kJ/mol --> eV
% 
% Cp = simplify(subs(diff(new_H,T),[x],[0]));
% 
% figure
% hold on
% fs5 = fsurf(Cp*96.487*1000,[200 700 0 0.5]);
% 
% %fs5.Color = 'r';
% 
% 
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' \partialH/\partialT (J/mol*K)')
% title ('C^O_{P}');
% grid on
% box on
% hold off

%%

% syms x y
% syms A [1 11]
% syms B [1 11]
% syms C [1 11]
% 
% 
% Ya_Sr = 1-x;
% Ya_Ba = x;
% %Yb_3 = 6-2*y;   % In terms of y where y = (3 - delta)
% Yb_3 = 2*y;     % In terms of delta
% %Yb_4 = 2*y-5;   % In terms of y where y = (3 - delta)
% Yb_4 = 1-2*y;   % In terms of delta
% %Yo_o = y/3;     % In terms of y where y = (3 - delta)
% Yo_o = 1-y/3;   % In terms of delta
% %Yo_Va = 1-y/3;  % In terms of y where y = (3 - delta)
% Yo_Va = y/3;    % In terms of delta
% 
% 
% nine = Yb_3*Yb_4*Ya_Sr*Yo_o*A9;
% ten = Yb_3*Yb_4*Ya_Sr*Yo_o*(Yb_3 - https://www.vasp.at/wiki/index.php/Band_decomposed_charge_densitiesYb_4)*A10;
% eleven = Yb_3*Yb_4*Ya_Ba*Yo_o*A11;
% 
% 
% all = nine + ten + eleven;
% 
% diff_all = simplify(diff(all,y)) % excess
% 
% diff_S = simplify(diff(S,y) - 6*m*r*T) % configurational entropy terms
% 
% diff_end = simplify(my_diffG + 6*m*r*T)



