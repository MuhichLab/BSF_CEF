clear;
clc;

%%

expansion = 2;

% Constants
r = 8.617333262E-5; % eV/K;


[G_soln,G_excess,G_end_mem] =CEF_for_plots_no_dft(expansion);


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


%% Expand L terms in excess to be G(T) terms
syms A [1 24] real
syms B [1 24] real
syms C [1 24] real
syms D [1 24] real

if expansion == 1

    sub_me = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13...
        L14 L15 L16 L17 L18 L19 L20 L22 L24];
    sub_value = [(A1 + B1*T) (A2 + B2*T) (A3 + B3*T) (A4 + B4*T)...
        (A5 + B5*T) (A6 + B6*T) (A7 + B7*T) (A8 + B8*T) (A9 + B9*T)...
        (A10 + B10*T) (A11 + B11*T) (A12 + B12*T) (A13 + B13*T)...
        (A14 + B14*T) (A15 + B15*T) (A16 + B16*T) (A17 + B17*T)...
        (A18 + B18*T) (A19 + B19*T) (A20 + B20*T) (A22 + B22*T) (A24 + B24*T)];
    sub_value_dft = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13...
            A14 A15 A16 A17 A18 A19 A20 A22 A24];
    
    G_excess = subs(G_excess,sub_me,sub_value);
    G_soln = subs(G_soln,sub_me,sub_value);


elseif expansion == 2

    sub_me = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13...
        L14 L15 L16 L17 L18 L19 L20 L22 L24];
    
    sub_value = ...
        [(A1 + B1*T + C1*T*log(T)) (A2 + B2*T + C2*T*log(T)) ...
        (A3 + B3*T + C3*T*log(T)) (A4 + B4*T + C4*T*log(T)) ...
        (A5 + B5*T + C5*T*log(T)) (A6 + B6*T + C6*T*log(T)) ...
        (A7 + B7*T + C7*T*log(T)) (A8 + B8*T + C8*T*log(T)) ...
        (A9 + B9*T + C9*T*log(T)) (A10 + B10*T + C10*T*log(T)) ...
        (A11 + B11*T + C11*T*log(T)) (A12 + B12*T + C12*T*log(T)) ...
        (A13 + B13*T + C13*T*log(T)) (A14 + B14*T + C14*T*log(T)) ...
        (A15 + B15*T + C15*T*log(T)) (A16 + B16*T + C16*T*log(T)) ...
        (A17 + B17*T + C17*T*log(T)) (A18 + B18*T + C18*T*log(T)) ...
        (A19 + B19*T + C19*T*log(T)) (A20 + B20*T + C20*T*log(T)) ...
        (A22 + B22*T + C22*T*log(T)) ...
        (A24 + B24*T + C24*T*log(T))];
    
    sub_value_dft = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13...
            A14 A15 A16 A17 A18 A19 A20 A22 A24];
    
    G_excess = subs(G_excess,sub_me,sub_value);
    G_soln = subs(G_soln,sub_me,sub_value);


elseif expansion == 3

    sub_me = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13...
        L14 L15 L16 L17 L18 L19 L20 L22 L24];
    
    sub_value = ...
        [(A1 + B1*T + C1*T*log(T) + D1*T^2/2) (A2 + B2*T + C2*T*log(T) + D2*T^2/2) ...
        (A3 + B3*T + C3*T*log(T) + D3*T^2/2) (A4 + B4*T + C4*T*log(T) + D4*T^2/2) ...
        (A5 + B5*T + C5*T*log(T) + D5*T^2/2) (A6 + B6*T + C6*T*log(T) + D6*T^2/2) ...
        (A7 + B7*T + C7*T*log(T) + D7*T^2/2) (A8 + B8*T + C8*T*log(T) + D8*T^2/2) ...
        (A9 + B9*T + C9*T*log(T) + D9*T^2/2) (A10 + B10*T + C10*T*log(T) + D10*T^2/2) ...
        (A11 + B11*T + C11*T*log(T) + D11*T^2/2) (A12 + B12*T + C12*T*log(T) + D12*T^2/2) ...
        (A13 + B13*T + C13*T*log(T) + D13*T^2/2) (A14 + B14*T + C14*T*log(T) + D14*T^2/2) ...
        (A15 + B15*T + C15*T*log(T) + D15*T^2/2) (A16 + B16*T + C16*T*log(T) + D16*T^2/2) ...
        (A17 + B17*T + C17*T*log(T) + D17*T^2/2) (A18 + B18*T + C18*T*log(T) + D18*T^2/2) ...
        (A19 + B19*T + C19*T*log(T) + D19*T^2/2) (A20 + B20*T + C20*T*log(T) + D20*T^2/2) ...
        (A22 + B22*T + C22*T*log(T) + D22*T^2/2) ...
        (A24 + B24*T + C24*T*log(T) + D24*T^2/2)];
    
    sub_value_dft = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13...
            A14 A15 A16 A17 A18 A19 A20 A22 A24];
        
    G_excess = subs(G_excess,sub_me,sub_value);
    G_soln = subs(G_soln,sub_me,sub_value);
end

% How many terms we want

load('model_guess0001_allfits.mat') % This is all_fits
load('model_guess0001_d0s.mat') % this is d0_progress  % only needed for check
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
b1          = all_fits(23,end-(fit_col)+1);
b2          = all_fits(24,end-(fit_col)+1);
b3          = all_fits(25,end-(fit_col)+1);
b4          = all_fits(26,end-(fit_col)+1);
b5          = all_fits(27,end-(fit_col)+1);
b6          = all_fits(28,end-(fit_col)+1);
b7          = all_fits(29,end-(fit_col)+1);
b8          = all_fits(30,end-(fit_col)+1);
b9          = all_fits(31,end-(fit_col)+1);
b10         = all_fits(32,end-(fit_col)+1);
b11         = all_fits(33,end-(fit_col)+1);
b12         = all_fits(34,end-(fit_col)+1);
b13         = all_fits(35,end-(fit_col)+1);
b14         = all_fits(36,end-(fit_col)+1);
b15         = all_fits(37,end-(fit_col)+1);
b16         = all_fits(38,end-(fit_col)+1);
b17         = all_fits(39,end-(fit_col)+1);
b18         = all_fits(40,end-(fit_col)+1);
b19         = all_fits(41,end-(fit_col)+1);
b20         = all_fits(42,end-(fit_col)+1);
b22         = all_fits(43,end-(fit_col)+1);
b24         = all_fits(44,end-(fit_col)+1);
c1          = all_fits(45,end-(fit_col)+1);
c2          = all_fits(46,end-(fit_col)+1);
c3          = all_fits(47,end-(fit_col)+1);
c4          = all_fits(48,end-(fit_col)+1);
c5          = all_fits(49,end-(fit_col)+1);
c6          = all_fits(50,end-(fit_col)+1);
c7          = all_fits(51,end-(fit_col)+1);
c8          = all_fits(52,end-(fit_col)+1);
c9          = all_fits(53,end-(fit_col)+1);
c10         = all_fits(54,end-(fit_col)+1);
c11         = all_fits(55,end-(fit_col)+1);
c12         = all_fits(56,end-(fit_col)+1);
c13         = all_fits(57,end-(fit_col)+1);
c14         = all_fits(58,end-(fit_col)+1);
c15         = all_fits(59,end-(fit_col)+1);
c16         = all_fits(60,end-(fit_col)+1);
c17         = all_fits(61,end-(fit_col)+1);
c18         = all_fits(62,end-(fit_col)+1);
c19         = all_fits(63,end-(fit_col)+1);
c20         = all_fits(64,end-(fit_col)+1);
c22         = all_fits(65,end-(fit_col)+1);
c24         = all_fits(66,end-(fit_col)+1);
gDiffa1     = all_fits(67,end-(fit_col)+1);
gDiffa2     = all_fits(68,end-(fit_col)+1);    
gDiffb1     = all_fits(69,end-(fit_col)+1);
gDiffb2     = all_fits(70,end-(fit_col)+1);
gDiffc1     = all_fits(71,end-(fit_col)+1);
gDiffc2     = all_fits(72,end-(fit_col)+1);


sub_me = [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                A15,A16,A17,A18,A19,A20,A22,A24,B1,B2,B3,B4,B5,B6,...
                B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,B18,B19,B20,B22,...
                B24,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,...
                C17,C18,C19,C20,C22,C24,gDiffA1,gDiffA2,gDiffB1,gDiffB2,gDiffC1,gDiffC2];
            
sub_value = [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,...
                a15,a16,a17,a18,a19,a20,a22,a24,b1,b2,b3,b4,b5,b6,...
                b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b22,...
                b24,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,...
                c17,c18,c19,c20,c22,c24,gDiffa1,gDiffa2,gDiffb1,gDiffb2,gDiffc1,gDiffc2];


final_Fulleqn = simplify(subs(G_soln,[sub_me],[sub_value]));
final_G_ex = simplify(subs(G_excess,[sub_me(1:66)],[sub_value(1:66)]));
final_G_end = simplify(subs(G_end_mem,[sub_me(67:72)],[sub_value(67:72)]));

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
zlabel(' S_end + S_excess (R)')



%% G

% figure
% fsurf(subs(G,[T],[temp1]),[0 1 0 0.5])
% title (['Gibbs Free Energy @ ' num2str(temp1) ' K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' G (eV/BSF)')

%% -dG/ds = muO

% figure
% fsurf(-subs(diff(G,y),[T],[temp1]),[0 1 0 0.5])
% title (['\partialG/\partial\delta = \mu^O @ ' num2str(temp1) ' K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' \mu^O')

%% S

% figure
% fsurf(subs(S/r,[T],[temp1]),[0 1 0 0.5])
% title (['Entropy @ ' num2str(temp1) ' K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' S (R)')

%% dS/ds

figure
fsurf(diff(subs(S/r,[T],[temp1]),y),[0 1 0.01 0.48])
title (['\partialS/\partial\delta @ ' num2str(temp1) ' K']);
xlabel('Ba mol Fraction (x)')
ylabel('Extent of Reduction (\delta)')
zlabel(' \partialS/\partial\delta (R)')

%% H where H = G + TS  where S = -dG/dT

% figure
% fsurf(subs(H,[T],[temp1]),[0 1 0 0.5])
% title (['Enthalpy @ ' num2str(temp1) ' K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' H (eV/BSF)')

%% Compare H to DFT data

% DFT data
dft_data = table2array(readtable('BSF_data_phase_shift.xlsx'));
dft_data_copy = dft_data;
rows = any(isnan(dft_data),2);
dft_data(rows,:) = [];

X_dft = dft_data(:,2); % mol fract Ba
Y_dft = dft_data(:,1); % delta
Z_dft = dft_data(:,3)/8; % E eV per ABO3

% figure
% hold on
% fsurf(subs(H,[T],[0]),[0 1 0 0.5])
% scatter3(X_dft,Y_dft,Z_dft,100,'or')
% title (['Enthalpy @ 0 K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' ~G @ 0 K (eV/ mol BSF)')

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
% legend('X = 0','X = 0.125', 'X = 0.5', 'X = 1.0', 'DFT Data')
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

x_frac = 0.999;

[Ho_1 So_1] = get_O2_thermo(400+273.15);
[Ho_2 So_2] = get_O2_thermo(575+273.15);
[Ho_3 So_3] = get_O2_thermo(750+273.15);
[Ho_4 So_4] = get_O2_thermo(925+273.15);
[Ho_5 So_5] = get_O2_thermo(1100+273.15);


H_red_400 = simplify(diff(subs(H_end,[x T],[x_frac (400+273.15)]),y)+ 0.5*Ho_1/98.4875);
H_red_575 = simplify(diff(subs(H_end,[x T],[x_frac (575+273.15)]),y)+ 0.5*Ho_2/98.4875);
H_red_750 = simplify(diff(subs(H_end,[x T],[x_frac (750+273.15)]),y)+ 0.5*Ho_3/98.4875);
H_red_925 = simplify(diff(subs(H_end,[x T],[x_frac (925+273.15)]),y)+ 0.5*Ho_4/98.4875);
H_red_1100 = simplify(diff(subs(H_end,[x T],[x_frac (1100+273.15)]),y)+ 0.5*Ho_5/98.4875);


figure
hold on
fs5 = fplot(H_red_1100*96.487*2,[0 0.5],'linewidth',2.0);
fs4 = fplot(H_red_925*96.487*2,[0 0.5],'linewidth',2.0);
fs3 = fplot(H_red_750*96.487*2,[0 0.5],'linewidth',2.0);
fs2 = fplot(H_red_575*96.487*2,[0 0.5],'linewidth',2.0);
fs1 = fplot(H_red_400*96.487*2,[0 0.5],'linewidth',2.0);

ylim([0 225])
%title (['Full Model SrFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
xlabel('Extent of Reduction (\delta)')
ylabel(' \partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1})')
legend('1100 ^\circC','925   ^\circC','750   ^\circC','575   ^\circC','400   ^\circC')
box on
hold off

%% look at T dependence fo dS/ds for x = 0

x_frac = 0.999;

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

ylim([0 500])
%title (['Full Model SrFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
xlabel('Extent of Reduction (\delta)')
ylabel(' \partialS/\partial\delta (J \bullet (mol O_{2} \bullet K)^{-1})')
%legend('1100 ^\circC','925   ^\circC','750   ^\circC','575   ^\circC','400   ^\circC')
box on
hold off


%% reduction energy compares

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
% 
% dH_0 = diff(subs(H,[x T],[0 0]),y);
% dH_0125 = diff(subs(H,[x T],[0.125 0]),y);
% dH_05 = diff(subs(H,[x T],[0.5 0]),y);
% dH_1 = diff(subs(H,[x T],[0.999 0]),y);
% 
% 
% double(int(dH_0,0,0.5));
% dH_ds0 = matlabFunction(dH_0);
% double(int(dH_05,0,0.25));
% dH_ds05 = matlabFunction(dH_05);
% 
% s = 0:0.01:0.5;
% 
% figure
% hold on
% %fplot(dH_0,[0 0.5],'linewidth',2.0);
% area(s,dH_ds0(s),'FaceColor','#8C1D40','Linewidth',4,'FaceAlpha',0.75)
% plot([0.0630 0.0630],[0 dH_ds0(0.0630)],'Color','#FFC627','Linewidth',3)
% plot([0.1250 0.1250],[0 dH_ds0(0.1250)],'Color','#FFC627','Linewidth',3)
% plot([0.1880 0.1880],[0 dH_ds0(0.1880)],'Color','#FFC627','Linewidth',3)
% plot([0.2500 0.2500],[0 dH_ds0(0.2500)],'Color','#FFC627','Linewidth',3)
% plot([0.5000 0.5000],[0 dH_ds0(0.5000)],'Color','#FFC627','Linewidth',3)
% title (['SrFeO_{3-\delta} @ 0 K']);
% xlabel('Extent of Reduction (\delta)')
% ylabel('H (eV \bullet (mol BSF)^{-1})')
% box on
% hold off
%%

%% reduction energy compares
% % DFT data
% dft_data = table2array(readtable('BSF_data_phase_shift.xlsx'));
% dft_data_copy = dft_data;
% rows = any(isnan(dft_data),2);
% dft_data(rows,:) = [];
% 
% X_dft = dft_data(:,2); % mol fract Ba
% Y_dft = dft_data(:,1); % delta
% Z_dft = dft_data(:,3)/8; % E eV per ABO3
% 
% 
% test_T = 0;
% 
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
% H_0 = subs(H,[x T],[0 test_T]);
% dH_0 = diff(subs(H,[x T],[0 test_T]),y);
% dH_0125 = diff(subs(H,[x T],[0.125 test_T]),y);
% dH_05 = diff(subs(H,[x T],[0.5 test_T]),y);
% dH_1 = diff(subs(H,[x T],[0.999 test_T]),y);
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
% scatter(0.0630,double(int(dH_0,0,0.0630)),75,'d','MarkerFaceColor','#FFC627')
% scatter(0.1250,double(int(dH_0,0,0.1250)),75,'d','MarkerFaceColor','#FFC627')
% scatter(0.1880,double(int(dH_0,0,0.1880)),75,'d','MarkerFaceColor','#FFC627')
% scatter(0.2500,double(int(dH_0,0,0.2500)),75,'d','MarkerFaceColor','#FFC627')
% scatter(0.5000,double(int(dH_0,0,0.5000)),75,'d','MarkerFaceColor','#FFC627')
% 
% plot([0.0630 0.5],[double(int(dH_0,0,0.0630)) double(int(dH_0,0,0.0630))],'--','Color','#FFC627','Linewidth',2)
% plot([0.1250 0.5],[double(int(dH_0,0,0.1250)) double(int(dH_0,0,0.1250))],'--','Color','#FFC627','Linewidth',2)
% plot([0.1880 0.5],[double(int(dH_0,0,0.1880)) double(int(dH_0,0,0.1880))],'--','Color','#FFC627','Linewidth',2)
% plot([0.2500 0.5],[double(int(dH_0,0,0.2500)) double(int(dH_0,0,0.2500))],'--','Color','#FFC627','Linewidth',2)
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


%% dH/ds end only x = 0
% 
% x_frac = 0;
% 
% [Ho_1 So_1] = get_O2_thermo(400+273.15);
% [Ho_2 So_2] = get_O2_thermo(575+273.15);
% [Ho_3 So_3] = get_O2_thermo(750+273.15);
% [Ho_4 So_4] = get_O2_thermo(925+273.15);
% [Ho_5 So_5] = get_O2_thermo(1100+273.15);
% 
% 
% H_red_400 = simplify(diff(subs(H_end,[x T],[x_frac (400+273.15)]),y)+ 0.5*Ho_1);
% H_red_575 = simplify(diff(subs(H_end,[x T],[x_frac (575+273.15)]),y)+ 0.5*Ho_2);
% H_red_750 = simplify(diff(subs(H_end,[x T],[x_frac (750+273.15)]),y)+ 0.5*Ho_3);
% H_red_925 = simplify(diff(subs(H_end,[x T],[x_frac (925+273.15)]),y)+ 0.5*Ho_4);
% H_red_1100 = simplify(diff(subs(H_end,[x T],[x_frac (1100+273.15)]),y)+ 0.5*Ho_5);
% 
% 
% figure
% hold on
% fs1 = fplot(H_red_400*96.487*2,[0 0.5],'linewidth',2.0);
% fs2 = fplot(H_red_575*96.487*2,[0 0.5],'linewidth',2.0);
% fs3 = fplot(H_red_750*96.487*2,[0 0.5],'linewidth',2.0);
% fs4 = fplot(H_red_925*96.487*2,[0 0.5],'linewidth',2.0);
% fs5 = fplot(H_red_1100*96.487*2,[0 0.5],'linewidth',2.0);
% 
% ylim([0 250])
% title (['End-Member Only SrFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
% xlabel('Extent of Reduction (\delta)')
% ylabel(' \partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1})')
% legend('400 C', '575 C', '750 C', '925 C','1100 C')
% box on
% hold off
% 
% %% dS/ds gend only x = 0
% 
% x_frac = 0.00000000000000000000000000001;
% 
% S_red_400 = simplify(diff(subs(S_end,[x T],[x_frac (400+273.15)]),y) + 0.5*So_1);
% S_red_575 = simplify(diff(subs(S_end,[x T],[x_frac (575+273.15)]),y) + 0.5*So_2);
% S_red_750 = simplify(diff(subs(S_end,[x T],[x_frac (750+273.15)]),y) + 0.5*So_3);
% S_red_925 = simplify(diff(subs(S_end,[x T],[x_frac (925+273.15)]),y) + 0.5*So_4);
% S_red_1100 = simplify(diff(subs(S_end,[x T],[x_frac (1100+273.15)]),y) + 0.5*So_5);
% 
% 
% figure
% hold on
% fs1 = fplot(S_red_400*96.487*2*1000,[0 0.5],'linewidth',2.0); %+205.152
% fs2 = fplot(S_red_575*96.487*2*1000,[0 0.5],'linewidth',2.0);
% fs3 = fplot(S_red_750*96.487*2*1000,[0 0.5],'linewidth',2.0);
% fs4 = fplot(S_red_925*96.487*2*1000,[0 0.5],'linewidth',2.0);
% fs5 = fplot(S_red_1100*96.487*2*1000,[0 0.5],'linewidth',2.0);
% 
% ylim([0 250])
% title (['End-Member Only SrFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
% xlabel('Extent of Reduction (\delta)')
% ylabel(' \partialS/\partial\delta (J \bullet (mol O_{2} \bullet K)^{-1})')
% legend('400 C', '575 C', '750 C', '925 C','1100 C')
% box on
% hold off
% 
% %% dH/ds end only x = 0
% 
% x_frac = 0;
% 
% [Ho_1 So_1] = get_O2_thermo(400+273.15);
% [Ho_2 So_2] = get_O2_thermo(575+273.15);
% [Ho_3 So_3] = get_O2_thermo(750+273.15);
% [Ho_4 So_4] = get_O2_thermo(925+273.15);
% [Ho_5 So_5] = get_O2_thermo(1100+273.15);
% 
% 
% H_red_400 = simplify(diff(subs(H_excess,[x T],[x_frac (400+273.15)]),y)+ 0.5*Ho_1);
% H_red_575 = simplify(diff(subs(H_excess,[x T],[x_frac (575+273.15)]),y)+ 0.5*Ho_2);
% H_red_750 = simplify(diff(subs(H_excess,[x T],[x_frac (750+273.15)]),y)+ 0.5*Ho_3);
% H_red_925 = simplify(diff(subs(H_excess,[x T],[x_frac (925+273.15)]),y)+ 0.5*Ho_4);
% H_red_1100 = simplify(diff(subs(H_excess,[x T],[x_frac (1100+273.15)]),y)+ 0.5*Ho_5);
% 
% 
% figure
% hold on
% fs1 = fplot(H_red_400*96.487*2,[0 0.5],'linewidth',2.0);
% fs2 = fplot(H_red_575*96.487*2,[0 0.5],'linewidth',2.0);
% fs3 = fplot(H_red_750*96.487*2,[0 0.5],'linewidth',2.0);
% fs4 = fplot(H_red_925*96.487*2,[0 0.5],'linewidth',2.0);
% fs5 = fplot(H_red_1100*96.487*2,[0 0.5],'linewidth',2.0);
% 
% ylim([0 250])
% title (['Excess Only SrFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
% xlabel('Extent of Reduction (\delta)')
% ylabel(' \partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1})')
% legend('400 C', '575 C', '750 C', '925 C','1100 C')
% box on
% hold off
% 
% %% dS/ds gend only x = 0
% 
% x_frac = 0.00000000000000000000000000001;
% 
% S_red_400 = simplify(diff(subs(S_excess,[x T],[x_frac (400+273.15)]),y) + 0.5*So_1);
% S_red_575 = simplify(diff(subs(S_excess,[x T],[x_frac (575+273.15)]),y) + 0.5*So_2);
% S_red_750 = simplify(diff(subs(S_excess,[x T],[x_frac (750+273.15)]),y) + 0.5*So_3);
% S_red_925 = simplify(diff(subs(S_excess,[x T],[x_frac (925+273.15)]),y) + 0.5*So_4);
% S_red_1100 = simplify(diff(subs(S_excess,[x T],[x_frac (1100+273.15)]),y) + 0.5*So_5);
% 
% 
% figure
% hold on
% fs1 = fplot(S_red_400*96.487*2*1000,[0 0.5],'linewidth',2.0); %+205.152
% fs2 = fplot(S_red_575*96.487*2*1000,[0 0.5],'linewidth',2.0);
% fs3 = fplot(S_red_750*96.487*2*1000,[0 0.5],'linewidth',2.0);
% fs4 = fplot(S_red_925*96.487*2*1000,[0 0.5],'linewidth',2.0);
% fs5 = fplot(S_red_1100*96.487*2*1000,[0 0.5],'linewidth',2.0);
% 
% ylim([0 250])
% title (['Excess Only SrFeO_{3-\delta} :  ' num2str(fit_col) ' L Term(s) fit']);
% xlabel('Extent of Reduction (\delta)')
% ylabel(' \partialS/\partial\delta (J \bullet (mol O_{2} \bullet K)^{-1})')
% legend('400 C', '575 C', '750 C', '925 C','1100 C')
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



