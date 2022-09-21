clear;
clc;

%% Setting Up

% Which Model are we loading? 
file_save_name = "model_0001";

load(file_save_name + "_allfits.mat") % This is all_fits
load(file_save_name + "_d0s.mat") % this is d0_progress

% Which # of Excess terms are we intrested in?
fit_col = 3; % 3 = 3 excess terms in G_soln

% Constants
r = 8.617333262E-5; % eV/K;

%% Load Data
EXP_file = "Transformed_exp_data.xlsx";
DFT_file = "BSF_data_phase_shift.xlsx";

[X,Y,Z,x_exp,Temp,Press,dd,muhg_o] = Load_Data(DFT_file,EXP_file);

%% Generate Symbolic CEF equtions
syms x y T real
syms G1 G2 G3 G4
syms L [1 24] real
syms gDiffA [1 2] real
syms gDiffB [1 2] real
syms gDiffC [1 2] real
syms A [1 24] real
syms B [1 24] real
syms C [1 24] real

%Generate base CEF
[G_soln,G_ex,G_o] = CEF();

% Expand Gibbs endmember terms to the for G(T) = A + B*T + C*T*ln(T)
[G_sol,G_dft,Go,Go_dft] = Temp_Expansion_endmembers(G_soln,G_o);

% Expand Gibbs excess terms to the for G(T) = A + B*T + C*T*ln(T)
[Gsol,Gdft,Gex] = Temp_Expansion_excess(G_sol,G_dft,G_ex);

%% Since we have DFT data we can "localize" the the enthalpic energy
% Note we still allow correction to this - explained in paper

[Gsol,Go,Gdft,Go_dft] = Localize_H(Gsol,Go,Gdft,Go_dft);

%% Assing Parameters based on number of excess terms

d0s = d0_progress(:,end-(fit_col)+1);

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


final_Fulleqn = simplify(expand(subs(Gsol,[sub_me],[sub_value])));
final_G_ex = simplify(expand(subs(Gex,[sub_me(1:66)],[sub_value(1:66)])));
final_G_end = simplify(expand(subs(Go,[sub_me(67:72)],[sub_value(67:72)])));

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

%% Temperature for most plots 

temp1 = 800 + 273.15; % K

%% NOTE: These commented out plots are good to have, but not necessary for quick evaluation

% %% Entropy comapres between G_end and G_excess
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
% zlabel(' S_{end} + S_{excess} (R)')
% 
% 
% %% G
% 
% figure
% fsurf(subs(G,[T],[temp1]),[0 1 0 0.5])
% title (['Gibbs Free Energy @ ' num2str(temp1) ' K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' G (eV/BSF)')
% 
% %% -dG/ds = muO
% 
% figure
% fsurf(-subs(diff(G,y),[T],[temp1]),[0 1 0 0.5])
% title (['\partialG/\partial\delta = \mu^O @ ' num2str(temp1) ' K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' \mu^O')
% 
% %% S constant T
% 
% figure
% fsurf(subs(S/r,[T],[temp1]),[0 1 0 0.5])
% title (['Entropy @ ' num2str(temp1) ' K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' S (R)')
% 
% %% S constant x
% 
% figure
% fsurf(subs(S/r,[x],[1e-6]),[100 1500 0 0.5])
% title (['Entropy @  x = 0']);
% xlabel('Temp (K)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' S (R)')
% 
% %% dS/ds
% 
% figure
% fsurf(diff(subs(S/r,[T],[temp1]),y),[0 1 0.01 0.48])
% title (['\partialS/\partial\delta @ ' num2str(temp1) ' K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' \partialS/\partial\delta (R)')
% 
% %% dS/ds @ x = 0
% 
% figure
% fsurf(diff(subs(S/r,[x],[1e-6]),y),[100 1500 0.02 0.48])
% title (['\partialS/\partial\delta @ x = 0']);
% xlabel('Temp (K)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' \partialS/\partial\delta (R)')
% %% H where H = G + TS  where S = -dG/dT
% 
% figure
% fsurf(subs(H,[T],[temp1]),[0 1 0 0.5])
% title (['Enthalpy @ ' num2str(temp1) ' K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' H (eV/BSF)')
% 
% %% Compare H to DFT data
% 
% figure
% hold on
% fsurf(subs(H,[T],[0]),[0 1 0 0.5])
% scatter3(X,Y,Z,100,'or')
% title (['Enthalpy @ 0 K']);
% xlabel('Ba mol Fraction (x)')
% ylabel('Extent of Reduction (\delta)')
% zlabel(' ~G @ 0 K (eV/ mol BSF)')

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
%cb.Limits = [80, 200];

%legend('1100 ^\circC','925   ^\circC','750   ^\circC','575   ^\circC','400   ^\circC')
box on
hold off


%% look at T dependence for dS/ds for x = 0

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

%legend('1100 ^\circC','925   ^\circC','750   ^\circC','575   ^\circC','400   ^\circC')
box on
hold off
%% reduction energy compares

% Need DFT data in table format
dft_data = table2array(readtable(DFT_file));
dft_data_copy = dft_data;
rows = any(isnan(dft_data),2);
dft_data(rows,:) = [];

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
