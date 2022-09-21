clear;
clc;
%% Pull in Data


d0_dirs = ["../Logic_Cross/model_guess0001_d0s.mat",...
    "../Logic_Corners/model_guess0001_d0s.mat",...
    "../Logic_No_DFT/model_guess0001_d0s.mat",...
    "../No_Logic_cross/model_guess0001_d0s.mat",...
    "../No_Logic_No_DFT/model_guess0001_d0s.mat"];

MAE_dirs = ["../Logic_Cross/model_guess0001_MAEs.mat",...
    "../Logic_Corners/model_guess0001_MAEs.mat",...
    "../Logic_No_DFT/model_guess0001_MAEs.mat",...
    "../No_Logic_cross/model_guess0001_MAEs.mat",...
    "../No_Logic_No_DFT/model_guess0001_MAEs.mat"];


fit_col = 3;
% d0s

d0s = [];
for x = 1:length(d0_dirs)
    load(d0_dirs(x));
    d0s = [d0s d0_progress(:,end-(fit_col)+1)]; %the d0_progress varible comes from the load
end

% MAEs and STDV
maes = [];
stdv = [];
for x = 1:length(MAE_dirs)
    load(MAE_dirs(x));
    maes = [maes; MAEs(1,:)]; %the d0_progress varible comes from the load
    stdv = [stdv; MAEs(2,:)];
end

%% Plot MAEs

font_s = 24; % text size for plots

Xaxes = cell(22,1);
Xaxes{1} = '1 Terms';
Xaxes{2} = '2 Terms';
Xaxes{3} = '3 Terms';
Xaxes{4} = '4 Terms';
Xaxes{5} = '5 Terms';
Xaxes{6} = '6 Terms';
Xaxes{7} = '7 Terms';
Xaxes{8} = '8 Terms';
Xaxes{9} = '9 Terms';
Xaxes{10} = '10 Terms';
Xaxes{11} = '11 Terms';
Xaxes{12} = '12 Terms';
Xaxes{13} = '13 Terms';
Xaxes{14} = '14 Terms';
Xaxes{15} = '15 Terms';
Xaxes{16} = '16 Terms';
Xaxes{17} = '17 Terms';
Xaxes{18} = '18 Terms';
Xaxes{19} = '19 Terms';
Xaxes{20} = '20 Terms';
Xaxes{21} = '21 Terms';
Xaxes{22} = '22 Terms';

Legend = cell(5,1);
Legend{1} = 'Cross Fit';
Legend{2} = 'Corners Only';
Legend{3} = 'Experimental Only';
Legend{4} = 'No Logcis Cross';
Legend{5} = 'No Logic No DFT';
% Legend{6} = 'Model 6';


figure()
plot([1:22],maes,'Linewidth',2)
%set(h, {'color'}, {[0 0 1]; [0 0 0]; [1 0 0]});
xticks([1:22])
xticklabels(flip(Xaxes))
legend(Legend)
title('MAE trend','FontSize',font_s)
xlabel('# of Terms','FontSize',font_s)
ylabel('MAE in \delta','FontSize',font_s)

%% Plot delta predictions at T,X,P

res = 100;
Temps = linspace(573.15,1573.15,res);
mol_fracs = linspace(0,1,res);  
Preses = linspace(0.01,1,res);

load('CrossFit_eqn.mat');
load('CornersFit_eqn.mat');
load('NoDFTFit_eqn.mat');
load('NoLogicCrossFit_eqn.mat');
load('NoLogicNoDFTFit_eqn.mat');

%%
% 
% 
% final_Fulleqn = [crossfit_3term_model, Corners_anchored_only, no_dft_model];
% 
% 
% syms y x T
% tic
% M = length(Temps);
% %my_x = [0,0.05,0.1,0.15,0.2];
% dd_data = [];
% parfor m = 1:3 % number of fits
%     count = 1;
%     v = [];
%     for n = 1:M % loop through temps
%         for j = 1:M % loop throuhg preses
%             for k = 1:M % loop through mol_fracs
%                 subbed_dG_soln = matlabFunction(simplify(expand(subs(diff(final_Fulleqn(m),y),[T x],[Temps(n) mol_fracs(k)]))))
%                 % solve using fzero
%                  eqn = @(y) -(subbed_dG_soln(y)) - get_mu_o(Temps(n),Preses(j));
% %                 dd_pred = fzero(eqn,.25);
%                 % solve using vpasolve
% %                 vpa_eqn = subs(diff(final_Fulleqn(m),y),[T x],[Temps(n) mol_fracs(k)])
% %                 dd_pred = vpasolve(vpa_eqn == get_mu_o(Temps(n),Preses(j)),y);
%                 % solve using fsolve
%                 dd_pred = real(fsolve(eqn,0.02,optimset('FunValCheck', 'off', 'Display', 'off')));
%                 v(count,:) = [mol_fracs(k), Temps(n), Preses(j), dd_pred]; 
%                 count = count + 1;
%             end
%         end
%     end
%     dd_data(:,:,m) = v;
% end
% toc
% %dd_data

%%

[xx,yy,zz] = meshgrid(linspace(min(Temps),max(Temps),res),linspace(min(Preses),max(Preses),res),linspace(min(mol_fracs),max(mol_fracs),res));
final_Fulleqn = [CrossFit, CornersFit, NoDFTFit];


syms y x T
tic
%my_x = [0,0.05,0.1,0.15,0.2];
dd_data = [];
for m = 1:3 % number of fits
    v = zeros(numel(xx),4);
    parfor n = 1:numel(xx) % loop through temps
        subbed_dG_soln = matlabFunction(simplify(expand(subs(diff(final_Fulleqn(m),y),[T x],[xx(n) zz(n)]))))
        eqn = @(y) -(subbed_dG_soln(y)) - get_mu_o(xx(n),yy(n));
        dd_pred = real(fsolve(eqn,0.02,optimset('FunValCheck', 'off', 'Display', 'off')));
        v(n,:) = [zz(n), xx(n), yy(n), dd_pred]; 
    end
    dd_data(:,:,m) = v;
end
toc
%% muh_o check

get_mu_o(300+273.15,0.9)
get_mu_o(800+273.15,1.0)

%% Plotting predictions

model = 1;
figure()
scatter3(dd_data(:,1,model),dd_data(:,4,model),dd_data(:,2,model),100,dd_data(:,3,model));
cb = colorbar;  
box on
xlabel('Ba Fraction')
ylabel('Predicted \delta')
zlabel('Temperature (K)')
colormap('parula')
cb.Label.String = 'P_{O2} (Bar)';



%% compare to exp

model = 1;

dref = d0s(1,model);
% read in data and create arrays
exp_data = readtable('Transformed_exp_data.xlsx');

Press = table2array(exp_data(:,2));  % partial pressure O2
Temp = table2array(exp_data(:,3))+273.15;  % K
x_exp = table2array(exp_data(:,4));  % mol fraciton Ba
dd = table2array(exp_data(:,5));  % delta of delta

%comapre for x = 0
remove = (x_exp(:,1) ~= 0);
Press(remove)=[];
Temp(remove) = [];
x_exp(remove )= [];
dd(remove) = [];


exp_data_x_0 = [x_exp, Temp, Press, dd+dref];

check_data = dd_data(:,:,model);
%% To compare to x = 0
remove = (check_data(:,1) ~= 0);
check_data(remove,:) = [];

figure()
hold on
scatter3(check_data(1:end,3),check_data(1:end,2),check_data(1:end,4),75,'ok');
scatter3(exp_data_x_0(:,3),exp_data_x_0(:,2),exp_data_x_0(:,4),75,'or');
box on
grid on
xlabel('P_{O2} (Bar)')
zlabel('\delta')
ylabel('Temperature (K)')
legend('Predicted','Experiemntal')
hold off


%%
Y = check_data(:,3);
Z = check_data(:,4);
X = check_data(:,2);


% Create 100x100 grid mesh (x,y) points
[xGrid,yGrid] = meshgrid(linspace(min(X),max(X)),linspace(min(Y),max(Y)));
% Interpolation
zGrid = griddata(X(:),Y(:),Z(:),xGrid(:),yGrid(:),'cubic');
zGrid = reshape(zGrid,size(xGrid));

figure
hold on
surf(xGrid,yGrid,zGrid,'FaceAlpha',0.75,'EdgeColor','none')%,'FaceColor','g',)
scatter3(exp_data_x_0(:,2),exp_data_x_0(:,3),exp_data_x_0(:,4),75,'ok');
cb = colorbar; 
colormap jet
box on
grid on
ylabel('P_{O2} (Bar)')
zlabel('\delta')
xlabel('Temperature (K)')
legend('Predicted','Experiemntal')

%% Comparing 3 models

model = 1;

dref = d0s(1,model);
% read in data and create arrays
exp_data = readtable('Transformed_exp_data.xlsx');

Press = table2array(exp_data(:,2));  % partial pressure O2
Temp = table2array(exp_data(:,3))+273.15;  % K
x_exp = table2array(exp_data(:,4));  % mol fraciton Ba
dd = table2array(exp_data(:,5));  % delta of delta

%To Remove specific pieces of data
remove = (x_exp(:,1) > 0);
Press(remove)=[];
Temp(remove) = [];
x_exp(remove )= [];
dd(remove) = [];


exp_data_x_0 = [x_exp, Temp, Press, dd];

check_data = dd_data(:,:,model);
remove = (check_data(:,1) > 0);
check_data(remove,:) = [];


Y1 = check_data(:,3);
Z1 = check_data(:,4) - dref;
X1 = check_data(:,2);


model = 2;

dref = d0s(1,model);
% read in data and create arrays
exp_data = readtable('Transformed_exp_data.xlsx');

Press = table2array(exp_data(:,2));  % partial pressure O2
Temp = table2array(exp_data(:,3))+273.15;  % K
x_exp = table2array(exp_data(:,4));  % mol fraciton Ba
dd = table2array(exp_data(:,5));  % delta of delta

%To Remove specific pieces of data
remove = (x_exp(:,1) > 0);
Press(remove)=[];
Temp(remove) = [];
x_exp(remove )= [];
dd(remove) = [];


exp_data_x_0 = [x_exp, Temp, Press, dd];

check_data = dd_data(:,:,model);
remove = (check_data(:,1) > 0);
check_data(remove,:) = [];

Y2 = check_data(:,3);
Z2 = check_data(:,4) - dref;
X2 = check_data(:,2);

model = 3;

dref = d0s(1,model);
% read in data and create arrays
exp_data = readtable('Transformed_exp_data.xlsx');

Press = table2array(exp_data(:,2));  % partial pressure O2
Temp = table2array(exp_data(:,3))+273.15;  % K
x_exp = table2array(exp_data(:,4));  % mol fraciton Ba
dd = table2array(exp_data(:,5));  % delta of delta

%To Remove specific pieces of data
remove = (x_exp(:,1) > 0);
Press(remove)=[];
Temp(remove) = [];
x_exp(remove )= [];
dd(remove) = [];


exp_data_x_0 = [x_exp, Temp, Press, dd];

check_data = dd_data(:,:,model);
remove = (check_data(:,1) > 0);
check_data(remove,:) = [];

Y3 = check_data(:,3);
Z3 = check_data(:,4) - dref;
X3 = check_data(:,2);


% Create 100x100 grid mesh (x,y) points
[xGrid1,yGrid1] = meshgrid(linspace(min(X1),max(X1)),linspace(min(Y1),max(Y1)));
% Interpolation
zGrid1 = griddata(X1(:),Y1(:),Z1(:),xGrid1(:),yGrid1(:),'cubic');
zGrid1 = reshape(zGrid1,size(xGrid1));

% Create 100x100 grid mesh (x,y) points
[xGrid2,yGrid2] = meshgrid(linspace(min(X2),max(X2)),linspace(min(Y2),max(Y2)));
% Interpolation
zGrid2 = griddata(X1(:),Y1(:),Z2(:),xGrid1(:),yGrid1(:),'cubic');
zGrid2 = reshape(zGrid2,size(xGrid2));

% Create 100x100 grid mesh (x,y) points
[xGrid3,yGrid3] = meshgrid(linspace(min(X3),max(X3)),linspace(min(Y3),max(Y3)));
% Interpolation
zGrid3 = griddata(X1(:),Y1(:),Z3(:),xGrid1(:),yGrid1(:),'cubic');%% Equations from Gibbs

syms T
G = simplify(final_Fulleqn);
% S = -dG/dT
S = simplify(-(diff(G,T)));
% H where H = G + TS  
H = simplify(G + T*S);
zGrid3 = reshape(zGrid3,size(xGrid3));

figure
hold on
surf(xGrid1,yGrid1,zGrid1,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','g')
surf(xGrid1,yGrid1,zGrid2,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','r')
surf(xGrid1,yGrid1,zGrid3,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','b')
scatter3(exp_data_x_0(:,2),exp_data_x_0(:,3),exp_data_x_0(:,4),75,'ok');
box on
grid on
title ('SrFeO_{3-\delta} (x = 0)');
ylabel('P_{O2} (Bar)')
zlabel('\Delta\delta')
xlabel('Temperature (K)')
%legend('Predicted \Delta\delta (Cross Fit)','Predicted \Delta\delta (End-Mem Anchors Only)','Predicted \Delta\delta (Experimental Only)','Observed \Delta\delta')
legend('Predicted \Delta\delta (Cross Fit)','Observed \Delta\delta')
%% this dependes on models loaded above

% compareDdplots(xGrid1,yGrid1,zGrid1,zGrid2,exp_data_x_0(:,2),exp_data_x_0(:,3),exp_data_x_0(:,4),75,'ok')

%%

figure
hold on
surf(xGrid1,yGrid1,zGrid1 - zGrid2,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','g')
surf(xGrid1,yGrid1,zGrid1 - zGrid3,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','r')
surf(xGrid1,yGrid1,zGrid2 - zGrid3,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','b')

box on
grid on
ylabel('P_{O2} (Bar)')
zlabel('\delta')
xlabel('Temperature (K)')
legend('Cross Fit - End-Mem Only Anchored ','Cross Fit - Experimental Only ','End-Mem Only Anchored - Experimental Only')

%%
%% Equations from Gibbs for each fit

syms T
G_cross = simplify(CrossFit);
G_corners = simplify(CornersFit);
G_nodft = simplify(NoDFTFit);
G_nologcross = simplify(NoLogicCrossFit);
G_nolognodft = simplify(NoLogicNoDFTFit);
% S = -dG/dT
S_cross = simplify(-(diff(G_cross,T)));
S_corners = simplify(-(diff(G_corners,T)));
S_nodft = simplify(-(diff(G_nodft,T)));
S_nologcross = simplify(-(diff(G_nologcross,T)));
S_nolognodft = simplify(-(diff(G_nolognodft,T)));
% H where H = G + TS  
H_cross = simplify(G_cross + T*S_cross);
H_corners = simplify(G_corners + T*S_corners);
H_nodft = simplify(G_nodft + T*S_nodft);
H_nologcross = simplify(G_nologcross + T*S_nologcross);
H_nolognodft = simplify(G_nolognodft + T*S_nolognodft);

%%  Compare dH/ds and dS/ds for x = 0 + oxygen enthalpy

x_frac = 0;

[Ho So] = get_O2_thermo(800+273.15);

H_red_cross = simplify(diff(subs(H_cross,[x T],[x_frac (800+273.15)]),y)+ 0.5*Ho/98.4875);
H_red_corners = simplify(diff(subs(H_corners,[x T],[x_frac (800+273.15)]),y)+ 0.5*Ho/98.4875);
H_red_nodft = simplify(diff(subs(H_nodft,[x T],[x_frac (800+273.15)]),y)+ 0.5*Ho/98.4875);
H_red_nologcross = simplify(diff(subs(H_nologcross,[x T],[x_frac (800+273.15)]),y)+ 0.5*Ho/98.4875);
H_red_nolognodft = simplify(diff(subs(H_nolognodft,[x T],[x_frac (800+273.15)]),y)+ 0.5*Ho/98.4875);

figure
hold on
fplot(H_red_cross*96.487*2,[0 0.5],'-k','linewidth',2.0);
fplot(H_red_corners*96.487*2,[0 0.5],'-r','linewidth',2.0);
fplot(H_red_nodft*96.487*2,[0 0.5],'-b','linewidth',2.0);
fplot(H_red_nologcross*96.487*2,[0 0.5],'--k','linewidth',2.0);
fplot(H_red_nolognodft*96.487*2,[0 0.5],'--b','linewidth',2.0);


ylim([0 250])
title ('SrFeO_{3-\delta} @ 800 ^\circC');
xlabel('Extent of Reduction (\delta)')
ylabel(' \partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1})')
legend('Logic Cross Fit','Logic Corners Fit','Logic No DFT Fit','No Logic Cross Fit','No Logic No DFT Fit')
box on
hold off
%%
x_frac = 1e-6;

S_red_cross = simplify(diff(subs(S_cross,[x T],[x_frac (800+273.15)]),y) + 0.5*So/98.4875);
S_red_corners = simplify(diff(subs(S_corners,[x T],[x_frac (800+273.15)]),y) + 0.5*So/98.4875);
S_red_nodft = simplify(diff(subs(S_nodft,[x T],[x_frac (800+273.15)]),y) + 0.5*So/98.4875);
S_red_nologcross = simplify(diff(subs(S_nologcross,[x T],[x_frac (800+273.15)]),y) + 0.5*So/98.4875);
S_red_nolognodft = simplify(diff(subs(S_nolognodft,[x T],[x_frac (800+273.15)]),y) + 0.5*So/98.4875);

figure
hold on

fplot(S_red_cross*96.487*2*1000,[0 0.5],'-k','linewidth',2.0);
fplot(S_red_corners*96.487*2*1000,[0 0.5],'-r','linewidth',2.0);
fplot(S_red_nodft*96.487*2*1000,[0 0.5],'-b','linewidth',2.0);
fplot(S_red_nologcross*96.487*2*1000,[0 0.5],':k','linewidth',2.0);
fplot(S_red_nolognodft*96.487*2*1000,[0 0.5],':b','linewidth',2.0);


ylim([0 250])
title ('SrFeO_{3-\delta} @ 800 ^\circC');
xlabel('Extent of Reduction (\delta)')
ylabel(' \partialS/\partial\delta (J \bullet (mol O_{2} \bullet K)^{-1})')
legend('Logic Cross Fit','Logic Corners Fit','Logic No DFT Fit','No Logic Cross Fit','No Logic No DFT Fit')
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


test_T = 800+273.15;

test_T2 = 0;

%x = 0
line0 = dft_data(dft_data(:,2)==0,:);
(line0(end,3) - line0(1,3))/8;
(line0(end-1,3) - line0(1,3))/8;
(line0(end-2,3) - line0(1,3))/8;
(line0(end-3,3) - line0(1,3))/8;
(line0(end-4,3) - line0(1,3))/8;
(line0(end-5,3) - line0(1,3))/8;
%x = 0.125
line0125 = dft_data(dft_data(:,2)==0.125,:);
%x = 0.5
line05 = dft_data(dft_data(:,2)==0.5,:);
% x = 1.0
line1 = dft_data(dft_data(:,2)==1,:);
(line1(end,3) - line1(1,3))/8;


H_0_cross = subs(H_cross,[x T],[0 test_T]); 
H_0_corners = subs(H_corners,[x T],[0 test_T]);
H_0_no_dft = subs(H_nodft,[x T],[0 test_T]);

dH_0_cross = diff(subs(H_cross,[x T],[0 test_T]),y); 
dH_0_corners = diff(subs(H_corners,[x T],[0 test_T]),y);
dH_0_no_dft = diff(subs(H_nodft,[x T],[0 test_T]),y);

dH_ds0_cross = matlabFunction(dH_0_cross);
matH_0_cross = matlabFunction(H_0_cross);
dH_ds0_corners = matlabFunction(dH_0_corners);
matH_0_corners = matlabFunction(H_0_corners);
dH_ds0_no_dft = matlabFunction(dH_0_no_dft);
matH_0_no_dft = matlabFunction(H_0_no_dft);

H_0_cross2 = subs(H_cross,[x T],[0 test_T2]); 
H_0_corners2 = subs(H_corners,[x T],[0 test_T2]);
H_0_no_dft2 = subs(H_nodft,[x T],[0 test_T2]);

dH_0_cross2 = diff(subs(H_cross,[x T],[0 test_T2]),y); 
dH_0_corners2 = diff(subs(H_corners,[x T],[0 test_T2]),y);
dH_0_no_dft2 = diff(subs(H_nodft,[x T],[0 test_T2]),y);

dH_ds0_cross2 = matlabFunction(dH_0_cross2);
matH_0_cross2 = matlabFunction(H_0_cross2);
dH_ds0_corners2 = matlabFunction(dH_0_corners2);
matH_0_corners2 = matlabFunction(H_0_corners2);
dH_ds0_no_dft2 = matlabFunction(dH_0_no_dft2);
matH_0_no_dft2 = matlabFunction(H_0_no_dft2);

s = 0:0.01:0.5;

figure
hold on
fplot(H_0_cross-subs(H_0_cross,y,0),[0 0.5],'-k','Linewidth',2)
fplot(H_0_corners-subs(H_0_corners,y,0),[0 0.5],'-r','Linewidth',2)
fplot(H_0_no_dft-subs(H_0_no_dft,y,0),[0 0.5],'-b','Linewidth',2)

plotX = [0.0630,0.1250,0.1880,0.5000];
plotY = [(line0(end-4,3) - line0(1,3))/8,(line0(end-3,3) - line0(1,3))/8,...
    (line0(end-2,3) - line0(1,3))/8,(line0(end  ,3) - line0(1,3))/8];

scatter(plotX,plotY,100,[140/255 29/255 64/255],'d')

fplot(H_0_cross2-subs(H_0_cross2,y,0),[0 0.5],'-.k','Linewidth',2)
fplot(H_0_corners2-subs(H_0_corners2,y,0),[0 0.5],'--r','Linewidth',2)
fplot(H_0_no_dft2-subs(H_0_no_dft2,y,0),[0 0.5],'--b','Linewidth',2)

ylim([0,0.5])
title (['T = 800 C^\circ (solid)', ' T = 0 K (dashed)' ]);
xlabel('Extent of Reduction (\delta)')
ylabel('H - H(\delta=0) (eV \bullet (mol BSF)^{-1})')
legend('Cross Fit Model','Corners Only Model','No DFT Model','DFT Data')
box on
hold off

%% comapre delta H reductio to lit
[Ho So] = mu_o_thermo(298);
[H1 S1] = mu_o_thermo(873);
[H2 S2] = mu_o_thermo(973);
[H3 S3] = mu_o_thermo(1073);
[H4 S4] = mu_o_thermo(1173);
[H5 S5] = mu_o_thermo(1273);

myH_cross = matlabFunction(subs(H_cross,x,0));

fsurf(subs(H_cross,x,0)*96.487,[800 1300 0 0.5])

dHo = (myH_cross(298,0))*96.487;
dH1 = (myH_cross(873,0.25))*96.487;
dH2 = (myH_cross(973,0.3))*96.487;
dH3 = (myH_cross(1073,0.325))*96.487;
dH4 = (myH_cross(1173,0.36))*96.487;
dH5 = (myH_cross(1273,0.41))*96.487;

red1 = dH1 - dHo
red2 = dH2 - dHo
red3 = dH3 - dHo
red4 = dH4 - dHo
red5 = dH5 - dHo


%%  Compare dH/ds and dS/ds across x + oxygen enthalpy
% 
% x_frac1 = 0;
% x_frac2 = 0.05;
% x_frac3 = 0.1;
% x_frac4 = 0.15;
% x_frac5 = 0.2;
% x_frac6 = 0.5;
% 
% 
% [Ho, So] = get_O2_thermo(800+273.15);
% 
% H_red_cross1 = simplify(diff(subs(H_cross,[x T],[x_frac1 (800+273.15)]),y)+ 0.5*Ho/98.4875);
% H_red_cross2 = simplify(diff(subs(H_cross,[x T],[x_frac2 (800+273.15)]),y)+ 0.5*Ho/98.4875);
% H_red_cross3 = simplify(diff(subs(H_cross,[x T],[x_frac3 (800+273.15)]),y)+ 0.5*Ho/98.4875);
% H_red_cross4 = simplify(diff(subs(H_cross,[x T],[x_frac4 (800+273.15)]),y)+ 0.5*Ho/98.4875);
% H_red_cross5 = simplify(diff(subs(H_cross,[x T],[x_frac5 (800+273.15)]),y)+ 0.5*Ho/98.4875);
% H_red_cross6 = simplify(diff(subs(H_cross,[x T],[x_frac6 (800+273.15)]),y)+ 0.5*Ho/98.4875);
% %H_red_cross7 = simplify(diff(subs(H_cross,[x T],[x_frac7 (800+273.15)]),y)+ 0.5*Ho/98.4875);
% 
% 
% figure
% hold on
% fplot(H_red_cross1*96.487*2,[0 0.5],'linewidth',2.0);
% fplot(H_red_cross2*96.487*2,[0 0.5],'linewidth',2.0);
% fplot(H_red_cross3*96.487*2,[0 0.5],'linewidth',2.0);
% fplot(H_red_cross4*96.487*2,[0 0.5],'linewidth',2.0);
% fplot(H_red_cross5*96.487*2,[0 0.5],'linewidth',2.0);
% fplot(H_red_cross6*96.487*2,[0 0.5],'linewidth',2.0);
% 
% 
% ylim([0 250])
% title ('Ba_{x}Sr_{1-x}FeO_{3-\delta} @ 800 ^\circC');
% xlabel('Extent of Reduction (\delta)')
% ylabel(' \partialH/\partial\delta (kJ \bullet (mol O_{2})^{-1})')
% legend('x=0','x=0.05','x=0.10','x=0.15','x=0.2','x=0.5')
% box on
% hold off
% 
% %%  Compare dH/ds and dS/ds across x + oxygen enthalpy
% 
% x_frac1 = 0.00000000000000000000000000000000000000000000000000001;
% x_frac2 = 0.05;
% x_frac3 = 0.1;
% x_frac4 = 0.15;
% x_frac5 = 0.2;
% x_frac6 = 0.5;
% 
% 
% [Ho, So] = get_O2_thermo(800+273.15);
% 
% S_red_cross1 = simplify(diff(subs(S_cross,[x T],[x_frac1 (800+273.15)]),y)+ 0.5*So/98.4875);
% S_red_cross2 = simplify(diff(subs(S_cross,[x T],[x_frac2 (800+273.15)]),y)+ 0.5*So/98.4875);
% S_red_cross3 = simplify(diff(subs(S_cross,[x T],[x_frac3 (800+273.15)]),y)+ 0.5*So/98.4875);
% S_red_cross4 = simplify(diff(subs(S_cross,[x T],[x_frac4 (800+273.15)]),y)+ 0.5*So/98.4875);
% S_red_cross5 = simplify(diff(subs(S_cross,[x T],[x_frac5 (800+273.15)]),y)+ 0.5*So/98.4875);
% S_red_cross6 = simplify(diff(subs(S_cross,[x T],[x_frac6 (800+273.15)]),y)+ 0.5*So/98.4875);
% S_red_cross7 = simplify(diff(subs(S_cross,[x T],[x_frac7 (800+273.15)]),y)+ 0.5*So/98.4875);
% 
% 
% figure
% hold on
% fplot(S_red_cross1*96.487*2*1000,[0 0.5],'linewidth',2.0);
% fplot(S_red_cross2*96.487*2*1000,[0 0.5],'linewidth',2.0);
% fplot(S_red_cross3*96.487*2*1000,[0 0.5],'linewidth',2.0);
% fplot(S_red_cross4*96.487*2*1000,[0 0.5],'linewidth',2.0);
% fplot(S_red_cross5*96.487*2*1000,[0 0.5],'linewidth',2.0);
% fplot(S_red_cross6*96.487*2*1000,[0 0.5],'linewidth',2.0);
% 
% 
% ylim([0 250])
% title ('Ba_{x}Sr_{1-x}FeO_{3-\delta} @ 800 ^\circC');
% xlabel('Extent of Reduction (\delta)')
% ylabel(' \partialS/\partial\delta (kJ \bullet (mol O_{2})^{-1})')
% legend('x=0','x=0.05','x=0.10','x=0.15','x=0.2','x=0.5')
% box on
% hold off


%%
save("res_100_compare");

