%function get_errors_fun(file_save_name)
%% How man terms we want
clear ;
clc;

file_save_name = "model_guess0001_no_intial_fits";

load(file_save_name + "_allfits.mat") % This is all_fits
load(file_save_name + "_d0s.mat") % this is d0_progress

fit_col = 3;

font_s = 24; % text size for plots

%% read in data and create arrays
exp_data = readtable('Transformed_exp_data.xlsx');

Press = table2array(exp_data(:,2));  % partial pressure O2
Temp = table2array(exp_data(:,3))+273.15;  % K
x = table2array(exp_data(:,4));  % mol fraciton Ba
dd = table2array(exp_data(:,5));  % delta of delta
x_exp = x;

% use remove logical on each array
% remove = (Temp(:) == 811.15);
% 
% Press(remove)=[];
% Temp(remove) = [];
% x(remove )= [];
% x_exp = x;
% dd(remove) = [];

%% DFT data
dft_data = table2array(readtable('BSF_data_phase_shift.xlsx'));
dft_data_copy = dft_data;
rows = any(isnan(dft_data),2);
dft_data(rows,:) = [];

X = dft_data(:,2); % mol fract Ba
Y = dft_data(:,1); % delta
Z = dft_data(:,3)/8; % E eV per ABO3


%% Grids for ploting

[xx_exp,yy_exp] = meshgrid(linspace(min(x),max(x)),linspace(min(dd),max(dd)));
[xx_DFT,yy_DFT] = meshgrid(linspace(min(X),max(X)),linspace(min(Y),max(Y)));



%% calcualte Z -1/RT dGsolid/dd = 1/2RT GO2
% O2 chemical potential 

Ho = [];
So = [];
for t = 1:length(Temp);
    [H,S] = get_O2_thermo(Temp(t));
    Ho(t) = H;
    So(t) = S;
end
conv = 96.487;  % 1 eV = 96.487 kJ/mol
Po = 21278.25; %Pa
Patm = 21278.25; %Pa
R_gas = 8.314462/1000/conv; %eV/K
Ho = transpose(Ho);
So = transpose(So);

muhg_o2 = Ho/conv - (Temp).*So/conv + R_gas.*(Temp).*log(Press);

muhg_o = 0.5*muhg_o2;


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


%% Expand L terms in excess to be G(T) terms
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


%% Table of terms
% final_tab = array2table(round(all_fits(:,end-(fit_col)+1),5), ...
%       'VariableNames',{'Coeff'});
% 
% final_tab.Properties.RowNames = ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9",...
%     "A10", "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19",...
%     "A20", "A22", "A24", "B1", "B2", "B3", "B4", "B5", "B6",...
%     "B7", "B8", "B9", "B10", "B11", "B12", "B13", "B14", "B15", "B16", "B17",...
%     "B18", "B19", "B20", "B22", "B24", "C1", "C2", "C3", "C4", "C5", "C6",...
%     "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17",...
%     "C18", "C19", "C20", "C22", "C24","gDiffA1","gDiffA2","gDiffB1",...
%     "gDiffB2","gDiffC1","gDiffC2"]
% 
% %% Table of all fits
% all_tab = array2table(round(all_fits,5));
% 
% all_tab.Properties.RowNames = ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9",...
%     "A10", "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19",...
%     "A20", "A22", "A24", "B1", "B2", "B3", "B4", "B5", "B6",...
%     "B7", "B8", "B9", "B10", "B11", "B12", "B13", "B14", "B15", "B16", "B17",...
%     "B18", "B19", "B20", "B22", "B24", "C1", "C2", "C3", "C4", "C5", "C6",...
%     "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17",...
%     "C18", "C19", "C20", "C22", "C24","gDiffA1","gDiffA2","gDiffB1",...
%     "gDiffB2","gDiffC1","gDiffC2"]
%% plot of new fits with dropped terms

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
    
G_soln_test = subs(sym(G_soln),[sub_me],[sub_value]); % eV @ 300 C
dG_soln_test = diff(G_soln_test,y);

G_soln_test = matlabFunction(G_soln_test);
dG_soln_test = matlabFunction(dG_soln_test);


% %az = 90;
% %el = 0;
% figure()
% hold on
% color = ['r','g','b','y'];
% count = 1;
% for Tspot = 400:400:1600
%     surf(xx_exp,yy_exp,-dG_subs(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
%                     A15,A16,A17,A18,A19,A20,A22,A24,B1,B2,B3,B4,B5,B6,...
%                     B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,B18,B19,B20,B22,...
%                     B24,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,...
%                     C17,C18,C19,C20,C22,C24,Tspot,dref,gDiffA1,gDiffA2,gDiffB1,...
%                     gDiffB2,gDiffC1,gDiffC2,xx_exp,yy_exp),'FaceColor',color(count))
%     count = count + 1;
% end
% ylabel(['\delta + ', num2str(dref)],'FontSize',24)
% xlabel('x','FontSize',24)
% zlabel('\mu^{solid}_o','FontSize',24)
% Legend=cell(4,1);%  two positions 
% Legend{1}=' T = 400 K' ;
% Legend{2}=' T = 800 K'; 
% Legend{3}=' T = 1200 K';
% Legend{4}=' T = 1600 K';
% legend(Legend);
% title('Dropped Terms')
% hold off
% 
% 
% figure
% hold on
% surf(xx_DFT,yy_DFT,G_subs(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
%                     A15,A16,A17,A18,A19,A20,A22,A24,dref,gDiffA1,gDiffA2,xx_DFT,yy_DFT))
% scatter3(X,Y,Z,'r')
% ylabel('\delta','FontSize',24)
% xlabel('x','FontSize',24)
% zlabel('E (0 K)','FontSize',24)
% title('Dropped Terms')
% hold off


%% Predcting

% dG_subs_pred = subs(sym(dG_dy),[L21 L23],[L9 L11]);
% 
% T_dG_subs_pred = subs(dG_subs_pred,[L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13 L14 L15 L16 L17 L18 L19 L20 L22 L24],...
%     [(A1 + B1*T + C1*T*log(T)) (A2 + B2*T + C2*T*log(T)) (A3 + B3*T + C3*T*log(T)) ...
%     (A4 + B4*T + C4*T*log(T)) (A5 + B5*T + C5*T*log(T)) (A6 + B6*T + C6*T*log(T)) (A7 + B7*T + C7*T*log(T)) (A8 + B8*T + C8*T*log(T)) ...
%     (A9 + B9*T + C9*T*log(T)) (A10 + B10*T + C10*T*log(T)) (A11 + B11*T + C11*T*log(T)) ...
%     (A12 + B12*T + C12*T*log(T)) (A13 + B13*T + C13*T*log(T)) (A14 + B14*T + C14*T*log(T)) ...
%     (A15 + B15*T + C15*T*log(T)) (A16 + B16*T + C16*T*log(T)) (A17 + B17*T + C17*T*log(T)) (A18 + B18*T + C18*T*log(T)) ...
%     (A19 + B19*T + C19*T*log(T)) (A20 + B20*T + C20*T*log(T)) ...
%     (A22 + B22*T + C22*T*log(T)) (A24 + B24*T + C24*T*log(T))]);
% 
% T_dG_subs_pred = matlabFunction(T_dG_subs_pred);
% 
% [xx_pre,yy_pre] = meshgrid(linspace(0,1.0),linspace(0,0.5));
% 
% P_ratio = 0.9;

% for T_pre = 400:200:1600
%     figure
%     surf(xx_pre,yy_pre,-T_dG_subs_pred(T_pre,gDiffA1,gDiffA2,gDiffB1,...
%                     gDiffB2,gDiffC1,gDiffC2,...
%         xx_pre,yy_pre),'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none')
%     hold on
%     surf(xx_pre,yy_pre,ones(size(xx_pre))*get_mu_o(T_pre,P_ratio),'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
%     
% %     surf(xx_exp,yy_exp,-dG_subs(A1,A2,A3,A4,A5,A6,A9,A10,A11,A12,A13,A14,...
% %                     A15,A16,B1,B2,B3,B4,B5,B6,B9,B10,B11,B12,B13,B14,B15,...
% %                     B16,C1,C2,C3,C4,C5,C6,C9,C10,C11,C12,C13,C14,C15,C16,...
% %                     Tspot,dref,gB1,gB2,gB3,gB4,gC1,gC2,gC3,gC4,xx_exp,yy_exp),'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none')
% %     hold on
% %     surf(xx_exp,yy_exp,ones(size(xx_exp))*get_mu_o(T_pre,P_ratio),'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
%     title("Ba_{x}Sr_{1-x}FeO_{3-\delta} @ T = "+T_pre+ " K; P = "+P_ratio+" Atm",'FontSize',font_s)
%     ylabel('\delta','FontSize',font_s)
%     xlabel('x','FontSize',font_s)
%     zlabel('\mu_{o}','FontSize',font_s)
%     legend("\mu^{s}", "\mu^{g}",'FontSize',font_s)
%     hold off
% end

%% Checking the abs Error of muhg_o

my_muh_error = zeros(1,length(dd));
for n = [1:length(x_exp)] 
    x_vals = [0,0.05,0.1,0.15,0.2];
    for k = 1:length(x_vals)
        if x_exp(n) == x_vals(k)
            dref = 0.02;
            break
        end
    end
    muhg_pred = -dG_soln_test(Temp(n),x_exp(n),dd(n)+dref);
    my_muh_error(n) = abs(muhg_pred - muhg_o(n));
end
my_muh_err_avg = mean(my_muh_error);
muh_err_std = std(my_muh_error);
fprintf('This is the error for calualting chemical potenital given T, X, and delta: %f +- %f eV (%f +- %f kJ/mol).\n',my_muh_err_avg,muh_err_std,my_muh_err_avg*96.487,muh_err_std*96.487); 

pd = fitdist(transpose(my_muh_error),'Normal');
mu = pd.mu;
sigma = pd.sigma;
y_norm = normpdf(sort(my_muh_error),mu,sigma);
plot(sort(my_muh_error),y_norm,'linewidth',2.0)

title('Normal Plot','FontSize',font_s)
xlabel('Abs Error in Chemical Potential (eV)','FontSize',font_s)
ylabel('Probability Density','FontSize',font_s)



%% Checking the abs Error by solving for dd

my_dd_error = zeros(1,length(dd));
m = length(x_exp);
heat_mp = [];
parfor n = [1:m]
    x_vals = [0,0.05,0.1,0.15,0.2];
    for k = 1:length(x_vals)
        if x_exp(n) == x_vals(k)
            dref = 0.02;
            break
        end
    end
    subbed_dG_soln = matlabFunction(subs(sym(dG_soln_test),[T x],[Temp(n) x_exp(n)]));
    eqn = @(y) -(subbed_dG_soln(y)) - get_mu_o(Temp(n),Press(n));
    dd_pred = fzero(eqn,dd(n)+dref);
    my_dd_error(n) = abs(dd_pred - (dd(n)+dref));
    my_true_dd_error(n) = (dd_pred - (dd(n)+dref));
    heat_mp = [heat_mp ; [x_exp(n) Temp(n) Press(n) my_true_dd_error(n)]];
end

%% Normal Plot absolute error

pd = fitdist(transpose(my_dd_error),'Normal');
mu = pd.mu;
sigma = pd.sigma;
y_norm = normpdf(sort(my_dd_error),mu,sigma);
fprintf('This is the mean absolute error for predicting delta given T, X, and P: %f +- %f.\n',mu,sigma);
figure
plot(sort(my_dd_error),y_norm,'linewidth',2.0)

title('Normal Plot of Absolute Error','FontSize',font_s)
xlabel('Error in \Delta\delta','FontSize',font_s)
ylabel('Probability Density','FontSize',font_s)
%legend('Normal Distribution','FontSize',font_s)


%% Normal Plot true error

pd_true = fitdist(transpose(my_true_dd_error),'Normal');
mu_true = pd_true.mu;
sigma_true = pd_true.sigma;
y_norm_true = normpdf(sort(my_true_dd_error),mu_true,sigma_true);
fprintf('This is the mean true error for predicting delta given T, X, and P: %f +- %f.\n',mu_true,sigma_true); 
figure
plot(sort(my_true_dd_error),y_norm_true,'linewidth',2.0)

title('Normal Plot of True Error','FontSize',font_s)
xlabel('Error in \Delta\delta','FontSize',font_s)
ylabel('Probability Density','FontSize',font_s)
%legend('Normal Distribution','FontSize',font_s)

%% Heat Map of dd Error

ht_x = heat_mp(:,1);  % x
ht_y = heat_mp(:,2);  % temp
ht_z = heat_mp(:,3);  % Press
cc = heat_mp(:,4);    % Error

figure
scatter3(ht_x,ht_y,cc,40,ht_z)
cb = colorbar;  
xlabel('Ba Fraction')
ylabel('Temperature (K)')
zlabel('Error in \Delta\delta')
cb.Label.String = 'Partial Pressure O_2 (Bar)';


%% Getting all MAEs

MAEs = zeros(2,22);
my_dd_error = zeros(22,length(x_exp));
my_true_dd_error = zeros(22,length(x_exp));
y_norm = zeros(22,length(x_exp));
for spt = 1:22
    fit_col = spt;
    
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

    G_soln_sub = subs(sym(G_soln),[sub_me],[sub_value]); 
    dG_soln_sub = simplify(expand(diff(G_soln_sub,y)));

    G_soln_sub = matlabFunction(G_soln_sub);
    dG_soln_sub = matlabFunction(dG_soln_sub);
    
    %my_dd_error = zeros(1,length(dd));
    m = length(x_exp);
    parfor n = [1:m]
        x_vals = [0,0.05,0.1,0.15,0.2];
        for k = 1:length(x_vals)
            if x_exp(n) == x_vals(k)
                dref = 0.02;
                break
            end
        end
        subbed_dG_soln = matlabFunction(subs(sym(dG_soln_sub),[T x],[Temp(n) x_exp(n)]))
        eqn = @(y) -(subbed_dG_soln(y)) - get_mu_o(Temp(n),Press(n));
        dd_pred = fzero(eqn,dd(n)+dref);
        my_dd_error(fit_col,n) = abs(dd_pred - (dd(n)+dref));
        my_true_dd_error(fit_col,n) = (dd_pred - (dd(n)+dref));
    end
    pd = fitdist(transpose(my_dd_error(fit_col,:)),'Normal');
    mu = pd.mu;
    sigma = pd.sigma;
    fprintf('The mean absolute error for %d term(s) predicting delta is %f +- %f.\n',spt,mu,sigma);
    y_norm(fit_col,:) = normpdf(sort(my_dd_error(fit_col,:)),mu,sigma);
    MAEs(1,23-spt) = mu;
    MAEs(2,23-spt) = sigma;
    
end
%% Normal Plots for All

Legend = cell(22,1);
Legend{1} = '1 Terms';
Legend{2} = '2 Terms';
Legend{3} = '3 Terms';
Legend{4} = '4 Terms';
Legend{5} = '5 Terms';
Legend{6} = '6 Terms';
Legend{7} = '7 Terms';
Legend{8} = '8 Terms';
Legend{9} = '9 Terms';
Legend{10} = '10 Terms';
Legend{11} = '11 Terms';
Legend{12} = '12 Terms';
Legend{13} = '13 Terms';
Legend{14} = '14 Terms';
Legend{15} = '15 Terms';
Legend{16} = '16 Terms';
Legend{17} = '17 Terms';
Legend{18} = '18 Terms';
Legend{19} = '19 Terms';
Legend{20} = '20 Terms';
Legend{21} = '21 Terms';
Legend{22} = '22 Terms';

L = 1;
N = 1;
M = 22;

figure()
hold on
for plt = L:N:M
    plot(sort(my_dd_error(plt,:)),y_norm(plt,:),'linewidth',2.0)
end
for plt = L:N:M
    plot([mean(my_dd_error(plt,:)) mean(my_dd_error(plt,:))]...
        ,[0 max(y_norm(plt,:))],'k:','linewidth',2.0)
end
title('Normal Plot of Absolute Error','FontSize',font_s)
xlabel('Error in \Delta\delta','FontSize',font_s)
ylabel('Probability Density','FontSize',font_s)
legend(Legend{L:N:M})
hold off

figure()
plot(MAEs(1,:))
xticks([1:22])
xticklabels(flip(Legend))
ylim([0.002 0.01])
title('MAE trend','FontSize',font_s)
xlabel('# of Terms','FontSize',font_s)
ylabel('MAE in \delta','FontSize',font_s)


%% d0 trends plot

Legend = cell(22,1);
Legend{1} = '1 Terms';
Legend{2} = '2 Terms';
Legend{3} = '3 Terms';
Legend{4} = '4 Terms';
Legend{5} = '5 Terms';
Legend{6} = '6 Terms';
Legend{7} = '7 Terms';
Legend{8} = '8 Terms';
Legend{9} = '9 Terms';
Legend{10} = '10 Terms';
Legend{11} = '11 Terms';
Legend{12} = '12 Terms';
Legend{13} = '13 Terms';
Legend{14} = '14 Terms';
Legend{15} = '15 Terms';
Legend{16} = '16 Terms';
Legend{17} = '17 Terms';
Legend{18} = '18 Terms';
Legend{19} = '19 Terms';
Legend{20} = '20 Terms';
Legend{21} = '21 Terms';
Legend{22} = '22 Terms';

L = 1;
N = 1;
M = 10;
d0_progress
d0_progress_flip = flip(d0_progress,2)
x_vals = [0,0.05,0.1,0.15,0.2];

figure()
hold on
for plt = L:N:M
    plot(x_vals,d0_progress_flip(:,plt),'linewidth',2.0)
end

%title('Normal Plot of Absolute Error','FontSize',font_s)
xlabel('Ba mol fraciton (x)','FontSize',font_s)
ylabel('\delta^{\circ}','FontSize',font_s)
legend(Legend{L:N:M})
hold off
%%
save(file_save_name + '_get_errors')
save(file_save_name + '_MAEs','MAEs')
%end