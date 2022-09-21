clear;
clc;

font_s = 24; % text size for plots
%% Load Calculated Errors workspace

error_file = "model_guess0001_get_errors";

load(error_file) % This is all_fits

%% Normal Plot of muh_o absolute error 

pd = fitdist(transpose(my_muh_error),'Normal');
mu = pd.mu;
sigma = pd.sigma;
y_norm = normpdf(sort(my_muh_error),mu,sigma);
plot(sort(my_muh_error),y_norm,'linewidth',2.0)

title('Normal Plot','FontSize',font_s)
xlabel('Abs Error in Chemical Potential (eV)','FontSize',font_s)
ylabel('Probability Density','FontSize',font_s)

%% Normal Plot absolute error of dd

pd = fitdist((my_dd_error),'Normal');
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


%% Normal Plot true error of dd

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

%% Normal Plots for All combinatiosn of excess terms

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


%% d0 trends plot during minimzation

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