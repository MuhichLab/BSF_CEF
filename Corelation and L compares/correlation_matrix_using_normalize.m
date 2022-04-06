 clear;
clc;

rng(42)

Ya_Sr = @(x) 1-x;
Ya_Ba = @(x) x;
Yb_3 = @(y) 2.*y;     
Yb_4 = @(y) 1-2.*y;   
Yo_o = @(y) 1-y./3;   
Yo_Va = @(y) y./3;    

Title = cell(26,1);
Title{1} = 'L1';
Title{2} = 'L2';
Title{3} = 'L3';
Title{4} = 'L4';
Title{5} = 'L5';
Title{6} = 'L6';
Title{7} = 'L7';
Title{8} = 'L8';
Title{9} = 'L9';
Title{10} = 'L10';
Title{11} = 'L11';
Title{12} = 'L12';
Title{13} = 'L13';
Title{14} = 'L14';
Title{15} = 'L15';
Title{16} = 'L16';
Title{17} = 'L17';
Title{18} = 'L18';
Title{19} = 'L19';
Title{20} = 'L20';
Title{21} = 'L21';
Title{22} = 'L22';
Title{23} = 'L23';
Title{24} = 'L24';
Title{25} = 'L9 + L21';
Title{26} = 'L11 + L23';


% exp_data = readtable('Transformed_exp_data.xlsx');
% data = [table2array(exp_data(:,4)) , table2array(exp_data(:,5))];
% ind=randperm(size(data,1));
% x_exp = table2array(exp_data(ind,4));  % mol fraciton Ba
% dd = table2array(exp_data(ind,5));  % delta of delta
% 
% % DFT data
% dft_data = table2array(readtable('BSF_data_phase_shift.xlsx'));
% dft_data_copy = dft_data;
% rows = any(isnan(dft_data),2);
% dft_data(rows,:) = [];
% ind2=randperm(size(dft_data,1));
% X = dft_data(ind2,2); % mol fract Ba
% d = dft_data(ind2,1); % delta

% X_dir = x_exp;
% d_dir = dd;



% X = rand(100,1);
% 
% a = 0;
% b = 0.5;
% d = (b-a).*rand(100,1) + a;
% 
% Xa = min(x);
% Xb = max(x);
% X_dir = (Xb-Xa).*rand(100,1) + Xa;
% 
% 
% dda = min(dd);
% ddb = max(dd);
% d_dir = (ddb-dda).*rand(100,1) + dda;


X = rand(1000,1);
%X = .25;

a = 0;
b = 0.5;
d = (b-a).*rand(1000,1) + a;

%d = 0.1;

X_dir = X;
d_dir = d;


%% Check plots

figure
hold on
scatter((Yo_o(d).*Yo_Va(d).*Ya_Sr(X).*Yb_3(d)),(Yo_o(d).*Yo_Va(d).*Ya_Sr(X).*Yb_3(d).*(Yo_o(d)-Yo_Va(d))))


%%
A = [(Ya_Sr(X).*Ya_Ba(X).*Yb_3(d).*Yo_o(d)) ...                      %L1
    (Ya_Sr(X).*Ya_Ba(X).*Yb_3(d).*Yo_o(d).*(Ya_Sr(X)-Ya_Ba(X))) ...  %L2
    (Ya_Sr(X).*Ya_Ba(X).*Yb_3(d).*Yo_Va(d)) ...                      %L3
    (Ya_Sr(X).*Ya_Ba(X).*Yb_3(d).*Yo_Va(d).*(Ya_Sr(X)-Ya_Ba(X)))...  %L4
    (Ya_Sr(X).*Ya_Ba(X).*Yb_4(d).*Yo_o(d)) ...                       %L5
    (Ya_Sr(X).*Ya_Ba(X).*Yb_4(d).*Yo_o(d).*(Ya_Sr(X)-Ya_Ba(X)))...   %L6
    (Ya_Sr(X).*Ya_Ba(X).*Yb_4(d).*Yo_Va(d)) ...                      %L7
    (Ya_Sr(X).*Ya_Ba(X).*Yb_4(d).*Yo_Va(d).*(Ya_Sr(X)-Ya_Ba(X)))...  %L8
    (Yb_3(d).*Yb_4(d).*Ya_Sr(X).*Yo_o(d)) ...                        %L9
    (Yb_3(d).*Yb_4(d).*Ya_Sr(X).*Yo_o(d).*(Yb_3(d)-Yb_4(d)))...      %L10
    (Yb_3(d).*Yb_4(d).*Ya_Ba(X).*Yo_o(d)) ...                        %L11
    (Yb_3(d).*Yb_4(d).*Ya_Ba(X).*Yo_o(d).*(Yb_3(d)-Yb_4(d)))...      %L12
    (Yb_3(d).*Yb_4(d).*Ya_Sr(X).*Yo_Va(d)) ...                       %L13
    (Yb_3(d).*Yb_4(d).*Ya_Sr(X).*Yo_Va(d).*(Yb_3(d)-Yb_4(d)))...     %L14
    (Yb_3(d).*Yb_4(d).*Ya_Ba(X).*Yo_Va(d)) ...                       %L15
    (Yb_3(d).*Yb_4(d).*Ya_Ba(X).*Yo_Va(d).*(Yb_3(d)-Yb_4(d)))...     %L16
    (Yo_o(d).*Yo_Va(d).*Ya_Sr(X).*Yb_3(d)) ...                       %L17
    (Yo_o(d).*Yo_Va(d).*Ya_Sr(X).*Yb_3(d).*(Yo_o(d)-Yo_Va(d)))...    %L18
    (Yo_o(d).*Yo_Va(d).*Ya_Ba(X).*Yb_3(d)) ...                       %L19
    (Yo_o(d).*Yo_Va(d).*Ya_Ba(X).*Yb_3(d).*(Yo_o(d)-Yo_Va(d)))...    %L20
    (Yo_o(d).*Yo_Va(d).*Ya_Sr(X).*Yb_4(d)) ...                       %L21
    (Yo_o(d).*Yo_Va(d).*Ya_Sr(X).*Yb_4(d).*(Yo_o(d)-Yo_Va(d)))...    %L22
    (Yo_o(d).*Yo_Va(d).*Ya_Ba(X).*Yb_4(d)) ...                       %L23
    (Yo_o(d).*Yo_Va(d).*Ya_Ba(X).*Yb_4(d).*(Yo_o(d)-Yo_Va(d)))];     %L24

names = {'L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11','L12'...
    ,'L13','L14','L15','L16','L17','L18','L19','L20','L21','L22','L23','L24'};

% [R,P] = corrcoef(A);
% % figure
% [Rfull,PValuefull] = corrplot(A,'varNames',names);%,'testR','on');


A = normalize(A,1);
%% regress instead

x1 = [X;X];
x2 = [d;d];
x = [ones(size(x1)) x1 x2 x1.*x2 x1.^2 x2.^2 x1.^2.*x2.^2 x1.^3 x2.^3 x2.^4 ...
    x1.*x2.^2 x1.*x2.^3 x1.*x2.^4 x1.^2.*x2 x1.^2.*x2.^3 x1.^2.*x2.^4 x1.^3.*x2 ...
    x1.^3.*x2.^2 x1.^3.*x2.^3 x1.^3.*x2.^4];
my_combos = combntns(1:size(A,2),2);
sts_dataA = zeros(size(my_combos,1),4);
for m = 1:size(my_combos,1)
    y = [A(:,my_combos(m,1));A(:,my_combos(m,2))];
    [b,~,~,~,stats] = regress(y,x);
    sts_dataA(m,:) = [my_combos(m,:), stats(1), stats(3)];
    
end


%%
n = 1;
m = 3;


y = [A(:,n);A(:,m)];
[b,~,~,~,stats] = regress(y,x);

figure
scatter3(x1(1:1000),x2(1:1000),y(1:1000),'r','filled')
hold on
scatter3(x1(1001:end),x2(1001:end),y(1001:end),'g','filled')
x1fit = min(x1):0.001:max(x1);
x2fit = min(x2):0.001:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT +...
    b(5)*X1FIT.^2 + b(6)*X2FIT.^2 + b(7)*X1FIT.^2.*X2FIT.^2 + b(8)*X1FIT.^3 +...
    b(9)*X2FIT.^3 + b(10)*X2FIT.^4 + b(11)*X1FIT.*X2FIT.^2 + ...
    b(12)*X1FIT.*X2FIT.^3 + b(13)*X1FIT.*X2FIT.^4 + b(14)*X1FIT.^2.*X2FIT + ...
    b(15)*X1FIT.^2.*X2FIT.^3 + b(16)*X1FIT.^2.*X2FIT.^4 + b(17)*X1FIT.^3.*X2FIT + ...
    b(18)*X1FIT.^3.*X2FIT.^2 + b(19)*X1FIT.^3.*X2FIT.^3 + b(20)*X1FIT.^3.*X2FIT.^4 ;
mesh(X1FIT,X2FIT,YFIT,'FaceColor','b','FaceAlpha',0.5, 'EdgeColor','none')
view(50,10)
hold off
legend(['L_',num2str(n)],['L_',num2str(m)],['Regression: ', newline, 'R^2 = ',num2str(round(stats(1),2))],'FontSize',20)
title([Title{n},' & ',Title{m},' Compare'])
xlabel('x','FontSize',20)
ylabel('\delta','FontSize',20)

%% L9 L21 and L11 L23 subs
B = [(Ya_Sr(X).*Ya_Ba(X).*Yb_3(d).*Yo_o(d)) ...                      %L1
    (Ya_Sr(X).*Ya_Ba(X).*Yb_3(d).*Yo_o(d).*(Ya_Sr(X)-Ya_Ba(X))) ...  %L2
    (Ya_Sr(X).*Ya_Ba(X).*Yb_3(d).*Yo_Va(d)) ...                      %L3
    (Ya_Sr(X).*Ya_Ba(X).*Yb_3(d).*Yo_Va(d).*(Ya_Sr(X)-Ya_Ba(X)))...  %L4
    (Ya_Sr(X).*Ya_Ba(X).*Yb_4(d).*Yo_o(d)) ...                       %L5
    (Ya_Sr(X).*Ya_Ba(X).*Yb_4(d).*Yo_o(d).*(Ya_Sr(X)-Ya_Ba(X)))...   %L6
    (Ya_Sr(X).*Ya_Ba(X).*Yb_4(d).*Yo_Va(d)) ...                      %L7
    (Ya_Sr(X).*Ya_Ba(X).*Yb_4(d).*Yo_Va(d).*(Ya_Sr(X)-Ya_Ba(X)))...  %L8
    (Yb_3(d).*Yb_4(d).*Ya_Sr(X).*Yo_o(d) + Yo_o(d).*Yo_Va(d).*Ya_Sr(X).*Yb_4(d)) ...  %L9 & L21
    (Yb_3(d).*Yb_4(d).*Ya_Sr(X).*Yo_o(d).*(Yb_3(d)-Yb_4(d)))...      %L10
    (Yb_3(d).*Yb_4(d).*Ya_Ba(X).*Yo_o(d) + Yo_o(d).*Yo_Va(d).*Ya_Ba(X).*Yb_4(d)) ...  %L11 & L23
    (Yb_3(d).*Yb_4(d).*Ya_Ba(X).*Yo_o(d).*(Yb_3(d)-Yb_4(d)))...      %L12
    (Yb_3(d).*Yb_4(d).*Ya_Sr(X).*Yo_Va(d)) ...                       %L13
    (Yb_3(d).*Yb_4(d).*Ya_Sr(X).*Yo_Va(d).*(Yb_3(d)-Yb_4(d)))...     %L14
    (Yb_3(d).*Yb_4(d).*Ya_Ba(X).*Yo_Va(d)) ...                       %L15
    (Yb_3(d).*Yb_4(d).*Ya_Ba(X).*Yo_Va(d).*(Yb_3(d)-Yb_4(d)))...     %L16
    (Yo_o(d).*Yo_Va(d).*Ya_Sr(X).*Yb_3(d)) ...                       %L17
    (Yo_o(d).*Yo_Va(d).*Ya_Sr(X).*Yb_3(d).*(Yo_o(d)-Yo_Va(d)))...    %L18
    (Yo_o(d).*Yo_Va(d).*Ya_Ba(X).*Yb_3(d)) ...                       %L19
    (Yo_o(d).*Yo_Va(d).*Ya_Ba(X).*Yb_3(d).*(Yo_o(d)-Yo_Va(d)))...    %L20
    (Yo_o(d).*Yo_Va(d).*Ya_Sr(X).*Yb_4(d).*(Yo_o(d)-Yo_Va(d)))...    %L22
    (Yo_o(d).*Yo_Va(d).*Ya_Ba(X).*Yb_4(d).*(Yo_o(d)-Yo_Va(d)))];     %L24

combo_names = {'L1','L2','L3','L4','L5','L6','L7','L8','L9_L21','L10','L11_L23','L12'...
    ,'L13','L14','L15','L16','L17','L18','L19','L20','L22','L24'};

B = normalize(B,1);

%%
x1 = [X;X];
x2 = [d;d];
x = [ones(size(x1)) x1 x2 x1.*x2 x1.^2 x2.^2 x1.^2.*x2.^2 x1.^3 x2.^3 x2.^4 ...
    x1.*x2.^2 x1.*x2.^3 x1.*x2.^4 x1.^2.*x2 x1.^2.*x2.^3 x1.^2.*x2.^4 x1.^3.*x2 ...
    x1.^3.*x2.^2 x1.^3.*x2.^3 x1.^3.*x2.^4];
my_combos = combntns(1:size(B,2),2);
sts_dataB = zeros(size(my_combos,1),4);
for m = 1:size(my_combos,1)
    y = [B(:,my_combos(m,1));B(:,my_combos(m,2))];
    [b,~,~,~,stats] = regress(y,x);
    sts_dataB(m,:) = [my_combos(m,:), stats(1), stats(3)];
    
end

%%
n = 1;
m = 3;


y = [B(:,n);B(:,m)];
[b,~,~,~,stats] = regress(y,x);

figure
scatter3(x1(1:1000),x2(1:1000),y(1:1000),'r','filled')
hold on
scatter3(x1(1001:end),x2(1001:end),y(1001:end),'g','filled')
x1fit = min(x1):0.001:max(x1);
x2fit = min(x2):0.001:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT +...
    b(5)*X1FIT.^2 + b(6)*X2FIT.^2 + b(7)*X1FIT.^2.*X2FIT.^2 + b(8)*X1FIT.^3 +...
    b(9)*X2FIT.^3 + b(10)*X2FIT.^4 + b(11)*X1FIT.*X2FIT.^2 + ...
    b(12)*X1FIT.*X2FIT.^3 + b(13)*X1FIT.*X2FIT.^4 + b(14)*X1FIT.^2.*X2FIT + ...
    b(15)*X1FIT.^2.*X2FIT.^3 + b(16)*X1FIT.^2.*X2FIT.^4 + b(17)*X1FIT.^3.*X2FIT + ...
    b(18)*X1FIT.^3.*X2FIT.^2 + b(19)*X1FIT.^3.*X2FIT.^3 + b(20)*X1FIT.^3.*X2FIT.^4 ;
mesh(X1FIT,X2FIT,YFIT,'FaceColor','b','FaceAlpha',0.5, 'EdgeColor','none')
view(50,10)
hold off
legend(['L_',num2str(n)],['L_',num2str(m)],['Regression: ', newline, 'R^2 = ',num2str(round(stats(1),2))],'FontSize',20)
title([Title{n},' & ',Title{m},' Compare'])
xlabel('x','FontSize',20)
ylabel('\delta','FontSize',20)

%% Non Derivative Plots
% fntsz = 18;
% 
% [xx,yy] = meshgrid((linspace(min(X),max(X))),(linspace(min(d),max(d))));
% 
% figure
% hold on
% surf(xx,yy,(Ya_Sr(xx).*Ya_Ba(xx).*Yb_4(yy).*Yo_o(yy).*(Ya_Sr(xx)-Ya_Ba(xx))),'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
% surf(xx,yy,(Ya_Sr(xx).*Ya_Ba(xx).*Yb_3(yy).*Yo_o(yy).*(Ya_Sr(xx)-Ya_Ba(xx)) + Ya_Sr(xx).*Ya_Ba(xx).*Yb_4(yy).*Yo_Va(yy).*(Ya_Sr(xx)-Ya_Ba(xx))),'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none')
% legend('L6','L2 & L8','FontSize',fntsz)
% title('(+) Correlation','FontSize',fntsz)
% xlabel('x','FontSize',fntsz)
% ylabel('\delta','FontSize',fntsz)


%% derivaitive of full term

syms x y

Ya_Sr = 1-x;
Ya_Ba =  x;
Yb_3 =  2.*y;     
Yb_4 =  1-2.*y;   
Yo_o =  1-y./3;   
Yo_Va = y./3; 

eqn1 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_3.*Yo_o),y));                     %L1
eqn2 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_3.*Yo_o.*(Ya_Sr-Ya_Ba)),y));  %L2
eqn3 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_3.*Yo_Va),y));                      %L3
eqn4 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_3.*Yo_Va.*(Ya_Sr-Ya_Ba)),y));  %L4
eqn5 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_4.*Yo_o),y));                        %L5
eqn6 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_4.*Yo_o.*(Ya_Sr-Ya_Ba)),y));   %L6
eqn7 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_4.*Yo_Va),y));                     %L7
eqn8 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_4.*Yo_Va.*(Ya_Sr-Ya_Ba)),y)); %L8
eqn9 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Sr.*Yo_o),y));                         %L9
eqn10 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Sr.*Yo_o.*(Yb_3-Yb_4)),y));      %L10
eqn11 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Ba.*Yo_o),y));                         %L11
eqn12 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Ba.*Yo_o.*(Yb_3-Yb_4)),y));      %L12
eqn13 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Sr.*Yo_Va),y));                        %L13
eqn14 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Sr.*Yo_Va.*(Yb_3-Yb_4)),y));     %L14
eqn15 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Ba.*Yo_Va),y));                        %L15
eqn16 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Ba.*Yo_Va.*(Yb_3-Yb_4)),y));     %L16
eqn17 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Sr.*Yb_3),y));                       %L17
eqn18 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Sr.*Yb_3.*(Yo_o-Yo_Va)),y));    %L18
eqn19 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Ba.*Yb_3),y));                     %L19
eqn20 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Ba.*Yb_3.*(Yo_o-Yo_Va)),y));    %L20
eqn21 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Sr.*Yb_4),y));                    %L21
eqn22 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Sr.*Yb_4.*(Yo_o-Yo_Va)),y));    %L22
eqn23 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Ba.*Yb_4),y));                    %L23
eqn24 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Ba.*Yb_4.*(Yo_o-Yo_Va)),y)); %L24

DA = [eqn1(X_dir,d_dir) eqn2(X_dir,d_dir) eqn3(X_dir,d_dir) eqn4(X_dir,d_dir)...
    eqn5(X_dir,d_dir) eqn6(X_dir,d_dir) eqn7(X_dir,d_dir) eqn8(X_dir,d_dir) ...
    eqn9(X_dir,d_dir) eqn10(X_dir,d_dir) eqn11(X_dir,d_dir) eqn12(X_dir,d_dir)...
    eqn13(X_dir,d_dir) eqn14(X_dir,d_dir) eqn15(X_dir,d_dir) eqn16(X_dir,d_dir) ...
    eqn17(X_dir,d_dir) eqn18(X_dir,d_dir) eqn19(X_dir,d_dir) eqn20(X_dir,d_dir)...
    eqn21(X_dir,d_dir) eqn22(X_dir,d_dir) eqn23(X_dir,d_dir) eqn24(X_dir,d_dir)];

DA = normalize(DA,1);
%%
% figure
% [RDfull,PValueDfull] = corrplot(DA,'varNames',names);%,'testR','on');

x1 = [X;X];
x2 = [d;d];
x = [ones(size(x1)) x1 x2 x1.*x2 x1.^2 x2.^2 x1.^2.*x2.^2 x1.^3 x2.^3 x2.^4 ...
    x1.*x2.^2 x1.*x2.^3 x1.*x2.^4 x1.^2.*x2 x1.^2.*x2.^3 x1.^2.*x2.^4 x1.^3.*x2 ...
    x1.^3.*x2.^2 x1.^3.*x2.^3 x1.^3.*x2.^4];
my_combos = combntns(1:size(DA,2),2);
sts_dataDA = zeros(size(my_combos,1),4);
for m = 1:size(my_combos,1)
    y = [DA(:,my_combos(m,1));DA(:,my_combos(m,2))];
    [b,~,~,~,stats] = regress(y,x);
    sts_dataDA(m,:) = [my_combos(m,:), stats(1), stats(3)];
    
end

%% L9 L21 and L11 L23 subs Derivatives
syms x y

Ya_Sr = 1-x;
Ya_Ba =  x;
Yb_3 =  2.*y;     
Yb_4 =  1-2.*y;   
Yo_o =  1-y./3;   
Yo_Va = y./3; 

eqn1 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_3.*Yo_o),y));                     %L1
eqn2 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_3.*Yo_o.*(Ya_Sr-Ya_Ba)),y));  %L2
eqn3 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_3.*Yo_Va),y));                      %L3
eqn4 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_3.*Yo_Va.*(Ya_Sr-Ya_Ba)),y));  %L4
eqn5 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_4.*Yo_o),y));                        %L5
eqn6 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_4.*Yo_o.*(Ya_Sr-Ya_Ba)),y));   %L6
eqn7 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_4.*Yo_Va),y));                     %L7
eqn8 = matlabFunction(diff((Ya_Sr.*Ya_Ba.*Yb_4.*Yo_Va.*(Ya_Sr-Ya_Ba)),y)); %L8
eqn9_21 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Sr.*Yo_o + Yo_o.*Yo_Va.*Ya_Sr.*Yb_4),y)); %L9 and L21
eqn10 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Sr.*Yo_o.*(Yb_3-Yb_4)),y));      %L10
eqn11_23 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Ba.*Yo_o + Yo_o.*Yo_Va.*Ya_Ba.*Yb_4),y)); %L11 and L23
eqn12 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Ba.*Yo_o.*(Yb_3-Yb_4)),y));      %L12
eqn13 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Sr.*Yo_Va),y));                        %L13
eqn14 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Sr.*Yo_Va.*(Yb_3-Yb_4)),y));     %L14
eqn15 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Ba.*Yo_Va),y));                        %L15
eqn16 = matlabFunction(diff((Yb_3.*Yb_4.*Ya_Ba.*Yo_Va.*(Yb_3-Yb_4)),y));     %L16
eqn17 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Sr.*Yb_3),y));                       %L17
eqn18 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Sr.*Yb_3.*(Yo_o-Yo_Va)),y));    %L18
eqn19 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Ba.*Yb_3),y));                     %L19
eqn20 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Ba.*Yb_3.*(Yo_o-Yo_Va)),y));    %L20
eqn22 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Sr.*Yb_4.*(Yo_o-Yo_Va)),y));    %L22
eqn24 = matlabFunction(diff((Yo_o.*Yo_Va.*Ya_Ba.*Yb_4.*(Yo_o-Yo_Va)),y)); %L24

DB = [eqn1(X_dir,d_dir) eqn2(X_dir,d_dir) eqn3(X_dir,d_dir) eqn4(X_dir,d_dir)...
    eqn5(X_dir,d_dir) eqn6(X_dir,d_dir) eqn7(X_dir,d_dir) eqn8(X_dir,d_dir) ...
    eqn9_21(X_dir,d_dir) eqn10(X_dir,d_dir) eqn11_23(X_dir,d_dir) eqn12(X_dir,d_dir)...
    eqn13(X_dir,d_dir) eqn14(X_dir,d_dir) eqn15(X_dir,d_dir) eqn16(X_dir,d_dir) ...
    eqn17(X_dir,d_dir) eqn18(X_dir,d_dir) eqn19(X_dir,d_dir) eqn20(X_dir,d_dir)...
    eqn22(X_dir,d_dir) eqn24(X_dir,d_dir)];

DB = normalize(DB,1);
%%
% figure
% [RDB,PValueDB,H] = corrplot(DB,'varNames',combo_names);
% haxes  = [H.Parent]; 
% for x = 1:length(haxes);
%     haxes(x).XTick = [];
%     haxes(x).YTick = [];
%     haxes(x).ZTick = [];
%     haxes(x).XTickLabel = [];
%     haxes(x).YTickLabel = [];
%     haxes(x).ZTickLabel = [];

x1 = [X;X];
x2 = [d;d];
x = [ones(size(x1)) x1 x2 x1.*x2 x1.^2 x2.^2 x1.^2.*x2.^2 x1.^3 x2.^3 x2.^4 ...
    x1.*x2.^2 x1.*x2.^3 x1.*x2.^4 x1.^2.*x2 x1.^2.*x2.^3 x1.^2.*x2.^4 x1.^3.*x2 ...
    x1.^3.*x2.^2 x1.^3.*x2.^3 x1.^3.*x2.^4];
my_combos = combntns(1:size(DB,2),2);
sts_dataDB = zeros(size(my_combos,1),4);
for m = 1:size(my_combos,1)
    y = [DB(:,my_combos(m,1));DB(:,my_combos(m,2))];
    [b,~,~,~,stats] = regress(y,x);
    sts_dataDB(m,:) = [my_combos(m,:), stats(1), stats(3)];
    
end
%%
n = 1;
m = 3;


y = [DB(:,n);DB(:,m)];
[b,~,~,~,stats] = regress(y,x);

figure
scatter3(x1(1:1000),x2(1:1000),y(1:1000),'r','filled')
hold on
scatter3(x1(1001:end),x2(1001:end),y(1001:end),'g','filled')
x1fit = min(x1):0.001:max(x1);
x2fit = min(x2):0.001:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT +...
    b(5)*X1FIT.^2 + b(6)*X2FIT.^2 + b(7)*X1FIT.^2.*X2FIT.^2 + b(8)*X1FIT.^3 +...
    b(9)*X2FIT.^3 + b(10)*X2FIT.^4 + b(11)*X1FIT.*X2FIT.^2 + ...
    b(12)*X1FIT.*X2FIT.^3 + b(13)*X1FIT.*X2FIT.^4 + b(14)*X1FIT.^2.*X2FIT + ...
    b(15)*X1FIT.^2.*X2FIT.^3 + b(16)*X1FIT.^2.*X2FIT.^4 + b(17)*X1FIT.^3.*X2FIT + ...
    b(18)*X1FIT.^3.*X2FIT.^2 + b(19)*X1FIT.^3.*X2FIT.^3 + b(20)*X1FIT.^3.*X2FIT.^4 ;
mesh(X1FIT,X2FIT,YFIT,'FaceColor','b','FaceAlpha',0.5, 'EdgeColor','none')
view(50,10)
hold off
legend(['L^\prime_',num2str(n)],['L^\prime_',num2str(m)],['Regression: ', newline, 'R^2 = ',num2str(round(stats(1),2))],'FontSize',20)
title([Title{n},' & ',Title{m},' Compare'])
xlabel('x','FontSize',20)
ylabel('\delta','FontSize',20)



%% Derivative Plots
% fntsz = 18;
% 
% syms x y
% 
% Ya_Sr = 1-x;
% Ya_Ba =  x;
% Yb_3 =  2.*y;     
% Yb_4 =  1-2.*y;   
% Yo_o =  1-y./3;   
% Yo_Va = y./3; 
% 
% [xx,yy] = meshgrid((linspace(min(X),max(X))),(linspace(min(d),max(d))));
% 
% figure
% hold on
% surf(xx,yy,term_6(xx,yy),'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
% surf(xx,yy,term_2_8(xx,yy),'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none')
% legend('L6','L2 & L8','FontSize',fntsz)
% title('Good (-) Correlation','FontSize',fntsz)
% xlabel('x','FontSize',fntsz)
% ylabel('\delta','FontSize',fntsz)
% 
% figure
% hold on
% plot(term_6(X,d),term_2_8(X,d),'ok')
% plot(term_6(X,d),term_10_22(X,d),'or')
% xlabel('L6','FontSize',fntsz)
% ylabel('Correlated L','FontSize',fntsz)
% legend('L2 & L8','L10 & L21','FontSize',fntsz)
% 
% figure
% hold on
% surf(xx,yy,term_6(xx,yy),'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
% surf(xx,yy,term_10_22(xx,yy),'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none')
% legend('L6','L10 & L21','FontSize',fntsz)
% title('Bad Correlation','FontSize',fntsz)
% xlabel('x','FontSize',fntsz)
% ylabel('\delta','FontSize',fntsz)
% 
% figure
% hold on
% surf(xx,yy,term_13_17(xx,yy),'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
% surf(xx,yy,term_1_7(xx,yy),'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none')
% legend('L13 & L17','L1 & L7','FontSize',fntsz)
% title('Good (+) Correlation','FontSize',fntsz)
% xlabel('x','FontSize',fntsz)
% ylabel('\delta','FontSize',fntsz)
% 
% figure
% hold on
% surf(xx,yy,term_13_17(xx,yy),'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
% surf(xx,yy,term_9_21(xx,yy),'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none')
% legend('L13 & L17','L9 & L21','FontSize',fntsz)
% title('Bad Correlation','FontSize',fntsz)
% xlabel('x','FontSize',fntsz)
% ylabel('\delta','FontSize',fntsz)

%% R and P compare


% compR = [];
% for i = 1:length(RDB)
%     for j = 1:length(RDB)
%         if (i ~= j) && (j < i)
%             compR = [compR ; [i j RDB(i,j) Rb(i,j)]];
%         end
%     end
% end
% compR = sortrows(compR,3,'descend','ComparisonMethod','abs');
% compP = [];
% for i = 1:length(PValueDB)
%     for j = 1:length(PValueDB)
%         if (i ~= j) && (j < i)
%             compP = [compP ; [i j PValueDB(i,j) PValueb(i,j)]];
%         end
%     end
% end
% compP = sortrows(compP,3,'ComparisonMethod','abs');
% %% R and P compare all terms
% 
% compR_all = [];
% for i = 1:length(RDfull)
%     for j = 1:length(RDfull)
%         if (i ~= j) && (j < i)
%             compR_all = [compR_all ; [i j RDfull(i,j) Rfull(i,j)]];
%         end
%     end
% end
% compR_all = sortrows(compR_all,3,'descend','ComparisonMethod','abs');
% compP_all = [];
% for i = 1:length(PValueDfull)
%     for j = 1:length(PValueDfull)
%         if (i ~= j) && (j < i)
%             compP_all = [compP_all ; [i j PValueDfull(i,j) PValuefull(i,j)]];
%         end
%     end
% end
% compP_all = sortrows(compP_all,3,'ComparisonMethod','abs');


%% compares

big_listA = [sts_dataA,sts_dataDA(:,3:4)];

big_listB = [sts_dataB,sts_dataDB(:,3:4)];



%% L Compares
n = 2;
m = 4;

Title = cell(26,1);
Title{1} = 'L1';
Title{2} = 'L2';
Title{3} = 'L3';
Title{4} = 'L4';
Title{5} = 'L5';
Title{6} = 'L6';
Title{7} = 'L7';
Title{8} = 'L8';
Title{9} = 'L9';
Title{10} = 'L10';
Title{11} = 'L11';
Title{12} = 'L12';
Title{13} = 'L13';
Title{14} = 'L14';
Title{15} = 'L15';
Title{16} = 'L16';
Title{17} = 'L17';
Title{18} = 'L18';
Title{19} = 'L19';
Title{20} = 'L20';
Title{21} = 'L21';
Title{22} = 'L22';
Title{23} = 'L23';
Title{24} = 'L24';
Title{25} = 'L9 + L21';
Title{26} = 'L11 + L23';

syms x y

Ya_Sr = 1-x;
Ya_Ba = x;
Yb_3 = 2*y;     
Yb_4 = 1-2*y;   
Yo_o = 1-y/3;   
Yo_Va = y/3;   

l1 = Ya_Sr.*Ya_Ba.*Yb_3.*Yo_o;
l2 = Ya_Sr.*Ya_Ba.*Yb_3.*Yo_o.*(Ya_Sr - Ya_Ba);
l3 = Ya_Sr.*Ya_Ba.*Yb_3.*Yo_Va;
l4 = Ya_Sr.*Ya_Ba.*Yb_3.*Yo_Va.*(Ya_Sr - Ya_Ba);
l5 = Ya_Sr.*Ya_Ba.*Yb_4.*Yo_o;
l6 = Ya_Sr.*Ya_Ba.*Yb_4.*Yo_o.*(Ya_Sr - Ya_Ba);
l7 = Ya_Sr.*Ya_Ba.*Yb_4.*Yo_Va;
l8 = Ya_Sr.*Ya_Ba.*Yb_4.*Yo_Va.*(Ya_Sr - Ya_Ba);
l9 = Yb_3.*Yb_4.*Ya_Sr.*Yo_o;
l10 = Yb_3.*Yb_4.*Ya_Sr.*Yo_o.*(Yb_3 - Yb_4);
l11 = Yb_3.*Yb_4.*Ya_Ba.*Yo_o;
l12 = Yb_3.*Yb_4.*Ya_Ba.*Yo_o.*(Yb_3 - Yb_4);
l13 = Yb_3.*Yb_4.*Ya_Sr.*Yo_Va;
l14 = Yb_3.*Yb_4.*Ya_Sr.*Yo_Va.*(Yb_3 - Yb_4);
l15 = Yb_3.*Yb_4.*Ya_Ba.*Yo_Va;
l16 = Yb_3.*Yb_4.*Ya_Ba.*Yo_Va.*(Yb_3 - Yb_4);
l17 = Yo_o.*Yo_Va.*Ya_Sr.*Yb_3;
l18 = Yo_o.*Yo_Va.*Ya_Sr.*Yb_3.*(Yo_o - Yo_Va);
l19 = Yo_o.*Yo_Va.*Ya_Ba.*Yb_3;
l20 = Yo_o.*Yo_Va.*Ya_Ba.*Yb_3.*(Yo_o - Yo_Va);
l21 = Yo_o.*Yo_Va.*Ya_Sr.*Yb_4;
l22 = Yo_o.*Yo_Va.*Ya_Sr.*Yb_4.*(Yo_o - Yo_Va);
l23 = Yo_o.*Yo_Va.*Ya_Ba.*Yb_4;
l24 = Yo_o.*Yo_Va.*Ya_Ba.*Yb_4.*(Yo_o - Yo_Va);
l9l21 = l9 + l21;
l11l23 = l11 + l23;

eqns = [l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18 l19 l20 l21 l22 l23 l24 l9l21 l11l23];

x1 = [X;X];
x2 = [d;d];
x_regres = [ones(size(x1)) x1 x2 x1.*x2 x1.^2 x2.^2 x1.^2.*x2.^2 x1.^3 x2.^3 x2.^4 ...
    x1.*x2.^2 x1.*x2.^3 x1.*x2.^4 x1.^2.*x2 x1.^2.*x2.^3 x1.^2.*x2.^4 x1.^3.*x2 ...
    x1.^3.*x2.^2 x1.^3.*x2.^3 x1.^3.*x2.^4];
y_regres = [A(:,n);A(:,m)];
[b,~,~,~,stats] = regress(y_regres,x_regres);


fntsz = 20;

[xx,yy] = meshgrid((linspace(0,1)),(linspace(0,0.5)));
[X1FIT,X2FIT] = meshgrid((linspace(0,1)),(linspace(0,0.5)));
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT +...
    b(5)*X1FIT.^2 + b(6)*X2FIT.^2 + b(7)*X1FIT.^2.*X2FIT.^2 + b(8)*X1FIT.^3 +...
    b(9)*X2FIT.^3 + b(10)*X2FIT.^4 + b(11)*X1FIT.*X2FIT.^2 + ...
    b(12)*X1FIT.*X2FIT.^3 + b(13)*X1FIT.*X2FIT.^4 + b(14)*X1FIT.^2.*X2FIT + ...
    b(15)*X1FIT.^2.*X2FIT.^3 + b(16)*X1FIT.^2.*X2FIT.^4 + b(17)*X1FIT.^3.*X2FIT + ...
    b(18)*X1FIT.^3.*X2FIT.^2 + b(19)*X1FIT.^3.*X2FIT.^3 + b(20)*X1FIT.^3.*X2FIT.^4 ;

eqn1 = matlabFunction(eqns(n));
eqn2 = matlabFunction(eqns(m));

figure
hold on
surf(xx,yy,eqn1(xx,yy),'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
surf(xx,yy,eqn2(xx,yy),'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none')
mesh(X1FIT,X2FIT,YFIT,'FaceColor','b','FaceAlpha',0.5, 'EdgeColor','none')
legend(['L',num2str(n)],['L',num2str(m)],'Linear Regression','FontSize',fntsz)
title([Title{n},' & ',Title{m},' Compare'])
xlabel('x','FontSize',fntsz)
ylabel('\delta','FontSize',fntsz)