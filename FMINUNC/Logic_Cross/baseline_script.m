clear;
clc
%%
file_save_name = "model_guess0001_no_intial_fits";

total_start = tic;

exp_data_name = "Transformed_exp_data.xlsx";
dft_data_name = "BSF_data_phase_shift.xlsx";
    % options for CEF build `phase_shift` , no_phase_shift`, `no_dft`
CEF_build = "phase_shift";
    %pick the expansion of Cp(T) 1 = no Cp basically 2 = 0th order 3 = 1st order
    % 1 is wierd because you arent using a  regular Cp expansion
expansion = 2;
    % X0 = is the instial guess for the d0 solvers
X0 = 0.02;
    % x_vals are the experimental mol fractions gathered
x_vals = [0,0.05,0.1,0.15,0.2];
    %d00_val is the energy value where d0 = 0
d_chem = round(get_mu_o(300+273.15,0.9),4); % d0 = 0 at 300 C and 0.9 Po2
    % ball_park is a no self consistent d0 fit that localizes the coeffcients
    % before entering main loop (yes/no)
ball_park = "no";
    % before_loop_fit is a setting to fit all coeff befor eentering the loop
    % that drops terms, essientially this just does the first fit of the loop
before_loop_fit = "no";


add_me = 0.0001; %how much to add for a non zero initial guess

%% Load Experimental data
exp_data = readtable(exp_data_name);

Press = table2array(exp_data(:,2));  % partial pressure O2 in bar
Temp = table2array(exp_data(:,3))+273.15;  % K
x_exp = table2array(exp_data(:,4));  % mol fraciton Ba
dd = table2array(exp_data(:,5));  % delta of delta

%%To Remove specific pieces of data
% remove = (Temp(:) == 811.15);
% Press(remove)=[];
% Temp(remove) = [];
% x(remove )= [];G
% x_exp = x;
% dd(remove) = [];

%% Load DFT data

dft_data = table2array(readtable(dft_data_name));
dft_data_copy = dft_data;
rows = any(isnan(dft_data),2);
dft_data(rows,:) = [];

X = dft_data(:,2); % mol fract Ba
Y = dft_data(:,1); % delta
Z = dft_data(:,3)/8; % E eV per ABO3

%% O2 chemical potential (convert to eV/K)

Ho = [];
So = [];
for t = 1:length(Temp);
    [H,S] = get_O2_thermo(Temp(t));
    Ho(t) = H;
    So(t) = S;
end
conv = 96.487;  % 1 eV = 96.487 kJ/mol
R_gas = 8.314462/1000/conv; %eV/K
Ho = transpose(Ho);
So = transpose(So);

muhg_o2 = Ho/conv - (Temp).*So/conv + R_gas.*(Temp).*log(Press);

muhg_o = 0.5*muhg_o2;

 
%% Fit equtions

% all equtions come in as symbolic equations

[G_dft, dG_dy, G_end_dft, G_end_diff, G_soln] =CEF_with_phase_shift(expansion);

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

 dG_dy = subs(dG_dy,[L21 L23],[L9 L11]);
 G_soln = subs(G_soln,[L21 L23],[L9 L11]);
 G_dft = subs(G_dft,[L21 L23],[L9 L11]);


%% Expand L terms in excess to be G(T) terms
syms A [1 24] real
syms B [1 24] real
syms C [1 24] real
syms D [1 24] real

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

dG_dy = subs(dG_dy,sub_me,sub_value);
G_soln = subs(G_soln,sub_me,sub_value);
G_dft = subs(G_dft,sub_me,sub_value_dft);


%% Make all equatiosn matlab functions and simplify

G_dft = matlabFunction(simplify(expand(G_dft)));
dG_dy = matlabFunction(simplify(expand(dG_dy)));
G_end_dft = matlabFunction(simplify(expand(G_end_dft)));
G_end_diff = matlabFunction(simplify(expand(G_end_diff)));
G_soln = matlabFunction(simplify(expand(G_soln)));



%% fmincon options

options = optimoptions('fminunc','MaxFunctionEvaluations',...
    1000000,'MaxIterations',1000000,'UseParallel',true,...
    'FunctionTolerance',1e-10,'StepTolerance',1.0000e-12,...
    'OptimalityTolerance',1e-12,'PlotFcn',{'optimplotfval'});

%% Endmebers only fit based off extent of Cp expansion 

% NOTE: I pass in -muhg_o because the fit is -dG/dy = muh

guess0 = zeros(6,1)+add_me;

[guess_end] = ...
fminunc(@(guess_end)min_me_end_only_2(guess_end,X,Y,Z,x_exp,dd,X0,x_vals,-muhg_o,...
Temp,G_end_dft,G_end_diff,d_chem), guess0,options);

d0_end_only = solve_dref_end_only_2(guess_end,G_end_diff,X0,x_vals,d_chem)
d0_end = d0_end_only;
    
%% Endmember values for in progress watching

final_tab_endMem = array2table(round(guess_end(1:6),5), ...
      'VariableNames',{'Coeff'});
final_tab_endMem.Properties.RowNames = ["gDiffA1","gDiffA2","gDiffB1",...
    "gDiffB2","gDiffC1","gDiffC2"]


%% Plugging in endmembers into all equations needed fot fits and simplify


G_dft = matlabFunction(subs(sym(G_dft),[gDiffA1 gDiffA2],...
    [guess_end(1) guess_end(2)]));
dG_dy = matlabFunction(subs(sym(dG_dy),[gDiffA1 gDiffA2 gDiffB1 gDiffB2 gDiffC1 gDiffC2],...
    [guess_end(1) guess_end(2) guess_end(3) guess_end(4) guess_end(5) guess_end(6)]));
G_soln = matlabFunction(subs(sym(G_soln),[gDiffA1 gDiffA2 gDiffB1 gDiffB2 gDiffC1 gDiffC2],...
    [guess_end(1) guess_end(2) guess_end(3) guess_end(4) guess_end(5) guess_end(6)]));

G_dft = matlabFunction(simplify(expand(sym(G_dft))));
dG_dy = matlabFunction(simplify(expand(sym(dG_dy))));
G_soln = matlabFunction(simplify(expand(sym(G_soln))));


%% Cross Fits intital ball park fit using end-member d0s

if ball_park == "yes"
        guess0 = zeros(66,1)+add_me;
        [guess,error]  = ...
        fminunc(@(guess)min_me_first_2(guess,X,Y,Z,x_exp,dd,-muhg_o,...
        Temp,d0_end_only,dG_dy,G_dft,x_vals),guess0,options);
else
    guess = zeros(66,1)+add_me;
end



%% Intial fit of all coeffs before entering loop, basically just the first step of the while loop. 

if before_loop_fit == "yes"
    dp_terms = ones(length(guess),1);
    guess0 = guess;
    [guess,error]  =...
        fminunc(@(guess)min_me_2(guess,X,Y,Z,x_exp,dd,-muhg_o,Temp,...
        G_dft,dG_dy,x_vals,d0_end_only,d_chem,dp_terms),guess0,options);

else
    guess = zeros(66,1)+add_me;
end

%% Dropping terms and refitting

% The number 22 in the shapses come from intially having 24 L terms but
% settign L21 = L9 and L23 = L11 droppping the total number of terms to 22.

dp_terms = ones(length(guess),1);
d0_progress = zeros(length(x_vals),22);
all_fits = zeros(length(guess)+length(guess_end),22);
spot_track = [];
count = 1;
num_terms = 25; % arbitrary value
guess_drop = guess;
limits = ones([length(guess),1]);

while num_terms > 1

    % keeps values for remaing coeffs from last fit
    guess0 = guess_drop(dp_terms>0);
    
% %     resets intial guess to 0
%     guess0 = zeros(sum(dp_terms),1)+add_me;

    [guess_drop,drop_error]  =...
        fminunc(@(guess_drop)min_me_2(guess_drop,X,Y,Z,x_exp,dd,-muhg_o,Temp,...
        G_dft,dG_dy,x_vals,d0_end_only,d_chem,dp_terms),guess0,options);
    
    % rebuild the full guess_drop array
    og_guess = guess_drop;
    guess_drop = zeros(length(dp_terms),1);
    n = 1;
    for k = (1:length(dp_terms))
        if dp_terms(k) == 1
            guess_drop(k) = og_guess(n);
            n = n + 1;
        end
    end  

    all_fits(:,count) = [guess_drop; guess_end];
    d0 = solve_dref_2(guess_drop,dG_dy,x_vals,d0_end_only,d_chem);
    d0_progress(:,count) = d0 
    
    original_error = drop_error;
    [drop_this, combos1, combos2, combos3] = drop_me_2(guess_drop,...
        X,Y,Z,x_exp,dd,-muhg_o,Temp,G_dft,dG_dy,x_vals,original_error,...
        d0,spot_track);



    drop_por1 = drop_this(1:22);%.*dp_terms(1:22);
    drop_por2 = drop_this(23:253);
    drop_por3 = drop_this(254:1793);

    
    if num_terms > 3

        min_val1 = min(drop_por1(drop_por1>0));
        min_val2 = min(drop_por2(drop_por2>0));
        min_val3 = min(drop_por3(drop_por3>0));

        spot1 = find(drop_por1 == min_val1);
        spot2 = find(drop_por2 == min_val2);
        spot3 = find(drop_por3 == min_val3);

        combos1(spot1);
        combos2(spot2,:);
        combos3(spot3,:);

        fprintf('Combo 1 are these terms and this Value: %d       : %.3e \n',combos1(spot1),min_val1);
        fprintf('Combo 2 are these terms and this Value: %d %d   : %.3e \n',combos2(spot2,:),min_val2);
        fprintf('Combo 3 are these terms and this Value: %d %d %d : %.3e \n',combos3(spot3,:),min_val3);


       if (any(combos3(spot3,:) == combos2(spot2,1))) && ...
                (any(combos3(spot3,:) == combos2(spot2,2))) && min_val2 < min_val3
            new_min = min(drop_por1(combos2(spot2,:)));
            spot = find(drop_por1 == new_min);
            dp_terms(spot,1) = 0;
            dp_terms(spot+22,1) = 0;
            dp_terms(spot+44,1) = 0;
 
       elseif (any(combos3(spot3,:) == combos2(spot2,1))) && ...
                (any(combos3(spot3,:) == combos2(spot2,2))) && min_val2 > min_val3
            spot = setdiff(combos3(spot3,:),combos2(spot2,:));
            dp_terms(spot,1) = 0;
            dp_terms(spot+22,1) = 0;
            dp_terms(spot+44,1) = 0;
  
       else

            if min_val2 < min_val3
                new_min = min(drop_por1(combos2(spot2,:)));
                spot = find(drop_por1 == new_min);
                dp_terms(spot,1) = 0;
                dp_terms(spot+22,1) = 0;
                dp_terms(spot+44,1) = 0;

            elseif min_val3 < min_val2
                new_min = min(drop_por1(combos3(spot3,:)));
                spot = find(drop_por1 == new_min);
                dp_terms(spot,1) = 0;
                dp_terms(spot+22,1) = 0;
                dp_terms(spot+44,1) = 0;

            end
       end


    elseif num_terms <= 3
        min_val1 = min(drop_por1(drop_por1>0));
        spot1 = find(drop_por1 == min_val1);
        spot = spot1;
        dp_terms(spot,1) = 0;
        dp_terms(spot+22,1) = 0;
        dp_terms(spot+44,1) = 0;

    end

    fprintf('This was the dropped L term(spot): %d \n',spot);
    fprintf('There are %d L terms left. \n',sum(dp_terms(1:22)));

    spot_track(count) = spot;
    count = count + 1;
    num_terms = sum(dp_terms(1:22))+1; % Nuumber L terms fitted
end

%% Table of all fits

all_tab = array2table(round(all_fits,5));

all_tab.Properties.RowNames = ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9",...
    "A10", "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19",...
    "A20", "A22", "A24", "B1", "B2", "B3", "B4", "B5", "B6",...
    "B7", "B8", "B9", "B10", "B11", "B12", "B13", "B14", "B15", "B16", "B17",...
    "B18", "B19", "B20", "B22", "B24", "C1", "C2", "C3", "C4", "C5", "C6",...
    "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17",...
    "C18", "C19", "C20", "C22", "C24","gDiffA1","gDiffA2","gDiffB1",...
    "gDiffB2","gDiffC1","gDiffC2"]

%% Saving 
total_stop = toc(total_start);

save(file_save_name + "_workspace");
save(file_save_name + "_allfits",'all_fits');
save(file_save_name + "_d0s",'d0_progress'); 
