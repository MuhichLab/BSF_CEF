clear;
clc;
Start = tic();
%%
file_save_name = 'model_0001';

EXP_file = "Transformed_exp_data.xlsx";
DFT_file = "BSF_DFT_20220808.xlsx";

    % X0 = is the instial guess for the d0 solvers
X0 = 0.02;
    % x_vals are the experimental mol fractions gathered
x_vals = [0,0.05,0.1,0.15,0.2];
    %d00_val is the energy value where d0 = 0
d_chem = round(get_mu_o(300+273.15,0.9),4); % d0 = 0 at 300 C and 0.9 Po2
    % ball_park is a non-self consistent d0 fit that localizes the excess
    % coeffcients before entering main loop (yes/no)
ball_park = "no";
    % before_loop_fit is a setting to fit all excess coeffs self consistantly
    % before entering the main loop that drops terms
    % essientially this just does the first fit of the loop
before_loop_fit = "no";

add_me = 0.0001; %how much to add for a non zero initial guess

%% Load Data

[X,Y,Z,x_exp,Temp,Press,dd,muhg_o] = Load_Data(DFT_file,EXP_file);

% % You can visulaize the loaded data by running this line
% Plot_Data(X,Y,Z,Press,Temp,x_exp,dd)
 
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

% Expand Gibbs endmember terms to the form G(T) = A + B*T + C*T*ln(T)
[G_sol,G_dft,Go,Go_dft] = Temp_Expansion_endmembers(G_soln,G_o);

% Expand Gibbs excess terms to the form G(T) = A + B*T + C*T*ln(T)
[Gsol,Gdft,Gex] = Temp_Expansion_excess(G_sol,G_dft,G_ex);

%% Since we have DFT data we can "localize" the the enthalpic energy
% Note we still allow correction to this - explained in paper

[Gsol,Go,Gdft,Go_dft] = Localize_H(Gsol,Go,Gdft,Go_dft);

%% Take necessary derivatives

dG_dy = simplify(expand(diff(Gsol,y)));
dGo_dy = simplify(expand(diff(Go,y)));

%% Make all equations matlab functions and simplify

Gdft = matlabFunction(simplify(expand(Gdft)));
dG_dy = matlabFunction(simplify(expand(dG_dy)));
Go_dft = matlabFunction(simplify(expand(Go_dft)));
dGo_dy = matlabFunction(simplify(expand(dGo_dy)));
Gsol = matlabFunction(simplify(expand(Gsol)));

%% fmincon options
% Uncomment the 'PlotFcn if you wish to watch the minimization while
% running on your local machine - BE CAREFUL of parentheses

options = optimoptions('fminunc','MaxFunctionEvaluations',...
    1000000,'MaxIterations',1000000,'UseParallel',true,...
    'FunctionTolerance',1e-10,'StepTolerance',1.0000e-12,...
    'OptimalityTolerance',1e-12);%,'PlotFcn',{'optimplotfval'});

%% Endmebers only fit - i.e Go parameter minimization 

% Endmember self intilization
guess0 = zeros(6,1)+add_me;

% NOTE: I pass in -muhg_o because the fit is -dG/dy = muh
[guess_end] = ...
fminunc(@(guess_end)min_me_end_only_2(guess_end,X,Y,Z,x_exp,dd,X0,x_vals,-muhg_o,...
Temp,Go_dft,dGo_dy,d_chem), guess0,options);

d0_end_only = solve_dref_end_only_2(guess_end,dGo_dy,X0,x_vals,d_chem)
    
%% Endmember values for in progress watching

final_tab_endMem = array2table(round(guess_end(1:6),5), ...
      'VariableNames',{'Coeff'});
final_tab_endMem.Properties.RowNames = ["gDiffA1","gDiffA2","gDiffB1",...
    "gDiffB2","gDiffC1","gDiffC2"]

%% Plugging in endmembers into all equations needed for fits and simplify

Gdft = subs(sym(Gdft),[gDiffA1 gDiffA2],...
    [guess_end(1) guess_end(2)]);
dG_dy = subs(sym(dG_dy),[gDiffA1 gDiffA2 gDiffB1 gDiffB2 gDiffC1 gDiffC2],...
    [guess_end(1) guess_end(2) guess_end(3) guess_end(4) guess_end(5) guess_end(6)]);
Gsol = subs(sym(Gsol),[gDiffA1 gDiffA2 gDiffB1 gDiffB2 gDiffC1 gDiffC2],...
    [guess_end(1) guess_end(2) guess_end(3) guess_end(4) guess_end(5) guess_end(6)]);

% These are the only equaitons that are used later dGo_dy, Go_dft, and Go
% no longer needed
Gdft = matlabFunction(simplify(expand(sym(Gdft))));
dG_dy = matlabFunction(simplify(expand(sym(dG_dy))));
Gsol = matlabFunction(simplify(expand(sym(Gsol))));


%% Cross Fits intital ball park fit using end-member d0s
% ball_park is a non-self consistent d0 fit that localizes the excess
% coeffcients before entering main loop. This was not used for the
% Logic_Cross fit for the paper, however, theoretically, if you trust your
% endmember fit and the d0s assosited with them, this would provide a good
% starting point for your excess terms. They will still very with changes
% to guess0 (this is a problem with the optimization technique), but the 
% fixed endember terms and non self-consistent d0s will help prevent the
% excess terms from moving out of the "energy well" localized by the
% endmember terms. 

if ball_park == "yes"
        guess0 = zeros(66,1)+add_me;
        [guess,error]  = ...
        fminunc(@(guess)min_me_first_2(guess,X,Y,Z,x_exp,dd,-muhg_o,...
        Temp,d0_end_only,dG_dy,Gdft,x_vals),guess0,options);
else
    guess = zeros(66,1)+add_me;
end

%% Intial fit of all coeffs before entering loop, basically just the first step of the while loop. 
% before_loop_fit is a setting to fit all excess coeffs self consistantly
% before entering the main loop that drops terms essientially this is 
% just does the first step of the following loop

if before_loop_fit == "yes"
    dp_terms = ones(length(guess),1);
    guess0 = guess;
    [guess,error]  =...
        fminunc(@(guess)min_me_2(guess,X,Y,Z,x_exp,dd,-muhg_o,Temp,...
        Gdft,dG_dy,x_vals,d0_end_only,d_chem,dp_terms),guess0,options);

else
    guess = zeros(66,1)+add_me;
end

%% Logic dropping of terms and refitting


dp_terms = ones(length(guess),1);
d0_progress = zeros(length(x_vals),22);
all_fits = zeros(length(guess)+length(guess_end),22);
spot_track = zeros(1,22);
count = 1;
num_terms = 25; % arbitrary value
guess_drop = guess;
limits = ones([length(guess),1]);

while num_terms > 1
% dp_terms controls which terms to set to 0 effectivily dropping them

    % keeps previous values and removes dropped terms
    guess0 = guess_drop(dp_terms>0);

    % Minimization of current fit paramters
    [guess_drop,drop_error]  =...
        fminunc(@(guess_drop)min_me_2(guess_drop,X,Y,Z,x_exp,dd,-muhg_o,Temp,...
        Gdft,dG_dy,x_vals,d0_end_only,d_chem,dp_terms),guess0,options);
    
    % rebuild the full guess_drop with zeros in place for dropped terms - this is necessary here!
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
    
    %[Data_P,Data_L,Data_x,Data_muh,Data_T] = Data_Point_Calc(guess_drop,X,Y,Z,G_dft,dG_dy,x_vals,d0);
    % Significane determination and Drop Logic
    [dp_terms,spot] = Significance_and_Drop_Logic...
        (guess_drop,dp_terms,num_terms,drop_error,X,Y,Z,x_exp,dd,muhg_o,Temp,Gdft,...
        dG_dy,x_vals,d0,spot_track);  

    fprintf('This was the dropped L term(spot): %d \n',spot);
    fprintf('There are %d L terms left. \n',sum(dp_terms(1:22)));

    spot_track(count) = spot;
    count = count + 1;
    num_terms = sum(dp_terms(1:22))+1; % Number L terms fitted
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
    "gDiffB2","gDiffC1","gDiffC2"];

%% Saving 
Stop = toc(Start);

save(file_save_name + "_workspace");
save(file_save_name + "_allfits",'all_fits');
save(file_save_name + "_d0s",'d0_progress'); 
