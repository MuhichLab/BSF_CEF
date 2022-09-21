clear ;
clc;

%% Which Model are we loading 
file_save_name = "model_0001";

load(file_save_name + "_allfits.mat") % This is all_fits
load(file_save_name + "_d0s.mat") % this is d0_progress

% For calculating the error of only one combination of terms
fit_col = 3; % 3 = 3 excess terms in G_soln

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

%[G_dft, dG_dy, G_end_dft, G_end_diff, G_soln] =CEF_with_phase_shift(expansion);

%Generate base CEF
[G_soln,G_ex,G_o] = CEF();

% Expand Gibbs endmember terms to the for G(T) = A + B*T + C*T*ln(T)
[G_sol,G_dft,Go,Go_dft] = Temp_Expansion_endmembers(G_soln,G_o);

% Expand Gibbs excess terms to the for G(T) = A + B*T + C*T*ln(T)
[Gsol,Gdft,Gex] = Temp_Expansion_excess(G_sol,G_dft,G_ex);

%% Since we have DFT data we can "localize" the the enthalpic energy
% Note we still allow correction to this - explained in paper

[Gsol,Go,Gdft,Go_dft] = Localize_H(Gsol,Go,Gdft,Go_dft);

%% Take necessary derivatives

dG_dy = simplify(expand(diff(Gsol,y)));
dGo_dy = simplify(expand(diff(Go,y)));

%% Assign parmaeters based off seleceted number of exces terms (fit_col)

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
    
G_soln_test = subs(sym(Gsol),[sub_me],[sub_value]); % eV @ 300 C
dG_soln_test = diff(G_soln_test,y);

G_soln_test = matlabFunction(G_soln_test);
dG_soln_test = matlabFunction(dG_soln_test);

%% Checking the abs Error of muhg_o

my_muh_error = zeros(1,length(dd));
for n = 1:length(x_exp) 
    x_vals = [0,0.05,0.1,0.15,0.2];
    for k = 1:length(x_vals)
        if x_exp(n) == x_vals(k)
            dref = d0_progress(k,end-(fit_col)+1);
            break
        end
    end
    muhg_pred = -dG_soln_test(Temp(n),x_exp(n),dd(n)+dref);
    my_muh_error(n) = abs(muhg_pred - muhg_o(n));
end
my_muh_err_avg = mean(my_muh_error);
muh_err_std = std(my_muh_error);
fprintf('This is the error for calualting chemical potenital given T, X, and delta: %f +- %f eV (%f +- %f kJ/mol).\n',my_muh_err_avg,muh_err_std,my_muh_err_avg*96.487,muh_err_std*96.487); 

%% Checking the abs and true Error of dd

dd_error = zeros(1,length(dd));
m = length(x_exp);
heat_mp = [];
parfor n = [1:m]
    x_vals = [0,0.05,0.1,0.15,0.2];
    for k = 1:length(x_vals)
        if x_exp(n) == x_vals(k)
            dref = d0_progress(k,end-(fit_col)+1);
            break
        end
    end
    subbed_dG_soln = matlabFunction(subs(sym(dG_soln_test),[T x],[Temp(n) x_exp(n)]));
    eqn = @(y) -(subbed_dG_soln(y)) - get_mu_o(Temp(n),Press(n));
    dd_pred = fzero(eqn,dd(n)+dref);
    dd_error(n) = abs(dd_pred - (dd(n)+dref));
    true_dd_error(n) = (dd_pred - (dd(n)+dref));
    heat_mp = [heat_mp ; [x_exp(n) Temp(n) Press(n) true_dd_error(n)]];
end

%% Getting all MAEs for all combinations of excess terms

MAEs = zeros(2,22);
all_dd_error = zeros(22,length(x_exp));
all_true_dd_error = zeros(22,length(x_exp));
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

    G_soln_sub = subs(sym(Gsol),[sub_me],[sub_value]); 
    dG_soln_sub = simplify(expand(diff(G_soln_sub,y)));

    G_soln_sub = matlabFunction(G_soln_sub);
    dG_soln_sub = matlabFunction(dG_soln_sub);
    
    m = length(x_exp);
    parfor n = [1:m]
        x_vals = [0,0.05,0.1,0.15,0.2];
        for k = 1:length(x_vals)
            if x_exp(n) == x_vals(k)
                dref = d0_progress(k,end-(fit_col)+1);
                break
            end
        end
        subbed_dG_soln = matlabFunction(subs(sym(dG_soln_sub),[T x],[Temp(n) x_exp(n)]))
        eqn = @(y) -(subbed_dG_soln(y)) - get_mu_o(Temp(n),Press(n));
        dd_pred = fzero(eqn,dd(n)+dref);
        all_dd_error(fit_col,n) = abs(dd_pred - (dd(n)+dref));
        all_true_dd_error(fit_col,n) = (dd_pred - (dd(n)+dref));
    end
    pd = fitdist(transpose(all_dd_error(fit_col,:)),'Normal');
    mu = pd.mu;
    sigma = pd.sigma;
    fprintf('The mean absolute error for %d term(s) predicting delta is %f +- %f.\n',spt,mu,sigma);
    y_norm(fit_col,:) = normpdf(sort(all_dd_error(fit_col,:)),mu,sigma);
    MAEs(1,23-spt) = mu;
    MAEs(2,23-spt) = sigma;
    
end

%%
save(file_save_name + '_get_errors')
save(file_save_name + '_MAEs','MAEs')
