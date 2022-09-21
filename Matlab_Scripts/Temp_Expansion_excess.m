function [Gsol,Gdft,Gex] = Temp_Expansion_excess(G_soln,G_dft,G_ex)


syms x y T real
syms G1 G2 G3 G4
syms L [1 24] real
syms gDiffA [1 2] real
syms gDiffB [1 2] real
syms gDiffC [1 2] real
syms A [1 24] real
syms B [1 24] real
syms C [1 24] real


%% First - We found that L9 = 6*L21 and L11 = 6*L23 thus we combine them into one term
% NOTE THIS IS SPECIFIC TO THE 3 SUBLATICE 2 COMPONENT SYSTEM OF ABO3

 Gsol = subs(G_soln,[L21 L23],[L9 L11]);
 Gdft = subs(G_dft,[L21 L23],[L9 L11]);
 Gex = subs(G_ex,[L21 L23],[L9 L11]);
 
 %% Now We Exampd the L terms to the form G(T) = A + B*T + C*T*ln(T)
 
 %% Expand L terms in excess to be G(T) terms
% We have two versions one for Gsol and one for the no T eqution Gdft

sub_me = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13...
    L14 L15 L16 L17 L18 L19 L20 L22 L24];

% Used for Gsol and Gex
sub_value_Gsol = ...
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

% Used for Gdft
sub_value_dft = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13...
        A14 A15 A16 A17 A18 A19 A20 A22 A24];

Gsol = subs(Gsol,sub_me,sub_value_Gsol);
Gex = subs(Gex,sub_me,sub_value_Gsol);
Gdft = subs(Gdft,sub_me,sub_value_dft);

end