function [Gsol,Go,Gdft,Go_dft] = Localize_H(Gsol,Go,Gdft,Go_dft)

syms x y T real
syms G1 G2 G3 G4
syms L [1 24] real
syms gDiffA [1 2] real
syms gDiffB [1 2] real
syms gDiffC [1 2] real
syms A [1 24] real
syms B [1 24] real
syms C [1 24] real


% These come from the DFT data - end point values per mol ABO3
% THESE MAY NEED TO CHANGE - YOU SHOULD VARIFY
g1 = -550.31/8;  %SrFeO3   (Fe 4+) (x = 0 y = 0)
g2 = -683.90/8;  %BaFeO3   (Fe 4+) (x = 1 y = 0)
g3 = -546.87/8;  %SrFeO2.5 (Fe 3+) (x = 0 y = 0.5)
g4 = -682.27/8;  %BaFeO2.5 (Fe 3+) (x = 1 y = 0.5)

syms Ho [1 2]

% Full Solution Model
Gsol = subs(Gsol,[gDiffA1 gDiffA2],...
    [(Ho1 + gDiffA1) (Ho2 + gDiffA2)]);

Gsol = simplify(expand(subs(Gsol,[Ho1 Ho2 G3 G4],...
    [(g3 - g1) (g4 - g2) g3 g4])));

% Go Solution Model
Go = subs(Go,[gDiffA1 gDiffA2],...
    [(Ho1 + gDiffA1) (Ho2 + gDiffA2)]);

Go = simplify(expand(subs(Go,[Ho1 Ho2 G3 G4],...
    [(g3 - g1) (g4 - g2) g3 g4])));

% Gdft Solution model  (No T terms)
Gdft = subs(Gdft,[gDiffA1 gDiffA2],...
    [(Ho1 + gDiffA1) (Ho2 + gDiffA2)]);

Gdft = simplify(expand(subs(Gdft,[Ho1 Ho2 G3 G4],...
    [(g3 - g1) (g4 - g2) g3 g4])));

% Go_dft Solution model  (No T terms)
Go_dft = subs(Go_dft,[gDiffA1 gDiffA2],...
    [(Ho1 + gDiffA1) (Ho2 + gDiffA2)]);

Go_dft = simplify(expand(subs(Go_dft,[Ho1 Ho2 G3 G4],...
    [(g3 - g1) (g4 - g2) g3 g4])));

end