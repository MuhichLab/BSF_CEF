function [G_sol,G_dft,Go,Go_dft] = Temp_Expansion_endmembers(G_soln,G_o)
% In this script we expand our gibbs energy models based on a temperature
% expansion out to the first three terms G = A + B*T + C*T*ln(T)

    syms T real
    syms gDiff [1 2]
    syms Ho [1 2]
    syms G1 G2 G3 G4
    syms gDiffA [1 2]
    syms gDiffB [1 2]
    syms gDiffC [1 2]
    
% Note here we also make the mathmatical simplifications of the endmember
% terms descirbed in the paper
    %% Full model subs 
    % First we build two full solution models that have the excess terms in
    % them, tho we do not do anyting with the excess terms.
    
    % Full Solution model
    G_sol = subs(G_soln,[G1 G2],...
        [(G3 - gDiff1) (G4 - gDiff2)]);

    G_sol = simplify(expand(subs(G_sol,[gDiff1 gDiff2],[(gDiffA1 + gDiffB1*T + gDiffC1*T*log(T))...
        (gDiffA2 + gDiffB2*T + gDiffC2*T*log(T))])));

    % Full Solution Model DFT Side (non-derivative no T dependece)
    % NOTE you can not sub 0 in for T in G_soln because you get NaN 
    % if C*T*ln(T) is present. So you must build from scratch

    G_dft = subs(G_soln,[G1 G2 T],[(G3 - gDiff1) (G4 - gDiff2) 0]);

    % Notice there is no T expression for the DFT solution model
    G_dft = simplify(expand(subs(G_dft,[gDiff1 gDiff2],[gDiffA1 gDiffA2])));

    %% Go model subs Go = G_end - T*S

    % Go subsitutions
    Go = subs(G_o,[G1 G2],...
        [(G3 - gDiff1) (G4 - gDiff2)]);%,'Vars',{[L],x,y});

    Go = simplify(expand(subs(Go,[gDiff1 gDiff2],[(gDiffA1 + gDiffB1*T + gDiffC1*T*log(T))...
        (gDiffA2 + gDiffB2*T + gDiffC2*T*log(T))])));

    % Go_dft
    Go_dft = subs(G_o,[G1 G2 T],...
        [(G3 - gDiff1) (G4 - gDiff2) 0]);

    Go_dft = simplify(expand(subs(Go_dft,[gDiff1 gDiff2],[gDiffA1 gDiffA2])));

end