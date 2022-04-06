function [G_dft, dG_dy, G_end_dft, G_end_diff, G_soln] = CEF_with_phase_shift(expansion)
    %end point values for pahse shift 
    g1 = -550.31/8;  %SrFeO3   (Fe 4+)
    g2 = -683.90/8;  %BaFeO3   (Fe 4+)
    g3 = -546.87/8;  %SrFeO2.5 (Fe 3+)
    g4 = -682.27/8;  %BaFeO2.5 (Fe 3+)

    %Constants
    m = 5/6*log(5/6) + 1/6*log(1/6);  % This can chaneg depenting on checmial makeup. 
    r = 8.617333262E-5; % eV/K;

    syms x y real

    Ya_Sr = 1-x;
    Ya_Ba = x;
    Yb_3 = 2*y;     % In terms of delta
    Yb_4 = 1-2*y;   % In terms of delta
    Yo_o = 1-y/3;   % In terms of delta
    Yo_Va = y/3;    % In terms of delta

    % Entropy
    syms R
    S = -R*(Ya_Sr*log(Ya_Sr) + Ya_Ba*log(Ya_Ba) + Yb_3*log(Yb_3) + Yb_4*log(Yb_4))...
        -3*R*(Yo_o*log(Yo_o) + Yo_Va*log(Yo_Va));

    % End-Member
    syms G1 G2 G3 G4 Go G_Sr3O G_Ba3O R T M real  % M = (5/6*log(5/6)+1/6*log(1/6))

    G_Sr4O = G1;
    G_Sr4Va = G1 - 3/2*Go;
    G_Ba4O = G2;
    G_Ba4Va = G2 - 3/2*Go;
    G_Sr3Va = G_Sr3O - G_Sr4O + G_Sr4Va;
    G_Sr3O = solve(((G3 - 1/6*G_Sr3Va - 3*R*T*M)/(5/6))==G_Sr3O,G_Sr3O);
    G_Ba3Va = G_Ba3O - G_Ba4O + G_Ba4Va;
    G_Ba3O = solve((G4 - 1/6*G_Ba3Va - 3*R*T*M)/(5/6)==G_Ba3O,G_Ba3O);

    G_Sr3Va = G_Sr3O - G_Sr4O + G_Sr4Va; % Redefine now that G_Sr3O is solved
    G_Ba3Va = G_Ba3O - G_Ba4O + G_Ba4Va; % Redefine now that G_Ba3O is solved

    G_end = simplify(Ya_Sr*Yb_4*Yo_o*G_Sr4O + Ya_Ba*Yb_4*Yo_o*G_Ba4O + Ya_Sr*Yb_4*Yo_Va*G_Sr4Va + Ya_Ba*Yb_4*Yo_Va*(G_Ba4Va) +...
        Ya_Sr*Yb_3*Yo_o*G_Sr3O + Ya_Ba*Yb_3*Yo_o*G_Ba3O + Ya_Sr*Yb_3*Yo_Va*G_Sr3Va + Ya_Ba*Yb_3*Yo_Va*G_Ba3Va);

    % Excess
    syms L [1 24] real

    G_ex = Ya_Sr*Ya_Ba*Yb_3*Yo_o*(L1 + (Ya_Sr - Ya_Ba)*L2) + Ya_Sr*Ya_Ba*Yb_3*Yo_Va*(L3 + (Ya_Sr - Ya_Ba)*L4) + ...
        Ya_Sr*Ya_Ba*Yb_4*Yo_o*(L5 + (Ya_Sr - Ya_Ba)*L6) + Ya_Sr*Ya_Ba*Yb_4*Yo_Va*(L7 + (Ya_Sr - Ya_Ba)*L8) + ...
        Yb_3*Yb_4*Ya_Sr*Yo_o*(L9 + (Yb_3 - Yb_4)*L10) + Yb_3*Yb_4*Ya_Ba*Yo_o*(L11 + (Yb_3 - Yb_4)*L12) + ...
        Yb_3*Yb_4*Ya_Sr*Yo_Va*(L13 + (Yb_3 - Yb_4)*L14) + Yb_3*Yb_4*Ya_Ba*Yo_Va*(L15 + (Yb_3 - Yb_4)*L16) + ...
        Yo_o*Yo_Va*Ya_Sr*Yb_3*(L17 + (Yo_o - Yo_Va)*L18) + Yo_o*Yo_Va*Ya_Ba*Yb_3*(L19 + (Yo_o - Yo_Va)*L20) + ...
        Yo_o*Yo_Va*Ya_Sr*Yb_4*(L21 + (Yo_o - Yo_Va)*L22) + Yo_o*Yo_Va*Ya_Ba*Yb_4*(L23 + (Yo_o - Yo_Va)*L24);
    

    %% Complete Solution Model 
    G_sol = G_end - T*S + G_ex;

    %% Make expansion of G -> G(T) 1 = B*T 2 = B*T + C*T*ln(T) 3 = B*T + C*T*ln(T) + T^2/2
    
    syms T real
    syms gDiff [1 2] real
    syms Ho [1 2] real
    syms gDiffA [1 2] real
    syms gDiffB [1 2] real
    syms gDiffC [1 2] real
    syms gDiffD [1 2] real
    
    if expansion == 1
        %% Full Solution model and derivative for exp fit 
        G_soln = subs(G_sol,[G1 G2 M R],...
            [(G3 - gDiff1) (G4 - gDiff2) m r]);

        G_soln = subs(G_soln,[gDiff1 gDiff2],[(Ho1 + gDiffB1*T)...
            (Ho2 + gDiffB2*T)]);

        G_soln = simplify(expand(subs(G_soln,[Ho1 Ho2 G3 G4],[(g3 - g1 + gDiffA1) (g4 - g2 + gDiffA2) g3 g4])));

        dG_dy = simplify(expand(diff(G_soln,y)));

        %% Full Solution Model DFT Side (non-derivative no T dependece) 
        % NOTE you can not sub 0 in for T in G_soln of G_sol because you get NaN if C*T*ln(T) is present 
        % So you must build from scratch
        
        G_dft = subs(G_sol,[G1 G2 M R T],[(G3 - gDiff1) (G4 - gDiff2) m r 0]);

        G_dft = subs(G_dft,[gDiff1 gDiff2],[(Ho1) (Ho2)]); %this is repetitive but just following the derivative format

        G_dft = simplify(expand(subs(G_dft,[Ho1 Ho2 G3 G4],[(g3 - g1 + gDiffA1) (g4 - g2 + gDiffA2) g3 g4])));

        %% End Mebers only Model
        % note you can not sub 0 in for T from G_soln because you get NaN 
        % so you must build two seperate models 

        % Build dG_end
        G_end_only = subs((G_end-T*S),[G1 G2 M R],...
            [(G3 - gDiff1) (G4 - gDiff2) m r]);%,'Vars',{[L],x,y});

        G_end_only = subs(G_end_only,[gDiff1 gDiff2],[(Ho1 + gDiffB1*T)...
            (Ho2 + gDiffB2*T)]);

        G_end_mem = simplify(expand(subs(G_end_only,[Ho1 Ho2 G3 G4],[(g3 - g1 + gDiffA1) (g4 - g2 + gDiffA2) g3 g4])));

        G_end_diff = simplify(expand(diff(sym(G_end_mem),y)));

        % Build G_end_dft
        G_end_only = subs(G_end-T*S,[G1 G2 M R T],...
            [(G3 - gDiff1) (G4 - gDiff2) m r 0]);%,'Vars',{[L],x,y});

        G_end_only = subs(G_end_only,[gDiff1 gDiff2],[Ho1 Ho2]);

        G_end_dft = simplify(expand(subs(G_end_only,[Ho1 Ho2 G3 G4],[(g3 - g1 + gDiffA1) (g4 - g2 + gDiffA2) g3 g4])));
        
    elseif expansion == 2
        %% Full Solution model and derivative for exp fit 
        G_soln = subs(G_sol,[G1 G2 M R],...
            [(G3 - gDiff1) (G4 - gDiff2) m r]);

        G_soln = subs(G_soln,[gDiff1 gDiff2],[(Ho1 + gDiffB1*T + gDiffC1*T*log(T))...
            (Ho2 + gDiffB2*T + gDiffC2*T*log(T))]);

        G_soln = simplify(expand(subs(G_soln,[Ho1 Ho2 G3 G4],[(g3 - g1 + gDiffA1) (g4 - g2 + gDiffA2) g3 g4])));

        dG_dy = simplify(expand(diff(G_soln,y)));

        %% Full Solution Model DFT Side (non-derivative no T dependece) 
        % NOTE you can not sub 0 in for T in G_soln of G_sol because you get NaN if C*T*ln(T) is present 
        % So you must build from scratch
        
        G_dft = subs(G_sol,[G1 G2 M R T],[(G3 - gDiff1) (G4 - gDiff2) m r 0]);

        G_dft = subs(G_dft,[gDiff1 gDiff2],[(Ho1) (Ho2)]); %this is repetitive but just following the derivative format

        G_dft = simplify(expand(subs(G_dft,[Ho1 Ho2 G3 G4],[(g3 - g1 + gDiffA1) (g4 - g2 + gDiffA2) g3 g4])));
        %% End Mebers only Model
        % note you can not sub 0 in for T from G_soln because you get NaN 
        % so you must build two seperate models 

        % Build dG_end
        G_end_only = subs((G_end-T*S),[G1 G2 M R],...
            [(G3 - gDiff1) (G4 - gDiff2) m r]);%,'Vars',{[L],x,y});

        G_end_only = subs(G_end_only,[gDiff1 gDiff2],[(Ho1 + gDiffB1*T + gDiffC1*T*log(T))...
            (Ho2 + gDiffB2*T + gDiffC2*T*log(T))]);

        G_end_mem = simplify(expand(subs(G_end_only,[Ho1 Ho2 G3 G4],[(g3 - g1 + gDiffA1) (g4 - g2 + gDiffA2) g3 g4])));

        G_end_diff = simplify(expand(diff(sym(G_end_mem),y)));

        % Build G_end_dft
        G_end_only = subs(G_end-T*S,[G1 G2 M R T],...
            [(G3 - gDiff1) (G4 - gDiff2) m r 0]);%,'Vars',{[L],x,y});

        G_end_only = subs(G_end_only,[gDiff1 gDiff2],[Ho1 Ho2]);

        G_end_dft = simplify(expand(subs(G_end_only,[Ho1 Ho2 G3 G4],[(g3 - g1 + gDiffA1) (g4 - g2 + gDiffA2) g3 g4])));

    elseif expansion == 3
        %% Full Solution model and derivative for exp fit 
        G_soln = subs(G_sol,[G1 G2 M R],...
            [(G3 - gDiff1) (G4 - gDiff2) m r]);

        G_soln = subs(G_soln,[gDiff1 gDiff2],[(Ho1 + gDiffB1*T + gDiffC1*T*log(T) + gDiffD1*T^2/2)...
            (Ho2 + gDiffB2*T + gDiffC2*T*log(T) + gDiffD2*T^2/2)]);

        G_soln = simplify(expand(subs(G_soln,[Ho1 Ho2 G3 G4],[(g3 - g1 + gDiffA1) (g4 - g2 + gDiffA2) g3 g4])));

        dG_dy = simplify(expand(diff(G_soln,y)));

        %% Full Solution Model DFT Side (non-derivative no T dependece) 
        % NOTE you can not sub 0 in for T in G_soln of G_sol because you get NaN if C*T*ln(T) is present 
        % So you must build from scratch

        G_dft = subs(G_sol,[G1 G2 M R T],[(G3 - gDiff1) (G4 - gDiff2) m r 0]);

        G_dft = subs(G_dft,[gDiff1 gDiff2],[(Ho1) (Ho2)]); %this is repetitive but just following the derivative format

        G_dft = simplify(expand(subs(G_dft,[Ho1 Ho2 G3 G4],[(g3 - g1 + gDiffA1) (g4 - g2 + gDiffA2) g3 g4])));
        %% End Mebers only Model
        % note you can not sub 0 in for T from G_soln because you get NaN 
        % so you must build two seperate models 

        % Build dG_end
        G_end_only = subs((G_end-T*S),[G1 G2 M R],...
            [(G3 - gDiff1) (G4 - gDiff2) m r]);%,'Vars',{[L],x,y});

        G_end_only = subs(G_end_only,[gDiff1 gDiff2],[(Ho1 + gDiffB1*T + gDiffC1*T*log(T) + gDiffD1*T^2/2)...
            (Ho2 + gDiffB2*T + gDiffC2*T*log(T) + gDiffD2*T^2/2)]);

        G_end_mem = simplify(expand(subs(G_end_only,[Ho1 Ho2 G3 G4],[(g3 - g1 + gDiffA1) (g4 - g2 + gDiffA2) g3 g4])));

        G_end_diff = simplify(expand(diff(sym(G_end_mem),y)));

        % Build G_end_dft
        G_end_only = subs(G_end-T*S,[G1 G2 M R T],...
            [(G3 - gDiff1) (G4 - gDiff2) m r 0]);%,'Vars',{[L],x,y});

        G_end_only = subs(G_end_only,[gDiff1 gDiff2],[Ho1 Ho2]);

        G_end_dft = simplify(expand(subs(G_end_only,[Ho1 Ho2 G3 G4],[(g3 - g1 + gDiffA1) (g4 - g2 + gDiffA2) g3 g4])));
    end
end

