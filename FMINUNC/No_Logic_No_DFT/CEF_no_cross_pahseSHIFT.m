function [dG_dy, G_end_diff, G] = CEF();

    %Constants
    m = 5/6*log(5/6) + 1/6*log(1/6);
    r = 8.617333262E-5; % eV/K;

    syms x y

    Ya_Sr = 1-x;
    Ya_Ba = x;
    %Yb_3 = 6-2*y;   % In terms of y where y = (3 - delta)
    Yb_3 = 2*y;     % In terms of delta
    %Yb_4 = 2*y-5;   % In terms of y where y = (3 - delta)
    Yb_4 = 1-2*y;   % In terms of delta
    %Yo_o = y/3;     % In terms of y where y = (3 - delta)
    Yo_o = 1-y/3;   % In terms of delta
    %Yo_Va = 1-y/3;  % In terms of y where y = (3 - delta)
    Yo_Va = y/3;    % In terms of delta

    % Entropy
    syms R
    S = -R*(Ya_Sr*log(Ya_Sr) + Ya_Ba*log(Ya_Ba) + Yb_3*log(Yb_3) + Yb_4*log(Yb_4))...
        -3*R*(Yo_o*log(Yo_o) + Yo_Va*log(Yo_Va));

    %End-Member
    syms G1 G2 G3 G4 Go G_Sr3O G_Ba3O R T M  % M = (5/6*log(5/6)+1/6*log(1/6))

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

    %Excess
    syms L [1 24]

%     G_excess = (Ya_Sr*Ya_Ba*Yb_3*Yo_o*(L1 + (Ya_Sr - Ya_Ba)*L2) + Ya_Sr*Ya_Ba*Yb_3*Yo_Va*(L3 + (Ya_Sr - Ya_Ba)*L4) + ...
%         Ya_Sr*Ya_Ba*Yb_4*Yo_o*(L5 + (Ya_Sr - Ya_Ba)*L6) + Ya_Sr*Ya_Ba*Yb_4*Yo_Va*(L7 + (Ya_Sr - Ya_Ba)*L8) + ...
%         Yb_3*Yb_4*Ya_Sr*Yo_o*(L9 + (Yb_3 - Yb_4)*L10) + Yb_3*Yb_4*Ya_Ba*Yo_o*(L11 + (Yb_3 - Yb_4)*L12) + ...
%         Yb_3*Yb_4*Ya_Sr*Yo_Va*(L13 + (Yb_3 - Yb_4)*L14) + Yb_3*Yb_4*Ya_Ba*Yo_Va*(L15 + (Yb_3 - Yb_4)*L16) + ...
%         Yo_o*Yo_Va*Ya_Sr*Yb_3*(L17 + (Yb_3 - Yb_4)*L18) + Yo_o*Yo_Va*Ya_Ba*Yb_3*(L19 + (Yb_3 - Yb_4)*L20) + ...
%         Yo_o*Yo_Va*Ya_Sr*Yb_4*(L21 + (Yb_3 - Yb_4)*L22) + Yo_o*Yo_Va*Ya_Ba*Yb_4*(L23 + (Yb_3 - Yb_4)*L24));
        
    G_excess = (Ya_Sr*Ya_Ba*Yb_3*Yo_o*(L1 + (Ya_Sr - Ya_Ba)*L2) + Ya_Sr*Ya_Ba*Yb_3*Yo_Va*(L3 + (Ya_Sr - Ya_Ba)*L4) + ...
        Ya_Sr*Ya_Ba*Yb_4*Yo_o*(L5 + (Ya_Sr - Ya_Ba)*L6) + Ya_Sr*Ya_Ba*Yb_4*Yo_Va*(L7 + (Ya_Sr - Ya_Ba)*L8) + ...
        Yb_3*Yb_4*Ya_Sr*Yo_o*(L9 + (Yb_3 - Yb_4)*L10) + Yb_3*Yb_4*Ya_Ba*Yo_o*(L11 + (Yb_3 - Yb_4)*L12) + ...
        Yb_3*Yb_4*Ya_Sr*Yo_Va*(L13 + (Yb_3 - Yb_4)*L14) + Yb_3*Yb_4*Ya_Ba*Yo_Va*(L15 + (Yb_3 - Yb_4)*L16) + ...
        Yo_o*Yo_Va*Ya_Sr*Yb_3*(L17 + (Yo_o - Yo_Va)*L18) + Yo_o*Yo_Va*Ya_Ba*Yb_3*(L19 + (Yo_o - Yo_Va)*L20) + ...
        Yo_o*Yo_Va*Ya_Sr*Yb_4*(L21 + (Yo_o - Yo_Va)*L22) + Yo_o*Yo_Va*Ya_Ba*Yb_4*(L23 + (Yo_o - Yo_Va)*L24));
    
    % 
    % G_excess = (Ya_Sr*Ya_Ba*Yb_3*Yo_o*(L1 + (Ya_Sr - Ya_Ba)*L2) + Ya_Sr*Ya_Ba*Yb_3*Yo_Va*(L1 + (Ya_Sr - Ya_Ba)*L2) + ...
    %     Ya_Sr*Ya_Ba*Yb_4*Yo_o*(L1 + (Ya_Sr - Ya_Ba)*L2) + Ya_Sr*Ya_Ba*Yb_4*Yo_Va*(L1 + (Ya_Sr - Ya_Ba)*L2) + ...
    %     Yb_3*Yb_4*Ya_Sr*Yo_o*(L1 + (Yb_3 - Yb_4)*L2) + Yb_3*Yb_4*Ya_Ba*Yo_o*(L1 + (Yb_3 - Yb_4)*L2) + ...
    %     Yb_3*Yb_4*Ya_Sr*Yo_Va*(L1 + (Yb_3 - Yb_4)*L2) + Yb_3*Yb_4*Ya_Ba*Yo_Va*(L1 + (Yb_3 - Yb_4)*L2) + ...
    %     Yo_o*Yo_Va*Ya_Sr*Yb_3*(L1 + (Yb_3 - Yb_4)*L2) + Yo_o*Yo_Va*Ya_Ba*Yb_3*(L1 + (Yb_3 - Yb_4)*L2) + ...
    %     Yo_o*Yo_Va*Ya_Sr*Yb_4*(L1 + (Yb_3 - Yb_4)*L2) + Yo_o*Yo_Va*Ya_Ba*Yb_4*(L1 + (Yb_3 - Yb_4)*L2));

    
    G_soln = simplify(expand(G_end - T*S + G_excess));
%     G = matlabFunction(simplify(expand(subs(G_soln,[G1 G2 G3 G4 M R],[g1 g2 g3 g4 m r]))));%,'Vars',{[L],x,y});
%     %Derivatives
% 
%     % w.r.t y(delta)
%     dG_dy = matlabFunction(simplify(expand(diff(subs(G_soln,[G1 G2 G3 G4 M R],[g1 g2 g3 g4 m r]),y))));%,'Vars',{[L],x,y});
% 
%     %w.r.t x
%     dG_dx = matlabFunction(simplify(expand(diff(subs(G_soln,[G1 G2 G3 G4 M R],[g1 g2 g3 g4 m r]),x))));%,'Vars',{[L],x,y});

% %    With G(T)
    
    syms T
    syms gDiff [1 2]
    syms Ho [1 2]
    syms gDiffA [1 2]
    syms gDiffB [1 2]
    syms gDiffC [1 2]
    
    G = simplify(expand(subs(G_soln,[M R],...
        [m r])));%,'Vars',{[L],x,y});
   
    %Derivatives

    % w.r.t y(delta)
    dG_dy = simplify(diff(subs(G_soln,[G1 G2 M R],...
        [(G3 - gDiff1) (G4 - gDiff2) m r]),y));%,'Vars',{[L],x,y});
    
    dG_dy = subs(dG_dy,[gDiff1 gDiff2],[(Ho1 + gDiffB1*T + gDiffC1*T*log(T))...
        (Ho2 + gDiffB2*T + gDiffC2*T*log(T))]);
    
    dG_dy = subs(dG_dy,[Ho1 Ho2],[(gDiffA1) (gDiffA2)]);
    
    G_end_diff = simplify(expand(diff((G_end - T*S),y)));
    
    G_end_diff= simplify(subs(G_end_diff,[G1 G2 M R],...
        [(G3 - gDiff1) (G4 - gDiff2) m r]));%,'Vars',{[L],x,y});
    
    G_end_diff = subs(G_end_diff,[gDiff1 gDiff2],[(Ho1 + gDiffB1*T + gDiffC1*T*log(T))...
        (Ho2 + gDiffB2*T + gDiffC2*T*log(T))]);
    
    G_end_diff = subs(G_end_diff,[Ho1 Ho2],[(gDiffA1) (gDiffA2)]);

%     %w.r.t x
%     dG_dx = simplify(diff(subs(G_soln,[G1 G2 M R],...
%         [(G3 - gDiff1) (G4 - gDiff2) m r]),x));%,'Vars',{[L],x,y});
%     
%     dG_dx = subs(dG_dx,[gDiff1 gDiff2],[(gDiffA1 + gDiffB1*T + gDiffC1*T*log(T))...
%     (gDiffA2 + gDiffB2*T + gDiffC2*T*log(T))]);
%     
%     dG_dx = matlabFunction(subs(dG_dx,[gDiffA1 gDiffA2],[(g3 - g1) (g4 - g2)]));

end

