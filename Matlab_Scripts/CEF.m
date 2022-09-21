function [G_soln,G_ex,G_o] = CEF()
% This function creates a symbolic matlab expression for the CEF expression
% of an ABO3 - 3 lattice 2 compostions on each lattice 

% The Function of Stat
%              G = H - TS
    
    %Constants
    m = 5/6*log(5/6) + 1/6*log(1/6); % 1/6 comes from the selction 0f 0.5 delta as an endmember 
    r = 8.617333262E-5; % eV/K;

    syms x y
    % Site fractions
    % these expression are in terms of mol frac x (x) and extent of
    % reduction delta (y).
    Ya_Sr = 1-x;
    Ya_Ba = x;
    Yb_3 = 2*y;
    Yb_4 = 1-2*y;
    Yo_o = 1-y/3;
    Yo_Va = y/3;    

    %% Entropy
    syms R
    S = -R*(Ya_Sr*log(Ya_Sr) + Ya_Ba*log(Ya_Ba) + Yb_3*log(Yb_3) + Yb_4*log(Yb_4))...
        -3*R*(Yo_o*log(Yo_o) + Yo_Va*log(Yo_Va));
    % sum r defines above allows easier conversion from ev to kJ if desired
    S = subs(S,R,r);

    %% End-Member
    syms G1 G2 G3 G4 Go G_Sr3O G_Ba3O T M  % Here M = (5/6*log(5/6)+1/6*log(1/6))
    
    % Here M = (5/6*log(5/6)+1/6*log(1/6)) is used to simplify the script
    % and m is plugged in for m later. Note this also allows one to esilt
    % change the netural enmber term from 0.5 to say 1 if needed

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
    
    % sub m and r
    G_end = subs(G_end,[M R], [m r]);

    %% Excess
    syms L [1 24] real

    G_ex = Ya_Sr*Ya_Ba*Yb_3*Yo_o*(L1 + (Ya_Sr - Ya_Ba)*L2) + Ya_Sr*Ya_Ba*Yb_3*Yo_Va*(L3 + (Ya_Sr - Ya_Ba)*L4) + ...
        Ya_Sr*Ya_Ba*Yb_4*Yo_o*(L5 + (Ya_Sr - Ya_Ba)*L6) + Ya_Sr*Ya_Ba*Yb_4*Yo_Va*(L7 + (Ya_Sr - Ya_Ba)*L8) + ...
        Yb_3*Yb_4*Ya_Sr*Yo_o*(L9 + (Yb_3 - Yb_4)*L10) + Yb_3*Yb_4*Ya_Ba*Yo_o*(L11 + (Yb_3 - Yb_4)*L12) + ...
        Yb_3*Yb_4*Ya_Sr*Yo_Va*(L13 + (Yb_3 - Yb_4)*L14) + Yb_3*Yb_4*Ya_Ba*Yo_Va*(L15 + (Yb_3 - Yb_4)*L16) + ...
        Yo_o*Yo_Va*Ya_Sr*Yb_3*(L17 + (Yo_o - Yo_Va)*L18) + Yo_o*Yo_Va*Ya_Ba*Yb_3*(L19 + (Yo_o - Yo_Va)*L20) + ...
        Yo_o*Yo_Va*Ya_Sr*Yb_4*(L21 + (Yo_o - Yo_Va)*L22) + Yo_o*Yo_Va*Ya_Ba*Yb_4*(L23 + (Yo_o - Yo_Va)*L24);
    
    %% Complete Solution Model 
    G_soln = G_end - T*S + G_ex;
    
    %% Model no excess terms we call this the zeroth order model G_o in the paper
    G_o = G_end - T*S;

end

