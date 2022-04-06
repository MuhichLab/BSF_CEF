function d0 = solve_dref_end_only_2(guess,G_end_diff,X0,x_vals,d_chem)
    syms x y T
    syms gDiffA [1 2]
    syms gDiffB [1 2]
    syms gDiffC [1 2]
    syms gDiffD [1 2]
    
    
    gDiffa1     = guess(1,1);
    gDiffa2     = guess(2,1);    
    gDiffb1     = guess(3,1);
    gDiffb2     = guess(4,1);
    gDiffc1     = guess(5,1);
    gDiffc2     = guess(6,1);
    
    sub_me = [gDiffA1,gDiffA2,gDiffB1,gDiffB2,gDiffC1,gDiffC2,T];
    sub_value = [gDiffa1,gDiffa2,gDiffb1,gDiffb2,gDiffc1,gDiffc2,573.15];
    
    dG_dd_og = subs(sym(G_end_diff),sub_me,sub_value); % eV @ 300 C
    
    d0 = zeros(length(x_vals),1);
    for i=1:length(x_vals)
        dG_dd = subs(sym(dG_dd_og),x,x_vals(i));
        dG_dd = matlabFunction(simplify(expand(dG_dd)));
        coords = @(d) d_chem + dG_dd(d); % d00_val is the condition which d0 = 0
        dref = fsolve(coords,X0,optimset('FunValCheck', 'off','Display', 'off'));
        if dref<0 || ~isreal(dref) 
            d0(i)=0.001;
        else
            d0(i) = dref;
        end
    end    
end


