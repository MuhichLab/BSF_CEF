function error = min_me_end_only_2(guess,X,Y,Z,x_exp,dd,X0,x_vals,muhg_o,...
    Temp,G_end_dft,G_end_diff,d_chem)
    % Don't forget muhg_o was passed in as -muhg_o
    
    lambda = 1; % how much to multiply DFT error by
    
    gDiffA1     = guess(1,1);
    gDiffA2     = guess(2,1);    
    gDiffB1     = guess(3,1);
    gDiffB2     = guess(4,1);
    gDiffC1     = guess(5,1);
    gDiffC2     = guess(6,1);


    
    
    d0_end_only = solve_dref_end_only_2(guess,G_end_diff,X0,x_vals,d_chem); % need to use G_soln because it does not have th edref variable
     
    
    error1 = 0;
    for i = 1:length(x_exp)
        for k = 1:length(x_vals)
            if x_exp(i) == x_vals(k)
                dref = d0_end_only(k);
                break
            end
        end
        value1 = G_end_diff(Temp(i),gDiffA1,gDiffA2,gDiffB1,...
                    gDiffB2,gDiffC1,gDiffC2,x_exp(i),dd(i)+dref);
        if  ~isreal(value1)
            value1 = 1e6;
        end
        error1 = error1 + (muhg_o(i) - value1)^2;
    end
    
    error2 = 0;
    for j = 1:length(X)
        value2 = G_end_dft(gDiffA1,gDiffA2,X(j),Y(j));
        error2 = error2 + lambda*(Z(j) - value2)^2;
    end
        
    error = error1 + error2;


end

