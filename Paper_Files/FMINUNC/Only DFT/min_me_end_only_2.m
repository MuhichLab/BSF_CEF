function error = min_me_end_only_2(guess,X,Y,Z,x_exp,dd,X0,x_vals,muhg_o,...
    Temp,G_end_dft,G_end_diff,d_chem)
    % Don't forget muhg_o was passed in as -muhg_o
    
    lambda = 1; % how much to multiply DFT error by
    
    gDiffA1     = guess(1,1);
    gDiffA2     = guess(2,1);    

    
    error2 = 0;
    for j = 1:length(X)
        value2 = G_end_dft(gDiffA1,gDiffA2,X(j),Y(j));
        error2 = error2 + lambda*(Z(j) - value2)^2;
    end
        
    error = error2;


end

