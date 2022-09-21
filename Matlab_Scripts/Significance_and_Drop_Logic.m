function [dp_terms,spot] = Significance_and_Drop_Logic...
    (guess_drop,dp_terms,num_terms,drop_error,X,Y,Z,x_exp,dd,muhg_o,Temp,G_dft,dG_dy,x_vals,...
        d0,spot_track)  
    
    
    % Determine Significance of 1,2, and 3 combinatiosn of terms
    original_error = drop_error;
    [drop_this, combos1, combos2, combos3] = drop_me_2(guess_drop,...
        X,Y,Z,x_exp,dd,-muhg_o,Temp,G_dft,dG_dy,x_vals,original_error,...
        d0,spot_track);

    drop_por1 = drop_this(1:22);
    drop_por2 = drop_this(23:253);
    drop_por3 = drop_this(254:1793);
 
    % Drop Logic
    if num_terms > 3
    
        % Find the least significant terms
        min_val1 = min(drop_por1(drop_por1>0));
        min_val2 = min(drop_por2(drop_por2>0));
        min_val3 = min(drop_por3(drop_por3>0));

        spot1 = find(drop_por1 == min_val1);
        spot2 = find(drop_por2 == min_val2);
        spot3 = find(drop_por3 == min_val3);

        combos1(spot1);
        combos2(spot2,:);
        combos3(spot3,:);

        fprintf('Combo 1 are these terms and this Value: %d       : %.3e \n',combos1(spot1),min_val1);
        fprintf('Combo 2 are these terms and this Value: %d %d   : %.3e \n',combos2(spot2,:),min_val2);
        fprintf('Combo 3 are these terms and this Value: %d %d %d : %.3e \n',combos3(spot3,:),min_val3);

       % Logic to compare htese terms against eachother
       if (any(combos3(spot3,:) == combos2(spot2,1))) && ...
                (any(combos3(spot3,:) == combos2(spot2,2))) && min_val2 < min_val3
            new_min = min(drop_por1(combos2(spot2,:)));
            spot = find(drop_por1 == new_min);
            dp_terms(spot,1) = 0;
            dp_terms(spot+22,1) = 0;
            dp_terms(spot+44,1) = 0;
 
       elseif (any(combos3(spot3,:) == combos2(spot2,1))) && ...
                (any(combos3(spot3,:) == combos2(spot2,2))) && min_val2 > min_val3
            spot = setdiff(combos3(spot3,:),combos2(spot2,:));
            dp_terms(spot,1) = 0;
            dp_terms(spot+22,1) = 0;
            dp_terms(spot+44,1) = 0;
  
       else

            if min_val2 < min_val3
                new_min = min(drop_por1(combos2(spot2,:)));
                spot = find(drop_por1 == new_min);
                dp_terms(spot,1) = 0;
                dp_terms(spot+22,1) = 0;
                dp_terms(spot+44,1) = 0;

            elseif min_val3 < min_val2
                new_min = min(drop_por1(combos3(spot3,:)));
                spot = find(drop_por1 == new_min);
                dp_terms(spot,1) = 0;
                dp_terms(spot+22,1) = 0;
                dp_terms(spot+44,1) = 0;

            end
       end

    elseif num_terms <= 3
        min_val1 = min(drop_por1(drop_por1>0));
        spot1 = find(drop_por1 == min_val1);
        spot = spot1;
        dp_terms(spot,1) = 0;
        dp_terms(spot+22,1) = 0;
        dp_terms(spot+44,1) = 0;

    end
end