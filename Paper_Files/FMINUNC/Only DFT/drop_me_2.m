function [drop, combos1, combos2, combos3] = drop_me_2(guess,X,Y,Z,x_exp,dd,...
    muhg_o,Temp,G_dft,dG_dy,x_vals,original_error,spot_track)
    orig = guess;
    combos1 = nchoosek(1:22,1);
    combos2 = nchoosek(1:22,2);
    combos3 = nchoosek(1:22,3);
    
    lambda = 1; % how much to multiply DFT error by
    count = 1;
    drop = zeros(length(combos1)+length(combos2)+length(combos3),1);
    % don't check terms that have already been droppped
    for k = 1:length(combos1)
        if ismember(combos1(k),spot_track)
            drop(count,1) = 0;
            count = count + 1;
            continue
        end
        guess([combos1(k)],1) = guess([combos1(k)],1)*(1.1);
        A1          = guess(1,1);
        A2          = guess(2,1);
        A3          = guess(3,1);
        A4          = guess(4,1);
        A5          = guess(5,1);
        A6          = guess(6,1);
        A7          = guess(7,1);
        A8          = guess(8,1);
        A9          = guess(9,1);
        A10         = guess(10,1);
        A11         = guess(11,1);
        A12         = guess(12,1);
        A13         = guess(13,1);
        A14         = guess(14,1);
        A15         = guess(15,1);
        A16         = guess(16,1);
        A17         = guess(17,1);
        A18         = guess(18,1);
        A19         = guess(19,1);
        A20         = guess(20,1);
        A22         = guess(21,1);
        A24         = guess(22,1);


        error2 = 0;
        for j = 1:length(X)
            value2 = G_dft(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                A15,A16,A17,A18,A19,A20,A22,A24,X(j),Y(j));
            error2 = error2 + lambda*(Z(j) - value2)^2;
        end

        drop(count,1) = abs(((error2) - original_error)/original_error)*100;
        guess = orig;
        count = count + 1;
    end
    
    for k = 1:length(combos2)
        if ismember(combos2(k,1),spot_track) || ismember(combos2(k,2),spot_track)
            drop(count,1) = 0;
            count = count + 1;
            continue
        end
        guess([combos2(k,:)],1) = guess([combos2(k,:)],1)*(1.1);
        A1          = guess(1,1);
        A2          = guess(2,1);
        A3          = guess(3,1);
        A4          = guess(4,1);
        A5          = guess(5,1);
        A6          = guess(6,1);
        A7          = guess(7,1);
        A8          = guess(8,1);
        A9          = guess(9,1);
        A10         = guess(10,1);
        A11         = guess(11,1);
        A12         = guess(12,1);
        A13         = guess(13,1);
        A14         = guess(14,1);
        A15         = guess(15,1);
        A16         = guess(16,1);
        A17         = guess(17,1);
        A18         = guess(18,1);
        A19         = guess(19,1);
        A20         = guess(20,1);
        A22         = guess(21,1);
        A24         = guess(22,1);

        error2 = 0;
        for j = 1:length(X);
            value2 = G_dft(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                A15,A16,A17,A18,A19,A20,A22,A24,X(j),Y(j));
            error2 = error2 + lambda*(Z(j) - value2)^2;
        end

        drop(count,1) = abs(((error2) - original_error)/original_error)*100;
        guess = orig;
        count = count + 1;
    end

    for k = 1:length(combos3)
        if ismember(combos3(k,1),spot_track) || ismember(combos3(k,2),spot_track) || ismember(combos3(k,3),spot_track)
            drop(count,1) = 0;
            count = count + 1;
            continue
        end  
        guess([combos3(k,:)],1) = guess([combos3(k,:)],1)*(1.1);
        A1          = guess(1,1);
        A2          = guess(2,1);
        A3          = guess(3,1);
        A4          = guess(4,1);
        A5          = guess(5,1);
        A6          = guess(6,1);
        A7          = guess(7,1);
        A8          = guess(8,1);
        A9          = guess(9,1);
        A10         = guess(10,1);
        A11         = guess(11,1);
        A12         = guess(12,1);
        A13         = guess(13,1);
        A14         = guess(14,1);
        A15         = guess(15,1);
        A16         = guess(16,1);
        A17         = guess(17,1);
        A18         = guess(18,1);
        A19         = guess(19,1);
        A20         = guess(20,1);
        A22         = guess(21,1);
        A24         = guess(22,1);

        error2 = 0;
        for j = 1:length(X);
            value2 = G_dft(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                A15,A16,A17,A18,A19,A20,A22,A24,X(j),Y(j));
            error2 = error2 + lambda*(Z(j) - value2)^2;
        end

        drop(count,1) = abs(((error2) - original_error)/original_error)*100;
        guess = orig;
        count = count + 1;
    end
    
end