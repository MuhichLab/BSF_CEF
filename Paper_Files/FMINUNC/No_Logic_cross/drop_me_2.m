function [drop, combos1, combos2, combos3] = drop_me_2(guess,X,Y,Z,x_exp,dd,...
    muhg_o,Temp,G_dft,dG_dy,x_vals,original_error,d0,spot_track)
    orig = guess;
    combos1 = nchoosek(1:22,1);
    combos2 = nchoosek(1:22,2);
    combos3 = nchoosek(1:22,3);
    
    
    %lambda = 0.5;
    count = 1;
    drop = zeros(length(combos1)+length(combos2)+length(combos3),1);
    % don't check terms that have already been droppped
    for k = 1:length(combos1)
        if ismember(combos1(k),spot_track)
            drop(count,1) = 0;
            count = count + 1;
            continue
        end
        guess([combos1(k) combos1(k)+22 combos1(k)+44],1) = guess([combos1(k) combos1(k)+22 combos1(k)+44],1)*(1.1);
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
        B1          = guess(23,1);
        B2          = guess(24,1);
        B3          = guess(25,1);
        B4          = guess(26,1);
        B5          = guess(27,1);
        B6          = guess(28,1);
        B7          = guess(29,1);
        B8          = guess(30,1);
        B9          = guess(31,1);
        B10         = guess(32,1);
        B11         = guess(33,1);
        B12         = guess(34,1);
        B13         = guess(35,1);
        B14         = guess(36,1);
        B15         = guess(37,1);
        B16         = guess(38,1);
        B17         = guess(39,1);
        B18         = guess(40,1);
        B19         = guess(41,1);
        B20         = guess(42,1);
        B22         = guess(43,1);
        B24         = guess(44,1);
        C1          = guess(45,1);
        C2          = guess(46,1);
        C3          = guess(47,1);
        C4          = guess(48,1);
        C5          = guess(49,1);
        C6          = guess(50,1);
        C7          = guess(51,1);
        C8          = guess(52,1);
        C9          = guess(53,1);
        C10         = guess(54,1);
        C11         = guess(55,1);
        C12         = guess(56,1);
        C13         = guess(57,1);
        C14         = guess(58,1);
        C15         = guess(59,1);
        C16         = guess(60,1);
        C17         = guess(61,1);
        C18         = guess(62,1);
        C19         = guess(63,1);
        C20         = guess(64,1);
        C22         = guess(65,1);
        C24         = guess(66,1);
        dref        = guess(67,1);

        error1 = 0;
        for i = 1:length(x_exp)

            value1 = dG_dy(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                A15,A16,A17,A18,A19,A20,A22,A24,B1,B2,B3,B4,B5,B6,...
                B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,B18,B19,B20,B22,...
                B24,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,...
                C17,C18,C19,C20,C22,C24,Temp(i),x_exp(i),(dd(i)+dref));
            error1 = error1 + (muhg_o(i) - value1)^2; % Note! (-)muhg_o is passed in
        end

        error2 = 0;
        for j = 1:length(X)
            value2 = G_dft(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                A15,A16,A17,A18,A19,A20,A22,A24,X(j),Y(j));
            error2 = error2 + (Z(j) - value2)^2;
        end

        drop(count,1) = abs(((error1 + error2) - original_error)/original_error)*100;
        guess = orig;
        count = count + 1;
    end
    
    for k = 1:length(combos2)
        if ismember(combos2(k,1),spot_track) || ismember(combos2(k,2),spot_track)
            drop(count,1) = 0;
            count = count + 1;
            continue
        end
        guess([combos2(k,:) combos2(k,:)+22 combos2(k,:)+44],1) = guess([combos2(k,:) combos2(k,:)+22 combos2(k,:)+44],1)*(1.1);
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
        B1          = guess(23,1);
        B2          = guess(24,1);
        B3          = guess(25,1);
        B4          = guess(26,1);
        B5          = guess(27,1);
        B6          = guess(28,1);
        B7          = guess(29,1);
        B8          = guess(30,1);
        B9          = guess(31,1);
        B10         = guess(32,1);
        B11         = guess(33,1);
        B12         = guess(34,1);
        B13         = guess(35,1);
        B14         = guess(36,1);
        B15         = guess(37,1);
        B16         = guess(38,1);
        B17         = guess(39,1);
        B18         = guess(40,1);
        B19         = guess(41,1);
        B20         = guess(42,1);
        B22         = guess(43,1);
        B24         = guess(44,1);
        C1          = guess(45,1);
        C2          = guess(46,1);
        C3          = guess(47,1);
        C4          = guess(48,1);
        C5          = guess(49,1);
        C6          = guess(50,1);
        C7          = guess(51,1);
        C8          = guess(52,1);
        C9          = guess(53,1);
        C10         = guess(54,1);
        C11         = guess(55,1);
        C12         = guess(56,1);
        C13         = guess(57,1);
        C14         = guess(58,1);
        C15         = guess(59,1);
        C16         = guess(60,1);
        C17         = guess(61,1);
        C18         = guess(62,1);
        C19         = guess(63,1);
        C20         = guess(64,1);
        C22         = guess(65,1);
        C24         = guess(66,1);
        dref        = guess(67,1);
        
        error1 = 0;
        for i = 1:length(x_exp);

            value1 = dG_dy(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                A15,A16,A17,A18,A19,A20,A22,A24,B1,B2,B3,B4,B5,B6,...
                B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,B18,B19,B20,B22,...
                B24,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,...
                C17,C18,C19,C20,C22,C24,Temp(i),x_exp(i),(dd(i)+dref));
            error1 = error1 + (muhg_o(i) - value1)^2; % Note! (-)muhg_o is passed in
        end

        error2 = 0;
        for j = 1:length(X);
            value2 = G_dft(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                A15,A16,A17,A18,A19,A20,A22,A24,X(j),Y(j));
            error2 = error2 + (Z(j) - value2)^2;
        end

        drop(count,1) = abs(((error1 + error2) - original_error)/original_error)*100;
        guess = orig;
        count = count + 1;
    end

    for k = 1:length(combos3)
        if ismember(combos3(k,1),spot_track) || ismember(combos3(k,2),spot_track) || ismember(combos3(k,3),spot_track)
            drop(count,1) = 0;
            count = count + 1;
            continue
        end  
        guess([combos3(k,:) combos3(k,:)+22 combos3(k,:)+44],1) = guess([combos3(k,:) combos3(k,:)+22 combos3(k,:)+44],1)*(1.1);
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
        B1          = guess(23,1);
        B2          = guess(24,1);
        B3          = guess(25,1);
        B4          = guess(26,1);
        B5          = guess(27,1);
        B6          = guess(28,1);
        B7          = guess(29,1);
        B8          = guess(30,1);
        B9          = guess(31,1);
        B10         = guess(32,1);
        B11         = guess(33,1);
        B12         = guess(34,1);
        B13         = guess(35,1);
        B14         = guess(36,1);
        B15         = guess(37,1);
        B16         = guess(38,1);
        B17         = guess(39,1);
        B18         = guess(40,1);
        B19         = guess(41,1);
        B20         = guess(42,1);
        B22         = guess(43,1);
        B24         = guess(44,1);
        C1          = guess(45,1);
        C2          = guess(46,1);
        C3          = guess(47,1);
        C4          = guess(48,1);
        C5          = guess(49,1);
        C6          = guess(50,1);
        C7          = guess(51,1);
        C8          = guess(52,1);
        C9          = guess(53,1);
        C10         = guess(54,1);
        C11         = guess(55,1);
        C12         = guess(56,1);
        C13         = guess(57,1);
        C14         = guess(58,1);
        C15         = guess(59,1);
        C16         = guess(60,1);
        C17         = guess(61,1);
        C18         = guess(62,1);
        C19         = guess(63,1);
        C20         = guess(64,1);
        C22         = guess(65,1);
        C24         = guess(66,1);
        dref        = guess(67,1);
        error1 = 0;
        for i = 1:length(x_exp);

            value1 = dG_dy(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                A15,A16,A17,A18,A19,A20,A22,A24,B1,B2,B3,B4,B5,B6,...
                B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,B18,B19,B20,B22,...
                B24,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,...
                C17,C18,C19,C20,C22,C24,Temp(i),x_exp(i),(dd(i)+dref));
            error1 = error1 + (muhg_o(i) - value1)^2; % Note! (-)muhg_o is passed in
        end

        error2 = 0;
        for j = 1:length(X);
            value2 = G_dft(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                A15,A16,A17,A18,A19,A20,A22,A24,X(j),Y(j));
            error2 = error2 + (Z(j) - value2)^2;
        end

        drop(count,1) = abs(((error1 + error2) - original_error)/original_error)*100;
        guess = orig;
        count = count + 1;
    end
    
end