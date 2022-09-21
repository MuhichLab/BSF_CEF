function error = min_me_end_bad(guess,x_exp,dd,muhg_o,Temp,dG_dy,x_vals,X0,d_chem)
    %lambda = 0.5;
    
    A1          = guess(1,1)*0;
    A2          = guess(2,1)*0;
    A3          = guess(3,1)*0;
    A4          = guess(4,1)*0;
    A5          = guess(5,1)*0;
    A6          = guess(6,1)*0;
    A7          = guess(7,1)*0;
    A8          = guess(8,1)*0;
    A9          = guess(9,1)*0;
    A10         = guess(10,1)*0;
    A11         = guess(11,1)*0;
    A12         = guess(12,1)*0;
    A13         = guess(13,1)*0;
    A14         = guess(14,1)*0;
    A15         = guess(15,1)*0;
    A16         = guess(16,1)*0;
    A17         = guess(17,1)*0;
    A18         = guess(18,1)*0;
    A19         = guess(19,1)*0;
    A20         = guess(20,1)*0;
    A22         = guess(21,1)*0;
    A24         = guess(22,1)*0;
    B1          = guess(23,1)*0;
    B2          = guess(24,1)*0;
    B3          = guess(25,1)*0;
    B4          = guess(26,1)*0;
    B5          = guess(27,1)*0;    
    B6          = guess(28,1)*0;
    B7          = guess(29,1)*0;
    B8          = guess(30,1)*0;
    B9          = guess(31,1)*0;
    B10         = guess(32,1)*0;
    B11         = guess(33,1)*0;
    B12         = guess(34,1)*0;
    B13         = guess(35,1)*0;
    B14         = guess(36,1)*0;
    B15         = guess(37,1)*0;
    B16         = guess(38,1)*0;
    B17         = guess(39,1)*0;
    B18         = guess(40,1)*0;
    B19         = guess(41,1)*0;
    B20         = guess(42,1)*0;
    B22         = guess(43,1)*0;
    B24         = guess(44,1)*0;
    C1          = guess(45,1)*0;
    C2          = guess(46,1)*0;
    C3          = guess(47,1)*0;
    C4          = guess(48,1)*0;
    C5          = guess(49,1)*0;
    C6          = guess(50,1)*0;
    C7          = guess(51,1)*0;
    C8          = guess(52,1)*0;
    C9          = guess(53,1)*0;
    C10         = guess(54,1)*0;
    C11         = guess(55,1)*0;
    C12         = guess(56,1)*0;
    C13         = guess(57,1)*0;
    C14         = guess(58,1)*0;
    C15         = guess(59,1)*0;
    C16         = guess(60,1)*0
    C17         = guess(61,1)*0;
    C18         = guess(62,1)*0;
    C19         = guess(63,1)*0;
    C20         = guess(64,1)*0;
    C22         = guess(65,1)*0;
    C24         = guess(66,1)*0;
    dref        = guess(67,1);
    gDiffA1     = guess(68,1);
    gDiffA2     = guess(69,1);    
    gDiffB1     = guess(70,1);
    gDiffB2     = guess(71,1);
    gDiffC1     = guess(72,1);
    gDiffC2     = guess(73,1);



    
    error1 = 0;
    for i = 1:length(x_exp);
        value1 = dG_dy(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                    A15,A16,A17,A18,A19,A20,A22,A24,B1,B2,B3,B4,B5,B6,...
                    B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,B18,B19,B20,B22,...
                    B24,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,...
                    C17,C18,C19,C20,C22,C24,Temp(i),dref,gDiffA1,gDiffA2,gDiffB1,...
                    gDiffB2,gDiffC1,gDiffC2,x_exp(i),dd(i));
                
        if  ~isreal(value1)
            value1 = 1e6;
        end
        error1 = error1 + (muhg_o(i) - value1)^2;
    end
    
        
    error = error1;
end

