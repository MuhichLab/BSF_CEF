function d0_plots(guess,G_soln)

    a1          = guess(1,1);
    a2          = guess(2,1);
    a3          = guess(3,1);
    a4          = guess(4,1);
    a5          = guess(5,1);
    a6          = guess(6,1);
    a7          = guess(7,1);
    a8          = guess(8,1);
    a9          = guess(9,1);
    a10         = guess(10,1);
    a11         = guess(11,1);
    a12         = guess(12,1);
    a13         = guess(13,1);
    a14         = guess(14,1);
    a15         = guess(15,1);
    a16         = guess(16,1);
    a17         = guess(17,1);
    a18         = guess(18,1);
    a19         = guess(19,1);
    a20         = guess(20,1);
    a22         = guess(21,1);
    a24         = guess(22,1);
    b1          = guess(23,1);
    b2          = guess(24,1);
    b3          = guess(25,1);
    b4          = guess(26,1);
    b5          = guess(27,1);
    b6          = guess(28,1);
    b7          = guess(29,1);
    b8          = guess(30,1);
    b9          = guess(31,1);
    b10         = guess(32,1);
    b11         = guess(33,1);
    b12         = guess(34,1);
    b13         = guess(35,1);
    b14         = guess(36,1);
    b15         = guess(37,1);
    b16         = guess(38,1);
    b17         = guess(39,1);
    b18         = guess(40,1);
    b19         = guess(41,1);
    b20         = guess(42,1);
    b22         = guess(43,1);
    b24         = guess(44,1);
    c1          = guess(45,1);
    c2          = guess(46,1);
    c3          = guess(47,1);
    c4          = guess(48,1);
    c5          = guess(49,1);
    c6          = guess(50,1);
    c7          = guess(51,1);
    c8          = guess(52,1);
    c9          = guess(53,1);
    c10         = guess(54,1);
    c11         = guess(55,1);
    c12         = guess(56,1);
    c13         = guess(57,1);
    c14         = guess(58,1);
    c15         = guess(59,1);
    c16         = guess(60,1);
    c17         = guess(61,1);
    c18         = guess(62,1);
    c19         = guess(63,1);
    c20         = guess(64,1);
    c22         = guess(65,1);
    c24         = guess(66,1);
  
    
    syms A [1 24] real
    syms B [1 24] real
    syms C [1 24] real
    syms y real 
    syms T real
    syms x real 
    sub_me = [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,...
                    A15,A16,A17,A18,A19,A20,A22,A24,B1,B2,B3,B4,B5,B6,...
                    B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,B18,B19,B20,B22,...
                    B24,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,...
                    C17,C18,C19,C20,C22,C24,T];
    sub_value = [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,...
                    a15,a16,a17,a18,a19,a20,a22,a24,b1,b2,b3,b4,b5,b6,...
                    b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b22,...
                    b24,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,...
                    c17,c18,c19,c20,c22,c24,573.15];
    
    guess = transpose(guess);
    
    
    dG_dd = diff(subs(G_soln,[sub_me],[sub_value]),y); % eV @ 300 C
    dG_dd = matlabFunction(dG_dd);
    
    figure
    hold on
    for x_t = [0 0.05 0.1 0.15 0.2];
        coords = @(d) -0.6274 + dG_dd(x_t,d);
        fplot(coords,[-0.1 0.6])
    end
    legend('x = 0','x = 0.05','x = 0.1','x = 0.15','x = 0.2')
    hold off
end