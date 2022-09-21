clear;
clc;

syms x y

Ya_Sr = 1-x;
Ya_Ba = x;
Yb_3 = 2*y;     
Yb_4 = 1-2*y;   
Yo_o = 1-y/3;   
Yo_Va = y/3;   


%% L Compares
l1 = Ya_Sr.*Ya_Ba.*Yb_3.*Yo_o;
l2 = Ya_Sr.*Ya_Ba.*Yb_3.*Yo_o.*(Ya_Sr - Ya_Ba);
l3 = Ya_Sr.*Ya_Ba.*Yb_3.*Yo_Va;
l4 = Ya_Sr.*Ya_Ba.*Yb_3.*Yo_Va.*(Ya_Sr - Ya_Ba);
l5 = Ya_Sr.*Ya_Ba.*Yb_4.*Yo_o;
l6 = Ya_Sr.*Ya_Ba.*Yb_4.*Yo_o.*(Ya_Sr - Ya_Ba);
l7 = Ya_Sr.*Ya_Ba.*Yb_4.*Yo_Va;
l8 = Ya_Sr.*Ya_Ba.*Yb_4.*Yo_Va.*(Ya_Sr - Ya_Ba);
l9 = Yb_3.*Yb_4.*Ya_Sr.*Yo_o;
l10 = Yb_3.*Yb_4.*Ya_Sr.*Yo_o.*(Yb_3 - Yb_4);
l11 = Yb_3.*Yb_4.*Ya_Ba.*Yo_o;
l12 = Yb_3.*Yb_4.*Ya_Ba.*Yo_o.*(Yb_3 - Yb_4);
l13 = Yb_3.*Yb_4.*Ya_Sr.*Yo_Va;
l14 = Yb_3.*Yb_4.*Ya_Sr.*Yo_Va.*(Yb_3 - Yb_4);
l15 = Yb_3.*Yb_4.*Ya_Ba.*Yo_Va;
l16 = Yb_3.*Yb_4.*Ya_Ba.*Yo_Va.*(Yb_3 - Yb_4);
l17 = Yo_o.*Yo_Va.*Ya_Sr.*Yb_3;
l18 = Yo_o.*Yo_Va.*Ya_Sr.*Yb_3.*(Yo_o - Yo_Va);
l19 = Yo_o.*Yo_Va.*Ya_Ba.*Yb_3;
l20 = Yo_o.*Yo_Va.*Ya_Ba.*Yb_3.*(Yo_o - Yo_Va);
l21 = Yo_o.*Yo_Va.*Ya_Sr.*Yb_4;
l22 = Yo_o.*Yo_Va.*Ya_Sr.*Yb_4.*(Yo_o - Yo_Va);
l23 = Yo_o.*Yo_Va.*Ya_Ba.*Yb_4;
l24 = Yo_o.*Yo_Va.*Ya_Ba.*Yb_4.*(Yo_o - Yo_Va);

l9l21 = l9 + l21;
l11l23 = l11 + l23;

eqns = [l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18 l19 l20 l21 l22 l23 l24 l9l21 l11l23];
 %%

 
Title = cell(26,1);
Title{1} = 'L1';
Title{2} = 'L2';
Title{3} = 'L3';
Title{4} = 'L4';
Title{5} = 'L5';
Title{6} = 'L6';
Title{7} = 'L7';
Title{8} = 'L8';
Title{9} = 'L9';
Title{10} = 'L10';
Title{11} = 'L11';
Title{12} = 'L12';
Title{13} = 'L13';
Title{14} = 'L14';
Title{15} = 'L15';
Title{16} = 'L16';
Title{17} = 'L17';
Title{18} = 'L18';
Title{19} = 'L19';
Title{20} = 'L20';
Title{21} = 'L21';
Title{22} = 'L22';
Title{23} = 'L23';
Title{24} = 'L24';
Title{25} = 'L9 + L21';
Title{26} = 'L11 + L23';

Legend = cell(6,1);
Legend{1} = 'x = 0.00';
Legend{2} = 'x = 0.20';
Legend{3} = 'x = 0.40';
Legend{4} = 'x = 0.60';
Legend{5} = 'x = 0.80';
Legend{6} = 'x = 1.00';

linS = {'-ok','-^k','-*k','-xk','-dk','-vk'};
figure
 for sub1 = 19:24
    eqn = eqns(sub1);
    subplot(3,2,sub1-18)
    count = 1;
    hold on
    for sub2 = 0:0.2:1
        eqn2 = subs(eqn,x,sub2);
        fplot(eqn2,linS{count});
%         if sub1 < 7;
%             set(gca,'xtick',[])
%         else
%             xlabel('\delta')
%         end
        count = count+1;
    end
    axis tight
    legend(Legend);
    title(Title{sub1})
    xlabel('\delta')
    xlim([0 0.5])
    %ylim([-.2 .2])
    set(gca,'FontSize',20)
    box on
    hold off
 end

 %%

n = 1;

count = 1;
eqn = eqns(n);
linS = {'-ok','-^k','-*k','-xk','-dk','-vk'};
figure
hold on
for sub2 = 0:0.2:1
    eqn2 = subs(eqn,x,sub2);
    fplot(eqn2,linS{count});
    count = count + 1;
end
axis tight
title(Title{n})
xlim([0 0.5])
%ylim([0 0.15])
xlabel('\delta')
set(gca,'FontSize',20)
box on
legend(Legend);
%savefig(Title{n})

%% Surface Plto Comapre
fntsz = 20;

[xx,yy] = meshgrid((linspace(0,1)),(linspace(0,0.5)));

n = 7;
m = 11;

eqn1 = matlabFunction(eqns(n));
eqn2 = matlabFunction(eqns(m));

figure
hold on
surf(xx,yy,eqn1(xx,yy),'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
surf(xx,yy,eqn2(xx,yy),'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none')
legend(['L',num2str(n)],['L',num2str(m)],'FontSize',fntsz)
title([Title{n},' & ',Title{m},' Compare'])
xlabel('x','FontSize',fntsz)
ylabel('\delta','FontSize',fntsz)


%% Looking at A + B*T + C*T*ln(T)

% n = 9;
% t = 800;
% sub_x = 0.2
% syms T A B C
% eqn = eqns(n).*(A + B.*T + C.*T.*log(T));
% 
% for a = 0:0.5:1
%     for b = 0:0.5:1
%         for c = 0:0.5:1
%             figure
%             fplot(subs(eqn,[x T A B C],[sub_x t a b c]))
%             xlim([0 0.5])
%             xlabel('\delta')
%             
%         end
%     end
% end

 %% What I was lookign at
% 
% n = 13;
% m = 17;
% 
% eqn1 = simplify((eqns(n) - eqns(m))/eqns(m));
% eqn2 = simplify((eqns(n) - eqns(m))/eqns(n));
% 
% figure
% hold on
% fplot(Yb_4/Yo_o)
% %fplot(Yo_o)
% xlim([0 0.5])
% ylim([-4 4])
% xlabel('\delta')
% hold off

%% derivaitives
Title_prime = cell(26,1);
Title_prime{1} = 'L1^\prime';
Title_prime{2} = 'L2^\prime';
Title_prime{3} = 'L3^\prime';
Title_prime{4} = 'L4^\prime';
Title_prime{5} = 'L5^\prime';
Title_prime{6} = 'L6^\prime';
Title_prime{7} = 'L7^\prime';
Title_prime{8} = 'L8^\prime';
Title_prime{9} = 'L9^\prime';
Title_prime{10} = 'L10^\prime';
Title_prime{11} = 'L11^\prime';
Title_prime{12} = 'L12^\prime';
Title_prime{13} = 'L13^\prime';
Title_prime{14} = 'L14^\prime';
Title_prime{15} = 'L15^\prime';
Title_prime{16} = 'L16^\prime';
Title_prime{17} = 'L17^\prime';
Title_prime{18} = 'L18^\prime';
Title_prime{19} = 'L19^\prime';
Title_prime{20} = 'L20^\prime';
Title_prime{21} = 'L21^\prime';
Title_prime{22} = 'L22^\prime';
Title_prime{23} = 'L23^\prime';
Title_prime{24} = 'L24^\prime';
Title_prime{25} = 'L9^\prime + L21^\prime';
Title_prime{26} = 'L11^\prime + L23^\prime';

%figure
%  for sub1 = 17:24
%     eqn = diff(eqns(sub1),y);
%     subplot(4,2,sub1-16)
%     count = 1;
%     hold on
%     for sub2 = 0:0.25:1
%         eqn2 = subs(eqn,x,sub2);
%         fplot(eqn2);
%         count = count+1;
%     end
%     %legend(Legend);
%     title(Title_prime{sub1})
%     xlim([0 0.5])
%     %ylim([-2 2])
%     xlabel('\delta')
%     hold off
%  end
Legend = cell(6,1);
Legend{1} = 'x = 0.00';
Legend{2} = 'x = 0.20';
Legend{3} = 'x = 0.40';
Legend{4} = 'x = 0.60';
Legend{5} = 'x = 0.80';
Legend{6} = 'x = 1.00';

linS = {'-ok','-^k','-*k','-xk','-dk','-vk'};
figure
 for sub1 = 13:18
    eqn = diff(eqns(sub1),y);
    subplot(3,2,sub1-12)
    count = 1;
    hold on
    for sub2 = 0:0.2:1
        eqn2 = subs(eqn,x,sub2);
        fplot(eqn2,linS{count});
%         if sub1 < 7;
%             set(gca,'xtick',[])
%         else
%             xlabel('\delta')
%         end
        count = count+1;
    end
    axis tight
    legend(Legend);
    title(Title_prime{sub1})
    xlabel('\delta')
    xlim([0 0.5])
    %ylim([-.2 .2])
    set(gca,'FontSize',20)
    box on
    hold off
 end
 %%

n = 26;

count = 1;
linS = {'-ok','-^k','-*k','-xk','-dk','-vk'};
eqn = diff(eqns(n),y)
figure
hold on
for sub2 = 0:0.2:1
    eqn2 = subs(eqn,x,sub2);
    fplot(eqn2,linS{count});
    count = count+1;
end

axis tight
title(Title_prime{n})
xlim([0 0.5])
%ylim([0 0.15])
xlabel('\delta')
set(gca,'FontSize',20)
box on
%legend(Legend);
savefig(Title_prime{n})

%% Surface Plto Comapre
fntsz = 20;

[xx,yy] = meshgrid((linspace(0,1)),(linspace(0,0.5)));

n = 21;
m = 22;

eqn1 = matlabFunction(diff(eqns(n),y));
eqn2 = matlabFunction(diff(eqns(m),y));

figure
hold on
surf(xx,yy,eqn1(xx,yy),'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
surf(xx,yy,eqn2(xx,yy),'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none')
legend(['L',num2str(n)],['L',num2str(m)],'FontSize',fntsz)
title([Title_prime{n},' & ',Title_prime{m},' Compare'])
xlabel('x','FontSize',fntsz)
ylabel('\delta','FontSize',fntsz)


