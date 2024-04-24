function [est_g7,I7] = My_est_g(y_w,MAX_y_k,GG,K_p,MSE,max_e1,K,G1,M,L_r,L_t,RR_W)
%y_w接收信号
%MAX_y_k从K个子载波内选出来的子载波数
%GG = G_t*G_r
%K_p选择子载波的个数
%MSE均方误差
%max_el = 1;
%K载波个数
%G1分组的个数
%M传输符号个数
%L_r接收端射频链路个数
%L_t反射端天线个数
%RR_W测量矩阵
    RR_WW1 = RR_W(1:M*L_r*L_t,1:GG);RR_WW2 = RR_W(1:M*L_r*L_t,1+GG:GG*2);
    RR_WW3 = RR_W(1:M*L_r*L_t,1+2*GG:GG*3);RR_WW4 = RR_W(1:M*L_r*L_t,1+3*GG:GG*4);
    RR_WW5 = RR_W(1:M*L_r*L_t,1+4*GG:GG*5);RR_WW6 = RR_W(1:M*L_r*L_t,1+5*GG:GG*6);
    RR_WW7 = RR_W(1:M*L_r*L_t,1+6*GG:GG*7);RR_WW8 = RR_W(1:M*L_r*L_t,1+7*GG:GG*8);
    RR_WW9 = RR_W(1:M*L_r*L_t,1+8*GG:GG*9);RR_WW10 = RR_W(1:M*L_r*L_t,1+9*GG:GG*10);
    RR_WW11 = RR_W(1:M*L_r*L_t,1+10*GG:GG*11);RR_WW12 = RR_W(1:M*L_r*L_t,1+11*GG:GG*12);
    RR_WW13 = RR_W(1:M*L_r*L_t,1+12*GG:GG*13);RR_WW14 = RR_W(1:M*L_r*L_t,1+13*GG:GG*14);
    RR_WW15 = RR_W(1:M*L_r*L_t,1+14*GG:GG*15);RR_WW16 = RR_W(1:M*L_r*L_t,1+15*GG:GG*16);
    Res16 = y_w;
    MAX_y_k5 = MAX_y_k;
    C5 = zeros(GG,K_p);
    I7 = [];
    P5 = zeros(GG,1);
    MSE2 = MSE;
    while(MSE2 > max_e1)
        for k_p = 1:K_p
            pp=ceil(MAX_y_k5(1,k_p)/(K/G1));
            RR_WW = RR_W(1:M*L_r*L_t,1+(pp-1)*GG:GG*pp);
            C5(:,k_p) =  RR_WW'* Res16(:,MAX_y_k5(1,k_p));
        end
        for m = 1:GG
            P5(m,1) = sum(abs(C5(m,:)));
        end
        pos7 = find(P5 == max(P5));%%找到最大的投影
        I7 = [I7 pos7(1)];%%%Gt*Gr; 
        est_gg7 = zeros(length(I7),K);
        res7 = zeros(1,K);
        pinv_RR_WW1 = pinv(RR_WW1(:,I7));pinv_RR_WW2 = pinv(RR_WW2(:,I7));
        pinv_RR_WW3 = pinv(RR_WW3(:,I7));pinv_RR_WW4 = pinv(RR_WW4(:,I7));
        pinv_RR_WW5 = pinv(RR_WW5(:,I7));pinv_RR_WW6 = pinv(RR_WW6(:,I7));
        pinv_RR_WW7 = pinv(RR_WW7(:,I7));pinv_RR_WW8 = pinv(RR_WW8(:,I7));
        pinv_RR_WW9 = pinv(RR_WW9(:,I7));pinv_RR_WW10 = pinv(RR_WW10(:,I7));
        pinv_RR_WW11 = pinv(RR_WW11(:,I7));pinv_RR_WW12 = pinv(RR_WW12(:,I7));
        pinv_RR_WW13 = pinv(RR_WW13(:,I7));pinv_RR_WW14 = pinv(RR_WW14(:,I7));
        pinv_RR_WW15 = pinv(RR_WW15(:,I7));pinv_RR_WW16 = pinv(RR_WW16(:,I7));
        for k = 1:K
            pp=ceil(k/(K/G1));
%             RR_WW = RR_W(1:M*L_r*L_t,1+(pp-1)*GG:GG*pp);
%             pinv_RR_WW = pinv(RR_WW(:,I7));
%             [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW,y_w(:,k),RR_WW(:,I7));
            if(pp==1)
                pinv_RR_WW22 = pinv_RR_WW1;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW1(:,I7));
            elseif(pp==2)
                pinv_RR_WW22 = pinv_RR_WW2;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW2(:,I7));
            elseif(pp==3)
                pinv_RR_WW22 = pinv_RR_WW3;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW3(:,I7));
            elseif(pp==4)
                pinv_RR_WW22 = pinv_RR_WW4;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW4(:,I7));
            elseif(pp==5)
                pinv_RR_WW22 = pinv_RR_WW5;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW5(:,I7));
            elseif(pp==6)
                pinv_RR_WW22 = pinv_RR_WW6;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW6(:,I7));
            elseif(pp==7)
                pinv_RR_WW22 = pinv_RR_WW7;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW7(:,I7));
            elseif(pp==8)
                pinv_RR_WW22 = pinv_RR_WW8;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW8(:,I7));
            elseif(pp==9)
                pinv_RR_WW22 = pinv_RR_WW9;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW9(:,I7));
            elseif(pp==10)
                pinv_RR_WW22 = pinv_RR_WW10;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW10(:,I7));
            elseif(pp==11)
                pinv_RR_WW22 = pinv_RR_WW11;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW11(:,I7));
            elseif(pp==12)
                pinv_RR_WW22 = pinv_RR_WW12;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW12(:,I7));
            elseif(pp==13)
                pinv_RR_WW22 = pinv_RR_WW13;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW13(:,I7));
            elseif(pp==14)
                pinv_RR_WW22 = pinv_RR_WW14;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW14(:,I7));
            elseif(pp==15)
                pinv_RR_WW22 = pinv_RR_WW15;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW15(:,I7));
            else
                pinv_RR_WW22 = pinv_RR_WW16;
                [Res16(:,k),  est_gg7(:,k)] = MyRes(pinv_RR_WW22,y_w(:,k),RR_WW16(:,I7));
            end
             res7(1,k) =  Res16(:,k)'*Res16(:,k);
        end
        MSE2 = sum(res7)/(K*L_r*M*L_t);
    end
    GG17 = zeros(1,length(I7));
    p_av7  = zeros(1,length(I7));
    for L1 = 1: length(I7)
        GG17(1,L1) = est_gg7(L1,:)* est_gg7(L1,:)'; 
        p_av7(1,L1) = GG17(1,L1)/K;        
    end
    PP7 = find(GG17 == max(GG17));
    EE7 = [];
    for L1 = 1:length(I7)
        if(p_av7(1,L1) >= 0.025*GG17(1,PP7)/K)
            EE7 = [EE7 L1];
        end
    end
    est_g7 = zeros(length(I7),K);
    if(MSE > max_e1)
        for k = 1:K
            est_g7(EE7,k) = est_gg7(EE7,k);
        end
    end
end