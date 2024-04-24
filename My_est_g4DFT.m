function [est_g14,I14] = My_est_g4DFT(y_w,MAX_y_k,GG,K_p,MSE,max_e1,K,G2,M,L_r,L_t,RR_W1)
    RR_W10 = RR_W1(1:M*L_r*L_t,1:GG);
    RR_W11 = RR_W1(1:M*L_r*L_t,1+GG:GG*2);
    RR_W12 = RR_W1(1:M*L_r*L_t,1+2*GG:GG*3);
    RR_W13 = RR_W1(1:M*L_r*L_t,1+3*GG:GG*4);
    Res14 = y_w;
    MAX_y_k14 = MAX_y_k;
    C14 = zeros(GG,K_p);
    I14 = [];
    P14 = zeros(GG,1);
    MSE14 = MSE;
   while(length(I14) <= 20)
        for k_p = 1:K_p
            pp = ceil(MAX_y_k14(1,k_p)/(K/G2));
            RR_WW14 = RR_W1(1:M*L_r*L_t,1+(pp-1)*GG:GG*pp);
            C14(:,k_p) =  RR_WW14'* Res14(:,MAX_y_k14(1,k_p));
        end
        for m = 1:GG
            P14(m,1) = sum(abs(C14(m,:)));
        end
        pos14 = find(P14 == max(P14));%%找到最大的投影
        I14 = [I14 pos14(1)];%%%Gt*Gr; 
        est_gg14 = zeros(length(I14),K);
        res14 = zeros(1,K);
        pinv_RR_WW10 = pinv(RR_W10(:,I14)); 
        pinv_RR_WW11 = pinv(RR_W11(:,I14)); 
        pinv_RR_WW12 = pinv(RR_W12(:,I14)); 
        pinv_RR_WW13 = pinv(RR_W13(:,I14)); 
        for k = 1:K
            pp=ceil(k/(K/G2));
            if(pp == 1)
               pinv_RR_WW14 = pinv_RR_WW10; 
               [Res14(:,k),  est_gg14(:,k)] = MyRes(pinv_RR_WW14,y_w(:,k),RR_W10(:,I14));
            elseif(pp==2)
               pinv_RR_WW14 = pinv_RR_WW11; 
               [Res14(:,k),  est_gg14(:,k)] = MyRes(pinv_RR_WW14,y_w(:,k),RR_W11(:,I14));
            elseif(pp==3)
               pinv_RR_WW14 = pinv_RR_WW12; 
               [Res14(:,k),  est_gg14(:,k)]= MyRes(pinv_RR_WW14,y_w(:,k),RR_W12(:,I14));
            else
               pinv_RR_WW14 = pinv_RR_WW13;
               [Res14(:,k),  est_gg14(:,k)] = MyRes(pinv_RR_WW14,y_w(:,k),RR_W13(:,I14));
            end
            res14(1,k) =  Res14(:,k)'*Res14(:,k);
        end
        MSE14 = sum(res14)/(K*L_r*M*L_t);
    end
    GG14 = zeros(1,length(I14));
    p_av14  = zeros(1,length(I14));
    for L1 = 1: length(I14)
        GG14(1,L1) = est_gg14(L1,:)* est_gg14(L1,:)'; 
        p_av14(1,L1) = GG14(1,L1)/K;        
    end
    PP14 = find(GG14 == max(GG14));
    EE14 = [];
    for L1 = 1:length(I14)
        if(p_av14(1,L1) >= 0.025*GG14(1,PP14)/K)
            EE14 = [EE14 L1];
        end
    end
    est_g14 = zeros(length(I14),K);
    if(MSE > max_e1)
        for k = 1:K
            est_g14(EE14,k) = est_gg14(EE14,k);
        end
    end
end