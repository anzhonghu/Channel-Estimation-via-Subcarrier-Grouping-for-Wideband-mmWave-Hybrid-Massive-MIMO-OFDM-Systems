function [est_g,I] = My_est_g1(y_w,MAX_y_k,GG,K_p,MSE,max_e1,K,M,L_r,R_W,L_t) 
    Res = y_w; 
    C = zeros(GG,K_p);
    P = zeros(GG,1);
    I = [];
    MSE1 = MSE;
    while(MSE1 > max_e1)
        for k_p = 1:K_p
            C(:,k_p) = R_W'* Res(:,MAX_y_k(1,k_p));
        end
        for m = 1:GG
            P(m,1) = sum(abs(C(m,:)));
        end
        pos = find(P == max(P));%%找到最大的投影
        I = [I pos(1)];%%%Gt*Gr; 
        est_gg = zeros(length(I),K);
        res = zeros(1,K);
        pinv_R_W = pinv(R_W(:,I));
        for k = 1:K
             est_gg(:,k) = pinv_R_W*y_w(:,k);
             Res(:,k) = y_w(:,k) - R_W(:,I)*est_gg(:,k);%更新res
             res(1,k) =  Res(:,k)'*Res(:,k);
        end
        MSE1 = sum(res)/(K*L_r*M*L_t);
    end
    GG1 = zeros(1,length(I));
    p_av  = zeros(1,length(I));
    for L1 = 1: length(I)
        GG1(1,L1) = est_gg(L1,:)* est_gg(L1,:)'; 
        p_av(1,L1) = GG1(1,L1)/K;        
    end
    PP = find(p_av == max(p_av));
    EE = [];
    for L1 = 1:length(I)
        if(p_av(1,L1) >= 0.025*p_av(1,PP))
            EE = [EE L1];
        end
    end
    est_g = zeros(length(I),K);
    if(MSE > max_e1)
        for k = 1:K
            est_g(EE,k) = est_gg(EE,k);
        end
    end
end