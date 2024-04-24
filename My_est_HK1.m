function [est_k,I1] = My_est_HK1(y_w,PSI_K,GG,L,NN,K,Fai,G_t,G_r,lamda_k)
    RES2 = y_w(:,1);
    MAX_C = [];
    PPsi_K11 = Fai*PSI_K(1:NN,1:GG);
     P1 = zeros(GG,1);
    while(1)
        if(length(MAX_C) >= L)
            break;
        end
        C1= PPsi_K11'* RES2;
        for m = 1:GG
            P1(m,1) = sum(abs(C1(m,:)));
        end
        pos3 = find(abs(P1) == max(abs(P1)));
        MAX_C = [MAX_C,pos3(1)];
        est_gg2 = zeros(length(MAX_C),1);
        est_gg2(:,1) = pinv(PPsi_K11(:,MAX_C))*y_w(:,1);
        RES2(:,1) = y_w(:,1) - PPsi_K11(:,MAX_C)*est_gg2(:,1);
    end
    PSi_D_index = ceil(MAX_C/G_t);
    PSi_A_index = MAX_C-(PSi_D_index-1)*G_r;
    psi_a = lamda_k(1,1)/lamda_k(1,1)*sin((PSi_A_index-1)/G_r*pi-pi/2);
    PSi_A2 = asin(psi_a);
    psi_d = lamda_k(1,1)/lamda_k(1,1)*sin((PSi_D_index-1)/G_r*pi-pi/2);
    PSi_D2 = asin(psi_d);
    PSi_A = PSi_A2;
    PSi_D = PSi_D2;
    XXX_A = zeros(K,L);
    XXX_D = zeros(K,L);
    for k = 1:K
        XXX_A(k,:) = 1/lamda_k(1,k)*sin(PSi_A);
        XXX_D(k,:) = 1/lamda_k(1,k)*sin(PSi_D);
    end
    min_d_XXX_A = zeros(K-1,L);
    min_d_XXX_D = zeros(K-1,L);
    I1 = zeros(K,L);
    for k = 2:K
        d_XXX_A = zeros(G_t,L);
        d_XXX_D = zeros(G_t,L); 
        for g = 1:G_t
            d_XXX_A(g,:) = XXX_A(k,:)-1/lamda_k(1,1)*sin(-pi/2+pi*(g-1)/G_t);
            d_XXX_D(g,:) = XXX_D(k,:)-1/lamda_k(1,1)*sin(-pi/2+pi*(g-1)/G_t);
        end 
        for ll = 1:L
            min_d_XXX_A(k-1,ll) = find(abs(d_XXX_A(:,ll)) == min(abs(d_XXX_A(:,ll))));
            min_d_XXX_D(k-1,ll) = find(abs(d_XXX_D(:,ll)) == min(abs(d_XXX_D(:,ll))));
        end
    end
    for k = 1:K
        if(k == 1)
            I1(k,:) = MAX_C;
        else
            I1(k,:) = (min_d_XXX_D(k-1,:)-1)*G_t + min_d_XXX_A(k-1,:);
        end
    end
    est_k = zeros(L,K);
    for k = 1:K
        PPsi_K12 = PPsi_K11;
        est_k(:,k) = pinv(PPsi_K12(:,I1(k,:)))*y_w(:,k);
    end
end
