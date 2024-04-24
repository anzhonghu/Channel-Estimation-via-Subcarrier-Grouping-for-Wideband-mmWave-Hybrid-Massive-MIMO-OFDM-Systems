function [ est_k223,I12] = My_est_HK4(y_w,PSI_K,GG,L,NN,K,Fai,G_t,G_r,lamda_k)
        RES21 = y_w(:,1);
        RES22 = y_w(:,9);
        RES23 = y_w(:,5);
        RES24 = y_w(:,13);
        MAX_C21 = [];
        MAX_C22 = [];
        MAX_C23 = [];
        MAX_C24 = [];
        %% 1
         P12 = zeros(GG,1);
         PPsi_K121 = Fai*PSI_K(1:NN,1:GG);
        while(1)
            if(length(MAX_C21) >= L)
                break;
            end
             C11= PPsi_K121'*  RES21;
            for m = 1:G_r*G_t
                P12(m,1) = sum(abs(C11(m,:)));
            end
            pos31 = find(abs(P12) == max(abs(P12)));
            MAX_C21 = [MAX_C21,pos31(1)];
            est_gg21 = zeros(length(MAX_C21),1);
            est_gg21(:,1) = pinv(PPsi_K121(:,MAX_C21))*y_w(:,1);
            RES21(:,1) = y_w(:,1) - PPsi_K121(:,MAX_C21)*est_gg21(:,1);
        end
        %% 9
         PPsi_K122 = Fai*PSI_K(1:NN,(K/2)*GG + 1:GG*(K/2 + 1));
%         PPsi_K122 = Fai*PSI_K(1:NN,2*GG + 1:GG*3);
        P22 = zeros(GG,1);
        while(1)
            if(length(MAX_C22) >= L)
                break;
            end
             C12= PPsi_K122'*  RES22;
            for m = 1:G_r*G_t
                P22(m,1) = sum(abs(C12(m,:)));
            end
            pos32 = find(abs(P22) == max(abs(P22)));
            MAX_C22 = [MAX_C22,pos32(1)];
            est_gg22 = zeros(length(MAX_C22),1);
            est_gg22(:,1) = pinv(PPsi_K122(:,MAX_C22))*y_w(:,9);
            RES22(:,1) = y_w(:,9) - PPsi_K122(:,MAX_C22)*est_gg22(:,1);
        end
         %% 5
         PPsi_K123 = Fai*PSI_K(1:NN,(K/4)*GG + 1:GG*(K/4 + 1));
%         PPsi_K123 = Fai*PSI_K(1:NN,GG + 1:GG*2);
        P23 = zeros(GG,1);
        while(1)
            if(length(MAX_C23) >= L)
                break;
            end
             C13= PPsi_K123'*  RES23;
            for m = 1:G_r*G_t
                P23(m,1) = sum(abs(C13(m,:)));
            end
            pos33 = find(abs(P23) == max(abs(P23)));
            MAX_C23 = [MAX_C23,pos33(1)];
            est_gg23 = zeros(length(MAX_C23),1);
            est_gg23(:,1) = pinv(PPsi_K123(:,MAX_C23))*y_w(:,5);
            RES23(:,1) = y_w(:,5) - PPsi_K123(:,MAX_C23)*est_gg23(:,1);
        end
        %% 13
         PPsi_K124 = Fai*PSI_K(1:NN,(3*K/4)*GG + 1:GG*(3*K/4+1));
%         PPsi_K124 = Fai*PSI_K(1:NN,3*GG + 1:GG*4);
        P24 = zeros(GG,1);
        while(1)
            if(length(MAX_C24) >= L)
                break;
            end
             C14= PPsi_K124'*  RES24;
            for m = 1:G_r*G_t
                P24(m,1) = sum(abs(C14(m,:)));
            end
            pos34 = find(abs(P24) == max(abs(P24)));
            MAX_C24 = [MAX_C24,pos34(1)];
            est_gg24 = zeros(length(MAX_C24),1);
            est_gg24(:,1) = pinv(PPsi_K124(:,MAX_C24))*y_w(:,13);
            RES24(:,1) = y_w(:,13) - PPsi_K124(:,MAX_C24)*est_gg24(:,1);
        end
         
     %% 1
        PSi_D_index21 = ceil(MAX_C21/G_t);
        PSi_A_index21 = MAX_C21-(PSi_D_index21-1)*G_r;%A_r
        PSi_A1 = (PSi_A_index21-1)/G_r*pi-pi/2;
        PSi_D1 = (PSi_D_index21-1)/G_r*pi-pi/2;
       %% 9
        PSi_D_index22 = ceil(MAX_C22/G_t);
        PSi_A_index22 = MAX_C22-(PSi_D_index22-1)*G_r;%A_r
        PSi_A2 = (PSi_A_index22-1)/G_r*pi-pi/2;
        PSi_D2 = (PSi_D_index22-1)/G_r*pi-pi/2;
        %% 5
        PSi_D_index23 = ceil(MAX_C23/G_t);
        PSi_A_index23 = MAX_C23-(PSi_D_index23-1)*G_r;%A_r
        PSi_A3 = (PSi_A_index23-1)/G_r*pi-pi/2;
        PSi_D3 = (PSi_D_index23-1)/G_r*pi-pi/2;
       %% 13
        PSi_D_index24 = ceil(MAX_C24/G_t);
        PSi_A_index24 = MAX_C24-(PSi_D_index24-1)*G_r;%A_r
        PSi_A4 = (PSi_A_index24-1)/G_r*pi-pi/2;
        PSi_D4 = (PSi_D_index24-1)/G_r*pi-pi/2;

        
        XXX_A1 = zeros(K,L);
        XXX_D1 = zeros(K,L);
        for k = 1:K
            if(k <= K/4 )
                XXX_A1(k,:) = 1/lamda_k(1,k)*sin(PSi_A1);
                XXX_D1(k,:) = 1/lamda_k(1,k)*sin(PSi_D1);
            elseif(k <= K/2 && k >K/4)
                XXX_A1(k,:) = 1/lamda_k(1,k)*sin(PSi_A3);
                XXX_D1(k,:) = 1/lamda_k(1,k)*sin(PSi_D3);
            elseif(k <= 3*K/4 && k > K/2)
                XXX_A1(k,:) = 1/lamda_k(1,k)*sin(PSi_A2);
                XXX_D1(k,:) = 1/lamda_k(1,k)*sin(PSi_D2);
            else
                XXX_A1(k,:) = 1/lamda_k(1,k)*sin(PSi_A4);
                XXX_D1(k,:) = 1/lamda_k(1,k)*sin(PSi_D4);
            end
        end
        
        
        %% 1
        min_d_XXX_A11 = zeros(K/4-1,L);
        min_d_XXX_D11 = zeros(K/4-1,L);
       
        for k = 2:K/4
            d_XXX_A11 = zeros(G_t,L);
            d_XXX_D11 = zeros(G_t,L);
            for g = 1:G_t
                d_XXX_A11(g,:) = XXX_A1(k,:)-1/lamda_k(1,1)*sin(-pi/2+pi*(g-1)/G_t);
                d_XXX_D11(g,:) = XXX_D1(k,:)-1/lamda_k(1,1)*sin(-pi/2+pi*(g-1)/G_t);
            end
            for ll = 1:L
                min_d_XXX_A11(k-1,ll) = find(abs(d_XXX_A11(:,ll)) == min(abs(d_XXX_A11(:,ll))));
                min_d_XXX_D11(k-1,ll) = find(abs(d_XXX_D11(:,ll)) == min(abs(d_XXX_D11(:,ll))));
            end
        end
        %% 5
        min_d_XXX_A13 = zeros(K/4-1,L);
        min_d_XXX_D13 = zeros(K/4-1,L);
        for k = K/4+2:K/2
            d_XXX_A13 = zeros(G_t,L);
            d_XXX_D13 = zeros(G_t,L);
            for g = 1:G_t
                d_XXX_A13(g,:) = XXX_A1(k,:)-1/lamda_k(1,5)*sin(-pi/2+pi*(g-1)/G_t);
                d_XXX_D13(g,:) = XXX_D1(k,:)-1/lamda_k(1,5)*sin(-pi/2+pi*(g-1)/G_t);
            end
            for ll = 1:L
                min_d_XXX_A13(k-1,ll) = find(abs(d_XXX_A13(:,ll)) == min(abs(d_XXX_A13(:,ll))));
                min_d_XXX_D13(k-1,ll) = find(abs(d_XXX_D13(:,ll)) == min(abs(d_XXX_D13(:,ll))));
            end
        end
        %% 9
        min_d_XXX_A12 = zeros(K/4-1,L);
        min_d_XXX_D12 = zeros(K/4-1,L);
        for k = K/2 + 2:3*K/4
            d_XXX_A12 = zeros(G_t,L);
            d_XXX_D12 = zeros(G_t,L);
            for g = 1:G_t
                d_XXX_A12(g,:) = XXX_A1(k,:)-1/lamda_k(1,9)*sin(-pi/2+pi*(g-1)/G_t);
                d_XXX_D12(g,:) = XXX_D1(k,:)-1/lamda_k(1,9)*sin(-pi/2+pi*(g-1)/G_t);
            end
            for ll = 1:L
                min_d_XXX_A12(k-1,ll) = find(abs(d_XXX_A12(:,ll)) == min(abs(d_XXX_A12(:,ll))));
                min_d_XXX_D12(k-1,ll) = find(abs(d_XXX_D12(:,ll)) == min(abs(d_XXX_D12(:,ll))));
            end
        end
        %% 13
        min_d_XXX_A14 = zeros(K/4-1,L);
        min_d_XXX_D14 = zeros(K/4-1,L);
        
        for k =3*K/4+2:K
            d_XXX_A14 = zeros(G_t,L);
            d_XXX_D14 = zeros(G_t,L);
            for g = 1:G_t
                 d_XXX_A14(g,:) = XXX_A1(k,:)-1/lamda_k(1,13)*sin(-pi/2+pi*(g-1)/G_t);
                d_XXX_D14(g,:) = XXX_D1(k,:)-1/lamda_k(1,13)*sin(-pi/2+pi*(g-1)/G_t);
            end
            for ll = 1:L
                min_d_XXX_A14(k-1,ll) = find(abs(d_XXX_A14(:,ll)) == min(abs(d_XXX_A14(:,ll))));
                min_d_XXX_D14(k-1,ll) = find(abs(d_XXX_D14(:,ll)) == min(abs(d_XXX_D14(:,ll))));
            end
        end
        
        I12 = zeros(K,L);
        for k = 1:K/4
            if(k == 1)
                I12(k,:) = MAX_C21;
            else
                I12(k,:) = (min_d_XXX_D11(k-1,:)-1)*G_t + min_d_XXX_A11(k-1,:);
            end
        end
        for k = K/4+1:K/2
            if(k == 5)
                I12(k,:) = MAX_C23;
            else
                I12(k,:) =  (min_d_XXX_D13(k-1,:)-1)*G_t + min_d_XXX_A13(k-1,:);
            end
        end
        for k = K/2 + 1:3*K/4
            if(k == 9)
                I12(k,:) = MAX_C22;
            else
                I12(k,:) = (min_d_XXX_D12(k-1,:)-1)*G_t + min_d_XXX_A12(k-1,:);
            end
        end
        for k = 3*K/4 + 1:K
            if(k == 13)
                I12(k,:) = MAX_C24;
            else
                I12(k,:) =(min_d_XXX_D14(k-1,:)-1)*G_t + min_d_XXX_A14(k-1,:);
            end
        end
        est_k223 = zeros(L,K);
        for k = 1:K
            if(k <= 4)
                PPsi_K12 = PPsi_K121;
            elseif(k <= K/2 && k > 4)
                PPsi_K12 =  PPsi_K123 ;
            elseif(k <= 12 && k > 8)
                PPsi_K12 = PPsi_K122;%9
            else
                PPsi_K12 = PPsi_K124;
            end
            est_k223(:,k) = pinv(PPsi_K12(:,I12(k,:)))*y_w(:,k);
        end
end