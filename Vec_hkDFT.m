function [n_H_K1,n_H_K4,n_H_K16,n_H_K3,n_H_K11,n_H_K14,n_H_K161]...
    =Vec_hkDFT(K,PSI_K,NN,GG,est_g,I,est_g16,I16,est_g14,I4,est_g3,I3,est_k16,I161,est_k,I1,est_k4,I14,MSE,max_e1)

    n_H_K1 = zeros(NN,K);
    n_H_K16 = zeros(NN,K);
    n_H_K4 = zeros(NN,K);
    n_H_K11 = zeros(NN,K);
    n_H_K161 = zeros(NN,K);
    n_H_K14 = zeros(NN,K);
    n_H_K3 = zeros(NN,K);
% 
    PPsi_K1 = PSI_K(1:NN,1:GG);
    PPsi_K21 = PSI_K(1:NN,1+4*GG:GG*5);
    PPsi_K3 = PSI_K(1:NN,1+8*GG:GG*9);
    PPsi_K4 = PSI_K(1:NN,1+12*GG:GG*13);
    for k = 1:K
        PPsi_K2 = PSI_K(1:NN,1+(k-1)*GG:GG*k);
        n_H_K1(:,k) = PPsi_K2(:,I)*est_g(:,k);%OMP分为1组
        n_H_K16(:,k) =PPsi_K2(:,I16)*est_g16(:,k);%OMP分为16组
        n_H_K4(:,k) = PPsi_K2(:,I4)*est_g14(:,k);%OMP分为4组
        n_H_K3(:,k) = PPsi_K2(:,I3(k,:))*est_g3(:,k);%下
        n_H_K11(:,k) = PPsi_K1(:,I1(k,:))*est_k(:,k);%新方法分为1组
        n_H_K161(:,k) = PPsi_K2(:,I161(k,:))*est_k16(:,k);%新方法分为16组
        if(k <= 4)
            PPsi_K12 = PPsi_K1;
        elseif(k <= K/2 && k > 4)
            PPsi_K12 = PPsi_K21 ;
        elseif(k <= 12 && k > 8)
            PPsi_K12 = PPsi_K3;%9
        else
            PPsi_K12 = PPsi_K4;
        end
        n_H_K14(:,k) = PPsi_K12(:,I14(k,:))*est_k4(:,k); %新方法分为4组
    end
end