function [norm_Fro_nHK,norm_Fro_nHK7,norm_Fro_nHK14,norm_Fro_nHK1,norm_Fro_nHK3,norm_Fro_nHK4,norm_Fro_nHK22]...
    =norm_Fro(K,PSI_K,NN,GG,est_g,I,est_g16,I16,est_g14,I4,est_g3,I3,est_k16,I161,est_k,I1,est_k4,I14,H_K,N_r,N_t,MSE,max_e1)
    norm_Fro_nHK1 = zeros(1,K);%�·�����Ϊ1��
    norm_Fro_nHK3 = zeros(1,K);%�½�
    norm_Fro_nHK4 = zeros(1,K);%�·�����Ϊ16��
    norm_Fro_nHK22 = zeros(1,K);%�·�����Ϊ4��
    norm_Fro_nHK = zeros(1,K);%OMP��Ϊ1��
    norm_Fro_nHK7 = zeros(1,K);%OMP��Ϊ16��
    norm_Fro_nHK14 = zeros(1,K);%OMP��Ϊ4��

    n_H_K = zeros(N_t,N_r);
    n_H_K7 = zeros(N_t,N_r);
    n_H_K14 = zeros(N_t,N_r);

    PPsi_K1 = PSI_K(1:NN,1:GG);
    PPsi_K21 = PSI_K(1:NN,1+4*GG:GG*5);
    PPsi_K3 = PSI_K(1:NN,1+8*GG:GG*9);
    PPsi_K4 = PSI_K(1:NN,1+12*GG:GG*13);
    for k = 1:K
        PPsi_K2 = PSI_K(1:NN,1+(k-1)*GG:GG*k);
       if(MSE > max_e1)
            n_H_K = reshape(PPsi_K2(:,I)*est_g(:,k),N_r,N_t);%OMP��Ϊ1��
            n_H_K7 = reshape(PPsi_K2(:,I16)*est_g16(:,k),N_r,N_t);%OMP��Ϊ16��
            n_H_K14 = reshape(PPsi_K2(:,I4)*est_g14(:,k),N_r,N_t);%OMP��Ϊ4��
       end           
        n_H_K1 = reshape(PPsi_K1(:,I1(k,:))*est_k(:,k),N_r,N_t);%�·�����Ϊ1��
        n_H_K3 = reshape(PPsi_K2(:,I3(k,:))*est_g3(:,k),N_r,N_t);%�½�
        n_H_K4 = reshape(PPsi_K2(:,I161(k,:))*est_k16(:,k),N_r,N_t);%�·�����Ϊ16��
        if(k <= 4)
            PPsi_K12 = PPsi_K1;
        elseif(k <= K/2 && k > 4)
            PPsi_K12 = PPsi_K21 ;
        elseif(k <= 12 && k > 8)
            PPsi_K12 = PPsi_K3;%9
        else
            PPsi_K12 = PPsi_K4;
        end
        n_H_K22 = reshape(PPsi_K12(:,I14(k,:))*est_k4(:,k),N_r,N_t); %�·�����Ϊ4��
        norm_Fro_nHK(1,k) = (norm(n_H_K - H_K(1:N_r,1+(k-1)*N_t:N_t*k),'fro'))^2;%OMP��Ϊ1��
        norm_Fro_nHK7(1,k) = (norm(n_H_K7 - H_K(1:N_r,1+(k-1)*N_t:N_t*k),'fro'))^2;%OMP��Ϊ16��
        norm_Fro_nHK14(1,k) = (norm(n_H_K14 - H_K(1:N_r,1+(k-1)*N_t:N_t*k),'fro'))^2;%OMP��Ϊ4��
        norm_Fro_nHK1(1,k) = (norm(n_H_K1-H_K(1:N_r,1+(k-1)*N_t:N_t*k),'fro'))^2;%�·�����Ϊ1��
        norm_Fro_nHK3(1,k) = (norm(n_H_K3-H_K(1:N_r,1+(k-1)*N_t:N_t*k),'fro'))^2;%�½�
        norm_Fro_nHK4(1,k) = (norm(n_H_K4-H_K(1:N_r,1+(k-1)*N_t:N_t*k),'fro'))^2;%�·�����Ϊ16��
        norm_Fro_nHK22(1,k) = (norm(n_H_K22-H_K(1:N_r,1+(k-1)*N_t:N_t*k),'fro'))^2;%�·�����Ϊ4��
    end
end