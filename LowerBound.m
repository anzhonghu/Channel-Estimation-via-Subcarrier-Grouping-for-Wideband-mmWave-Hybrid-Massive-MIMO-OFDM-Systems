function [est_g5,I5] = LowerBound(HK,PSI_K,GG,L1,NN,K)
%�½�
%HK ���ŵ�
%PSI_K ��A_t[k]��A_r[k]�Ŀ����ڿ˻�
%GG ��G_t*G_r
%NN ��N_t*N_r
%K �����ز���
    L = L1;
    I5 = zeros(K,L);
    Res4 = HK;
    CCCC = zeros(GG,K);
    est_g5 = zeros(L,K);
    for k = 1:K
        PPsi_K = PSI_K(1:NN,1+(k-1)*GG:GG*k);
        I4  = [];
        for ll = 1:L
            if(k <= K/2)
                 CCCC(:,k) = PPsi_K'* Res4(:,k); 
            else
                 CCCC(:,k) = PPsi_K'* Res4(:,k);
            end
            pos4 = find(abs(CCCC(:,k)) == max(abs(CCCC(:,k))));
            I4 =  [I4,pos4(1)];
            est_g4 = pinv(PPsi_K(:,I4))* HK(:,k);
            Res4(:,k) = HK(:,k)-PPsi_K(:,I4)*est_g4;
        end
        I5(k,:)= I4;
        est_g5(:,k) = est_g4;
    end
end