function [est_g6,I6] = My_est_HK16(y_w,PSI_K,GG,L,NN,K,Fai)
    I6 = zeros(K,L);
    CCCC11 = zeros(GG,K);
    Res44 = y_w;  
    est_g6 = zeros(L,K);
    for k = 1:K
         PPsi_K = Fai*PSI_K(1:NN,1+(k-1)*GG:GG*k);
         II = [];
         for ll = 1:L
              CCCC11(:,k) = PPsi_K'* Res44(:,k);
              pos55 = find(abs(CCCC11(:,k)) == max(abs(CCCC11(:,k))));
              II =  [II,pos55(1)];
              est_g44 = pinv(PPsi_K(:,II))* y_w(:,k);
              Res44(:,k) =  y_w(:,k) - PPsi_K(:,II)*est_g44;
         end
         I6(k,:)= II;
         est_g6(:,k) = est_g44;
    end
end