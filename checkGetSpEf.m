function [RR] = checkGetSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,S_Lr,S_Lt, NN_H_K,HK16)
        z_k_l = zeros(1,L_t1);
        d_k_l = zeros(1,L_t1);

%% 信道矩阵的伪逆
        A_K = sqrt(a_snr* P_Rc)* A_t_k(:,S_Lr)'*HK16*A_t_k(:,S_Lt)*F_BB;
       % inv_A_K = pinv(A_K);
%         H_e_k1 = NN_H_K-HK16;
%         A_e = sqrt(a_snr* P_Rc)* A_t_k(:,S_Lr)'*H_e_k1* A_t_k(:,S_Lt)*F_BB;%信道估计误差
        A_W = A_K*A_K'+A_t_k(:,S_Lr)'*A_t_k(:,S_Lr);
        WH_BB = A_K'*inv(A_W);
        AK = sqrt(a_snr* P_Rc)*A_t_k(:,S_Lr)'*NN_H_K*A_t_k(:,S_Lt)*F_BB;
%% 信号部分
        a_k =1;
%% 噪声部分
        for n1 = 1:L_t1
            z_k_l(1,n1) = WH_BB(n1,:)* A_t_k(:,S_Lr)'* A_t_k(:,S_Lr)* WH_BB(n1,:)';
        end
%% 估计误差部分
        for n1 = 1:L_t1
            d_k_l(1,n1)= (WH_BB(n1,:)* AK - F_BB(n1,:))*(WH_BB(n1,:)* AK - F_BB(n1,:))';
        end
        RR = sum(log2(1+a_k./(d_k_l+z_k_l)));
end