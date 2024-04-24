function [W_tr,C_W,Fai1,F_tr] = TrainFWDFT(L_t,L_r,N_r,M,N_t,NN,A_t_k1)
    N = sqrt(M*L_r);
    f_tr = zeros(N_t,L_t);%棰勭紪鐮佸櫒锛?
    w_tr = zeros(N_r,L_r);%缁勫悎鍣紱
    F_tr = zeros(N_t,L_t*M);
    W_tr = zeros(N_r,L_r*M);%缁勫悎鍣紱
    Fai1 = zeros(M*L_r*L_t,NN);%%娴嬮噺鐭╅樀
    C_W = zeros(L_r*M,L_r*M);
    for m = 1:M %璁＄畻Ftr
        f_tr = A_t_k1(:,ceil(m/(N/L_r)));
        for n = 1:L_r
           cc = (m-(ceil(m/(N/L_r))-1)*(N/L_r)-1)*L_r+n;
           w_tr(:,n) = A_t_k1(:,cc);
        end
          F_tr(1:N_r,1+(m-1)*L_t:L_t*m)=f_tr;
          W_tr(1:N_r,1+(m-1)*L_r:L_r*m)=w_tr;
          C_W(1+(m-1)*L_r:L_r*m,1+(m-1)*L_r:L_r*m) = w_tr'*w_tr;
          Fai1(1 + (m-1)*L_r*L_t:L_r*L_t*m,1:N_t*N_r)= kron(f_tr.',w_tr');
    end
end