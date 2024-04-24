function [W_tr,C_W,Fai1,F_tr] = TrainFW(L_t,L_r,N_r,M,N_t,NN)
    f_tr = zeros(N_t,L_t);%预编码器�?
    w_tr = zeros(N_r,L_r);%组合器；
    F_tr = zeros(N_t,L_t*M);
    W_tr = zeros(N_r,L_r*M);%组合器；
    Fai1 = zeros(M*L_r*L_t,NN);%%测量矩阵
    for m = 1:M %计算Ftr
        x = rand(N_t,L_t);
        f_tr(find((x<=0.25))) = -1/sqrt(N_t);
        f_tr(find((x>0.25))) = 1/sqrt(N_t);
        f_tr(find((x>0.5))) = -1j/sqrt(N_t);
        f_tr(find((x>0.75))) = 1j/sqrt(N_t);
        %计算Wtr
        x = rand(N_r,L_r);
        w_tr(find((x<=0.25))) = -1/sqrt(N_r);
        w_tr(find((x>0.25))) = 1/sqrt(N_r);
        w_tr(find((x>0.5))) = -1j/sqrt(N_r);
        w_tr(find((x> 0.75))) = 1j/sqrt(N_r);
        F_tr(1:N_r,1+(m-1)*L_t:L_t*m)=f_tr;
        W_tr(1:N_r,1+(m-1)*L_r:L_r*m)=w_tr;
        C_W(1+(m-1)*L_r:L_r*m,1+(m-1)*L_r:L_r*m) = w_tr'*w_tr;
        Fai1(1 + (m-1)*L_r*L_t:L_r*L_t*m,1:N_t*N_r)= kron(f_tr.',w_tr');
    end
end