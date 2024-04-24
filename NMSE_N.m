clear all 
close all;
%% ****************************************************训练帧的个数**************************************�?
M = 256;   
%% ************************************************接收天线和发射天线的个数*******************************�?
NN_t = [24,32,40,48,56];
NN_r = [24,32,40,48,56];
% NN_t = 24;
% NN_r = 24;
%% *******************************************************字典的大�?***********************************
G_t = 64;
G_r = 64;
GG = G_t*G_r;
g_t = 0:1:G_t-1;
g_r = 0:1:G_r-1;
%% *****************************************************子载波的个数*************************************�?
K = 16;
K_p = K/4;   %找到K_p个使得y[k]的LLF函数�?��的子载波�?
G = 1;%将K个子载波分成两组
G1 = 16;
G2 = 4;
%% *******************************************************带宽大小**************************************
f_c = 60*10^9;%载波的频率；
B = f_c*0.2;
c = 3.0*10^8;%速度
k = 0:1:K-1;
f = (k-K/2)*B/(K+1);
lamda_k = c./(f+f_c);%第K个子载波的波长；
lamda_g = zeros(1,G);
%%������ÿ�����ز��Ĳ���
for g = 1:G
    lamda_g(1,g) = G*sum(lamda_k(1+(g-1)*K/G:K*g/G))/K;
end
lamda_g1 = zeros(1,G1);
for g = 1:G1
    lamda_g1(1,g) = G1*sum(lamda_k(1+(g-1)*K/G1:K*g/G1))/K;
end
lamda_g2 = zeros(1,G2);
for g = 1:G2
    lamda_g2(1,g) = G2*sum(lamda_k(1+(g-1)*K/G2:K*g/G2))/K;
end
% max_e1 = 1;% max_e1 = norm(trace(inv_D_W1'*C_W*inv_D_W1)*N0/L_r/M);
d = c/(2*f_c);%%天线之间的距�?
L = 4;%路径�?
N0  = 1;
Nc = 4;
%% ******************************************************信道增益的计算以及实际信道的表示**************************************
CC = 200;
D = 0:1:K-1;
Nmse_nHK = zeros(1,length(NN_t));
Nmse_nHK1= zeros(1,length(NN_t));
Nmse_nHK14 = zeros(1,length(NN_t));
Nmse_nHK3 = zeros(1,length(NN_t));
Nmse_nHK4 = zeros(1,length(NN_t));
Nmse_nHK7 = zeros(1,length(NN_t));
Nmse_nHK22 = zeros(1,length(NN_t));

V_crb = zeros(1,length(NN_t));
P_Rc = 32.5 + 20*log10(200) + randn(1,1);
%% �����ŵ���F�����ľ�ֵ
Ts = (K+1)/K/B;
Tb = Ts*3;
sum_HK1 = 0;
EEE = exp(1j*2*pi.*D*K/2/K);
for d1 = 0:K-1
   KKK1 = exp(-1j*2*pi*d1.*D/K);
   E_K = EEE.*KKK1;
   SS1 = sinc(pi*(D.*Ts-Tb)./Tb);
   DD1 = (cos(pi*0.8*(D.*Ts-Tb)./Tb))./(1-4*0.64*((D.*Ts-Tb)./Tb).^2);
   WW1 = SS1.*DD1;
   DFT_WW1 = WW1*E_K.'; %b[k]
   DFT_sum = DFT_WW1* DFT_WW1';
   sum_HK1 = sum_HK1 + DFT_sum;
end
HHH = K/sum_HK1
%% ��������
KKKK = 0;
SNR = 20;
snr = 10^(SNR/10);
L_t = 1;%发�?端RF的个数；
L_r = 4;%接收端RF的个数；
max_e1 = 1;
A1 = sin(-pi/2+pi*g_t/G_t);
for s = 1:length(NN_t)
    N_t = NN_t(1,s);
    N_r = NN_r(1,s);
    NN = N_t*N_r;

    KKKK = KKKK + 1
    % A_t_k1_g = zeros(N_t,G_t);%离开�?

    %% ******************************************************将子载波分为2组获得的AOA/AOD字典矩阵**********************************************
    PSI_g = zeros(NN,G*GG);
    A_t_k1_g=zeros(N_t,G_t);
    N_tt = 0:N_t-1 ;
    n_t_A = N_tt.'*A1;
    for g = 1:G
        A_t_k1_g = exp(-1j*(2*pi*d)/lamda_g(1,g)*n_t_A)/sqrt(N_t);
        PSI_g(1:N_t*N_r,1+(g-1)*G_r*G_t:G_r*G_t*g) = kron(conj(A_t_k1_g),A_t_k1_g);
    end
    %% ******************************************************将子载波分为16组获得的AOA/AOD字典矩阵**********************************************
    PSI_G1 = zeros(NN,G1*GG);
    A_t_k2_g = zeros(N_t,G_t);
    for g = 1:G1
        A_t_k2_g = exp(-1j*(2*pi*d)/lamda_g1(1,g)*n_t_A)/sqrt(N_t);
        PSI_G1(1:N_t*N_r,1+(g-1)*G_r*G_t:G_r*G_t*g) = kron(conj(A_t_k2_g),A_t_k2_g);
    end
    %% ******************************************************将子载波分为4组获得的AOA/AOD字典矩阵**********************************************
    PSI_G2 = zeros(NN,G2*GG);
    A_t_k3_g = zeros(N_t,G_t);
    for g = 1:G2
        A_t_k3_g = exp(-1j*(2*pi*d)/lamda_g2(1,g)*n_t_A)/sqrt(N_t);
        PSI_G2(1:N_t*N_r,1+(g-1)*G_r*G_t:G_r*G_t*g) = kron(conj(A_t_k3_g),A_t_k3_g);
    end
    %% *************************************************************测量矩阵******************************************

    A_t_k1_k = zeros(N_t,G_t);
    PSI_K  = PSI_G1;

    NMSE_nHK = zeros(1,CC);
    NMSE_nHK1 = zeros(1,CC);
    NMSE_nHK3 = zeros(1,CC);
    NMSE_nHK4 = zeros(1,CC);
    NMSE_nHK7 = zeros(1,CC);
    NMSE_nHK22 = zeros(1,CC);
    NMSE_nHK14 = zeros(1,CC);
      V_CRB = zeros(1,CC); 
    
    nMSE_nHK1 = zeros(1,CC);
    nMSE_nHK = zeros(1,CC);
    nMSE_nHK3 = zeros(1,CC);
    nMSE_nHK4 = zeros(1,CC);
    nMSE_nHK7 = zeros(1,CC);
    nMSE_nHK22 = zeros(1,CC);
    nMSE_nHK14 = zeros(1,CC);
      v_crb = zeros(1,CC);
    for a = 1:CC    
        N_Q = 2;%相移器的量化�?
        HK = zeros(NN,K);
        [W_tr,C_W,Fai1] = TrainFW(L_t,L_r,N_r,M,N_t,NN);%����Ftr��Wtr
        Fai = sqrt(snr* P_Rc)*Fai1;
        R_W = zeros(M*L_r*L_t,G*GG);
        RR_W = zeros(M*L_r*L_t,G1*GG);
        RR_W1 = zeros(M*L_r*L_t,G2*GG);
        for g = 1:G
            R_W(1:M*L_r*L_t,1+(g-1)*GG:GG*g)= Fai* PSI_g(1:NN,1+(g-1)*GG:GG*g);
        end
        for g = 1:G1
            RR_W(1:M*L_r*L_t,1+(g-1)*GG:GG*g)=Fai* PSI_G1(1:NN,1+(g-1)*GG:GG*g);
        end
        for g = 1:G2
            RR_W1(1:M*L_r*L_t,1+(g-1)*GG:GG*g)= Fai* PSI_G2(1:NN,1+(g-1)*GG:GG*g);
        end
        noise2 = zeros(M*L_r,1);
        for m = 1:M
            noise1 = (randn(N_r,1) + 1j*randn(N_r,1))/sqrt(2);
            wtr = W_tr(1:N_r,1+(m-1)*L_r:L_r*m);
            noise2(1 + (m-1)*L_r:L_r*m,1)= wtr'* noise1;
        end
        noise = noise2;
        H_K = zeros(N_r,N_t*K);%K个子载波的信道矩�?
        G_d_k = zeros(L,L*K);%路径增益与Prc
        G_k = zeros(L,L);%对角矩阵对角线上的元素是每条路径上的增益
        Y = zeros(M*L_r*L_t,K);%Y为接收信号；
        y_w = zeros(M*L_r*L_t,K);
        y_K = zeros(1,K);%每个子载波的2范数�?
        %% H_K
        delay_L1= rand(1,L)*(Nc-1)*Ts; 
        alphal = sqrt(HHH)*(randn(1,L)+1j*randn(1,L))/sqrt(2);
        AOA = -pi/2 + rand(1,L) * pi;%到达�?
        AOD = -pi/2 + rand(1,L) * pi;%离开�?
        n_t_S_AOD = N_tt.'*sin(AOD);
        n_t_S_AOA = N_tt.'*sin(AOA);
        for l=1:L
            SS = sinc(pi*(D.*Ts-delay_L1(1,l)-Tb)./Tb);
            DD = (cos(pi*0.8*(D.*Ts-delay_L1(1,l)-Tb)./Tb))./(1-4*0.64*((D.*Ts-delay_L1(1,l)-Tb)./Tb).^2);
            WW = SS.*DD;
            G_d_k(l,1+(l-1)*K:K*l) = alphal(1,l)*WW*sqrt(N_t*N_r/L/P_Rc);
        end
        norm_Fro_HK = zeros(1,K);
        EEE = exp(1j*2*pi.*D*K/2/K);
        for k = 1:K
            KKK = exp(-1j*2*pi*(k-1).*D/K);
            E_K = EEE.*KKK;
            A_t = exp(-1j*(2*pi*d)/lamda_k(1,k)*n_t_S_AOD )/sqrt(N_t);
            A_r= exp(-1j*(2*pi*d)/lamda_k(1,k)*n_t_S_AOA)/sqrt(N_r);
            for l=1:L
                 G_k(l,l) =  G_d_k(l,1+(l-1)*K:K*l)*E_K.';%G_d_k傅里叶变�?
            end
            H_K(1:N_r,1+(k-1)*N_t:N_t*k) = A_r*G_k* A_t';
            HK(:,k)=reshape(H_K(1:N_r,1+(k-1)*N_t:N_t*k),[],1);
            Y(:,k) = Fai * HK(:,k) + noise;
            norm_Fro_HK(1,k)=(norm(H_K(1:N_r,1+(k-1)*N_t:N_t*k),'fro'))^2;
        end
        y_w = Y;
        %% ************************************************************计算合�?的E*************************************************
        MSE1 = trace(y_w'*y_w)/(K*L_r*M);
        MSE = MSE1;%          max_e1 = norm(trace(inv_D_W1'*C_W*inv_D_W1)*N0/L_r/M)
%% ********************************************************OMP
        for k = 1:K
            y_K(1,k) = Y(:,k)'*Y(:,k);
        end
        [A,Location] = sort(y_K);
        MAX_y_k = Location(K-K_p+1:K);
        L1 = 14;
        MAX_y_k = Location(K-K_p+1:K);
        [est_g,I] = My_est_g1(y_w,MAX_y_k,GG,K_p,MSE,max_e1,K,M,L_r,R_W,L_t); %%  分为1组的信道估计
        [est_g16,I16] = My_est_g(y_w,MAX_y_k,GG,K_p,MSE,max_e1,K,G1,M,L_r,L_t,RR_W);%%  分为16组的信道估计
        [est_g14,I4] = My_est_g4(y_w,MAX_y_k,GG,K_p,MSE,max_e1,K,G2,M,L_r,L_t,RR_W1);%%  分为4组的
        [est_g3,I3] = LowerBound(HK,PSI_K,GG,L1,NN,K);%�½�%  
        [est_k16,I161] = My_est_HK16(y_w,PSI_K,GG,L1,NN,K,Fai);%�·�������ʮ����
        [est_k,I1] = My_est_HK1(y_w,PSI_K,GG,L1,NN,K,Fai,G_t,G_r,lamda_k);%�·�������1��
        [est_k4,I14] = My_est_HK4(y_w,PSI_K,GG,L1,NN,K,Fai,G_t,G_r,lamda_k);%�·�������4��
        
        [norm_Fro_nHK,norm_Fro_nHK7,norm_Fro_nHK14,norm_Fro_nHK1,norm_Fro_nHK3,norm_Fro_nHK4,norm_Fro_nHK22]...
    =norm_Fro(K,PSI_K,NN,GG,est_g,I,est_g16,I16,est_g14,I4,est_g3,I3,est_k16,I161,est_k,I1,est_k4,I14,H_K,N_r,N_t,MSE,max_e1);
 %% CRB
        %S1 = 10;
        S1 = 14;
        Sum_vv_k = 0;
        R1 = real(C_W);
        RC = imag(C_W);
        C8 = zeros(2*M*L_r,2*M*L_r);
        C8(1:M*L_r,1:M*L_r)=R1;
        C8(1:M*L_r,1+M*L_r:2*M*L_r) = -RC;
        C8(1+M*L_r:2*M*L_r,1:M*L_r) = RC;
        C8(1+M*L_r:2*M*L_r,1+M*L_r:2*M*L_r) = R1;
        C9 = inv(C8/2);
        JJ = zeros(S1,S1);
        for k = 1:K
            PSI_KK = PSI_K(1:NN,(k-1)*GG+1:GG*k);
            A_k = Fai*PSI_KK(:,I3(k,:));
            R_A_k = real(A_k);
            I_A_k = imag(A_k);
            AA_K = zeros(2*M*L_r,2*S1);
            AA_K(1:M*L_r,1:S1) =  R_A_k ;
            AA_K(1:M*L_r,1+S1:2*S1)= -I_A_k;
            AA_K(1+M*L_r:2*M*L_r,1:S1)= I_A_k;
            AA_K(1+M*L_r:2*M*L_r,1+S1:2*S1) =  R_A_k ;
            J_k = AA_K.'*C9*AA_K;
            inv_J_k = inv(J_k);
            for s1 = 1:S1
                JJ(s1,s1) = inv_J_k(s1,s1) + inv_J_k(s1 + S1,s1 + S1);
            end
            vv_k = trace(PSI_KK(:,I3(k,:))*JJ*PSI_KK(:,I3(k,:))');
            Sum_vv_k = Sum_vv_k + vv_k;
        end

        NMSE_nHK(1,a) = sum(norm_Fro_nHK)/sum(norm_Fro_HK);%OMP��Ϊ1��
        NMSE_nHK7(1,a) = sum(norm_Fro_nHK7)/sum(norm_Fro_HK);%OMP��Ϊ16��?
        NMSE_nHK14(1,a) = sum(norm_Fro_nHK14)/sum(norm_Fro_HK);%OMP��Ϊ4��?
        NMSE_nHK1(1,a) = sum(norm_Fro_nHK1)/sum(norm_Fro_HK);%�·�����Ϊ1��
        NMSE_nHK4(1,a) = sum(norm_Fro_nHK4)/sum(norm_Fro_HK);%�·�����Ϊ16��
        NMSE_nHK22(1,a) = sum(norm_Fro_nHK22)/sum(norm_Fro_HK);%�·�����Ϊ4��
        NMSE_nHK3(1,a) = sum(norm_Fro_nHK3)/sum(norm_Fro_HK);%�½�
          V_CRB(1,a) = Sum_vv_k/sum(norm_Fro_HK);
     
        nMSE_nHK(1,a) = 10*log10(NMSE_nHK(1,a));%OMP��Ϊ1��
        nMSE_nHK7(1,a) = 10*log10(NMSE_nHK7(1,a));%OMP��Ϊ16��
        nMSE_nHK14(1,a) = 10*log10(NMSE_nHK14(1,a));%OMP��Ϊ4��? 
        nMSE_nHK1(1,a) = 10*log10(NMSE_nHK1(1,a));%�·�����Ϊ1��
        nMSE_nHK4(1,a) = 10*log10(NMSE_nHK4(1,a));%�·�����Ϊ16��
        nMSE_nHK22(1,a) = 10*log10(NMSE_nHK22(1,a));%�·�����Ϊ4��
        nMSE_nHK3(1,a) = 10*log10(NMSE_nHK3(1,a));%�½�
          v_crb(1,a) = 10*log10(V_CRB(1,a));
    end
    Nmse_nHK(1,s) = sum(nMSE_nHK)/CC;%OMP��Ϊ1��
    Nmse_nHK7(1,s) = sum(nMSE_nHK7)/CC;%OMP��Ϊ16��
    Nmse_nHK14(1,s) = sum(nMSE_nHK14)/CC;%OMP��Ϊ4��?
    Nmse_nHK1(1,s) = sum(nMSE_nHK1)/CC;%�·�����Ϊ1��
    Nmse_nHK4(1,s) = sum(nMSE_nHK4)/CC;%�·�����Ϊ16��
    Nmse_nHK22(1,s) = sum(nMSE_nHK22)/CC;%�·�����Ϊ4��
    Nmse_nHK3(1,s) = sum(nMSE_nHK3)/CC;%�½�
    V_crb(1,s) = sum(v_crb)/CC;
end
plot(NN_t,Nmse_nHK,'g-o','Markersize',7,'Linewidth',1);%SS-SW-OMP-Th,Th=0.025,P=1
hold on;
plot(NN_t,Nmse_nHK1,'r-o','Markersize',7,'Linewidth',1);%NewMethod,P=1
hold on;
plot(NN_t,Nmse_nHK3,'b-*','Markersize',7,'Linewidth',1);%LowerBound
hold on;
plot(NN_t,Nmse_nHK4,'r->','Markersize',7,'Linewidth',1);%NewMethod,P=16
hold on;
plot(NN_t,Nmse_nHK7,'g->','Markersize',7,'Linewidth',1);%SS-SW-OMP-Th,Th=0.025,P=16
hold on;
plot(NN_t,Nmse_nHK22,'r-<','Markersize',7,'Linewidth',1);%NewMethod,P=4
hold on;
plot(NN_t,Nmse_nHK14,'g-<','Markersize',7,'Linewidth',1);%SS-SW-OMP-Th,Th=0.025,P=4
hold on;
plot(NN_t,V_crb,'r-*','Markersize',7,'Linewidth',1);
xlabel('The number of antenna'),ylabel('NMSE[dB]');
% set(gca,'XLim',[32,64]);
legend('SS-SW-OMP+Th��P=1','NewMethod��P=1','LowerBound','NewMethod��P=16','SS-SW-OMP+Th��P=16','NewMethod��P=4','SS-SW-OMP+Th��P=4','CRB');



