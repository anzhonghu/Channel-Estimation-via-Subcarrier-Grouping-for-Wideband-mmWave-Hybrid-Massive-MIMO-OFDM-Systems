clear all; 
close all;
par.runId = 1;
rng(par.runId);
% ****************************************************传输符号个数**************************************
M = 256;               
% *******************************************************网格的大小***********************************
G_t = 64;%网格数
G_r = 64;
K = 16;%子载波个数
K_p = K/4;   %从K个子载波中选出K/4个子载波
G = 1;%不分组
G1 =16;%分为16组
G2 = 4;%分为4组

% ************************************************天线个数******************************
N_t = 32;
N_r = 32;
NN = N_t*N_r;
f_c = 60*10^9;%载波频率
B = f_c*0.2;
c = 3.0*10^8;%电磁波速度
k = 0:1:K-1;
f = (k-K/2)*B/(K+1);
lamda_k = c./(f+f_c);%每个子载波的波长

lamda_g = zeros(1,G);
for g = 1:G%分为1组的时候的该组的波长
    lamda_g(1,g) = G*sum(lamda_k(1+(g-1)*K/G:K*g/G))/K;
end
lamda_g1 = zeros(1,G1);
for g = 1:G1%分为16组的时候的每组的波长
    lamda_g1(1,g) = G1*sum(lamda_k(1+(g-1)*K/G1:K*g/G1))/K;
end
lamda_g2 = zeros(1,G2);
for g = 1:G2%分为G2组的时候的每组的波长
    lamda_g2(1,g) = G2*sum(lamda_k(1+(g-1)*K/G2:K*g/G2))/K;
end
d = c/(2*f_c);%%天线间的距离
L = 4;%路径个数
N0  = 1;
Nc = 4;
Ts = (K+1)/K/B;
TL = K/B;

% ******************************************************淇￠亾澧炵泭鐨勮绠椾互鍙婂疄闄呬俊閬撶殑琛ㄧず**************************************
CC = 1;
D = 0:1:K-1;
SNR = 20;
snr = 10.^(SNR/10);
SNR1 = -20:5:50;
snr1 = 10.^(SNR1/10);

max_e1 = 1;%   
P_Rc = 32.5 + 20*log10(200) + randn(1,1);
% 信道的均值           

Tb = Ts*3;
sum_HK1 = 0;
for d1 = 0:K-1
   KKK1 = exp(-1j*2*pi*d1.*D/K);
   SS1 = sinc(pi*(D.*Ts-Tb)./Tb);
   DD1 = (cos(pi*0.8*(D.*Ts-Tb)./Tb))./(1-4*0.64*((D.*Ts-Tb)./Tb).^2);
   WW1 = SS1.*DD1;
   DFT_WW1 = WW1*KKK1.'; %b[k]
   DFT_sum = DFT_WW1* DFT_WW1';
   sum_HK1 = sum_HK1 + DFT_sum;
end
HHH = K/sum_HK1
KKKK = 0;

L_t = 1;%鍙戦?绔疪F鐨勪釜鏁帮紱
L_r = 4;%鎺ユ敹绔疪F鐨勪釜鏁帮紱

N_Q = 2;%
N_tt = 0:N_t-1 ;
GG = G_t*G_r;
A_t_k1_g = zeros(N_t,G_t);%绂诲紑瑙?
A_r_k1_g = zeros(N_r,G_r);%鍒拌揪瑙?
A_r_K = zeros(N_r,G*G_r);
A_t_K = zeros(N_t,G*G_t);
g_t = 0:1:G_t-1;
g_r = 0:1:G_r-1;
% ******************************************************璁捐鐨勫埌AOA/AOD瀛楀吀鐭╅樀**********************************************

PSI_g = zeros(NN,G*GG);
A1 = sin(-pi/2+pi*g_t/G_t);
n_t_A = N_tt.'*A1;
for g = 1:G
    A_t_k1_g = exp(-1j*(2*pi*d)/lamda_g(1,g)*n_t_A)/sqrt(N_t);
    PSI_g(1:N_t*N_r,1+(g-1)*G_r*G_t:G_r*G_t*g) = kron(conj(A_t_k1_g),A_t_k1_g);
end
% ******************************************************灏嗗瓙杞芥尝鍒嗕负16缁勮幏寰楃殑AOA/AOD瀛楀吀鐭╅樀**********************************************
PSI_G1 = zeros(NN,G1*GG);
A_t_k2_g = zeros(N_t,G_t);
A_t_k0 = exp(-1j*(2*pi*d)/lamda_g1(1,1)*n_t_A)/sqrt(N_t);
for g = 1:G1
    A_t_k2_g = exp(-1j*(2*pi*d)/lamda_g1(1,g)*n_t_A)/sqrt(N_t);
    PSI_G1(1:N_t*N_r,1+(g-1)*G_r*G_t:G_r*G_t*g) = kron(conj(A_t_k2_g),A_t_k2_g);
end
% ******************************************************灏嗗瓙杞芥尝鍒嗕负4缁勮幏寰楃殑AOA/AOD瀛楀吀鐭╅樀**********************************************
PSI_G2 = zeros(NN,G2*GG);
A_t_k3_g = zeros(N_t,G_t);
for g = 1:G2
    A_t_k3_g = exp(-1j*(2*pi*d)/lamda_g2(1,g)*n_t_A)/sqrt(N_t);
    PSI_G2(1:N_t*N_r,1+(g-1)*G_r*G_t:G_r*G_t*g) = kron(conj(A_t_k3_g),A_t_k3_g);
end
PSI_K  = PSI_G1;
 max_e1 = 1;
KKKK = 0;
Rk = zeros(length(SNR1),K);
R1k = zeros(length(SNR1),K);
R4k = zeros(length(SNR1),K);
R16k = zeros(length(SNR1),K);
R3k = zeros(length(SNR1),K);
R11k = zeros(length(SNR1),K);
R14k = zeros(length(SNR1),K);
R161k = zeros(length(SNR1),K);
for a = 1:CC
    KKKK = KKKK + 1
    delay_L1= rand(1,L)*(Nc-1)*Ts; %均匀随机分布
    alphal = sqrt(HHH)*(randn(1,L)+1j*randn(1,L))/sqrt(2);
    AOA = -pi/2 + rand(1,L) * pi;%鍒拌揪瑙?
    AOD = -pi/2 + rand(1,L) * pi;%绂诲紑瑙?
    [W_tr,C_W,Fai1,F_tr] = TrainFW(L_t,L_r,N_r,M,N_t,NN);%生产Ftr和Wtr
    noise2 = zeros(M*L_r,1);
    for m = 1:M
        noise1 = (randn(N_r,1) + 1j*randn(N_r,1))/sqrt(2);
        wtr = W_tr(1:N_r,1+(m-1)*L_r:L_r*m);
        noise2(1 + (m-1)*L_r:L_r*m,1)= wtr'* noise1;
    end
    noise = noise2;
  
    for s = 1:length(SNR)
        L1= 14;
        L_t1 = L_r;
        Fai = sqrt(snr(1,s)* P_Rc)*Fai1; 
        R_W = zeros(M*L_r*L_t,G*GG);
        RR_W = zeros(M*L_r*L_t,G1*GG);
        RR_W1 = zeros(M*L_r*L_t,G2*GG);
        for g = 1:G
            R_W(1:M*L_r*L_t,1+(g-1)*GG:GG*g)= Fai* PSI_g(1:NN,1+(g-1)*GG:GG*g);
        end
        for g = 1:G1
            RR_W(1:M*L_r*L_t,1+(g-1)*GG:GG*g)= Fai* PSI_G1(1:NN,1+(g-1)*GG:GG*g);
        end
        for g = 1:G2
            RR_W1(1:M*L_r*L_t,1+(g-1)*GG:GG*g)= Fai* PSI_G2(1:NN,1+(g-1)*GG:GG*g);
        end
        HK = zeros(NN,K);
        H_K = zeros(N_r,N_t*K);%K涓瓙杞芥尝鐨勪俊閬撶煩闃?
        A_t = zeros(N_t,L);
        A_r = zeros(N_t,L);
        G_d_k = zeros(L,L*K);%璺緞澧炵泭涓嶱rc
        G_k = zeros(L,L);%瀵硅鐭╅樀瀵硅绾夸笂鐨勫厓绱犳槸姣忔潯璺緞涓婄殑澧炵泭
        Y = zeros(M*L_r*L_t,K);%Y涓烘帴鏀朵俊鍙凤
        y_K = zeros(1,K);%姣忎釜瀛愯浇娉㈢殑2鑼冩暟锛?
       % H_K

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
            A_t = exp(-1j*(2*pi*d)/lamda_k(1,k)*n_t_S_AOD)/sqrt(N_t);
            A_r = exp(-1j*(2*pi*d)/lamda_k(1,k)*n_t_S_AOA)/sqrt(N_r);
            for l=1:L
                G_k(l,l) =  G_d_k(l,1+(l-1)*K:K*l)*E_K.';%G_d_k鍌呴噷鍙跺彉鍖?
            end
            H_K(1:N_r,1+(k-1)*N_t:N_t*k) = A_r*G_k* A_t';
            HK(:,k)=reshape(H_K(1:N_r,1+(k-1)*N_t:N_t*k),[],1);
            Y(:,k) = Fai * HK(:,k) + noise;
            norm_Fro_HK(1,k)=(norm(H_K(1:N_r,1+(k-1)*N_t:N_t*k),'fro'))^2;
        end
        y_w = Y;
        %% ************************************************************璁＄畻鍚堥?鐨凟*************************************************
        MSE1 = trace(y_w'*y_w)/(K*L_r*M);
        MSE = MSE1;
       %% OMP
        for k = 1:K
            y_K(1,k) = Y(:,k)'*Y(:,k);
        end
        [A,Location] = sort(y_K);
        MAX_y_k = Location(K-K_p+1:K);
        [est_g,I] = My_est_g1(y_w,MAX_y_k,GG,K_p,MSE,max_e1,K,M,L_r,R_W,L_t); %%  鍒嗕负1缁勭殑淇￠亾浼拌
        [est_g16,I16] = My_est_g(y_w,MAX_y_k,GG,K_p,MSE,max_e1,K,G1,M,L_r,L_t,RR_W);%%  鍒嗕负16缁勭殑淇￠亾浼拌
        [est_g14,I4] = My_est_g4(y_w,MAX_y_k,GG,K_p,MSE,max_e1,K,G2,M,L_r,L_t,RR_W1);%%  鍒嗕负4缁勭殑
        [est_g3,I3] = LowerBound(HK,PSI_K,GG,L1,NN,K);%下界%  
        [est_k16,I161] = My_est_HK16(y_w,PSI_K,GG,L1,NN,K,Fai);%新方法：分十六组
        [est_k,I1] = My_est_HK1(y_w,PSI_K,GG,L1,NN,K,Fai,G_t,G_r,lamda_k);%新方法：分1组
        [est_k4,I14] = My_est_HK4(y_w,PSI_K,GG,L1,NN,K,Fai,G_t,G_r,lamda_k);%新方法：分4组
        [n_H_K1,n_H_K4,n_H_K16,n_H_K3,n_H_K11,n_H_K14,n_H_K161]...
    =Vec_hk(K,PSI_K,NN,GG,est_g,I,est_g16,I16,est_g14,I4,est_g3,I3,est_k16,I161,est_k,I1,est_k4,I14,MSE,max_e1);
    end
%%  选择链路
   A_t_k0 = exp(-1j*(2*pi*d)/lamda_g1(1,8)*n_t_A)/sqrt(N_t);  
   NN_H_K =  H_K(1:N_r,1+7*N_t:N_t*8);
   AA =  A_t_k0'*NN_H_K*A_t_k0;
   [S_Lr,S_Lt] = getRFChain(AA,L_r,L_t1,G_t,G_r);
   L_t1 = L_r;  
   AA1 = A_t_k0'*reshape(n_H_K1(:,8),N_r,N_t)*A_t_k0;
   AA4 = A_t_k0'*reshape(n_H_K4(:,8),N_r,N_t)*A_t_k0;
   AA16 = A_t_k0'*reshape(n_H_K16(:,8),N_r,N_t)*A_t_k0;
   AA3 = A_t_k0'*reshape(n_H_K3(:,8),N_r,N_t)*A_t_k0;
   NAA1 = A_t_k0'*reshape(n_H_K11(:,8),N_r,N_t)*A_t_k0;
   NAA4 = A_t_k0'*reshape(n_H_K14(:,8),N_r,N_t)*A_t_k0;
   NAA16 = A_t_k0'*reshape(n_H_K161(:,8),N_r,N_t)*A_t_k0;
  [S_Lr1,S_Lt1] = getRFChain(AA1,L_r,L_t1,G_t,G_r);
  [S_Lr4,S_Lt4] = getRFChain(AA4,L_r,L_t1,G_t,G_r);
  [S_Lr16,S_Lt16] = getRFChain(AA16,L_r,L_t1,G_t,G_r);
  [S_Lr3,S_Lt3] = getRFChain(AA3,L_r,L_t1,G_t,G_r);
  [N_S_Lr1,N_S_Lt1] = getRFChain(NAA1,L_r,L_t1,G_t,G_r);
  [N_S_Lr4,N_S_Lt4] = getRFChain(NAA4,L_r,L_t1,G_t,G_r);
  [N_S_Lr16,N_S_Lt16] = getRFChain(NAA16,L_r,L_t1,G_t,G_r);

 F_BB = eye(L_t1);
 A_t_k =  A_t_k0; 
 for ss = 1:length(SNR1)
     a_snr = snr1(1,ss);
     for k = 1:K
        HK1 = reshape(n_H_K1(:,k),N_r,N_t);
        HK4 = reshape(n_H_K4(:,k),N_r,N_t);
        HK16 = reshape(n_H_K16(:,k),N_r,N_t);
        HK3 = reshape(n_H_K3(:,k),N_r,N_t);

        NHK1 = reshape(n_H_K11(:,k),N_r,N_t);
        NHK4 = reshape(n_H_K14(:,k),N_r,N_t);
        NHK16 = reshape(n_H_K161(:,k),N_r,N_t);
        NN_H_K = H_K(1:N_r,1+(k-1)*N_t:N_t*k);
      %% 论文修改代码

        [WW] = checkGetSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,S_Lr,S_Lt,NN_H_K,NN_H_K);
        [WW1] = checkGetSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,S_Lr1,S_Lt1,NN_H_K,HK1);
        [WW4] = checkGetSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,S_Lr4,S_Lt4,NN_H_K,HK4);
        [WW16] = checkGetSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,S_Lr16,S_Lt16,NN_H_K,HK16);
        [WW3] =checkGetSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,S_Lr3,S_Lt3,NN_H_K,HK3);
        [WW11] = checkGetSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,N_S_Lr1,N_S_Lt1,NN_H_K,NHK1);
        [WW14] = checkGetSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,N_S_Lr4,N_S_Lt4,NN_H_K,NHK4);
        [WW161] = checkGetSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,N_S_Lr16,N_S_Lt16,NN_H_K,NHK16);
       %% 原论文代码
%             [WW] = getSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,S_Lr,S_Lt,NN_H_K,NN_H_K,k);
%             [WW1] = getSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,S_Lr1,S_Lt1,NN_H_K,HK1,k);
%             [WW4] = getSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,S_Lr4,S_Lt4,NN_H_K,HK4,k);
%             [WW16] = getSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,S_Lr16,S_Lt16,NN_H_K,HK16,k);
%             [WW3] = getSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,S_Lr3,S_Lt3,NN_H_K,HK3,k);
%             [WW11] = getSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,N_S_Lr1,N_S_Lt1,NN_H_K,NHK1,k);
%             [WW14] = getSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,N_S_Lr4,N_S_Lt4,NN_H_K,NHK4,k);
%             [WW161] = getSpEf(A_t_k,L_t1,F_BB,P_Rc,a_snr,N_S_Lr16,N_S_Lt16,NN_H_K,NHK16,k);
        Rk(ss,k) = WW  +  Rk(ss,k);
        R1k(ss,k) = WW1 + R1k(ss,k) ;
        R4k(ss,k) = WW4 +  R4k(ss,k) ;
        R16k(ss,k) = WW16 + R16k(ss,k);
        R3k(ss,k) = WW3 +  R3k(ss,k) ;
        R11k(ss,k) = WW11 + R11k(ss,k);
        R14k(ss,k) = WW14 + R14k(ss,k) ;
        R161k(ss,k) = WW161 + R161k(ss,k);
        %% 记录每个子载波
        RR(ss,k) = WW ;
        RR1(ss,k) = WW1 ;
        RR4(ss,k) = WW4 ;
        RR16(ss,k) = WW16;
        RR3(ss,k) = WW3 ;
        RR11(ss,k) = WW11 ;
        RR14(ss,k) = WW14 ;
        RR161(ss,k) = WW161 ;
     end
     ORR(a,ss) =sum(RR(ss,:))/(K+1);
     ORR1(a,ss) =sum(RR1(ss,:))/(K+1);
     ORR4(a,ss) =sum(RR4(ss,:))/(K+1);
     ORR16(a,ss) =sum(RR16(ss,:))/(K+1);
     LRR3(a,ss) =sum(RR3(ss,:))/(K+1);
     NRR11(a,ss) =sum(RR11(ss,:))/(K+1);
     NRR14(a,ss) =sum(RR14(ss,:))/(K+1);
     NRR16(a,ss) =sum(RR161(ss,:))/(K+1);
 end
end
Rk =Rk/CC;
R1k = R1k/CC ;
R4k =  R4k /CC;
R16k = R16k/CC;
R3k = R3k /CC;
R11k = R11k/CC;
R14k = R14k/CC ;
R161k = R161k/CC;

for s1 = 1 : length(SNR1)
    S_R(1,s1) = norm(sum(ORR(:,s1))/CC);
    S_R1(1,s1) = norm(sum(ORR1(:,s1))/CC);
    S_R4(1,s1) = norm(sum(ORR4(:,s1))/CC);
    S_R16(1,s1) = norm(sum(ORR16(:,s1))/CC);
    S_R3(1,s1) = norm(sum(LRR3(:,s1))/CC);
    S_R11(1,s1) = norm(sum(NRR11(:,s1))/CC);
    S_R14(1,s1) = norm(sum(NRR14(:,s1))/CC);
    S_R161(1,s1) = norm(sum(NRR16(:,s1))/CC);
    
%     Nmse_nHK(1,s1) = sum(nMSE_nHK(:,s1))/CC;%OMP分为1组
%     Nmse_nHK7(1,s1) = sum(nMSE_nHK7(:,s1))/CC;%OMP分为16组
%     Nmse_nHK14(1,s1) = sum(nMSE_nHK14(:,s1))/CC;%OMP分为4组?
%     Nmse_nHK1(1,s1) = sum(nMSE_nHK1(:,s1))/CC;%新方法分为1组
%     Nmse_nHK4(1,s1) = sum(nMSE_nHK4(:,s1))/CC;%新方法分为16组
%     Nmse_nHK22(1,s1) = sum(nMSE_nHK22(:,s1))/CC;%新方法分为4组
%     Nmse_nHK3(1,s1) = sum(nMSE_nHK3(:,s1))/CC;%下界
%    v(1,s1) = sum(Vc_rb(:,s1))/CC;
end
SS_RR = zeros(length(SNR1),7);
SS_RR(:,1) = S_R1;
SS_RR(:,2) = S_R4;
SS_RR(:,3) = S_R16;
SS_RR(:,4) = S_R;
SS_RR(:,5) = S_R11;
SS_RR(:,6) = S_R14;
SS_RR(:,7) = S_R161;
S_R_k = zeros(length(SNR1),3);
S_R_k(:,1) = abs(R11k(:,1));
S_R_k(:,2) = abs(R11k(:,8));
S_R_k(:,3) = abs(R11k(:,16));
figure(1)
plot(SNR1,S_R,'g-o','Markersize',7,'Linewidth',1);%SS-SW-OMP+Th，P=1
hold on;
plot(SNR1,S_R1,'g-*','Markersize',7,'Linewidth',1);%SS-SW-OMP+Th，P=1
hold on;
plot(SNR1,S_R4,'g->','Markersize',7,'Linewidth',1);%SS-SW-OMP+Th，P=1
hold on;
plot(SNR1,S_R16,'g-<','Markersize',7,'Linewidth',1);%SS-SW-OMP+Th，P=1
hold on;
plot(SNR1,S_R3,'b-*','Markersize',7,'Linewidth',1);%SS-SW-OMP+Th，P=1
hold on;
plot(SNR1,S_R11,'r-o','Markersize',7,'Linewidth',1);%SS-SW-OMP+Th，P=1
hold on;
plot(SNR1,S_R14,'r->','Markersize',7,'Linewidth',1);%SS-SW-OMP+Th，P=1
hold on;
plot(SNR1,S_R161,'r-<','Markersize',7,'Linewidth',1);%SS-SW-OMP+Th，P=1

 figure(2)
 plot(SNR1,R11k(:,1),'r-<','Markersize',7,'Linewidth',1);%SS-SW-OMP+Th，P=1
 hold on;
 plot(SNR1,R11k(:,8),'b-<','Markersize',7,'Linewidth',1);
 hold on;
 plot(SNR1,R11k(:,16),'g-<','Markersize',7,'Linewidth',1);
% plot(SNR,Nmse_nHK,'g-o','Markersize',7,'Linewidth',3);%SS-SW-OMP+Th，P=1
% hold on;
% plot(SNR,Nmse_nHK1,'r-o','Markersize',7,'Linewidth',3);%NewMethod，P=1
% hold on;
% plot(SNR,Nmse_nHK3,'b-*','Markersize',7,'Linewidth',3);%LowerBound
% hold on;
% plot(SNR,Nmse_nHK4,'r->','Markersize',7,'Linewidth',3);%NewMethod，P=16
% hold on;
% plot(SNR,Nmse_nHK7,'g->','Markersize',7,'Linewidth',3);%SS-SW-OMP+Th，P=16
% hold on;
% plot(SNR,Nmse_nHK22,'r-<','Markersize',7,'Linewidth',3);%NewMethod，P=4
% hold on;
% plot(SNR,Nmse_nHK14,'g-<','Markersize',7,'Linewidth',3);%SS-SW-OMP+Th，P=4
% hold on;
% plot(SNR,v,'r-*','Markersize',7,'Linewidth',1);
xlabel('SNR'),ylabel('NMSE[dB]');
% set(gca,'XLim',[32,96]);
legend('SS-SW-OMP+Th，P=1','NewMethod，P=1','LowerBound','NewMethod，P=16','SS-SW-OMP+Th，P=16','NewMethod，P=4','SS-SW-OMP+Th，P=4');

% 
