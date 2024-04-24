clear all 
close all;
par.runId = 1;
rng(par.runId);
%% ****************************************************璁粌甯х殑涓暟**************************************�??

%% ************************************************鎺ユ敹澶╃嚎鍜屽彂灏勫ぉ绾跨殑涓暟*******************************�??
NN_t = [40];
NN_r = [40];
% NN_t = 24;
% NN_r = 24;
%% *******************************************************瀛楀�?鐨勫ぇ灏?***********************************
G_t = 64;
G_r = 64;
GG = G_t*G_r;
g_t = 0:1:G_t-1;
g_r = 0:1:G_r-1;
%% *****************************************************瀛愯浇娉㈢殑涓�?*************************************�??
K = 16;
K_p = K/4;   %鎵惧埌K_p涓娇寰梱[k]鐨凩LF鍑芥暟鏈?ぇ鐨勫瓙杞芥尝锛?
G = 1;%灏咾涓瓙杞芥尝鍒嗘垚涓ょ粍
G1 = 16;
G2 = 4;
%% *******************************************************甯﹀澶у�?**************************************
f_c = 60*10^9;%杞芥尝鐨勯鐜囷�?
B = f_c*0.2;
c = 3.0*10^8;%閫熷�?
k = 0:1:K-1;
f = (k-K/2)*B/(K+1);
lamda_k = c./(f+f_c);%绗琄涓瓙杞芥尝鐨勬尝闀匡紱
lamda_g = zeros(1,G);
%%分组后的每组子载波的波长
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
d = c/(2*f_c);%%澶╃嚎涔嬮棿鐨勮窛绂?
L = 4;%璺緞鏁?
N0  = 1;
Nc = 4;
%% ******************************************************淇￠亾澧炵泭鐨勮绠椾互鍙婂疄闄呬俊閬撶殑琛ㄧず**************************************
CC = 2;
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
%% 计算信道的F范数的均�?
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
%% 参数设置
KKKK = 0;
SNR = 20;
snr = 10^(SNR/10);
L_t = 1;%鍙戦?绔疪F鐨勪釜鏁帮紱
L_r = 4;%鎺ユ敹绔疪F鐨勪釜鏁帮紱
max_e1 = 1;
A1 = sin(-pi/2+pi*g_t/G_t);
for s = 1:length(NN_t)
    N_t = NN_t(1,s);
    N_r = NN_r(1,s);
    NN = N_t*N_r;
    M = NN / L_r;   
    n_t = 0:1:N_t-1;
    n_r = 0:1:N_r-1;
    A_t_k1 = zeros(N_t,N_t);
    for n = 1:N_t
        A_t_k1(n,:) = exp(-1j*2*pi*(n-1)*n_t/N_t)/sqrt(N_t);
    end
    [W_tr,C_W,Fai1,F_tr] = TrainFWDFT(L_t,L_r,N_r,M,N_t,NN,A_t_k1);%生产Ftr和Wtr
    KKKK = KKKK + 1
    % A_t_k1_g = zeros(N_t,G_t);%绂诲紑瑙?

    %% ******************************************************灏嗗瓙杞芥尝鍒嗕�?2缁勮幏寰楃殑AOA/AOD瀛楀�?鐭╅�?**********************************************
    PSI_g = zeros(NN,G*GG);
    A_t_k1_g=zeros(N_t,G_t);
    N_tt = 0:N_t-1 ;
    n_t_A = N_tt.'*A1;
    for g = 1:G
        A_t_k1_g = exp(-1j*(2*pi*d)/lamda_g(1,g)*n_t_A)/sqrt(N_t);
        PSI_g(1:N_t*N_r,1+(g-1)*G_r*G_t:G_r*G_t*g) = kron(conj(A_t_k1_g),A_t_k1_g);
    end
    %% ******************************************************灏嗗瓙杞芥尝鍒嗕�?16缁勮幏寰楃殑AOA/AOD瀛楀�?鐭╅�?**********************************************
    PSI_G1 = zeros(NN,G1*GG);
    A_t_k2_g = zeros(N_t,G_t);
    for g = 1:G1
        A_t_k2_g = exp(-1j*(2*pi*d)/lamda_g1(1,g)*n_t_A)/sqrt(N_t);
        PSI_G1(1:N_t*N_r,1+(g-1)*G_r*G_t:G_r*G_t*g) = kron(conj(A_t_k2_g),A_t_k2_g);
    end
    %% ******************************************************灏嗗瓙杞芥尝鍒嗕�?4缁勮幏寰楃殑AOA/AOD瀛楀�?鐭╅�?**********************************************
    PSI_G2 = zeros(NN,G2*GG);
    A_t_k3_g = zeros(N_t,G_t);
    for g = 1:G2
        A_t_k3_g = exp(-1j*(2*pi*d)/lamda_g2(1,g)*n_t_A)/sqrt(N_t);
        PSI_G2(1:N_t*N_r,1+(g-1)*G_r*G_t:G_r*G_t*g) = kron(conj(A_t_k3_g),A_t_k3_g);
    end
    %% *************************************************************娴嬮噺鐭╅樀******************************************

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
        N_Q = 2;%鐩哥Щ鍣ㄧ殑閲忓寲�??
        HK = zeros(NN,K);
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
        H_K = zeros(N_r,N_t*K);%K涓瓙杞芥尝鐨勪俊閬撶煩�??
        G_d_k = zeros(L,L*K);%璺緞澧炵泭涓嶱rc
        G_k = zeros(L,L);%瀵硅鐭╅樀瀵硅绾夸笂鐨勫厓绱犳槸姣忔潯璺緞涓婄殑澧炵泭
        Y = zeros(M*L_r*L_t,K);%Y涓烘帴鏀朵俊鍙凤�?
        y_w = zeros(M*L_r*L_t,K);
        y_K = zeros(1,K);%姣忎釜�?�愯浇娉㈢殑2鑼冩暟锛?
        %% H_K
        delay_L1= rand(1,L)*(Nc-1)*Ts; 
        alphal = sqrt(HHH)*(randn(1,L)+1j*randn(1,L))/sqrt(2);
        AOA = -pi/2 + rand(1,L) * pi;%鍒拌揪瑙?
        AOD = -pi/2 + rand(1,L) * pi;%绂诲紑瑙?
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
                 G_k(l,l) =  G_d_k(l,1+(l-1)*K:K*l)*E_K.';%G_d_k鍌呴噷鍙跺彉�??
            end
            H_K(1:N_r,1+(k-1)*N_t:N_t*k) = A_r*G_k* A_t';
            HK(:,k)=reshape(H_K(1:N_r,1+(k-1)*N_t:N_t*k),[],1);
            Y(:,k) = Fai * HK(:,k) + noise;
            norm_Fro_HK(1,k)=(norm(H_K(1:N_r,1+(k-1)*N_t:N_t*k),'fro'))^2;
        end
        y_w = Y;
        %% ************************************************************璁＄畻鍚�??鐨凟*************************************************
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
        [est_g,I] = My_est_g1DFT(y_w,MAX_y_k,GG,K_p,MSE,max_e1,K,M,L_r,R_W,L_t); %%  鍒嗕�?1缁勭殑淇￠亾浼拌�?
        [est_g16,I16] = My_est_gDFT(y_w,MAX_y_k,GG,K_p,MSE,max_e1,K,G1,M,L_r,L_t,RR_W);%%  鍒嗕�?16缁勭殑淇￠亾浼拌�?
        [est_g14,I4] = My_est_g4DFT(y_w,MAX_y_k,GG,K_p,MSE,max_e1,K,G2,M,L_r,L_t,RR_W1);%%  鍒嗕�?4缁勭�?
        [est_g3,I3] = LowerBound(HK,PSI_K,GG,L1,NN,K);%下界%  
        [est_k16,I161] = My_est_HK16(y_w,PSI_K,GG,L1,NN,K,Fai);%新方法：分十六组
        [est_k,I1] = My_est_HK1(y_w,PSI_K,GG,L1,NN,K,Fai,G_t,G_r,lamda_k);%新方法：�?1�?
        [est_k4,I14] = My_est_HK4(y_w,PSI_K,GG,L1,NN,K,Fai,G_t,G_r,lamda_k);%新方法：�?4�?
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
            AA_K(1+M*L_r:2*M*L_r,1+S1:2*S1) =  R_A_k;
            J_k = AA_K.'*C9*AA_K;
            inv_J_k = inv(J_k);
            for s1 = 1:S1
                JJ(s1,s1) = inv_J_k(s1,s1) + inv_J_k(s1 + S1,s1 + S1);
            end
            vv_k = trace(PSI_KK(:,I3(k,:))*JJ*PSI_KK(:,I3(k,:))');
            Sum_vv_k = Sum_vv_k + vv_k;
        end


          V_CRB(1,a) = Sum_vv_k/sum(norm_Fro_HK);
     
         v_crb(1,a) = 10*log10(V_CRB(1,a));
    end
    V_crb(1,s) = sum(v_crb)/CC;
end
SS_RR = zeros(length(NN_t),7);
SS_RR(:,4) = real(V_crb).';

