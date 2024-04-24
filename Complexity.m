clear all; 
close all;
N = 40; 
L = 4;
M_t = 32;
M_r = 32;
G_t = 64;
G_r = 64;
D = 20;
K = 16;
S = 14;
Ps = [1; 2; 4; 8; 16];
Comp = zeros(5, 4);
for ii = 1 : 5
    P = Ps(ii, 1);
Comp(ii, 1) = P * N * L * M_r * M_t * G_r * G_t + D * K * N * L * G_r * G_t + P * D^3 * (N * L + D) + 2 * K * N * L * D^2 + K * M_r * M_t;
Comp(ii, 2) =  D * K * N * L * G_r * G_t + P * D^3 * (N * L + D) + 2 * K * N * L * D^2 + K * M_r * M_t;
Comp(ii, 3) = P * N * L * M_r * M_t * G_r * G_t + P * S * N * L * G_r * G_t + S^2 * (N * L + S) * (P * S + K) + N * L * S * (P * S + K) + P * N * L * S^2 + K * M_r * M_t * S;
Comp(ii, 4) = P * S * N * L * G_r * G_t + S^2 * (N * L + S) * (P * S + K) + N * L * S * (P * S + K) + P * N * L * S^2 + K * M_r * M_t * S;
end
semilogy(Comp(:, 1),'r-')
hold on
semilogy(Comp(:, 2),'r--o')
semilogy(Comp(:, 3),'b-')
semilogy(Comp(:, 4),'b--o')
x = [Ps, Comp];