clear; close all; clc;
N = 2000;     
sigma = 2;     
b = sigma * randn(1, N);
% Autocorrélation  

[Rb, ~] = xcorr(b, 'biased');
[Rub, ~] = xcorr(b, 'unbiased');
figure;
plot(Rb, 'b'); 
hold on;
plot(Rub, 'r');
xline(0, '--k');
legend('Autocorr. biaisée','Autocorr. non biaisée');
xlabel('Décalage k');
ylabel('R_b(k)');
title('Autocorrélation estimée du bruit blanc');
grid on;

% Spectre de puissance
B = fft(b);
S_est = (abs(B).^2) / N;
fe = 10000;
f = (0:N-1) * (fe / N);
figure;
plot(f, 10*log10(S_est));     
hold on;
yline(10*log10(sigma^2), 'r'); 
xlabel('Fréquence (Hz)');
ylabel('DSP (dB)');
title('Spectre du bruit blanc gaussien');
legend('Spectre estimé', 'DSP théorique');
grid on;


% Périodogramme de Daniell 
M = 20;
S_daniell = movmean_simple(S_est, 2*M+1);

% Périodogramme de Bartlett 
L = 500;                     
K = floor(N / L);            % nombre de segments
S_bartlett = zeros(1, L);

for i = 1:K
    x_segment = b((i-1)*L + 1 : i*L);
    X_seg = fft(x_segment);
    P_seg = (abs(X_seg).^2) / L;
    S_bartlett = S_bartlett + P_seg;
end
S_bartlett = S_bartlett / K;
f_bartlett = (0:L-1) * (fe / L);

% Périodogramme de Welch 
[S_welch, f_welch] = Mon_Welch(b, 512, fe);

% Corrélogramme 
S_corr = abs(fft(Rb));                
f_corr = (0:length(S_corr)-1) * (fe / length(S_corr));


figure;

f_daniell = (0:length(S_daniell)-1) * (fe / length(S_daniell));
plot(f_daniell, 10*log10(S_daniell), 'm', 'LineWidth', 1.2); hold on;yline(10*log10(sigma^2), 'r--', 'DSP théorique');
xlabel('Fréquence (Hz)');
ylabel('DSP (dB)');
title('Périodogramme de Daniell');
legend('Daniell', 'DSP théorique', 'Location', 'southoutside');
grid on;

figure;

plot(f_bartlett, 10*log10(S_bartlett), 'g', 'LineWidth', 1.2); hold on;
yline(10*log10(sigma^2), 'r--', 'DSP théorique');
xlabel('Fréquence (Hz)');
ylabel('DSP (dB)');
title('Périodogramme de Bartlett');
legend('Bartlett', 'DSP théorique', 'Location', 'southoutside');
grid on;

figure;


plot(f_welch, 10*log10(S_welch), 'b', 'LineWidth', 1.2); hold on;
yline(10*log10(sigma^2), 'r--', 'DSP théorique');
xlabel('Fréquence (Hz)');
ylabel('DSP (dB)');
title('Périodogramme de Welch');
legend('Welch', 'DSP théorique', 'Location', 'southoutside');
grid on;

figure;

plot(f_corr, 10*log10(S_corr), 'c', 'LineWidth', 2.2); hold on;
yline(10*log10(sigma^2), 'r', 'DSP théorique');
xlabel('Fréquence (Hz)');
ylabel('DSP (dB)');
title('Corrélogramme');
legend('Corrélogramme', 'DSP théorique', 'Location', 'southoutside');
grid on;

R = 100;          % nombre de réalisations
p = 2;            

SF_vals = zeros(1, R);        
SF_welch_vals = zeros(1, R);  
SF_p_vals = zeros(1, R);      

for r = 1:R
    %  Génération du bruit blanc gaussien 
    b = sigma * randn(1, N);
    B = fft(b);
    S_est = (abs(B).^2) / N;

    % Platitude spectrale classique 
    geo_mean = exp(mean(log(S_est + eps))); 
    arith_mean = mean(S_est);
    SF_vals(r) = geo_mean / arith_mean;

    %  Platitude spectrale avec Mon_Welch 
    [S_welch, f] = Mon_Welch(b, 512, fe);
    geo_mean_welch = exp(mean(log(S_welch + eps)));
    arith_mean_welch = mean(S_welch);
    SF_welch_vals(r) = geo_mean_welch / arith_mean_welch;

    %  Platitude avec moyenne d’ordre p 
    mean_p = (mean(S_est.^(1/p)))^p;
    SF_p_vals(r) = mean_p / arith_mean;
end

%  Moyennes et écarts-types 
SF_mean = mean(SF_vals);
SF_std = std(SF_vals);
SFw_mean = mean(SF_welch_vals);
SFw_std = std(SF_welch_vals);
SFp_mean = mean(SF_p_vals);
SFp_std = std(SF_p_vals);
% Autocorrélation

[Rb, lags]  = xcorr(b, 'biased');
[Rub, ~]    = xcorr(b, 'unbiased');

% Autocorrélation théorique du bruit blanc gaussien N(0, sigma^2)
R_theo = sigma^2 * (lags == 0);
% Autocorrélations : biaisée, non biaisée et théorique sur une même figure

figure;

% 1) Autocorrélation biaisée
subplot(3,1,1);
plot(lags, Rb, 'b');
xline(0, '--k');
legend('Autocorr. biaisée', 'Location', 'best');
xlabel('Décalage k');
ylabel('R_b(k)');
title('Autocorrélation estimée du bruit blanc (biaisée)');
grid on;

% 2) Autocorrélation non biaisée
subplot(3,1,2);
plot(lags, Rub, 'r');
xline(0, '--k');
legend('Autocorr. non biaisée', 'Location', 'best');
xlabel('Décalage k');
ylabel('R_b(k)');
title('Autocorrélation estimée du bruit blanc (non biaisée)');
grid on;

% 3) Autocorrélation théorique
subplot(3,1,3);
plot(lags, R_theo, 'k');
xline(0, '--k');
legend('Autocorr. théorique', 'Location', 'best');
xlabel('Décalage k');
ylabel('R_{\rm theo}(k)');
title('Autocorrélation théorique du bruit blanc');
grid on;

% Affichage des résultats de platitude spectrale 

fprintf('\n=== Platitude spectrale ===\n');

fprintf('Périodogramme brut :moyenne = %.4f, écart-type = %.4f\n', SF_mean, SF_std);

fprintf('Méthode de Welch :moyenne = %.4f, écart-type = %.4f\n', SFw_mean, SFw_std);

fprintf('Platitude généralisée (p=%g) : moyenne = %.4f, écart-type = %.4f\n', p, SFp_mean, SFp_std);




