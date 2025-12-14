clear; close all; clc;

fs = 8000;
data = load('fcno04fz.mat');  
s = data.fcno04fz;  
s = s(:);
len = length(s);


% Tramage
N   = 256;                      % taille de trame
hop = N/2;                      % recouvrement 50 %
w   = hamming(N,'periodic');    % fenêtre de Hamming
nfft = 512;

           

% RSB cibles et nombre de réalisations
RSB_targets = [5, 10, 15]; 
nRealz      = 20; 

% Stockage résultats
gain_RSB_moy = zeros(length(RSB_targets),1);
RSB_before   = zeros(length(RSB_targets),nRealz);
RSB_after    = zeros(length(RSB_targets),nRealz);

rng(0);  % pour reproductibilité

for idxRSB = 1:length(RSB_targets)
    targetRSB = RSB_targets(idxRSB);

    for r = 1:nRealz

        % Génération bruit blanc au bon RSB
        b = randn(len,1);

        P_s = mean(s.^2);
        desired_Pb = P_s / (10^(targetRSB/10)); 
        b = b * sqrt(desired_Pb / mean(b.^2));  
        
        y = s + b;  
        
        % Estimation de la puissance du bruit (bruit blanc)
        P_b_est = mean(b.^2);                
        
        % Analyse + soustraction spectrale trame par trame
        frame_idx = 1:hop:(len-N+1);        
        out_buf = zeros(len + N,1);
        win_sum = zeros(len + N,1);
        
        for i = 1:length(frame_idx)
            n0 = frame_idx(i);
            frame = y(n0:n0+N-1) .* w;
            
            Y = fft(frame, nfft);
            magY   = abs(Y);
            phaseY = angle(Y);
            
            % Spectre de puissance du signal bruité
            P_Y = magY.^2;
            
            % Estimation de la puissance du bruit dans le domaine fréquentiel
            P_B = (P_b_est * sum(w.^2)) * ones(size(P_Y));
            
            % Soustraction spectrale 
            P_S_hat = P_Y - P_B;
            P_S_hat = max(P_S_hat, 0); 
            
            % Module estimé du spectre de la parole
            magS_hat = sqrt(P_S_hat);

            % Reconstruction fréquentielle 
            S_hat_full = magS_hat .* exp(1j * phaseY);
            
            % Retour temporel
            s_frame_rec = real(ifft(S_hat_full, nfft));
            s_frame_rec = s_frame_rec(1:N) .* w; 
            
            out_buf(n0:n0+N-1)   = out_buf(n0:n0+N-1)   + s_frame_rec;
            win_sum(n0:n0+N-1)   = win_sum(n0:n0+N-1)   + (w.^2);
            
            % Affichage d'une trame représentative 
            k_affiche = 10; 
            if i == k_affiche && r == 1 && idxRSB == 1
                tvec = (0:N-1)/fs;
                figure(1); clf;
                
                % Temps
                subplot(2,1,1);
                plot(tvec, frame, 'r'); hold on;
                plot(tvec, s(n0:n0+N-1).*w, 'g');
                plot(tvec, s_frame_rec, 'b');
                legend('Bruité','Propre','Rehaussé');
                xlabel('Temps (s)');
                ylabel('Amplitude');
                title(['Trame ' num2str(i) ' (temporel)']);
                grid on;
                
                % Spectre
                subplot(2,1,2);
                f = (0:(nfft/2))*(fs/nfft);                 
                idx_half = 1:(nfft/2 + 1);
                
                plot(f, 20*log10(abs(magY(idx_half))+eps), 'r'); hold on;
                
                S_clean = fft(s(n0:n0+N-1).*w, nfft);
                plot(f, 20*log10(abs(S_clean(idx_half)) + eps), 'g');
                
                plot(f, 20*log10(magS_hat(idx_half)+eps), 'b');
                
                legend('|Y|','|S|','|S_{est}|');
                xlabel('Fréquence (Hz)');
                ylabel('Amplitude (dB)');
                title('Spectres (unilatéraux)');
                grid on;
                drawnow;
            end
        end
        
        % Normalisation Overlap-Add
        idx_nonzero = win_sum > 1e-8;
        out_buf(idx_nonzero) = out_buf(idx_nonzero) ./ win_sum(idx_nonzero);
        y_enh = out_buf(1:len);
        
        % Calcul des RSB avant / après
        noise_before = y - s;
        SNR_before = 10*log10( sum(s.^2) / sum(noise_before.^2) );
        
        error_after = y_enh - s;
        SNR_after  = 10*log10( sum(s.^2) / sum(error_after.^2) );
        
        RSB_before(idxRSB, r) = SNR_before;
        RSB_after(idxRSB,  r) = SNR_after;
    end

    gain_RSB_moy(idxRSB) = mean(RSB_after(idxRSB,:) - RSB_before(idxRSB,:));
end
soundsc(y_enh, fs)
% Affichage des résultats RSB
fprintf('\nTableau des résultats :\n');
fprintf('Cible (dB) | RSB Sortie (dB) | Gain (dB)\n');
for k = 1:length(RSB_targets)
    fprintf('   %3d     |     %6.3f      |  %6.3f\n', ...
        RSB_targets(k), mean(RSB_after(k,:)), gain_RSB_moy(k));
end

t_global = (0:len-1)/fs;

figure;
subplot(2,1,1);
plot(t_global, y_enh);
xlabel('Temps (s)');
ylabel('Amplitude');
title('Signal de parole rehaussé (temporel)');
grid on;

subplot(2,1,2);
spectrogram(y_enh, hamming(256), 128, 512, fs, 'yaxis');
title('Spectrogramme du signal de parole rehaussé');
xlabel('Temps (s)');
ylabel('Fréquence (kHz)');
