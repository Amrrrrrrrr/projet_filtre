clear; close all; clc;

fs = 8000;
t = 0:1/fs:2; 
data = load('fcno04fz.mat');  
s = data.fcno04fz;  

s = s(:);
len = length(s);

% Tram
N = 256;                
hop = N/2;             
w = hamming(N,'periodic');
nfft = 512;

% Paramètres soustraction spectrale
alpha = 1.0;            
spectral_floor = 0.002; 
beta = 0.02;            

% RSB (SNR) cibles et nombre réalisations
RSB_targets = [5, 10, 15]; 
nRealz = 3; % Réduit pour le test rapide (remettre à 30)

% Stockage résultats
gain_RSB_moy = zeros(length(RSB_targets),1);
RSB_before = zeros(length(RSB_targets),nRealz);
RSB_after  = zeros(length(RSB_targets),nRealz);

rng(0); 

for idxRSB = 1:length(RSB_targets)
    targetRSB = RSB_targets(idxRSB);
    for r = 1:nRealz
        % Générer bruit blanc
        b = randn(len,1);
        
        % Ajustement RSB
        P_s = mean(s.^2);
        desired_Pb = P_s / (10^(targetRSB/10)); 
        b = b * sqrt(desired_Pb / mean(b.^2)); 
        
        y = s + b;
        
        %% Estimation du spectre du bruit (Approximation puissance FFT)
        % Note : Comparer une puissance temporelle moyenne à des bins FFT
        % nécessite normalement un facteur d'échelle, mais on garde votre logique.
        P_b_est = mean(b.^2);                
        
        %% Analyse trame par trame
        frame_idx = 1:hop:(len-N+1);        
        out_buf = zeros(len + N,1);
        win_sum = zeros(len + N,1);
        
        for i = 1:length(frame_idx)
            n0 = frame_idx(i);
            frame = y(n0:n0+N-1) .* w;
            
            Y = fft(frame, nfft);
            magY = abs(Y);
            phaseY = angle(Y);
            
            % Soustraction (en puissance)
            P_Y = magY.^2;
            
            % Correction d'échelle simple pour l'estimation de bruit
            % (Le bruit estimé doit être à l'échelle de la FFT)
            P_B = (P_b_est * N) * ones(size(P_Y)); % Facteur N approximatif ajouté
            
            P_S_hat = P_Y - alpha * P_B;
            P_S_hat = max(P_S_hat, 0); % Half-wave rectification
            
            magS_hat = sqrt(P_S_hat);
            
            % --- CORRECTION RECONSTRUCTION ---
            % Puisque magS_hat et phaseY sont déjà sur 512 points (spectre complet)
            % on recombine et on fait l'IFFT directement.
            S_hat_full = magS_hat .* exp(1j * phaseY);
            
            s_frame_rec = real(ifft(S_hat_full, nfft));
            s_frame_rec = s_frame_rec(1:N) .* w; 
            
            % Overlap-add
            out_buf(n0:n0+N-1) = out_buf(n0:n0+N-1) + s_frame_rec;
            win_sum(n0:n0+N-1) = win_sum(n0:n0+N-1) + (w.^2);
            
            % --- CORRECTION PLOT ---
            k_affiche = 10; 
            if i == k_affiche && r == 1 % Affiche seulement pour la 1ère réal
                tvec = (0:N-1)/fs;
                figure(1); clf;
                
                subplot(2,1,1);
                plot(tvec, frame, 'r'); hold on;
                plot(tvec, s(n0:n0+N-1).*w, 'g');
                plot(tvec, s_frame_rec, 'b');
                legend('Bruité','Propre','Rehaussé');
                title(['Trame ' num2str(i) ' (Temporel)']);
                
                subplot(2,1,2);
                f = (0:(nfft/2))*(fs/nfft); % Vecteur de taille 257
                
                % On limite les données à nfft/2 + 1 pour correspondre à f
                idx_half = 1:(nfft/2 + 1);
                
                plot(f, 20*log10(abs(magY(idx_half))+eps), 'r'); hold on;
                
                S_clean = fft(s(n0:n0+N-1).*w, nfft);
                plot(f, 20*log10(abs(S_clean(idx_half)) + eps), 'g');
                
                plot(f, 20*log10(magS_hat(idx_half)+eps), 'b');
                
                legend('|Y|','|S|','|S_{est}|');
                title('Spectre (Unilatéral)');
                xlabel('Fréquence (Hz)');
                drawnow;
            end
        end
        
        % Normalisation overlap-add
        idx_nonzero = win_sum > 1e-8;
        out_buf(idx_nonzero) = out_buf(idx_nonzero) ./ win_sum(idx_nonzero);
        y_enh = out_buf(1:len);
        
        %% Calcul RSB
        noise_before = y - s;
        SNR_before = 10*log10( sum(s.^2) / sum(noise_before.^2) );
        
        error_after = y_enh - s;
        SNR_after = 10*log10( sum(s.^2) / sum(error_after.^2) );
        
        RSB_before(idxRSB, r) = SNR_before;
        RSB_after(idxRSB,  r) = SNR_after;
    end
    gain_RSB_moy(idxRSB) = mean(RSB_after(idxRSB,:) - RSB_before(idxRSB,:));
end

% Affichage résultats
fprintf('\nTableau des résultats :\n');
fprintf('Cible (dB) | RSB Sortie (dB) | Gain (dB)\n');
for k=1:length(RSB_targets)
    fprintf('   %3d     |     %6.3f      |  %6.3f\n', ...
        RSB_targets(k), mean(RSB_after(k,:)), gain_RSB_moy(k));
end

figure;
spectrogram(y_enh, hamming(256), 128, 512, fs, 'yaxis');
title('Spectrogramme Final');