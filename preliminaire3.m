clear; close all; clc;

load('fcno01fz.mat');
x = fcno01fz;     
Fe = 8000;        
N = length(x);
t = (0:N-1)/Fe;   

k0_list = [1 10 100];

for i = 1:length(k0_list)
    k0 = k0_list(i);

    h = [1 zeros(1, k0-1) 1];   

    y = filter(h, 1, x);

    [H, w] = freqz(h, 1, 1024);  
    f = w * 8000 / (2*pi);

    figure;
    stem(0:length(h)-1, h, 'filled');
    xlabel('k (échantillons)');
    ylabel('h(k)');
    title(['Réponse impulsionnelle du filtre (k0 = ' num2str(k0) ')']);
    grid on;

    figure;
    subplot(2,1,1);
    plot(t, y);
    xlabel('Temps (s)');
    ylabel('Amplitude');
    title(['Signal filtré - h(k) = δ(k) + δ(k - ' num2str(k0) ')']);
    grid on;

    subplot(2,1,2);
    spectrogram(y, 256, 200, 256, Fe, 'yaxis');
    title(['Spectrogramme du signal filtré (k0 = ' num2str(k0) ')']);
    ylim([0 4]);
    
    figure;
    zplane(h, 1);
    title(['Pôles et zéros du filtre (k0 = ' num2str(k0) ')']);
    
    figure;
    subplot(3,1,1);
    plot(f, abs(H));
    xlabel('Fréquence (Hz)');
    ylabel('|H(f)|');
    title(['Module de la réponse fréquentielle (k0 = ' num2str(k0) ')']);
    grid on;
    
    subplot(3,1,2);
    plot(f, angle(H));
    xlabel('Fréquence (Hz)');
    ylabel('Phase (rad)');
    title(['Phase du filtre (k0 = ' num2str(k0) ')']);
    grid on;

end

figure;
subplot(2,1,1);
plot(t, x);
xlabel('Temps (s)');
ylabel('Amplitude');
title('Signal de parole original');
grid on;

subplot(2,1,2);
spectrogram(x, 256, 200, 256, Fe, 'yaxis');
title('Spectrogramme du signal original');
ylim([0 4]);

