clear; close all; clc;

load('fcno01fz.mat');   
x = fcno01fz; 
Fe = 8000;
N = length(x);        
t = (0:N-1)/Fe;      
                    
RSB_list = [5 10 15];

for i = 1:length(RSB_list)
    RSB_dB = RSB_list(i);
    [x_bruite, bruit, RSB_reel] = bruiter_signal(x, RSB_dB);
    
    fprintf('RSB visé = %d dB, RSB réel obtenu = %.2f dB\n', RSB_dB, RSB_reel);
    
    figure;
    subplot(2,1,1);
    plot(t, x_bruite);
    xlabel('Temps (s)');
    ylabel('Amplitude');
    title(['Signal bruité - RSB = ' num2str(RSB_dB) ' dB']);
    grid on;

    subplot(2,1,2);
    spectrogram(x_bruite, 256, 200, 256, Fe, 'yaxis');
    title(['Spectrogramme du signal bruité - RSB = ' num2str(RSB_dB) ' dB']);
    ylim([0 8]);
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
ylim([0 8]); 

function [x_bruite, bruit, RSB_reel] = bruiter_signal(x, RSB_dB)
    Ps = mean(x.^2);                         
    Pb = Ps / (10^(RSB_dB / 10));            
    bruit = sqrt(Pb) * randn(size(x));      
    x_bruite = x + bruit;                   
    RSB_reel = 10*log10(mean(x.^2)/mean(bruit.^2));  
end