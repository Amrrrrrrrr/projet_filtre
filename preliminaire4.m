clear; close all; clc;

data = load("fcno04fz.mat");  
nom_champ = fieldnames(data);
x = data.(nom_champ{1}); 

Fe = 8000;             
Lx = length(x);         
t = (0:Lx-1)'/Fe;      

N = 256;                
R = N/2;                
w = hanning(N);         

% Zero-Padding
nb_zeros_debut = R; 
x_pad = [zeros(nb_zeros_debut, 1); x; zeros(N, 1)];
L_pad = length(x_pad);
y_reconstruit_pad = zeros(L_pad, 1);

% Boucle de Traitement 
n_trames = floor((L_pad - N) / R) + 1;

for k = 1:n_trames
    idx_debut = (k-1)*R + 1;
    idx_fin = idx_debut + N - 1;
    
    segment = x_pad(idx_debut:idx_fin);
    segment_fenetre = segment .* w; 
    
    segment_traite = segment_fenetre; 
    
    % Addition-Recouvrement
    y_reconstruit_pad(idx_debut:idx_fin) = y_reconstruit_pad(idx_debut:idx_fin) + segment_traite;
end

% Calcul de la fenêtre de normalisation globale
W_norm = zeros(L_pad, 1);
for k = 1:n_trames
    idx_debut = (k-1)*R + 1;
    idx_fin = idx_debut + N - 1;
    W_norm(idx_debut:idx_fin) = W_norm(idx_debut:idx_fin) + w;
end
W_norm(W_norm < 1e-10) = 1; 

% Normalisation
y_reconstruit_pad = y_reconstruit_pad ./ W_norm;

% Suppression du padding 
y_final = y_reconstruit_pad(nb_zeros_debut+1 : nb_zeros_debut+Lx);

erreur = x - y_final;

% Affichage du signal original
figure;
plot(t, x, 'b');
title('Signal original');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;

% Affichage du signal reconstruit
figure;
plot(t, y_final, 'r');
title('Signal reconstruit (après addition-recouvrement)');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;

% Affichage du signal d’erreur
figure;
plot(t, erreur);
title('Erreur de reconstruction (x - y_{final})');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;


disp(['Erreur maximale de reconstruction : ' num2str(max(abs(erreur)))]);