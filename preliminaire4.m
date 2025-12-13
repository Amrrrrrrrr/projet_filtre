clear; close all; clc;

%% 1. Chargement et Paramètres
data = load("fcno04fz.mat");  
nom_champ = fieldnames(data);
x = data.(nom_champ{1}); % Signal de parole chargé

Fe = 8000;              % Fréquence d'échantillonnage (8kHz)

% --- CORRECTION ICI ---
% Au lieu de fixer T=1, on calcule la durée réelle du signal
Lx = length(x);         % Nombre d'échantillons réel (57344)
t = (0:Lx-1)'/Fe;       % Vecteur temps adapté à la longueur de x
% ----------------------

% Paramètres d'analyse
N = 256;                % Taille de la fenêtre (ex: 32ms)
R = N/2;                % Recouvrement de 50% (Hop size = 128)
w = hanning(N);         % Fenêtre de Hanning

%% 2. Préparation (Zero-Padding)
nb_zeros_debut = R; 
x_pad = [zeros(nb_zeros_debut, 1); x; zeros(N, 1)];
L_pad = length(x_pad);
y_reconstruit_pad = zeros(L_pad, 1);

%% 3. Boucle de Traitement (Analyse -> Synthèse)
n_trames = floor((L_pad - N) / R) + 1;

for k = 1:n_trames
    % --- Analyse ---
    idx_debut = (k-1)*R + 1;
    idx_fin = idx_debut + N - 1;
    
    segment = x_pad(idx_debut:idx_fin);
    segment_fenetre = segment .* w; 
    
    % Traitement identité (pas de modification pour le moment)
    segment_traite = segment_fenetre; 
    
    % --- Synthèse (Addition-Recouvrement) ---
    y_reconstruit_pad(idx_debut:idx_fin) = y_reconstruit_pad(idx_debut:idx_fin) + segment_traite;
end

%% 4. Post-traitement et Vérification
% Calcul de la fenêtre de normalisation globale
W_norm = zeros(L_pad, 1);
for k = 1:n_trames
    idx_debut = (k-1)*R + 1;
    idx_fin = idx_debut + N - 1;
    W_norm(idx_debut:idx_fin) = W_norm(idx_debut:idx_fin) + w;
end

% Astuce : pour éviter la division par zéro dans les zones de padding pur 
% (où W_norm peut être 0), on remplace les 0 par des 1 juste pour la division.
% De toute façon, ces zones seront coupées à l'étape suivante.
W_norm(W_norm < 1e-10) = 1; 

% Normalisation
y_reconstruit_pad = y_reconstruit_pad ./ W_norm;

% Suppression du padding pour retrouver la taille originale exacte
y_final = y_reconstruit_pad(nb_zeros_debut+1 : nb_zeros_debut+Lx);

% --- Visualisation de l'erreur ---
erreur = x - y_final;

%% Affichage du signal original
figure;
plot(t, x, 'b');
title('Signal original');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;

%% Affichage du signal reconstruit
figure;
plot(t, y_final, 'r');
title('Signal reconstruit (après addition-recouvrement)');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;

%% Affichage du signal d’erreur
figure;
plot(t, erreur);
title('Erreur de reconstruction (x - y_{final})');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;


disp(['Erreur maximale de reconstruction : ' num2str(max(abs(erreur)))]);