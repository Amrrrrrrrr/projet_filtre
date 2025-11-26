clear; close all; clc;

fs = 8000;               
t = 0:1/fs:1;             
data = load("fcno04fz.mat");  
nom_champ = fieldnames(data);
x = data.(nom_champ{1});

N = 512;                  % taille de trame (en échantillons)
R = N/2;                  % pas de trame = 50% de recouvrement
w = hanning(N, 'symmetric'); 

%% Découpage du signal en trames
frames = buffer(x, N, N-R, 'nodelay'); % chaque colonne = 1 trame
nb_trames = size(frames, 2);

%% Application de la fenêtre sur chaque trame
for m = 1:nb_trames 
    frames(:,m) = frames(:,m) .* w;
end

%% Reconstruction du signal par addition-recouvrement
x_rec = zeros(length(x), 1);

for m = 1:nb_trames
    start_idx = (m-1)*R + 1;
    end_idx = start_idx + N - 1;

    % Éviter de dépasser la taille du signal
    if end_idx > length(x)
        break;
    end
    
    x_rec(start_idx:end_idx) = x_rec(start_idx:end_idx) + frames(:,m);
end

%% Normalisation pour corriger le fenêtrage 
% Calcul du facteur de recouvrement total pour chaque échantillon
win_sum = zeros(length(x),1);
for m = 1:nb_trames
    start_idx = (m-1)*R + 1;
    end_idx = start_idx + N - 1;
    if end_idx > length(x)
        break;
    end
    win_sum(start_idx:end_idx) = win_sum(start_idx:end_idx) + w;
end
x_rec = x_rec ./ (win_sum + eps); % éviter division par zéro

%% Vérification de la reconstruction
erreur = norm(x - x_rec) / norm(x);
fprintf('Erreur relative de reconstruction : %.2e\n', erreur);

%% Affichage des signaux
figure;
plot(x, 'b'); hold on;
plot(x_rec, 'r--');
legend('Signal original', 'Signal reconstruit');
title('Vérification de la reconstruction (addition-recouvrement)');
xlabel('Échantillons'); ylabel('Amplitude');
grid on;

%% Affichage du signal d’erreur
figure;
plot(x - x_rec);
title('Différence entre signal original et reconstruit');
xlabel('Échantillons'); ylabel('Erreur');
grid on;

