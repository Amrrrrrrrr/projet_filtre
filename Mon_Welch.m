function [S_welch, f_welch] = Mon_Welch(x, NFFT, Fe)

    x = x(:).';                  
    L = length(x);              

    %  Paramètres
    R = floor(NFFT / 2);
    window = hamming(NFFT).';
    U = mean(window.^2);          

    % Indices de début des segments 
    start_indices = 1 : (NFFT - R) : (L - NFFT + 1);
    M = length(start_indices);    

    S_sum = zeros(1, NFFT);

    for m = 1:M
        idx = start_indices(m) : start_indices(m) + NFFT - 1; 
        x_seg = x(idx) .* window;                              
        X = fft(x_seg, NFFT);                                  
        Pxx = (abs(X).^2) / (Fe * NFFT * U);                  
        S_sum = S_sum + Pxx;                                  
    end

    S_welch = S_sum / M;

    f_welch = (0:NFFT-1) * (Fe / NFFT);

end
