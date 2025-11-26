function y = movmean_simple(x, k)
    N = length(x);
    y = zeros(1, N);
    half = floor(k/2);
    
    for i = 1:N
        start_i = max(1, i - half);
        end_i = min(N, i + half);
        y(i) = mean(x(start_i:end_i));
    end
end
