function x_reconstruit = overlap_add(frames, Nf, overlap)
    step = round(Nf * (1 - overlap));
    [~, nbFrames] = size(frames);
    L = Nf + (nbFrames - 1)*step;
    x_reconstruit = zeros(1, L);
    win = hamming(Nf);

    for k = 1:nbFrames
        start = (k-1)*step + 1;
        x_reconstruit(start:start+Nf-1) = x_reconstruit(start:start+Nf-1) + (frames(:,k).*win)';
    end

    norm_factor = zeros(1, L);
    for k = 1:nbFrames
        start = (k-1)*step + 1;
        norm_factor(start:start+Nf-1) = norm_factor(start:start+Nf-1) + (win').^2;
    end
    x_reconstruit = x_reconstruit ./ (norm_factor + eps);
end
