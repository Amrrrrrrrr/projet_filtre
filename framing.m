function frames = framing(x, Nf, overlap)
    step = round(Nf * (1 - overlap));
    L = length(x);
    nbFrames = floor((L - Nf) / step) + 1;
    frames = zeros(Nf, nbFrames);
    for k = 1:nbFrames
        start = (k-1)*step + 1;
        frames(:,k) = x(start:start+Nf-1);
    end
end
