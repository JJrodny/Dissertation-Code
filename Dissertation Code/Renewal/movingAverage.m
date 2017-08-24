function y = movingAverage(x, w)
    k = ones(1, w) / w
    y = conv(x, k, 'same');
    A = x;
    B = w;
    L = filter(ones(B+1,1)/B+1,1,[A(:) zeros(floor(B/2),1)]);
    out = L(B/2+1:end);
end