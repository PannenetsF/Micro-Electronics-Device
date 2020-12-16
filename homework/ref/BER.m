function ret = BER(x)
    ret = x ./ (e^x - 1 + 1e-10);
end