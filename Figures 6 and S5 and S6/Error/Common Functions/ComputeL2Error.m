function L2_norm = ComputeL2Error(y, t, r)
    x = y - r;
    if length(x) ~= length(t)
        error('Signal and time vectors must have the same length');
    end
    dt = diff(t);
    L2_norm = sqrt(sum((x(1:end-1).^2) .* dt));
end
