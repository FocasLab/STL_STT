function [Lx, Ly, Ux, Uy] = evaluate_piecewise(t)
    % Load piecewise coefficients from CSV
    Cpiece = readmatrix('Robot.csv');

    % If the coefficients are stored as a column vector, reshape
    if size(Cpiece, 2) == 1 && mod(length(Cpiece), 4) == 0
        Cpiece = reshape(Cpiece, [], 4)';  % Reshape to 4 x N
    end

    % Piecewise breakpoints and polynomial degree
    breaks = [0 4 7 10 13 18];
    num_segments = length(breaks) - 1;
    degree = 5;
    coefs_per_seg = degree + 1;

    % Ensure t is a column vector
    t = t(:);

    % Initialize outputs
    Lx = zeros(size(t));
    Ly = zeros(size(t));
    Ux = zeros(size(t));
    Uy = zeros(size(t));

    for k = 1:length(t)
        tk = t(k);

        % Clamp t to valid range
        if tk < breaks(1)
            tk = breaks(1);
        elseif tk >= breaks(end)
            tk = breaks(end) - 1e-8;
        end

        % Find segment
        seg_idx = find(breaks <= tk, 1, 'last');
        if seg_idx >= num_segments + 1
            seg_idx = num_segments;
        end

        t_local = tk - breaks(seg_idx);

        % Evaluate each curve
        for curve = 1:4
            coefs = Cpiece(curve, (seg_idx-1)*coefs_per_seg + (1:coefs_per_seg));
            val = polyval(flip(coefs), t_local);
            switch curve
                case 1
                    Lx(k) = val;
                case 2
                    Ux(k) = val;
                case 3
                    Ly(k) = val;
                case 4
                    Uy(k) = val;
            end
        end
    end
end
