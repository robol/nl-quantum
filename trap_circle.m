function [Y, maxnorm] = trap_circle(f, c, r)

    ntrap = 2;
    tol = 1e-4;
    maxnorm = 0;
    
    z = c + r * exp(2i * pi * (0 : ntrap - 1 ) ./ ntrap);

    h = 1 / ntrap;
    fz = f(z(1)); maxnorm = max(maxnorm, norm(fz, 'fro'));
    Y = fz * r * h;
    for j = 2 : ntrap
        fz = f(z(j)); maxnorm = max(maxnorm, norm(fz, 'fro'));
        Y = Y + fz * r * h * exp(2i * pi * (j-1) / ntrap);
    end    

    while ntrap < 4096
        ntrap = 2 * ntrap;
        z = c + r * exp(2i * pi * (0 : ntrap - 1 ) ./ ntrap);
    
        h = 1 / ntrap;
        Yold = Y;
        Y = Y / 2;
        for j = 2 : 2 : ntrap
            fz = f(z(j)); maxnorm = max(maxnorm, norm(fz, 'fro'));                   
            Y = Y + fz * r * h * exp(2i * pi * (j-1) / ntrap);
        end

        err_est = norm(Y - Yold, 'fro') / maxnorm;
        fprintf('TRAP: %d points, err est = %e\n', ntrap, err_est);

        if err_est < tol
            break;
        end
    end
end