function [x, n, delta, error_flag] = eqsolver(stoichiometry, K, mass_conservation, total_masses, x, iter, w, tolerance)

    if min(total_masses) < 0
        error('One or more concentrations are negative, please check total_masses vector')
    end
    if min(total_masses) == 0
        warning('One or more concentrations are equal to zero. These have been replaced by an arbitrary small number.')
        total_masses(find(total_masses <= 0)) = 2.2e-16;
    end

    if (nargin<5) || isempty(x)
      x = ones(1,size(stoichiometry,2));      
    end

    if (nargin<6) || isempty(iter)
      iter = 20;      
    end

    if (nargin<7) || isempty(w)
      w = 0;      
    end

    if (nargin<8) || isempty(tolerance)
      tolerance = 1e2;      
    end

    for n = 1:iter
        [x_, delta_] = solve(stoichiometry, K, mass_conservation, total_masses, x);
        delta(n) = delta_;
        x = (x.*w+x_)./(w+1);
        if delta(n) <= tolerance * norm(eps(x));
            break
        end
    end
    if delta(n) > tolerance && n == iter
        warning('Tolerance not reached in the last step of solver, increase iterations or manually check if the solution is appropriate')
    end

    if ~isfinite(x) | isnan(x)
        error_flag = true;
        warning('The algorithm did not converge to a solution, increase w or change starting point')
    else
        error_flag = false;
    end
    x = exp(x);
end

function [x, delta_] = solve(stoichiometry, K, mass_conservation, total_masses, x)
    Cx = mass_conservation.*exp(x);
    W = Cx./sum(Cx, 2);
    W_ = prod((W./mass_conservation).^W, 2, "omitnan").*total_masses;

    X = [stoichiometry; W];
    Y = log([K; W_]);
    delta_ = norm(X*x'-Y);

    %if size(stoichiometry,2) == length(K) + length(total_masses)
    %    x = (X\Y)';
    %else
        x = (pinv(X)*Y)';
   %end
end
