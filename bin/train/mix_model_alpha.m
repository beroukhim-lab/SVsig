function [opt_alpha, f_bic] = mix_model_alpha(mfull, model1, model2, annot_tiles, num_param)
%MIX_MODEL_ALPHA Optimize mixture alpha values for complex weighted model.
%   MATLAB-native replacement for legacy mix_model_alpha.R workflow.
%   Returns:
%     opt_alpha - best alpha vector in [0,1]^num_param
%     f_bic     - best BIC objective value
%
%   This function is only called from complex+weights mode in
%   mix_model_param.m, so the simple model path is unaffected.

nume = sum(mfull(:));
if nume <= 0
    opt_alpha = 0.5 * ones(num_param, 1);
    f_bic = inf;
    return;
end

idx_nnz = (mfull > 0);
x = double(mfull(idx_nnz));
x = max(0, min(x, nume));

options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'off', ...
    'TolFun', 1e-6, ...
    'TolX', 1e-6, ...
    'DiffMinChange', 1e-6);

n_starts = 20;
f_bic = inf;
opt_alpha = 0.5 * ones(num_param, 1);

for s = 1:n_starts
    alpha0 = rand(num_param, 1);
    [alpha_try, bic_try] = fmincon(@complex_bic_fun, alpha0, [], [], [], [], ...
                                   zeros(num_param, 1), ones(num_param, 1), [], options);
    if bic_try < f_bic
        f_bic = bic_try;
        opt_alpha = alpha_try;
    end
end

    function bic = complex_bic_fun(alpha)
        mix_tmp = zeros(size(mfull));
        for cc = 1:num_param
            mix_tmp(annot_tiles(:,:,cc)) = mix_tmp(annot_tiles(:,:,cc)) + ...
                alpha(cc) * model1(annot_tiles(:,:,cc)) + ...
                (1 - alpha(cc)) * model2(annot_tiles(:,:,cc));
        end

        denom = sum(mix_tmp(:));
        if denom <= 0
            bic = inf;
            return;
        end
        mix_tmp = mix_tmp / denom;

        p = double(mix_tmp(idx_nnz));
        p = min(max(p, 1e-12), 1 - 1e-12);

        % Stable binomial log-likelihood via gammaln.
        ll = sum(gammaln(nume + 1) - gammaln(x + 1) - gammaln(nume - x + 1) + ...
                 x .* log(p) + (nume - x) .* log(1 - p));

        bic = -2 * ll + log(nume) * num_param;
    end
end
