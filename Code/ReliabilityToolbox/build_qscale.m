function Q = build_qscale(M_, T)
%BUILD_QSCALE  Rebuild heteroskedastic_shocks.Qscale from CURRENT M_.params.
%
%   Q = BUILD_QSCALE(M_, T)   n_shocks x T matrix of shock VARIANCE scales.
%
%   Starting point (CHANGED 2026-09-26): the STORED M_.heteroskedastic_shocks.Qscale
%   when the model carries one (NaN -> 1, cut or padded with ones to T columns),
%   else ones. Until now the function started from ones, which wiped every
%   scale that is not a kappa - the missing-data devices of the SPF model
%   (e_spfobs pre-1992) and of the Bauer-Rudebusch ESE model (ZLB yields, no
%   kappa parameters at all) - whenever a posterior-predictive routine called
%   it. On the baseline HLW model the two starting points coincide (its stored
%   scales are ones outside the kappa rows and periods), so nothing changes
%   there.
%
%   Then the model-specific part, hardcoded: the HLW kappas, when the model has
%   them, overwrite their rows and periods with the CURRENT M_.params values
%   (Dynare freezes the kappas at their calibrated values in the stored Qscale,
%   see model_assessment.m; this is what makes the rebuilt scales follow a draw).

if isfield(M_, 'heteroskedastic_shocks') && isfield(M_.heteroskedastic_shocks, 'Qscale') && ~isempty(M_.heteroskedastic_shocks.Qscale)
    Q = M_.heteroskedastic_shocks.Qscale;
    Q(isnan(Q)) = 1;
    if size(Q, 2) >= T
        Q = Q(:, 1:T);
    else
        Q = [Q, ones(size(Q, 1), T - size(Q, 2))];
    end
else
    % homoskedastic (uc3 and any model without a heteroskedastic block)
    Q = ones(M_.exo_nbr, T);
end

if ismember('kappa2020', M_.param_names) && isfield(M_, 'heteroskedastic_shocks')
    % Periods/blocks hard-coded to match HLW_2023_US_TD_ME_small.mod.
    k20 = M_.params(strcmp(M_.param_names,'kappa2020'));
    k21 = M_.params(strcmp(M_.param_names,'kappa2021'));
    k22 = M_.params(strcmp(M_.param_names,'kappa2022'));
    rows = find(ismember(M_.exo_names, {'e_pi','e_ygap'}));
    Q(rows, 242:244) = k20^2;
    Q(rows, 245:248) = k21^2;
    Q(rows, 249:252) = k22^2;
end

end
