function [L, chol_success] = robust_cholesky(KF_V)
    % ROBUST_CHOLESKY Compute Cholesky decomposition with fallback strategies
    %
    % RENAMED (2026-07): formerly f_robust_chol.m in Miscellaneous/. All call
    % sites inside ReliabilityToolbox/ use the new names; the full mapping
    % is in README.md ("Suggested file names").
    %
    % Inputs:
    %   KF_V - Square matrix for Cholesky decomposition
    %
    % Outputs:
    %   L - Lower triangular Cholesky factor (empty if all methods fail)
    %   chol_success - Logical flag indicating success (true) or failure (false)
    %
    % Method hierarchy:
    %   1. Direct Cholesky decomposition
    %   2. Eigenvalue regularization with symmetrization
    %   3. Diagonal loading with adaptive regularization
    %   4. Return failure flag if all methods fail
    
    % Initialize output
    L = [];
    chol_success = false;
    
    % First attempt: Direct Cholesky
    try
        L = chol(KF_V, 'lower');
        chol_success = true;
        return;
    catch
        % Continue to fallback methods
    end
    
    % First fallback: Eigenvalue regularization
    try
        KF_V_sym = (KF_V + KF_V') / 2; % Ensure symmetry
        [V, D] = eig(KF_V_sym);
        eigenvals = diag(D);
        
        % More conservative threshold
        min_eigenval = 1e-10;
        if any(eigenvals <= 0)
            eigenvals = max(eigenvals, min_eigenval);
            KF_V_psd = V * diag(eigenvals) * V';
        else
            KF_V_psd = KF_V_sym;
        end
        
        L = chol(KF_V_psd, 'lower');
        chol_success = true;
        return;
        
    catch
        % Continue to next fallback
    end
    
    % Second fallback: Regularized diagonal loading
    try
        reg_param = 1e-8 * trace(KF_V) / size(KF_V, 1); % Adaptive regularization
        KF_V_reg = KF_V + reg_param * eye(size(KF_V));
        L = chol(KF_V_reg, 'lower');
        chol_success = true;
        
    catch
        % Final fallback: return failure
        chol_success = false;
    end
    
end