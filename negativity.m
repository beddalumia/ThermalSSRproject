%% Compute negativity for a two-orbital RDM in usual EDIpack ordering
%
%   >> [E,F,N,M] = negativity(RDMij)
%
%      - E: Logarithmic negativity (standard "bosonic" partial transpose)
%      - F: Logarithmic negativity (twisted fermionic  partial transpose)
%      - N: "Original negativity"  (standard "bosonic" partial transpose)
%      - M: "Original negativity"  (twisted fermionic  partial transpose)
%
%  © Gabriele Bellomia and Frederic Bippus, 2025
function [E,F,N,M] = negativity(RDMij)
    if sum(sum(imag(RDMij)))>1e-8
    disp(imag(RDMij))
    error
    end
    % Preprocess the RDM
    [FiRDM, TiRDM,tindex] = partial_transpose(RDMij);
    % Get the eigenvalues
    p = eig(TiRDM);
    % Test for negativity
    N_ = -sum(p(p<0));
    N  = (sum(svd(TiRDM))-1)/2;
    if(abs(N-N_)>1e-8);
        warning('Negativity definitions mismatch: %g vs %g', N, N_);
    end
    p = eig(FiRDM);
    % Test for twisted vs untwisted negativity
    M_ = -sum(p(p<0)) + (trace(FiRDM) - 1)/2;
    M  = (sum(svd(FiRDM))-1)/2;
        if(abs(M-M_)>1e-8);
        warning('Negativity definitions (Fermionic) mismatch: %g vs %g', M, M_);
    end
    % Measure the entanglement
    E = log2(2*N+1);
    F = log2(2*M_+1);
 end