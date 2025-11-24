%% Compute negativity for a two-orbital RDM in usual EDIPACK ordering
%
%   >> [E,F,N,M] = negativity(RDMij)
%
%      - E: Logarithmic negativity (standard "bosonic" partial transpose)
%      - F: Logarithmic negativity (twisted fermionic  partial transpose)
%      - N: "Original negativity"  (standard "bosonic" partial transpose)
%      - M: "Original negativity"  (twisted fermionic  partial transpose)
%
%  © Gabriele, Bellomia and Frederic Bippus, 2025
function [E,F,N,M] = negativity(RDMij)
    if sum(sum(imag(RDMij)))>1e-8
    disp(imag(RDMij))
    error
    end
    % Preprocess the RDM
    [FiRDM, TiRDM,tindex] = partial_transpose(RDMij);
    [whatever , T] = partial_transpose_frederic(RDMij);
    % Get the eigenvalues
    p = eig(TiRDM);
    % Test for negativity
    N = -sum(p(p<0));
    N_= (trace(sqrtm(TiRDM'*TiRDM))-1)/2;
    if(abs(N-N_)>1e-8);
        warning('Negativity definitions mismatch: %g vs %g', N, N_);
    end
    p = eig(FiRDM);
    % Test for twisted vs untwisted negativity
    M = -sum(p(p<0)) + (trace(FiRDM) - 1)/2;
    %M_ = (trace(sqrtm(FiRDM'*FiRDM))-1)/2;
    M_ = (sum(svd(FiRDM))-1)/2;
        if(abs(M-M_)>1e-8);
        warning('Negativity definitions (Fermionic) mismatch: %g vs %g', M, M_);
    end
    % Measure the entanglement
    E = log2(2*N+1);
    F = log2(2*M_+1);
 end

 %% Partial transpose on the "left" single site
 function [Ti_RDM_ij, Fi_RDM_ij] = partial_transpose_frederic(RDM_ij)
    % HARDCODED MAGIC NUMBERS (we have a simple dimer)
    Nlat = 2;
    Norb = 1;
    % Let's assert the RDM has the right dimensions
    Nrdm = size(RDM_ij,1);
    Nlso = Nlat*Norb*2;
    Nrdm = 2^Nlso;
    assert(all(size(RDM_ij)==[Nrdm,Nrdm]))
    % Init the partial-transposed RDM to size(RDM)
    Ti_RDM_ij = zeros(Nrdm,Nrdm);
    %% Rotate the RDM_ij to a comfortable basis: |i_up i_dw j_up j_dw>
    Ti_RDM_ij(1,1) = RDM_ij(1,1);
    Ti_RDM_ij(2,2) = RDM_ij(9,9);
    Ti_RDM_ij(3,3) = RDM_ij(3,3);
    Ti_RDM_ij(4,4) = RDM_ij(5,5);
    Ti_RDM_ij(5,5) = RDM_ij(2,2);
    Ti_RDM_ij(6,6) = RDM_ij(11,11);
    Ti_RDM_ij(7,7) = RDM_ij(13,13);
    Ti_RDM_ij(8,8) = RDM_ij(10,10);
    Ti_RDM_ij(9,9) = RDM_ij(7,7);
    Ti_RDM_ij(10,10) = RDM_ij(4,4);
    Ti_RDM_ij(11,11) = RDM_ij(6,6);
    Ti_RDM_ij(12,12) = RDM_ij(15,15);
    Ti_RDM_ij(13,13) = RDM_ij(12,12);
    Ti_RDM_ij(14,14) = RDM_ij(14,14);
    Ti_RDM_ij(15,15) = RDM_ij(8,8);
    Ti_RDM_ij(16,16) = RDM_ij(16,16);
    %> Off-diagonals with kept sign
    Ti_RDM_ij(5,3) = RDM_ij(2,3);
    Ti_RDM_ij(3,5) = RDM_ij(3,2);
    Ti_RDM_ij(4,2) = RDM_ij(5,9);
    Ti_RDM_ij(2,4) = RDM_ij(9,5);
    Ti_RDM_ij(11,9) = RDM_ij(6,7);
    Ti_RDM_ij(9,11) = RDM_ij(7,6);
    Ti_RDM_ij(11,8) = RDM_ij(6,10);
    Ti_RDM_ij(8,11) = RDM_ij(10,6) ;
    Ti_RDM_ij(11,6) = RDM_ij(6,11);
    Ti_RDM_ij(6,11) = RDM_ij(11,6);
    Ti_RDM_ij(9,8) = RDM_ij(7,10);
    Ti_RDM_ij(8,9) = RDM_ij(10,7);
    Ti_RDM_ij(9,6) = RDM_ij(7,11);
    Ti_RDM_ij(6,9) = RDM_ij(11,7);
    Ti_RDM_ij(15,13) = RDM_ij(8,12);
    Ti_RDM_ij(13,15) = RDM_ij(12,8);
    Ti_RDM_ij(8,6) = RDM_ij(10,11);
    Ti_RDM_ij(6,8) = RDM_ij(11,10);
    Ti_RDM_ij(14,12) = RDM_ij(14,15);
    Ti_RDM_ij(12,14) = RDM_ij(15,14);
    %> FERMIONIC SIGNS
    Ti_RDM_ij(11,9) = -Ti_RDM_ij(11,9);
    Ti_RDM_ij(9,11) = -Ti_RDM_ij(9,11);
    Ti_RDM_ij(8,9) = -Ti_RDM_ij(8,9);
    Ti_RDM_ij(9,8) = -Ti_RDM_ij(9,8);
    Ti_RDM_ij(6,9) = -Ti_RDM_ij(6,9);
    Ti_RDM_ij(9,6) = -Ti_RDM_ij(9,6);
    Ti_RDM_ij(13,15) = -Ti_RDM_ij(13,15);
    Ti_RDM_ij(15,13) = -Ti_RDM_ij(15,13);
    Ti_RDM_ij(12,14) = -Ti_RDM_ij(12,14);
    Ti_RDM_ij(14,12) = -Ti_RDM_ij(14,12);
    %% "Fermionic" partial transpose (© Frederic Bippus)
    Fi_RDM_ij = [Ti_RDM_ij(1,1) 0 0 0 0 0 1j*Ti_RDM_ij(2,4) 0 0 1j*Ti_RDM_ij(3,5) 0 0 0 0 0 Ti_RDM_ij(11,6);...
        0 Ti_RDM_ij(2,2) 0 0 0 0 0 0 0 0 0 0 1j*Ti_RDM_ij(8,6) 0 0 0;...
        0 0 Ti_RDM_ij(3,3) 0 0 0 0 0 0 0 0 1j*Ti_RDM_ij(9,6) 0 0 0 0;...
        0 0 0 Ti_RDM_ij(4,4) 0 0 0 0 0 0 0 0 0 0 1j*Ti_RDM_ij(9,11) 0;...
        0 0 0 0 Ti_RDM_ij(5,5) 0 0 0 0 0 0 0 0 1j*Ti_RDM_ij(8,11) 0 0;...
        0 0 0 0 0 Ti_RDM_ij(6,6) 0 0 0 0 0 0 0 0 0 0;...
        1j*Ti_RDM_ij(2,4) 0 0 0 0 0 Ti_RDM_ij(7,7) 0 0 Ti_RDM_ij(8,9) 0 0 0 0 0 1j*Ti_RDM_ij(12,14);...
        0 0 0 0 0 0 0 Ti_RDM_ij(8,8) 0 0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0 0 Ti_RDM_ij(9,9) 0 0 0 0 0 0 0;...
        1j*Ti_RDM_ij(3,5) 0 0 0 0 0 Ti_RDM_ij(8,9) 0 0 Ti_RDM_ij(10,10) 0 0 0 0 0 1j*Ti_RDM_ij(13,15);...
        0 0 0 0 0 0 0 0 0 0 Ti_RDM_ij(11,11) 0 0 0 0 0;...
        0 0 1j*Ti_RDM_ij(9,6) 0 0 0 0 0 0 0 0 Ti_RDM_ij(12,12) 0 0 0 0;...
        0 1j*Ti_RDM_ij(8,6) 0 0 0 0 0 0 0 0 0 0 Ti_RDM_ij(13,13) 0 0 0;...
        0 0 0 0 1j*Ti_RDM_ij(8,11) 0 0 0 0 0 0 0 0 Ti_RDM_ij(14,14) 0 0;...
        0 0 0 1j*Ti_RDM_ij(9,11) 0 0 0 0 0 0 0 0 0 0 Ti_RDM_ij(15,15) 0;...
        Ti_RDM_ij(11,6) 0 0 0 0 0 1j*Ti_RDM_ij(12,14) 0 0 1j*Ti_RDM_ij(13,15) 0 0 0 0 0 Ti_RDM_ij(16,16)];
    %% "Bosonic" partial transpose (© Frederic Bippus)
    Ti_RDM_ij = [Ti_RDM_ij(1,1) 0 0 0 0 0 Ti_RDM_ij(2,4) 0 0 Ti_RDM_ij(3,5) 0 0 0 0 0 Ti_RDM_ij(11,6);...
        0 Ti_RDM_ij(2,2) 0 0 0 0 0 0 0 0 0 0 Ti_RDM_ij(8,6) 0 0 0;...
        0 0 Ti_RDM_ij(3,3) 0 0 0 0 0 0 0 0 Ti_RDM_ij(9,6) 0 0 0 0;...
        0 0 0 Ti_RDM_ij(4,4) 0 0 0 0 0 0 0 0 0 0 Ti_RDM_ij(9,11) 0;...
        0 0 0 0 Ti_RDM_ij(5,5) 0 0 0 0 0 0 0 0 Ti_RDM_ij(8,11) 0 0;...
        0 0 0 0 0 Ti_RDM_ij(6,6) 0 0 0 0 0 0 0 0 0 0;...
        Ti_RDM_ij(2,4) 0 0 0 0 0 Ti_RDM_ij(7,7) 0 0 Ti_RDM_ij(8,9) 0 0 0 0 0 Ti_RDM_ij(12,14);...
        0 0 0 0 0 0 0 Ti_RDM_ij(8,8) 0 0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0 0 Ti_RDM_ij(9,9) 0 0 0 0 0 0 0;...
        Ti_RDM_ij(3,5) 0 0 0 0 0 Ti_RDM_ij(8,9) 0 0 Ti_RDM_ij(10,10) 0 0 0 0 0 Ti_RDM_ij(13,15);...
        0 0 0 0 0 0 0 0 0 0 Ti_RDM_ij(11,11) 0 0 0 0 0;...
        0 0 Ti_RDM_ij(9,6) 0 0 0 0 0 0 0 0 Ti_RDM_ij(12,12) 0 0 0 0;...
        0 Ti_RDM_ij(8,6) 0 0 0 0 0 0 0 0 0 0 Ti_RDM_ij(13,13) 0 0 0;...
        0 0 0 0 Ti_RDM_ij(8,11) 0 0 0 0 0 0 0 0 Ti_RDM_ij(14,14) 0 0;...
        0 0 0 Ti_RDM_ij(9,11) 0 0 0 0 0 0 0 0 0 0 Ti_RDM_ij(15,15) 0;...
        Ti_RDM_ij(11,6) 0 0 0 0 0 Ti_RDM_ij(12,14) 0 0 Ti_RDM_ij(13,15) 0 0 0 0 0 Ti_RDM_ij(16,16)];
 end



 
