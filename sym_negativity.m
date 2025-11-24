%% Compute imbalance-resolved negativity for a two-orbital RDM in usual EDIPACK ordering
%
%   >> [N0,N1,N2] = sym_negativity(RDMij)
%
%  © Gabriele Bellomia, 2025

function [N0,N1,N2] = sym_negativity(RDMij)
    
% Build charge-imbalance blocks RDM_q (q: charge imbalance)

    % q = 0
    RDM_0 = zeros(4,4);
    RDM_0(1,1) = RDMij(4,4);
    RDM_0(2,3) = RDMij(7,10);
    RDM_0(3,2) = RDMij(10,7);
    RDM_0(4,4) = RDMij(13,13);
    
    % q = 1
    % > unimplemented (todo...) 
    
    % q = 2
    RDM_2 = zeros(4,4);
    RDM_2(1,1) = RDMij(1,1);
    RDM_2(2,3) = RDMij(6,11);
    RDM_2(3,2) = RDMij(11,6);
    RDM_2(4,4) = RDMij(16,16);
    
% Build partial transposes (no need to enforce twisted phases under P-SSR)
    
    % q = 0
    RDM_0(1,4) = -RDM_0(2,3); % Fermionic sign in the
    RDM_0(4,1) = -RDM_0(3,2); % partial transpose! :O
    RDM_0(3,2) = 0;
    RDM_0(2,3) = 0;
    
    % q = 1 [here we need to take care of twisted phases]
    % TODO: build this directly instead of using the messy code by Frederic

    
    % q = 2
    RDM_2(1,4) = RDM_2(2,3);
    RDM_2(4,1) = RDM_2(3,2);
    RDM_2(3,2) = 0;
    RDM_2(2,3) = 0;
    
    % Evaluate transposed spectra...
    p0 = eig(RDM_0);
    p2 = eig(RDM_2);
    % ...and their negativity :)
    N0 = -sum(p0(p0<0));
    N1 = NaN; % Unimplemented (TODO)
    N2 = -sum(p2(p2<0));
    
    % > Let's leverage the additivity of negativities!
    [FT, BT] = partial_transpose(RDMij);
    Ntot = 0.5*(trace(sqrtm(FT'*FT))-1);
    Ntot_bosonic = 0.5*(trace(sqrtm(BT'*BT))-1);
    if abs(Ntot_bosonic - Ntot) < 1e-12
       warning("You might want to investigate why bosonic and fermionic PPT coincide");
    end
    N1 = Ntot - N2 - N0;
    
end
 




 %% Full partial transpose on the "left" single site
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

 
