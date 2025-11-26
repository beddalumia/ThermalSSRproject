function RDM = RDM0(OBDM,N,SCOT)
%% N-mode reduced density matrices from N-mode correlation matrices
%
%  >> RDM = RDM0(OBDM,N,SCOT)
%
%     - OBDM is a valid N-site correlation matrix [2N,2N]
%     - N is the number of given fermionic modes (orbitals, sites)
%     - SCOT is a precomputed Slater-Condon overlap tensor [2N,2N,4^N,4^N]
%
% See also slater_condon
%
% © Gabriele Bellomia, 2025

        % Input assertions
        assert(isequal(size(OBDM),[2*N,2*N]))
        assert(isequal(size(SCOT),[2*N,2*N,4^N,4^N]))


        % Entanglement Hamiltonian (Peschel-Chung-Eisler theorem)
        [V,G] = eig(OBDM,"vector");
        D = log((1-G)./G);
        H = V*diag(D)*V^(-1); H = H';

        % Many-Body representation of the Entanglement Hamiltonian (EH)
        mbH = zeros(4^N,4^N);
        % ∑_s ∑_ij <istate| H_ij cdg_is c_js |jstate>
        for ispin=1:2
            for ilat=1:N
                for jlat=1:N
                    % <istate| H_ij cdg_is c_js |jstate> = H_ij * <istate| cdg_is c_js |jstate>
                    mbH(:,:) = H(ilat+(ispin-1)*N,jlat+(ispin-1)*N) *...
                               squeeze(SCOT(ilat+(ispin-1)*N,jlat+(ispin-1)*N,:,:)) + mbH(:,:);
                end
            end
        end

        % Reduced density matrix from many-body EH
        RDM = expm(-mbH); RDM = RDM ./ trace(RDM); 

        % Assert Slater-Condon rules (never too sure...)
        test_1bdm = zeros(2*N,2*N);
        for n = 1:2*N
            for m = 1:2*N
                test_1bdm(n,m) = trace(RDM * squeeze(SCOT(n,m,:,:)));
            end
        end
        assert(norm(test_1bdm-OBDM)<1d-10);

end