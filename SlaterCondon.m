function smatrix = SlaterCondon(nmodes)
    % SLATERCONDON : Implementation of Slater-Condon rules for fermions
    %                It pre-computes all the ❬istate| cdg_is c_js |jstate❭
    %                matrix elements, storing them in a 4D array.
    %                This can be used to build many-body representations of
    %                one-body operators. NB: it assumes Nup,Ndw conservation.
    %
    %  >> smatrix = SlaterCondon(nmodes :: number of single-fermion modes)
    %     smatrix :: 4D array [2*nmodes,2*nmodes,4^nmodes,4^nmodes]
    %     
    % ! This can be made way faster by implementing Slater-Condon
    %   rules in an efficient machine-tuned way, as discussed in
    %   https://arxiv.org/abs/1311.6244 (hal-01539072)
    %
    N = 4^nmodes;
    smatrix = zeros(2*nmodes,2*nmodes,N,N);
    for istate = 0:1:N-1 % nmode-orbital states
        for jstate = 0:1:N-1 % nmode-orbital states
            % ∑_s ∑_ij ❬istate| cdg_is c_js |jstate❭
            for ispin=1:2
                for imode=1:nmodes
                    for jmode=1:nmodes
                        % Apply cdg_is to ❬istate|
                        ibra = build_ket(istate,nmodes);
                        if ibra(imode+(ispin-1)*nmodes)==0
                            continue
                        end
                        ibra(imode+(ispin-1)*nmodes) = ibra(imode+(ispin-1)*nmodes)-1;
                        isign = (-1)^(sum(ibra(1:imode+(ispin-1)*nmodes)));
                        % Apply c_js to |jstate❭
                        jket = build_ket(jstate,nmodes);
                        if jket(jmode+(ispin-1)*nmodes)==0
                            continue
                        end
                        jket(jmode+(ispin-1)*nmodes) = jket(jmode+(ispin-1)*nmodes)-1;
                        jsign = (-1)^(sum(jket(1:jmode+(ispin-1)*nmodes)));
                        % Overlap ❬istate| cdg_is c_js |jstate❭
                        if isequal(ibra,jket)
                            smatrix(imode+(ispin-1)*nmodes,jmode+(ispin-1)*nmodes,istate+1,jstate+1) = isign * jsign;
                        end
                    end
                end
            end
        end
    end
end