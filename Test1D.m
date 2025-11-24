% Cross-check a numerical build of the 2-orbital RDM against the analyical
% results in Quantum Sci. Technol. 9 (2024) 015005 (Lexin Ding et al.)
% This script reproduces figure 2 in the paper, as well as it replicates 
% the results with our notion of symmetry-resolved negativity :)

L = 200; % Nsites
k = (1-mod(L,2)-round((L-1)/2):1:round((L-1)/2)) * 2*pi/L;
E = -cos(k); [~,indices] = sort(E); sorted_k = k(indices);

% Filling of the band η = N/(2*L)
global eta
fillings = [logspace(-4,-1,100),0.1:0.01:0.5]; %[1e-4,1e-3,1e-2:1e-2:0.09,0.1:0.01:0.5];
E_nSSR = zeros(4,length(fillings)); 
E_pSSR = zeros(4,length(fillings));
E_PPT_0 = zeros(4,length(fillings));
E_PPT = zeros(4,length(fillings));
for d=1:length(fillings)
eta = fillings(d) 
N = round(2*L*eta);
occupied = zeros(round(N/2),2);
for n = 1:round(N/2)
    for s = 1:2
        occupied(n,s) = sorted_k(n);
    end
end

%% Numerical 1-body density matrix!
OBDM = zeros(2*L,2*L); % [1,2,...,Nsites]_up[1,2,...,Nsites]_dw
for i = 1:L
    OBDM(i,i) = eta;
    OBDM(i+L,i+L) = eta;
    for j = i+1:L
        for n = 1:length(occupied(:,1))
            % Spin-up
            %OBDM(i,j) = OBDM(i,j) + exp(1i*(i-j)*occupied(n,1)) / L;
        end
        % Infinite chain limit:
        OBDM(i,j) = sin(pi*(i-j)*eta)/(pi*(i-j));
        for n = 1:length(occupied(:,2))
            % Spin-dw
            %OBDM(i+L,j+L) = OBDM(i+L,j+L) + exp(1i*(i-j)*occupied(n,2)) / L;
        end
        % Infinite chain limit:
        OBDM(i+L,j+L) = sin(pi*(i-j)*eta)/(pi*(i-j));
        if mod(length(occupied),2)==1
            analytic=sin(0.5*pi*(i-j)*N/L)/sin(pi*(i-j)/L) / L;
            %assert(abs(OBDM(i,j)-analytic)<1d-1);
        end
        if not(i==j)
            OBDM(j,i) = conj(OBDM(i,j));
            OBDM(j+L,i+L) = conj(OBDM(i+L,j+L));
        end
    end
end

% pre-compute Slater-Condon rules for 2 orbitals
sc_matrix = SlaterCondon(2); % 4D array [2*nmodes,2*nmodes,4^nmodes,4^nmodes]

%% Build two-orbital RDMs
OBDM_ij = zeros(4,4);
counter = 0;
for i = 1:1%L
    for j = [2,3,11,101]

        % Restrict the one-body density matrix to two orbitals
        % -> spin up
        OBDM_ij(1,1) = OBDM(i,i); 
        OBDM_ij(1,2) = OBDM(i,j); 
        OBDM_ij(2,1) = OBDM(j,i);
        OBDM_ij(2,2) = OBDM(j,j);
        % -> spin down
        OBDM_ij(3,3) = OBDM(i+L,i+L); 
        OBDM_ij(3,4) = OBDM(i+L,j+L); 
        OBDM_ij(4,3) = OBDM(j+L,i+L);
        OBDM_ij(4,4) = OBDM(j+L,j+L);
        % To check how close to 1 the eigvals are
        %  disp(eig(OBDM_ij))
        %  disp(sum(eig(OBDM_ij)))
        %  disp ********************

        % Entanglement Hamiltonian (Peshel-Cheong theorem)
        % H = (logm((1-OBDM_ij)/OBDM_ij))'; % Native matrix log not really working
        [V,G] = eig(OBDM_ij,"vector");
        D = log((1-G)./G);
        H = V*diag(D)*V^(-1); H = H';

        % Many-Body representation of the Entanglement Hamiltonian
        mbH = zeros(16,16);
        % ∑_s ∑_ij <istate| H_ij cdg_is c_js |jstate>
        for ispin=1:2
            for ilat=1:2
                for jlat=1:2
                    % <istate| H_ij cdg_is c_js |jstate> = H_ij * <istate| cdg_is c_js |jstate>
                    mbH(:,:) = H(ilat+(ispin-1)*2,jlat+(ispin-1)*2) * squeeze(sc_matrix(ilat+(ispin-1)*2,jlat+(ispin-1)*2,:,:)) + mbH(:,:);
                end
            end
        end

        % Reduced density matrix from maby-body EH
        RDM = expm(-mbH); RDM = RDM ./ trace(RDM); 

        % Assert Slater-Condon rules
        test_1bdm = zeros(4,4);

        for n = 1:4
            for m = 1:4
                test_1bdm(n,m) = trace(RDM * squeeze(sc_matrix(n,m,:,:)));
            end
        end

        assert(norm(test_1bdm-OBDM_ij)<1d-10);
        
        % Apply SSR
        counter = counter + 1;
        [N0,~,N2] = sym_negativity(RDM);
        E_PPT_0(counter,d) = log2(2*N0+1);
        E_PPT(counter,d) = log2(2*(N0+N2)+1);
        [E_pSSR(counter,d),E_nSSR(counter,d)] = build_SSR(RDM);


    end
end

end

plot(fillings,E_pSSR(1,:)*log(2),'r-'); hold on
plot(fillings,E_nSSR(1,:)*log(2),'b-');
plot(fillings,E_PPT_0(1,:)*log(2),'c-');
plot(fillings,E_PPT(1,:)*log(2),'m-');

plot(fillings,E_pSSR(2,:)*log(2),'r--');
plot(fillings,E_nSSR(2,:)*log(2),'b--');
plot(fillings,E_PPT_0(2,:)*log(2),'c--');
plot(fillings,E_PPT(2,:)*log(2),'m--');

plot(fillings,E_pSSR(3,:)*log(2),'r-.');
plot(fillings,E_nSSR(3,:)*log(2),'b-.'); 
plot(fillings,E_PPT_0(3,:)*log(2),'c-.');
plot(fillings,E_PPT(3,:)*log(2),'m-.');

plot(fillings,E_pSSR(4,:)*log(2),'r:');
plot(fillings,E_nSSR(4,:)*log(2),'b:'); 
plot(fillings,E_PPT_0(4,:)*log(2),'c:');
plot(fillings,E_PPT(4,:)*log(2),'m:');

set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlim([1e-4,0.5]);
xlabel("$\eta$",'Interpreter','latex')
ylabel("$E$",'Interpreter','latex')

print_basis


%% contains

    % FROM ED_SETUP:
    % |imp_up>|bath_up> * |imp_dw>|bath_dw>        <- 2*Nlat*Norb bits
    % |imp_sigma> = | (1…Norb)_1...(1…Norb)_Nlat > <--- Nlat*Norb bits
    % lso indices are: io = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
    function [vec,ket] = build_ket(state)
      % BUILD_KET : Puts together a bit representation of a Slater determinant
      %             (and, if asked, a pretty string to print it to screen...)
      %  
      %  >> [vec,ket] = build_ket(state)
      %
      %  state :: integer representation of a basis state (bits are occupation numbers)
      %
      %  This depends entirely on the Fock basis conventions choosen in all
      %  ED-based solvers from QcmPlab (and many other codes, actually)
      %
      %  Copyright 2022 Gabriele Bellomia
      %
      Nlat = 2;
      Norb = 1;
      for ilat = 1:Nlat
          for ispin = 1:2
              shift = (ilat-1)*Norb + (ispin-1)*Norb*Nlat;
              index = shift+(1:Norb);
              vec(index) = bitget(state,index);
          end
      end
      if nargout>1
          kup = num2str(vec(1:Norb*Nlat));
          kdw = num2str(vec(Norb*Nlat+1:end));
          ket = ['| ',strrep(kup,'1','↑'),' 〉⊗ ',...
              '| ',strrep(kdw,'1','↓'), ' 〉'];
          ket = strrep(ket,'0','•');
      end
  end
  %
  function print_basis()
   for state = 0:1:15
       [~,label] = build_ket(state);
       fprintf('%d\t',state+1)
       disp(label)
   end
  end
  %
  function smatrix = SlaterCondon(nmodes)
    % SLATERCONDON : Implementation of Slater-Condon rules for fermions
    %                It pre-computes all the <istate| cdg_is c_js |jstate>
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
    smatrix = zeros(nmodes,nmodes,N,N);
    for istate = 0:1:N-1 % nmode-orbital states
        for jstate = 0:1:N-1 % nmode-orbital states
            % ∑_s ∑_ij <istate| cdg_is c_js |jstate>
            for ispin=1:2
                for imode=1:nmodes
                    for jmode=1:nmodes
                        % Apply cdg_is to <istate|
                        ibra = build_ket(istate);
                        if ibra(imode+(ispin-1)*2)==0
                            continue
                        end
                        ibra(imode+(ispin-1)*2) = ibra(imode+(ispin-1)*2)-1;
                        isign = (-1)^(sum(ibra(1:imode+(ispin-1)*2)));
                        % Apply c_js to |jstate>
                        jket = build_ket(jstate);
                        if jket(jmode+(ispin-1)*2)==0
                            continue
                        end
                        jket(jmode+(ispin-1)*2) = jket(jmode+(ispin-1)*2)-1;
                        jsign = (-1)^(sum(jket(1:jmode+(ispin-1)*2)));
                        % Overlap <istate| cdg_is c_js |jstate>
                        if isequal(ibra,jket)
                            smatrix(imode+(ispin-1)*2,jmode+(ispin-1)*2,istate+1,jstate+1) = isign * jsign;
                        end
                    end
                end
            end
        end
    end
  end
 
