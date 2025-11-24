%% Extension of our test to finite temperature
% For T=0 this script reproduces Fig.2 in [1], as well as it replicates 
% the results with our notion of symmetry-resolved negativity :) 
% [1] Lexin Ding et al., Quantum Science and Technology 9 (2024) 015005

t = 1;   % Hopping
L = 200; % Nsites
T = 0;   % Temperature in units of t
k = (1-mod(L,2)-round((L-1)/2):1:round((L-1)/2)) * 2*pi/L;
E = -2*t*cos(k); [E,indices] = sort(E); sorted_k = k(indices);

% Filling of the band η = N/(2*L)
global eta
fillings = [5e-3,1e-2:1e-2:0.09,0.1:0.01:0.5];
E_nSSR = zeros(4,length(fillings)); 
E_pSSR = zeros(4,length(fillings));
E_PPT_0 = zeros(4,length(fillings));
E_PPT = zeros(4,length(fillings));
for d=1:length(fillings)
eta = fillings(d) 
N = round(2*L*eta);
spin_k = zeros(2,round(N/2));
for n = 1:L
    for s = 1:2
        spin_k(s,n) = sorted_k(n);
        if n == round(N/2)
           EF = E(n);
        end
    end
end
for n = 1:L
    dE = E(n)-EF;
    if abs(dE) < 1e-12
        weight(n) = 1/2; % To manage NaNs close to the Fermi level
    else
        weight(n) = 1/(1+exp(dE/(T*t))); % Temperature in units of t
    end
end 


%% Numerical 1-body density matrix!
OBDM = zeros(2*L,2*L); % [1,2,...,Nsites]_up[1,2,...,Nsites]_dw
for i = 1:L
    OBDM(i,i) = eta;
    OBDM(i+L,i+L) = eta;
    for j = i+1:L
        % Spin-up
        OBDM(i,j) = sum(weight.*exp(1i*(i-j)*spin_k(1,:)))/L;
        % Spin-dw
        OBDM(i+L,j+L) = sum(weight.*exp(1i*(i-j)*spin_k(2,:)))/L;
        if not(i==j)
            OBDM(j,i) = conj(OBDM(i,j)); % up
            OBDM(j+L,i+L) = conj(OBDM(i+L,j+L)); % dw
        end
    end
end

% pre-compute Slater-Condon rules for 2 orbitals
sc_matrix = SlaterCondon(2); % 4D array [2*nmodes,2*nmodes,4^nmodes,4^nmodes]

%% Build two-orbital RDMs
OBDM_ij = zeros(4,4);
counter = 0;
for i = 1:1%L
    for j = [2,3,11] %101

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
        [N0,N1,N2] = sym_negativity(RDM);
        [~,~,~,Ntot] = negativity(RDM);
        assert(abs(N0+N1+N2-Ntot)<1e-12);
        E_PPT_0(counter,d) = N0;%log2(2*N0+1);
        E_PPT(counter,d) = N0+N2;%log2(2*(N0+N2)+1);
        [E_pSSR(counter,d),E_nSSR(counter,d)] = build_SSR(RDM);

    end
end

end

tiledlayout(1, 2);

nexttile

plot(fillings,E_pSSR(1,:)*log(2),'r-'); hold on
plot(fillings,E_nSSR(1,:)*log(2),'b-');

plot(fillings,E_pSSR(2,:)*log(2),'r--');
plot(fillings,E_nSSR(2,:)*log(2),'b--');

plot(fillings,E_pSSR(3,:)*log(2),'r-.');
plot(fillings,E_nSSR(3,:)*log(2),'b-.'); 

%plot(fillings,E_pSSR(4,:)*log(2),'r:');
%plot(fillings,E_nSSR(4,:)*log(2),'b:'); 

%set(gca,'Xscale','log')
%set(gca,'Yscale','log')
xlim([1e-4,0.5]); ylim([0,0.11]);
xlabel("$\eta$",'Interpreter','latex')
ylabel("$\mathcal{E}$",'Interpreter','latex')

nexttile

plot(fillings,E_PPT_0(1,:)*log(2),'c-'); hold on
plot(fillings,E_PPT(1,:)*log(2),'m-');

plot(fillings,E_PPT_0(2,:)*log(2),'c--');
plot(fillings,E_PPT(2,:)*log(2),'m--');

plot(fillings,E_PPT_0(3,:)*log(2),'c-.');
plot(fillings,E_PPT(3,:)*log(2),'m-.');

%plot(fillings,E_PPT_0(4,:)*log(2),'c:');
%plot(fillings,E_PPT(4,:)*log(2),'m:');

%set(gca,'Xscale','log')
%set(gca,'Yscale','log')
xlim([1e-4,0.5]); ylim([0,0.11]);
xlabel("$\eta$",'Interpreter','latex')
ylabel("$\mathcal{N}$",'Interpreter','latex')

print_basis