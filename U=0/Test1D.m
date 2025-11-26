% Cross-check a numerical build of the 2-orbital RDM against the analyical
% results in Quantum Sci. Technol. 9 (2024) 015005 (Lexin Ding et al.)
% This script reproduces figure 2 in the paper, as well as it replicates 
% the results with our notion of symmetry-resolved negativity :)

t = 1;   % Hopping
L = 200; % Nsites
k = (1-mod(L,2)-round((L-1)/2):1:round((L-1)/2)) * 2*pi/L;
E = -2*t*cos(k); [E,indices] = sort(E); sorted_k = k(indices);

% Filling of the band η = N/(2*L)
global eta
fillings = [logspace(-4,-1,100),0.1:0.01:0.5];
E_nSSR = zeros(4,length(fillings)); 
E_pSSR = zeros(4,length(fillings));
E_PPT_0 = zeros(4,length(fillings));
E_PPT_2 = zeros(4,length(fillings));
for d=1:length(fillings)
eta = fillings(d) 
N = round(2*L*eta);
occupied = zeros(2,round(N/2));
for n = 1:round(N/2)
    for s = 1:2
        occupied(s,n) = sorted_k(n);
    end
end

%% Numerical 1-body density matrix!
OBDM = zeros(2*L,2*L); % [1,2,...,Nsites]_up[1,2,...,Nsites]_dw
for i = 1:L
    OBDM(i,i) = eta;
    OBDM(i+L,i+L) = eta;
    for j = i+1:L
        % Numerical computation
        % OBDM(i,j) = sum(weight.*exp(1i*(i-j)*spin_k(1,:)))/L; % up
        % OBDM(i+L,j+L) = sum(weight.*exp(1i*(i-j)*spin_k(2,:)))/L; % dw
        % Infinite chain limit:
        OBDM(i,j) = sin(pi*(i-j)*eta)/(pi*(i-j)); % up
        OBDM(i+L,j+L) = sin(pi*(i-j)*eta)/(pi*(i-j)); % dw
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

        RDM = RDM0(OBDM_ij,2,sc_matrix);
        
        %% Apply SSR and measure entanglement
        counter = counter + 1;
        [N0,~,N2] = sym_negativity(RDM);
        E_PPT_0(counter,d) = N0;%log2(2*N0+1);
        E_PPT_2(counter,d) = N0+N2;%log2(2*(N0+N2)+1);
        [E_pSSR(counter,d),E_nSSR(counter,d)] = build_SSR(RDM);


    end
end

end

figure("Position",[100 100 800 400])

tiledlayout(1, 2);

nexttile

plot(fillings,E_pSSR(1,:)*log(2),'r-'); hold on
plot(fillings,E_nSSR(1,:)*log(2),'b-');

plot(fillings,E_pSSR(2,:)*log(2),'r--');
plot(fillings,E_nSSR(2,:)*log(2),'b--');

% plot(fillings,E_pSSR(3,:)*log(2),'r-.');
% plot(fillings,E_nSSR(3,:)*log(2),'b-.'); 

% plot(fillings,E_pSSR(4,:)*log(2),'r:');
% plot(fillings,E_nSSR(4,:)*log(2),'b:'); 

%set(gca,'Xscale','log')
%set(gca,'Yscale','log')
xlim([0,0.5]); ylim([0,0.16]);
xlabel("$\eta$",'Interpreter','latex')
ylabel("$\mathcal{E}$",'Interpreter','latex')
legend(["P-SSR, NN","N-SSR, NN","P-SSR, NNN","N-SSR, NNN"],'Location','northwest')

nexttile

plot(fillings,E_PPT_2(1,:),'m-'); hold on
plot(fillings,E_PPT_0(1,:),'c-');

plot(fillings,E_PPT_2(2,:),'m--');
plot(fillings,E_PPT_0(2,:),'c--');

% plot(fillings,E_PPT_2(3,:),'m-.');
% plot(fillings,E_PPT_0(3,:),'c-.');

% plot(fillings,E_PPT_2(4,:),'m:');
% plot(fillings,E_PPT_0(4,:),'c:');

%set(gca,'Xscale','log')
%set(gca,'Yscale','log')
xlim([0,0.5]); ylim([0,0.16]);
xlabel("$\eta$",'Interpreter','latex')
ylabel("$\mathcal{N}$",'Interpreter','latex')
legend(["P-SSR, NN","N-SSR, NN","P-SSR, NNN","N-SSR, NNN"],'Location','northwest')

print_basis
