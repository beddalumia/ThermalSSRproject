%% Extension of our test to finite temperature
% For T=0 this script reproduces Fig.2 in [1], as well as it replicates 
% the results with our notion of symmetry-resolved negativity :) 
% [1] Lexin Ding et al., Quantum Science and Technology 9 (2024) 015005

t = 1;   % Hopping
L = 200; % Nsites
k = (1-mod(L,2)-round((L-1)/2):1:round((L-1)/2)) * 2*pi/L;
E = -2*t*cos(k); [E,indices] = sort(E); sorted_k = k(indices);

% Filling of the band η = N/(2*L)
global eta
fillings = [0.5,0.15];%[5e-3,1e-2:1e-2:0.09,0.1:0.01:0.5];
temperatures = linspace(0,4,100);
E_nSSR = zeros(4,length(fillings),length(temperatures)); 
E_pSSR = zeros(4,length(fillings),length(temperatures));
E_PPT_0 = zeros(4,length(fillings),length(temperatures));
E_PPT_2 = zeros(4,length(fillings),length(temperatures));
E_PPT = zeros(4,length(fillings),length(temperatures));

for d=1:length(fillings)

    %% Determine Fermi level
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

    for b = 1:length(temperatures)

        %% Determine thermal occupations
        T = temperatures(b); % Temperature in units of t 
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
        sc_matrix = slater_condon(2); % 4D array [2*nmodes,2*nmodes,4^nmodes,4^nmodes]

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

                RDM = RDM0(OBDM_ij,2,sc_matrix);
                
                %% Apply SSR and measure entanglement
                counter = counter + 1;
                [N0,N1,N2] = sym_negativity(RDM);
                [~,~,~,Ntot] = negativity(RDM);
                assert(abs(N0+N1+N2-Ntot)<1e-12);
                E_PPT_0(counter,d,b) = N0;%log2(2*N0+1);
                E_PPT_2(counter,d,b) = N0+N2;%log2(2*(N0+N2)+1);
                E_PPT(counter,d,b) = Ntot;
                [E_pSSR(counter,d,b),E_nSSR(counter,d,b)] = ssr_ree(RDM);

            end
        end
    end
end

figure("Name","Ground state SSR entanglement vs filling","Position",[100 100 800 400]);

tiledlayout(1, 2);

nexttile

plot(fillings,squeeze(E_pSSR(1,:,1))*log(2),'r-'); hold on
plot(fillings,squeeze(E_nSSR(1,:,1))*log(2),'b-');

plot(fillings,squeeze(E_pSSR(2,:,1))*log(2),'r--');
plot(fillings,squeeze(E_nSSR(2,:,1))*log(2),'b--');

plot(fillings,squeeze(E_pSSR(3,:,1))*log(2),'r-.');
plot(fillings,squeeze(E_nSSR(3,:,1))*log(2),'b-.'); 

%plot(fillings,E_pSSR(4,:,1)*log(2),'r:');
%plot(fillings,E_nSSR(4,:,1)*log(2),'b:'); 

%set(gca,'Xscale','log')
%set(gca,'Yscale','log')
xlim([1e-4,0.5]); ylim([0,0.11]);
xlabel("$\eta$",'Interpreter','latex')
ylabel("$\mathcal{E}$",'Interpreter','latex')

nexttile

plot(fillings,squeeze(E_PPT_2(1,:,1)),'m-'); hold on
plot(fillings,squeeze(E_PPT_0(1,:,1)),'c-'); 

plot(fillings,squeeze(E_PPT_2(2,:,1)),'m--');
plot(fillings,squeeze(E_PPT_0(2,:,1)),'c--');

plot(fillings,squeeze(E_PPT_2(3,:,1)),'m-.');
plot(fillings,squeeze(E_PPT_0(3,:,1)),'c-.');

%plot(fillings,squeeze(E_PPT_2(4,:,1))*log(2),'m:');
%plot(fillings,squeeze(E_PPT_0(4,:,1))*log(2),'c:');

%set(gca,'Xscale','log')
%set(gca,'Yscale','log')
xlim([1e-4,0.5]); ylim([0,0.11]);
xlabel("$\eta$",'Interpreter','latex')
ylabel("$\mathcal{N}^\mathrm{F}$",'Interpreter','latex')

figure("Name","Thermal SSR entanglement at half filling","Position",[100 100 800 400]);

tiledlayout(1, 2);

nexttile

plot(temperatures,squeeze(E_pSSR(1,1,:))*log(2),'r-'); hold on
plot(temperatures,squeeze(E_nSSR(1,1,:))*log(2),'b-');

%set(gca,'Xscale','log')
%set(gca,'Yscale','log')
%xlim([1e-4,0.5]); ylim([0,0.11]);
xlabel("$T/t$",'Interpreter','latex')
ylabel("$\mathcal{E}$",'Interpreter','latex')
legend(["P-SSR, NN","N-SSR, NN"],'Location','northeast')
ylim([0,0.5]) 

nexttile

plot(temperatures,squeeze(E_PPT(1,1,:)),'k-'); hold on
plot(temperatures,squeeze(E_PPT_2(1,1,:)),'m-'); hold on
plot(temperatures,squeeze(E_PPT_0(1,1,:)),'c-'); hold on


%set(gca,'Xscale','log')
%set(gca,'Yscale','log')
%xlim([1e-4,0.5]); ylim([0,0.11]);
xlabel("$T/t$",'Interpreter','latex')
ylabel("$\mathcal{N}^\mathrm{F}$",'Interpreter','latex')
legend(["Full, NN","P-SSR, NN","N-SSR, NN"],'Location','northeast')
ylim([0,0.5]) 

figure("Name","Thermal SSR entanglement at 0.15 filling","Position",[100 100 800 400]);

tiledlayout(1, 2);

nexttile

plot(temperatures,squeeze(E_pSSR(1,2,:))*log(2),'r-'); hold on
plot(temperatures,squeeze(E_nSSR(1,2,:))*log(2),'b-');

plot(temperatures,squeeze(E_pSSR(2,2,:))*log(2),'r--');
plot(temperatures,squeeze(E_nSSR(2,2,:))*log(2),'b--');

%set(gca,'Xscale','log')
%set(gca,'Yscale','log')
%xlim([1e-4,0.5]); ylim([0,0.11]);
xlabel("$T/t$",'Interpreter','latex')
ylabel("$\mathcal{E}$",'Interpreter','latex')
ylim([0,0.042]) 
legend(["P-SSR, NN","N-SSR, NN","P-SSR, NNN","N-SSR, NNN"],'Location','northeast')

nexttile

plot(temperatures,squeeze(E_PPT(1,2,:))*log(2),'k-'); hold on
plot(temperatures,squeeze(E_PPT_2(1,2,:))*log(2),'m-');
plot(temperatures,squeeze(E_PPT_0(1,2,:))*log(2),'c-');

plot(temperatures,squeeze(E_PPT(2,2,:))*log(2),'k--'); hold on
plot(temperatures,squeeze(E_PPT_2(2,2,:))*log(2),'m--');
plot(temperatures,squeeze(E_PPT_0(2,2,:))*log(2),'c--');

% plot(temperatures,squeeze(E_PPT(3,2,:))*log(2),'ko'); hold on
% plot(temperatures,squeeze(E_PPT_2(3,2,:))*log(2),'m-.');
% plot(temperatures,squeeze(E_PPT_0(3,2,:))*log(2),'c-.');

%set(gca,'Xscale','log')
%set(gca,'Yscale','log')
%xlim([1e-4,0.5]); ylim([0,0.11]);
xlabel("$T/t$",'Interpreter','latex')
ylabel("$\mathcal{N}^\mathrm{F}$",'Interpreter','latex')
legend(["Full, NN","P-SSR, NN","N-SSR, NN","Full, NNN","P-SSR, NNN","N-SSR, NNN"],'Location','northeast')
ylim([0,0.042]) 

print_basis