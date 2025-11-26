%% Extension of our test to the 2D square lattice
clear all
L = 20; % Nsites = L^2;
t = 1/4; % hopping amplitude
kspan = (1-mod(L,2)-round((L-1)/2):1:round((L-1)/2)) * 2*pi/L;
for ix = 1:L
    for iy = 1:L
        k{(ix-1)*L+iy} = [kspan(ix),kspan(iy)];
        E((ix-1)*L+iy) = -2*t*(cos(kspan(ix))+cos(kspan(iy)));
        Esurf(ix,iy) = E((ix-1)*L+iy);
    end
end
% imagesc(abs(Esurf)<1e-12)
[~,indices] = sort(E); sorted_k = k(indices);

% pre-compute Slater-Condon rules for 2 orbitals
sc_matrix = slater_condon(2); % 4D array [2*nmodes,2*nmodes,4^nmodes,4^nmodes]

% Filling of the band η = N/(2*L)
global eta
fillings = [0.5]%[1e-4,1e-3,1e-2:1e-2:0.09,0.1:0.01:0.5]; %[logspace(-4,-1,100),0.1:0.01:0.5]; 
E_nSSR = zeros(4,length(fillings)); 
E_pSSR = zeros(4,length(fillings));
E_PPT = zeros(4,length(fillings));
for d=1:length(fillings)
    eta = fillings(d) 
    N = round(2*L^2*eta);
    kx = zeros(round(N/2),2);
    ky = zeros(round(N/2),2);
    for n = 1:round(N/2)
        for s = 1:2
            kx(n,s) = sorted_k{n}(1);
            ky(n,s) = sorted_k{n}(2);
        end
    end

    %% Numerical 1-body density matrix
    OBDM = zeros(2*L^2,2*L^2); % [1,2,...,Nsites]_up[1,2,...,Nsites]_dw
    for ix = 1:L
        for iy = 1:L
            i = (ix-1)*L + iy;
            for jx = 1:L
                for jy = 1:L
                    j = (jx-1)*L + jy;
                    for n = 1:length(kx(:,1))
                        % Spin-up
                        OBDM(i,j) = OBDM(i,j) + exp(1i*((ix-jx)*kx(n,1)+(iy-jy)*ky(n,1))) / L^2;
                    end
                    for n = 1:length(kx(:,2))
                        % Spin-dw
                        OBDM(i+L^2,j+L^2) = OBDM(i+L^2,j+L^2) + exp(1i*((ix-jx)*kx(n,2)+(iy-jy)*ky(n,2))) / L^2;
                    end
                    % TODO: harder to vectorize in 2D, but we gotta do it :(
                end
            end
        end
    end

    %% Build two-orbital RDMs
    OBDM_ij = zeros(4,4);
    counter = 0;
    for i = 2%L
        for j = 1%[2,3,11,100]

            % Restrict the one-body density matrix to two orbitals
            % -> spin up
            OBDM_ij(1,1) = OBDM(i,i); 
            OBDM_ij(1,2) = OBDM(i,j); 
            OBDM_ij(2,1) = OBDM(j,i);
            OBDM_ij(2,2) = OBDM(j,j);
            % -> spin down
            OBDM_ij(3,3) = OBDM(i+L^2,i+L^2); 
            OBDM_ij(3,4) = OBDM(i+L^2,j+L^2); 
            OBDM_ij(4,3) = OBDM(j+L^2,i+L^2);
            OBDM_ij(4,4) = OBDM(j+L^2,j+L^2);
            % To check how close to 1 the eigvals are
            %  disp(eig(OBDM_ij))
            %  disp(sum(eig(OBDM_ij)))
            %  disp ********************

            RDM = RDM0(OBDM_ij,2,sc_matrix);
            
            %% Apply SSR and measure entanglement
            counter = counter + 1;
            E_PPT(counter,d) = negativity(RDM);
            [E_pSSR(counter,d),~,E_nSSR(counter,d)] = sym_negativity(RDM);


        end
    end

end

plot(fillings,E_pSSR(1,:)*log(2),'r-'); hold on
plot(fillings,E_nSSR(1,:)*log(2),'b-');
%plot(fillings,E_PPT(1,:)*log(2),'g-');

plot(fillings,E_pSSR(2,:)*log(2),'r--');
plot(fillings,E_nSSR(2,:)*log(2),'b--');
%plot(fillings,E_PPT(2,:)*log(2),'g--');

plot(fillings,E_pSSR(3,:)*log(2),'r-.');
plot(fillings,E_nSSR(3,:)*log(2),'b-.'); 
%plot(fillings,E_PPT(3,:)*log(2),'g-.');

plot(fillings,E_pSSR(4,:)*log(2),'r:');
plot(fillings,E_nSSR(4,:)*log(2),'b:'); 
%plot(fillings,E_PPT(4,:)*log(2),'g:');

set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlim([1e-4,0.5]);
xlabel("$\eta$",'Interpreter','latex')
ylabel("$E$",'Interpreter','latex')

print_basis
