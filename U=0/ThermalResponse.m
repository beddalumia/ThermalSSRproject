%% Local Dynamical Susceptibility from LDOS at finite T

clear all; clc;

k = 1;             % Boltzmann constant (8.617e-5 eV/K)
t = 1;             % NN hopping amplitude (our unit of energy)
T = 2;             % Temperature (units of the hopping)
mu = 0;            % Chemical potential (units of the hopping)

% Energy grid for LDOS
E = linspace(-5, 5, 2001); 
dE = E(2) - E(1);

% Dummy LDOS: A simple semi-circle band of width 2
rho = real(sqrt(2 - E.^2)) / (2*pi); 

% 1D chain LDOS (analytical, to take good care of divergences)
gloc = gloc_chain(E+1i*0.01,2*t);
%rho = -imag(gloc)/pi;

% Frequency grid for susceptibility (positive frequencies only)
omega = linspace(0.01, 20, 600);

%% Calculate the Imaginary Part (Absorption)
chi_imag = zeros(size(omega));

% Fermi-Dirac distribution
f_E = 1 ./ (exp((E - mu) / (k * T)) + 1);

% Compute the convolution (local bubble)
for i = 1:length(omega)
    w = omega(i);
    
    % Shift arrays and fill out-of-bounds with 0
    rho_shift = interp1(E, rho, E + w, 'linear', 0);
    f_shift   = interp1(E, f_E, E + w, 'linear', 0);
    
    integrand = rho .* rho_shift .* (f_E - f_shift);
    chi_imag(i) = pi * trapz(E, integrand);
end

%%%% >>> Careful, this requires the Signal Processing Toolbox <<<
%
% %% Calculate Real Part via Optimized Fast Kramers-Kronig
% % Step A: Symmetrize the imaginary part (odd function: chi''(-w) = -chi''(w))
% % We skip the first element of the flipped array to avoid duplicating w=0
% chi_imag_full = [-flip(chi_imag(2:end)), chi_imag];
% 
% % Step B: Apply the optimized fkkt function to the full symmetric array
% chi_real_full = fkkt(chi_imag_full);
% 
% % Step C: Extract only the positive frequencies
% % chi_imag_full has length 2*N - 1. The positive half starts exactly at index N.
% N_pos = length(chi_imag);
% chi_real = chi_real_full(N_pos:end);

%% Calculate the Real Part via folded Kramers-Kronig
chi_real = zeros(size(omega));

% Broadening parameter to regularize the pole at w' = w
% Make this roughly the size of your omega step (dw)
eta = omega(2) - omega(1); 

for i = 1:length(omega)
    w = omega(i);
    
    % Folded KK integrand: (w' * chi''(w')) / (w'^2 - w^2)
    % We add eta^2 to the denominator to prevent division by zero at the pole
    numerator = omega .* chi_imag;
    denominator = (omega.^2 - w^2);
    
    % Regularized integrand
    kk_integrand = (numerator .* denominator) ./ (denominator.^2 + eta^2);
    
    chi_real(i) = (2/pi) * trapz(omega, kk_integrand);
end

%% Plot the Results
figure;
plot(omega, chi_real, 'b-', 'LineWidth', 2); hold on;
plot(omega, chi_imag, 'r-', 'LineWidth', 2);
grid on;
xlabel('\omega', 'FontSize', 12);
ylabel('\chi_{loc}(\omega)', 'FontSize', 12);
%title('Dynamical Local Susceptibility at Finite T', 'FontSize', 14);
legend('\chi''', '\chi''''', 'FontSize', 12);


%% Utilities

function gloc = gloc_chain(zeta,D)
%% Local Green's function of a 1D tight-binding chain
%
%                   +π
%                   ⌠
%                   ⎮  dϕ          1                    1         
%           G(z) =  ⎮ ──── ⋅ ──────────────  =  ─────────────────
%                   ⎮  2π     z - D⋅cos(ϕ)      z⋅sqrt(1-(D/z)^2)
%                   ⌡
%                  -π
%
% As proven in https://doi.org/10.1007/3-540-28841-4_5 (Economou)
%
% Note: D is the half-bandwidth, i.e. twice the NN hopping amplitude.
    invz = 1./zeta;
    fact = D.*invz;
    gloc = invz./sqrt(1-fact.^2);    
end

function K = fkkt(F)
%% Fast Kramers-Kronig transform, exploiting NEXTPOW2 and built-in FFTW wrappers
%
%       K[F(..)](:) = FKKT(F), built through the chain:
%
%                   >   S[F] = hilbert(F) [Analytic-Signal]
%
%                   >   H[S] = imag(S)    [Hilbert-Transform]
%
%                   >   K[H] = -H         [Kramers-Kronig]
%
% See also NEXTPOW2, HILBERT, FFT
%
% Theoretical Background at:
%
%  https://en.wikipedia.org/wiki/Kramers–Kronig_relations
%
%  https://en.wikipedia.org/wiki/Hilbert_transform
%
% BSD 3-Clause License
%
%  Copyright (c) 2022, Gabriele Bellomia
%  All rights reserved.

    N = length(F)*4;        % Nyquist condition for the relevant FFTs
    P = pow2(nextpow2(N));  % FFT OPTIMIZATION TRICK: run < doc nextpow2 >
    S = hilbert(F,P);       % Wrapper of the FFTWs:   run < open hilbert >
    K = imag(-S(1:N/4));    % H = Im(AnalyticSignal): run < help hilbert >
    
end