%% Local Dynamical Susceptibility from LDOS at finite T

clear all; clc;

k = 1;             % Boltzmann constant (8.617e-5 eV/K)
t = 1;             % NN hopping amplitude (our unit of energy)
T = 2;             % Temperature (units of the hopping)
mu = 0;            % Chemical potential (units of the hopping)
lattice="Bethe";% "1D-chain" or "Bethe" or "2D-square"

% Single-shot computation to check the dynamical structure of chi
[chi_real, chi_imag, omega] = response(t, T, mu, lattice);
figure;
plot(omega, chi_real, 'b-', 'LineWidth', 2); hold on;
plot(omega, chi_imag, 'r-', 'LineWidth', 2);
grid on;
xlabel('\omega', 'FontSize', 12);
ylabel('\chi_{loc}(\omega)', 'FontSize', 12);
title(sprintf('Local dynamical response at $T=%g$, $\\mu=%g$ [%s]',...
    T,mu,lattice), 'Interpreter', 'latex');
legend('\chi''', '\chi''''', 'FontSize', 12);

% Loop over T and \mu values to see if something changes at small fillings
temperatures = [0.1:0.05:5]*t;
chemicalpots = [-4:0.05:0]*t;
equal_time_absorption = zeros(length(temperatures),length(chemicalpots));
for iT = 1:length(temperatures)
    T = temperatures(iT);
    for iN = 1:length(chemicalpots)
        mu = chemicalpots(iN); 
        fprintf('T = %g\n',T)
        fprintf('mu = %g\n',mu)
        [chi_real, chi_imag, omega] = response(t, T, mu, lattice);
        equal_time_absorption(iT,iN) = trapz(omega,chi_imag);
    end
end

% Contour plot
figure
imagesc(temperatures,chemicalpots,equal_time_absorption)
set(gca,'YDir','normal')
colorbar
xlabel('$T$','Interpreter','latex')
ylabel('$\mu$','Interpreter','latex')
title(sprintf("Lattice: %s",lattice))


function [chi_real,chi_imag,omega] = response(hopping,temperature,chempot,l)

    if nargin < 4
       l = "1D-chain";
    end

    k = 1;             % Boltzmann constant (8.617e-5 eV/K)
    t = hopping;       % NN hopping amplitude (our unit of energy)
    T = temperature;   % Temperature (units of the hopping)
    mu = chempot;      % Chemical potential (units of the hopping)

    % Energy grid for LDOS
    E = linspace(-5, 5, 2001); 
    dE = E(2) - E(1);
    
    %% Select LDOS

    switch(l)
    
        case("Bethe")
        % Dummy LDOS: A simple semi-circle band of width 4t
        rho = real(sqrt(4 - E.^2)) / (2*pi); 
    
        case("1D-chain")
        % 1D chain LDOS (analytical, to take good care of divergences)
        gloc = gloc_chain(E+1i*0.01,2*t);
        rho = -imag(gloc)/pi;

        case("2D-square")
        % 2D square LDOS (analytical, to take good care of divergences)
        gloc = gloc_square(E+1i*0.01,4*t);
        rho = -imag(gloc)/pi;

        case("2D-toy")
        % 2D toy-model LDOS (analytical, proposed by Economou)
        gloc = gloc_2dtoy(E+1i*0.01,4*t);
        rho = -imag(gloc)/pi;

        otherwise
        error("Invalid lattice!")

    end
    
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
    
end

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
% As proven in https://doi.org/10.1007/3-540-28841-4_5 [economou]
%
% Note: D is the half-bandwidth, i.e. twice the NN hopping amplitude.
    invz = 1./zeta;
    fact = D.*invz;
    gloc = invz./sqrt(1-fact.^2);    
end

function gloc = gloc_square(zeta,D)
%% Local Green's function of a 2D tight-binding square lattice
%
%  The 2d and 3d hypercubic lattices can be recasted as suitable elliptic 
%  integrals(†), for more info see [economou], [delves], [morita], [kogan].
%
% References:
%  [economou]   https://doi.org/10.1007/3-540-28841-4_5
%  [delves]     https://doi.org/10.1006/aphy.2001.6148
%  [kogan]      https://doi.org/10.4236/graphene.2021.101001
%  [morita]     https://doi.org/10.1063/1.1665693
%
% †) See also: elliptic
%
% Note: D is the half-bandwidth, i.e. four times the NN hopping amplitude.
    invz = 1./zeta;
    ellk = zellke(D^2*invz.^2);
    gloc = 2/pi.*invz.*ellk;
end

function gloc = gloc_2dtoy(zeta,D)
% Toy model for a generic 2D lattice, for when one wants to avoid elliptic
% integrals (but we are providing a fast implementation...).
% D is the half-bandwidth (=4t, if you want to compare to a square lattice)
% Ref. to https://doi.org/10.1007/3-540-28841-4_5 [economou]
    znum = zeta + D;
    zden = zeta - D;
    gloc = 1/(2*D).*log(znum./zden);
end

function [k,e] = zellke(m,tol)
%% ZELLKE Complete elliptic integrals. Allowing complex input.
%   [K,E] = ZELLKE(M) returns the value of the complete elliptic
%   integrals of the first and second kinds, evaluated for each
%   element of M. Domain: M ∈ ℂ (see [2] for details).
%   
%   [K,E] = ZELLKE(M,TOL) computes the complete elliptic integrals to
%   the accuracy TOL instead of the default TOL = EPS(CLASS(M)).
%
%   Some definitions of the complete elliptic integrals use the modulus
%   k instead of the parameter M. They are related by M = k^2.
%
%   For clarity:
%
%                 π
%                 ─
%                 2
%                 ⌠
%                 ⎮         dϕ
%         K(M) =  ⎮ ────────────────── , 
%                 ⎮    _______________
%                 ⎮   ╱          2
%                 ⎮ ╲╱  1 - M⋅sin (ϕ)
%                 ⌡
%                 0
%
%                 π
%                 ─
%                 2
%                 ⌠
%                 ⎮       _______________
%                 ⎮      ╱          2
%         E(M) =  ⎮ dϕ ╲╱  1 - M⋅sin (ϕ)  .
%                 ⎮
%                 ⌡
%                 0
%
%   Class support for input M:
%      float: double, single [MATLAB]
%      float: double [GNU Octave]
%
%   ZELLKE extends the native implementation provided by ELLIPKE,
%   generalizing the method of the arithmetic-geometric mean [1],
%   so to allow complex values of M, as described in [2].
%
%   The algorithm has been tested against ELLIPTICK and ELLIPTICE,
%   included in the SYMBOLIC MATH TOOLBOX: the returned values seem  
%   to match to machine precision, with an average ≈300x speedup.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%       Functions" Dover Publications, 1965, 17.6.
%   [2] B.C. Carlson, "Numerical computation of real or complex 
%       elliptic integrals" Numerical Algorithms volume 10, 
%       pages 13–26 (1995) [arXiv:math/9409227]
%
%   See also ELLIPKE, ELLIPTICK, ELLIPTICE.
%
%% BSD 3-Clause License
%
%  Copyright (c) 2022, Gabriele Bellomia
%  All rights reserved.

  if nargin<1
    error('Not enough input arguments.'); 
  end

  try
    classin = superiorfloat(m); % GNU Octave compatibility
  catch
    classin = 'double';
  end

  if nargin<2, tol = eps(classin); end

  if isempty(m), k = zeros(size(m),classin); e = k; return, end

  if ~isreal(tol)
      error('Second argument TOL must be real.');
  end

  if ~isscalar(tol) || tol < 0 || ~isfinite(tol)
    error('Second argument TOL must be a finite nonnegative scalar.');
  end

  a0 = 1;
  b0 = sqrt(1-m);
  c0 = NaN;
  s0 = m;
  i1 = 0; 
  mm = Inf;

  while mm > tol
        a1 = (a0+b0)/2;
        b1 = sqrt(a0.*b0);
        c1 = (a0-b0)/2;
        i1 = i1 + 1;
        w1 = 2^i1*norm(c1).^2;
        w2 = 2^i1*c1.^2;
        mm = max(w1(:));
        % Test for stagnation (may happen for TOL < machine precision)
        if isequal(c0, c1)
            error('ZELLKE did not converge. Consider increasing TOL.');
        end
        s0 = s0 + w2;  
        a0 = a1;  b0 = b1;  c0 = c1;
  end

  k = pi./(2*a1);
  e = k.*(1-s0/2);
  im = find(m==1);
  if ~isempty(im)
      e(im) = ones(length(im),1);
      k(im) = inf;
  end

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

