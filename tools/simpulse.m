function [M, t, B1, G] = simpulse(varargin)
% Function to simulate isochromat magnetization in response to umvsasl prep
%   pulses using bloch equation with RK4 integration method
%
% Run this function from a "pulse" directory in aslprep, containing the
%   files rho.txt, theta.txt, and grad.txt
%
% by David Frey
%
% Required: MIRT (git@github.com:JeffFessler/mirt.git)
%
% Arguments:
%   - B1max: peak B1 (mG) amplitude corresponding to uint16 max
%   - Gmax: peak gradient (G/cm) amplitude corresponding to uint16 max
%   - x0: initial displacement of isochromat (cm)
%   - v: isochromat velocity (cm/s)
%   - T1: isochromat longitudinal magnetization relaxation time constant (s)
%   - T2: isochromat transverse magnetization relaxation time constant (s)
%   - dt: sampling interval (s)
%

    % set defaults
    defaults.B1max = 117;
    defaults.Gmax = 1.5;
    defaults.x0 = 0;
    defaults.v = 0;
    defaults.T1 = Inf;
    defaults.T2 = Inf;
    defaults.dt = 4e-6;
    
    % parse input parameters
    args = vararg_pair(defaults, varargin);
    
    % load in and scale the pulse waveforms
    rho = load('rho.txt');
    theta = load('theta.txt');
    grad = load('grad.txt');
    
    B1 = args.B1max / 32766 * rho .* exp(1i * pi / 32766 * theta);
    G = 1e3*args.Gmax / 32766 * grad; % mG/cm
    t = args.dt * (0:length(G)-1);
    
    % ignore 2nd column (control case)
    B1 = B1(:,1);
    G = G(:,1);
    
    % set initial value
    M0 = [0; 0; 1];
    
    % initialize magnetization
    M = zeros(3, length(t));
    M(:,1) = M0;
    
    % define bloch equation
    function dMdt = bloch(M, B)
        gam = 4.258 * 2 * pi; % gyromagnetic ratio (radians/s/mG)
        dMdt = cross(M, gam * B) - ... % precession
               [M(1) / args.T2; M(2) / args.T2; 0] - ... % transverse relaxation
               [0; 0; (M(3) - M0(3)) / args.T1]; % longitudinal relaxation
    end

    % loop through time points & compute RK4 integration
    for i = 2:length(t)
        
        % compute the effective magnetic field
        B = [real(B1(i)); imag(B1(i)); 0] + [0; 0; G(i) * (args.x0 + args.v * t(i))];
        
        % compute k1, k2, k3, k4
        k1 = bloch(M(:,i-1), B) * args.dt;
        k2 = bloch(M(:,i-1) + k1 / 2, B) * args.dt;
        k3 = bloch(M(:,i-1) + k2 / 2, B) * args.dt;
        k4 = bloch(M(:,i-1) + k3, B) * args.dt;
        
        % update M
        M(:,i) = M(:,i-1) + (k1 + 2*k2 + 2*k3 + k4) / 6;
        
    end

end
