function [M, t, B1, G] = simpulse(varargin)

    defaults.B1max = 117; % peak B1 amplitude (mG)
    defaults.Gmax = 1.5; % peak gradient amplitude (G/cm)
    defaults.x0 = 0; % initial position (cm)
    defaults.v = 0; % velocity (cm/s)
    defaults.T1 = Inf; % longitudinal magnetization relaxation constant (s)
    defaults.T2 = Inf; % transverse magnetization relaxation constant (s)
    defaults.dt = 4e-6; % sampling interval (s)
    
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
    
    % loop through time points
    M = zeros(3, length(t));
    M(:,1) = M0; % initialize the first magnetization vector
    
    % define bloch equation
    function dMdt = bloch(M, B)
        gam = 4.258 * 2 * pi; % gyromagnetic ratio (radians/s/mG)
        dMdt = cross(M, gam * B) - ... % precession
               [M(1) / args.T2; M(2) / args.T2; 0] - ... % transverse relaxation
               [0; 0; (M(3) - M0(3)) / args.T1]; % longitudinal relaxation
    end

    % RK4 integration loop
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
