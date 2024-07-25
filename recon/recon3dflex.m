function x = recon3dflex(varargin)
% Function for performing CG-SENSE NUFFT reconstruction for umvsasl data
%
% by David Frey
% 
% Usage:
% Specify a data directory using the 'pfile' argument, or simply run this
%   function from the directory.
% The data directory must include the following:
%   - raw data pfile: (P*.7)
%   - kviews file (kviews*.txt)
%   - ktraj file (ktraj*.txt)
%
% Required paths:
%   - MIRT (git@github.com:JeffFessler/mirt.git)
%
% Arguments:
%   - pfile: pfile name search string, leave empty to use first P*.7 file
%       in current working directory
%   - smap: sensitivity map (must be [image size x ncoils]), leave empty
%       to compress coils
%   - niter: number of iterations for CG reconstruction
%   - coilwise: option to rearrange data for coil-wise recon of 1st frame,
%       this is useful for creating SENSE maps from fully-sampled data
%   - resfac: image space resolution upsampling factor
%   - ccfac: coil compression factor (i.e. ccfac = 4 will compressed data
%       to 1/4 of the channels)
%   - frames: frame indicies to reconstruct (default is all frames,
%       reconned sequentially)
%

    % check that mirt is set up
    aslrec.check4mirt();

    % set defaults
    defaults.pfile = [];
    defaults.smap = [];
    defaults.niter = 0;
    defaults.coilwise = 0;
    defaults.resfac = 1;
    defaults.ccfac = 1;
    defaults.frames = [];
    
    % parse input parameters
    args = vararg_pair(defaults,varargin);

    % get data from pfile
    [kdata,klocs,N,fov] = aslrec.read_data(args.pfile);
    if args.coilwise % rearrange for coil-wise reconstruction of frame 1 (for making SENSE maps)
        kdata = permute(kdata(:,:,1,:),[1,2,4,3]);
    end
    N = ceil(N*args.resfac); % upsample N
    
    % cut off first 50 pts of acquisition (sometimes gets corrupted)
    kdata(1:50,:,:,:) = [];
    klocs(1:50,:,:) = [];
    
    % get sizes
    nframes = size(kdata,3); % number of frames
    ncoils = size(kdata,4); % number of coils
    if isempty(args.frames)
        args.frames = 1:nframes; % default - use all frames
    end
    
    % do coil compression
    if isempty(args.smap) && (ncoils > 1)
        ncoils = 1;
        warning('sense map is empty, compressing data to 1 coil...');
        kdata = ir_mri_coil_compress(kdata,'ncoil',ncoils);
    elseif (args.ccfac > 1) && (size(args.smap,4) == ncoils)
        ncoils = ceil(ncoils/args.ccfac);
        fprintf('compressing data & SENSE map to %d coils...\n', ncoils);
        kdata = ir_mri_coil_compress(kdata,'ncoil',ncoils);
        args.smap = ir_mri_coil_compress(args.smap,'ncoil',ncoils);
    elseif size(args.smap,4) < ncoils
        ncoils = size(args.smap,4);
        warning('compressing data down to %d coils to match SENSE map...', ncoils);
        kdata = ir_mri_coil_compress(kdata,'ncoil',ncoils);
    end
    
    % set nufft arguments
    nufft_args = {N, 6*ones(1,3), 2*N, N/2, 'table', 2^10, 'minmax:kb'};

    % initialize x
    x = zeros([N(:)',length(args.frames)]);
    
    % calculate a new system operator
    omega = 2*pi*fov(:)'./N(:)'.*reshape(klocs,[],3);
    omega_msk = vecnorm(omega,2,2) < pi;
    omega = omega(omega_msk,:);
    A = Gnufft(true(N),[omega,nufft_args]); % NUFFT
    w = aslrec.pipedcf(A,3); % calculate density compensation
    if ncoils > 1 % sensitivity encoding
        A = Asense(A,args.smap);
    end
    
    % loop through frames and recon
    for i = 1:length(args.frames)
        framen = args.frames(i);

        % get data for current frame
        b = reshape(kdata(:,:,framen,:),[],ncoils);
        b = b(omega_msk,:);
        
        % initialize with density compensated adjoint solution
        fprintf("frame %d/%d: initializing solution x0 = A'*(w.*b)\n", i, length(args.frames))
        x0 = reshape( A' * (w.*b), N );
        x0 = ir_wls_init_scale(A, b, x0);
        
        % solve with CG
        x(:,:,:,i) = cg_solve(x0, A, b, args.niter, ...
            sprintf('frame %d/%d: ', i, length(args.frames))); % prefix the output message with frame number
        
    end
    
end

function x_star = cg_solve(x0, A, b, niter, msg_pfx)

    % set default message prefix
    if nargin < 5
        msg_pfx = '';
    end

    % loop through iterations of conjugate gradient descent
    x_set = zeros([size(x0),niter+1]);
    x_star = x0;
    x_set(:,:,:,1) = x_star;
    r = reshape(A'*(b - A*x_star), size(x_star));
    p = r;
    rsold = r(:)' * r(:);
    for n = 1:niter
        fprintf('%sCG iteration %d/%d, res: %.3g\n', msg_pfx, n, niter, rsold);
        
        % calculate the gradient descent step
        AtAp = reshape(A'*(A*p), size(x_star));
        alpha = rsold / (p(:)' * AtAp(:));
        x_star = x_star + alpha * p;
        x_set(:,:,:,n+1) = x_star;

        % calculate new residual
        r = r - alpha * AtAp;
        rsnew = r(:)' * r(:);
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;

        if exist('exitcg','var')
            break % set a variable called "exitcg" to exit at current iteration when debugging
        end

    end
    
end