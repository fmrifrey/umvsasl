function x = recon3dflex(varargin)
% Function for performing CG-SENSE NUFFT reconstruction for umvsasl data
%
% by David Frey
%
% Required: MIRT (git@github.com:JeffFessler/mirt.git)
%
% Arguments:
%   - pfile: pfile name search string, leave empty to use first P*.7 file
%       in current working directory
%   - smap: sensitivity map (must be [image size x ncoils]), leave empty
%       to compress coils
%   - niter: number of iterations for CG reconstruction
%   - recycle: option to recycle solution of first frame as initialization
%       for all subsequent frames
%   - coilwise: option to rearrange data for coil-wise recon of 1st frame,
%       this is useful for creating SENSE maps from fully-sampled data
%   - ccfac: coil compression factor (i.e. ccfac = 4 will compressed data
%       to 1/4 of the channels)
%

    % check that mirt is set up
    aslrec.check4mirt();

    % set defaults
    defaults.pfile = [];
    defaults.smap = [];
    defaults.niter = 0;
    defaults.recycle = 0;
    defaults.coilwise = 0;
    defaults.ccfac = 1;
    
    % parse input parameters
    args = vararg_pair(defaults,varargin);

    % get data from pfile
    [kdata,klocs,N,fov] = aslrec.read_data(args.pfile);
    if args.coilwise % rearrange for coil-wise reconstruction of frame 1 (for making SENSE maps)
        kdata = permute(kdata(:,:,1,:),[1,2,4,3]);
    end
    
    % get sizes
    nframes = size(kdata,3); % number of frames
    ncoils = size(kdata,4); % number of coils
    
    % check SENSE map
    if isempty(args.smap) && ncoils > 1
        warning('sense map is empty, compressing data to 1 coil...');
        ncoils = 1;
        kdata = ir_mri_coil_compress(kdata,'ncoil',ncoils);
    elseif (args.ccfac > 1)
        ncoils = ceil(ncoils/args.ccfac);
        kdata = ir_mri_coil_compress(kdata,'ncoil',ncoils);
        args.smap = ir_mri_coil_compress(args.smap,'ncoil',ncoils);
    end
    
    % set nufft arguments
    nufft_args = {N, 6*ones(1,3), 2*N, N/2, 'table', 2^10, 'minmax:kb'};

    % initialize x
    x = zeros([N(:)',nframes]);
    
    % calculate a new system operator
    omega = 2*pi*fov(:)'./N(:)' .* reshape(klocs, [], 3);
    A = Gnufft(true(N),cat(2,{omega},nufft_args)); % NUFFT
    w = aslrec.pipedcf(A); % calculate density compensation
    if ncoils > 1 % sensitivity encoding
        A = Asense(A,args.smap);
    end
    
    % loop through frames and recon
    for framen = 1:nframes
        
        % get data for current frame
        b = reshape(kdata(:,:,framen,:),[],ncoils);
        
        % solve reconstruction
        fprintf('Reconstructing frame %d/%d\n', framen, nframes);
        if (framen > 1) && (args.niter > 0) && args.recycle
            x0 = x(:,:,:,1); % recycle solution to frame 1 as initialization
            args.niter = min(10,ceil(args.niter/4));
        else
            x0 = reshape(A' * (w.*b), N);
            x0 = ir_wls_init_scale(A,b,x0);
        end
        
        % loop through iterations of conjugate gradient descent
        x_star = x0;
        r = reshape(A'*b - A'*(A*x_star), size(x_star));
        p = r;
        rsold = r(:)' * r(:);
        x_set = zeros([N(:)',args.niter]);
        for n = 1:args.niter
            % calculate the gradient descent step
            AtAp = reshape(A'*(A*p), size(x_star));
            alpha = rsold / (p(:)' * AtAp(:));
            x_star = x_star + alpha * p;
            x_set(:,:,:,n) = x_star;
            
            % calculate new residual
            r = r - alpha * AtAp;
            rsnew = r(:)' * r(:);
            p = r + (rsnew / rsold) * p;
            rsold = rsnew;
        end
        x(:,:,:,framen) = x_star;
        
    end
    
end
