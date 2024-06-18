function x = recon3dflex(varargin)
% Function for performing simple NUFFT-SENSE (0-iteration) reconstruction
% for asl3dflex data
%
% by David Frey
%
% Required: MIRT (git@github.com:JeffFessler/mirt.git)
%

    % check that mirt is set up
    aslrec.check4mirt();

    % set defaults
    defaults.pfile = []; % pfile search string
    defaults.smap = []; % sensitvity map ([im size x ncoils])
    defaults.recycle = 1; % option to recycle system operators from frame to frame
    defaults.coilwise = 0; % option to perform coil-wise reconstruction of first frame
    
    % parse input parameters
    args = vararg_pair(defaults,varargin);

    % get data from pfile
    [kdata,klocs,N,fov] = aslrec.get_pdata(args.pfile);
    if args.coilwise % rearrange for coil-wise reconstruction of frame 1 (for making SENSE maps)
        kdata = permute(kdata(:,:,1,:),[1,2,4,3]);
        klocs = repmat(klocs(:,:,:,1),[1,1,1,size(kdata,3)]);
    end
    
    % get sizes
    nframes = size(klocs,4); % number of frames
    ncoils = size(kdata,3); % number of coils

    % check SENSE map
    if isempty(args.smap) && ncoils > 1
        warning('sense map is empty, compressing data to 1 coil...');
        kdata = ir_mri_coil_compress(kdata,'ncoil',1);
        ncoils = 1;
    end
    
    % set nufft arguments
    nufft_args = {N, 6*ones(1,3), 2*N, N/2, 'table', 2^10, 'minmax:kb'};

    % initialize x
    x = zeros([N(:)',nframes]);
    
    % loop through frames and recon
    for framen = 1:nframes
        
        if (framen == 1) || ~args.recycle
            
            % calculate a new system operator
            omega = 2*pi*fov(:)'./N(:)' .* reshape(klocs(:,:,:,framen), [], 3);
            A = Gnufft(true(N),cat(2,{omega},nufft_args)); % NUFFT
            
            w = aslrec.pipedcf(A); % calculate density weighting
            
            if ncoils > 1 % sensitivity encoding
                A = Asense(A,args.smap);
            end
            
        end
        
        % get data for current frame
        y = reshape(kdata(:,:,framen,:),[],ncoils);
        
        % solve reconstruction
        x(:,:,:,framen) = reshape(A' * (w.*y), N);
        
    end
    
end
