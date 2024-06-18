function x = recon3dflex(varargin)
% Function for performing simple NUFFT-SENSE (0-iteration) reconstruction
% for asl3dflex data
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
%   - recycle: option to recycle NUFFT system operators from frame to frame
%   - coilwise: option to rearrange data for coil-wise recon of 1st frame,
%       this is useful for creating SENSE maps from fully-sampled data
%   - denscomp: option to use pipe-menon density compensation
%

    % check that mirt is set up
    aslrec.check4mirt();

    % set defaults
    defaults.pfile = [];
    defaults.smap = [];
    defaults.recycle = 1;
    defaults.coilwise = 0;
    defaults.denscomp = 1;
    
    % parse input parameters
    args = vararg_pair(defaults,varargin);

    % get data from pfile
    [kdata,klocs,N,fov] = aslrec.read_data(args.pfile);
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
            
            if args.denscomp
                w = aslrec.pipedcf(A); % calculate density weighting
            else
                w = ones(A.odim(1),1);
            end
            
            if ncoils > 1 % sensitivity encoding
                A = Asense(A,args.smap);
            end
            
        elseif any(abs(klocs(:,:,:,framen) - klocs(:,:,:,1)) > 0)
            % warn if kspace is not same on this frame
            warning('klocs{framen} != klocs{1}, recycling the system operator may be incorrect!\n');
        end
        
        % get data for current frame
        y = reshape(kdata(:,:,framen,:),[],ncoils);
        
        % solve reconstruction
        x(:,:,:,framen) = reshape(A' * (w.*y), N);
        
    end
    
end
