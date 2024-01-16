classdef seq
    
    properties
        kspace % raw complex data from pfile with kspace trajectory
        vols % data parsed into inidividual timeseries volumes
        
        % Data acquisition properties
        aqdims
        ncoils
        smap
        tr
        te
        
        % Recon image space properties
        dim
        fov
    end
    
    methods
        
        function obj = seq(varargin)
            
            % Assign defaults
            defaults = struct( ...
                'pfile', './P*.7', ...
                'kviews', './kviews*.txt', ...
                'ktraj', './ktraj*.txt', ...
                'ccfac', 1.0 ...
                );
            
            % Parse through arguments
            args = vararginparser(defaults, varargin{:});
            
            % Read in pfile
            pfile = dir(args.pfile);
            pfile = [pfile(1).folder,'/',pfile(1).name];
            [raw,phdr] = readpfile(pfile);
            
            % Read in kviews file
            kviews = dir(args.kviews);
            kviews = [kviews(1).folder,'/',kviews(1).name];
            kviews = load(kviews);
            
            % Read in ktraj file
            ktraj = dir(args.ktraj);
            ktraj = [ktraj(1).folder,'/',ktraj(1).name];
            ktraj = load(ktraj);
            
            % Save important header info
            obj.aqdims = [phdr.rdb.frame_size, phdr.rdb.user1, ...
                phdr.rdb.user2, phdr.rdb.user3]; % ndat, frames, shots, echoes
            obj.ncoils = phdr.rdb.dab(2) - phdr.rdb.dab(1) + 1;
            obj.dim = phdr.image.dim_X;
            obj.fov = phdr.image.dfov/10;
            obj.tr = phdr.image.tr*1e-3;
            obj.te = phdr.image.te*1e-3;
            
            % Set up the kspace matrix
            obj.kspace = repmat(struct('data',[],'locs',[]), obj.aqdims(2:4));
            
            % Compress coils
            if (args.ccfac < 1 && obj.ncoils > 1)
                nvcoils = ceil(args.ccfac*obj.ncoils);
                raw = ir_mri_coil_compress(raw, 'ncoil', nvcoils);
                obj.ncoils = nvcoils;
            end
            
            % Loop through acquisition dimensions
            Ni = obj.aqdims(2);
            Nj = obj.aqdims(3);
            Nk = obj.aqdims(4);
            for i = 1:Ni
                for j = 1:Nj
                    for k = 1:Nk
    
                        % Get transformation matrix for current index
                        matidx = (i-1)*Nj*Nk + (j-1)*Nk + k;
                        T = reshape(kviews(matidx,end-8:end)',3,3)';
                        
                        % Calculate kspace trajectory for current index
                        klocs = ktraj*T';
                        
                        % Assign data & locs to entry of kspace matrix
                        obj.kspace(i,j,k) = struct( ...
                            'data', squeeze(raw(:,(i-1)*Nj + j,k,1,:)), ...
                            'locs', klocs ...
                            );
                    
                    end
                end
            end
            
            % Initialize volumes array and sensitivity map
            obj.vols = [];
            obj.smap = [];
            
        end
        
        function obj = setup_vol(obj,aqidx1,aqidx2,aqidx3)
            
            vidx = length(obj.vols) + 1;
            fprintf('Setting up volume %d...\n', vidx);
            
            % Set 'all' case
            if strcmpi(aqidx1,'all')
                aqidx1 = 1:obj.aqdims(2);
            end
            
            % Set 'all' case
            if strcmpi(aqidx2,'all')
                aqidx2 = 1:obj.aqdims(3);
            end
            
            % Set 'all' case
            if strcmpi(aqidx3,'all')
                aqidx3 = 1:obj.aqdims(4);
            end
            
            % Initialize kspace locations/data
            klocs = [];
            kdat = [];
            
            % Loop through and concat acquisition data for desired volume
            for i = aqidx1
                for j = aqidx2
                    for k = aqidx3
                        klocs = [klocs; obj.kspace(i,j,k).locs];
                        kdat = [kdat; squeeze(obj.kspace(i,j,k).data)];
                    end
                end
            end
            
            % Check if we can copy the Gmri object / dcf from another vol
            twinidx = [];
            for i = 1:vidx-1
                if isequal(obj.vols(i).klocs, klocs)
                    twinidx = i;
                end
            end
            
            % Construct the volume
            if isempty(twinidx) % if no Gmri object already exists
                obj.vols = [obj.vols;
                    vol( klocs, kdat, obj.dim, obj.fov, [], [] )
                    ];
            else % if we can copy a Gmri object
                fprintf('Copying Gmri object and dcf from vol %d\n', twinidx);
                obj.vols = [obj.vols;
                    vol(klocs, kdat, obj.dim, obj.fov, obj.vols(twinidx).Gm, obj.vols(twinidx).dcf)
                    ];
            end
            
        end
        
        function obj = make_sense(obj)
            
            % Reconstruct the coil x coil image
            [~,imc] = obj.vols(1).recon();
            
            % Estimate sensitivity map using eSPIRIT
            obj.smap = bart('ecalib -b0 -m1', fftc(imc,1:3));
                
        end
        
        function im = recon(obj,varargin)
            
            % Set default arguments
            defaults = struct( ...
                'niter', 0, ...
                'smap', obj.smap, ...
                'parallelize', 1 ...
                );
            
            % Parse through default arguments
            args = vararginparser(defaults, varargin{:});
            
            % Initialize image matrix
            im = zeros(obj.dim,obj.dim,obj.dim,length(obj.vols));
            
            % Reconstruct
            if args.parallelize % using parallel pool
                % Copy variables to avoid broadcasting
                pp_vols = obj.vols;
                pp_smap = args.smap;
                niter = args.niter;
                
                % Loop through volumes
                parfor i = 1:length(obj.vols)
                    im(:,:,:,i) = pp_vols(i).recon( ...
                        'smap', pp_smap, 'niter', niter);
                end
            else % without using parallel pool
                % Loop through volumes
                for i = 1:length(obj.vols)
                    im(:,:,:,i) = obj.vols.recon( ...
                        'smap', args.smap, 'niter', args.niter);
                end
            end
                
        end
        
    end
end

