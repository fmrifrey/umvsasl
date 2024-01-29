classdef asl3df_frame
    
    properties
        klocs % kspace trajectory sampling locations for kspace volume
        kdata % kspace data at sampling locations
        Gm % Gmri object
        dcf % density compensation matrix
        fov % field of view
        dim % matrix size
        ncoils % number of coils
    end
    
    methods
        
        function obj = asl3df_frame(klocs, kdata, dim, fov, Gm, dcf)
            
            % Assign kspace locations and data
            obj.klocs = klocs;
            obj.kdata = kdata;
            
            % Assign image space properties
            obj.dim = dim;
            obj.fov = fov;
            obj.ncoils = size(kdata,2);
            
            % If no Gmri object is passed, create it
            if isempty(Gm) || nargin < 5
                nufft_args = {dim,...
                    6*ones(size(dim)),...
                    2*dim,...
                    dim/2,...
                    'table',...
                    2^10,...
                    'minmax:kb'};
                Gm = Gmri(klocs(:,1:length(dim)), true(dim), ...
                    'fov', fov, 'basis', {'rect'}, 'nufft', nufft_args(:)');
            end
            
            % If no density compensation function is passed, create it
            if isempty(dcf) || nargin < 6
                dcf = pipe_menon_dcf(Gm.Gnufft);
            end
            
            % Assign Gmri object and density compensation
            obj.Gm = Gm;
            obj.dcf = dcf;
            
        end
        
        function plotkspace(obj)
            if length(obj.dim) == 2 % 2D
                plot(obj.klocs(:,1),obj.klocs(:,2));
            else % 3D
                plot3(obj.klocs(:,1),obj.klocs(:,2),obj.klocs(:,3));
            end
        end
        
        function showraw(obj,type,coil)
            % Set default type
            if nargin < 2 || isempty(type)
                type = 'abs';
            end
            
            % Set coil default
            if nargin < 3 || isempty(coil)
                coil = 'mean';
            end
            
            % Get raw data
            raw = obj.kdata;
            
            % Determine coil combination to show
            if ischar(coil) && strcmpi(coil,'mean') % show coil-wise mean
                raw = mean(raw,2);
            elseif ischar(coil) && strcmpi(coil,'all') % show all coils
                raw = raw;
            else % show specific coil data
                raw = raw(:,coil);
            end
            
            switch lower(type)
                case 'abs' % plot magnitude
                    plot(abs(raw));
                case 'angle' % plot phase
                    plot(angle(raw));
                case 'real' % plot real
                    plot(real(raw));
                case 'imag' % plot imaginary
                    plot(imag(raw));
                case 'log' % plot log scale
                    raw(abs(raw(:)) < eps()) = eps();
                    plot(log(abs(raw)));
            end
            
        end
        
        function [im,imc] = recon(obj,varargin)
            
            % Set default arguments
            defaults = struct( ...
                'niter', 0, ...
                'smap', [] ...
                );
            
            % Parse through arguments
            args = vararginparser(defaults, varargin{:});
            
            % Convert dcf to matrix
            W = Gdiag(obj.dcf(:)./obj.Gm.arg.basis.transform);
            
            % Coil-wise recon with RMS coil combo
            if isempty(args.smap) || nargout == 2
                % Reconstruct the coil-wise images
                imc = reshape(obj.Gm'*(W*obj.kdata), [obj.dim,obj.ncoils]);

                if obj.ncoils == 1 % Save single coil data
                    im = imc;
                else % RMS coil combo
                    im = sqrt(mean(imc.^2,ndims(imc)));
                end
            end
            
            % PCG/CP-SENSE recon
            if ~isempty(args.smap)
                
                % Make regularizer
                R = Reg1(ones(obj.dim), 'beta', 2^-12 * numel(obj.kdata)/3, ...
                    'mask', true(obj.dim));
                C = R.C;
                
                % Incorporate sensitivity encoding into system matrix
                A = Asense(obj.Gm,args.smap);
                
                % Reshape density weighting matrix
                W = Gdiag(repmat(obj.dcf(:)./obj.Gm.arg.basis.transform,1,obj.ncoils));
                
                % Recon preconditioner using conjugate-phase
                im_cp = A' * reshape(W * obj.kdata, [], 1);
                im_cp = ir_wls_init_scale(A, obj.kdata(:), im_cp);
                im_cp = embed(im_cp,true(obj.dim));
                
                % Recon using preconditioned conjugate gradient (iterative)
                if args.niter > 0
                    im_pcg = qpwls_pcg1(im_cp(true(obj.dim)), A, 1, obj.kdata(:), 0, ...
                        'niter', nargs.iter);
                    im = embed(im_pcg,true(obj.dim));
                    
                else % ...or save image with CP recon
                    im = im_cp;
                end
                
            end
            
        end
        
    end
end

function Wi = pipe_menon_dcf(G,itrmax)

    % Set default for itrmax
    if nargin < 2 || isempty(itrmax)
        itrmax = 15;
    end
    
    % If G is a Gmri object, use its Gnufft object
    if isfield(G.arg,'Gnufft')
        G = G.Gnufft;
    end
    
    % Initialize weights to 1 (psf)
    Wi = ones(size(G,1),1);
    
    % Loop through iterations
    for itr = 1:itrmax
        
        % Pipe algorithm: W_{i+1} = W_{i} / (G * (G' * W_{i}))
        d = real( G.arg.st.interp_table(G.arg.st, ...
            G.arg.st.interp_table_adj(G.arg.st, Wi) ) );
        Wi = Wi ./ d;
        
    end
    
    % Normalize weights
    Wi = Wi / sum(abs(Wi));
    
end