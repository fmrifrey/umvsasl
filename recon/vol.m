classdef vol
    
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
        
        function obj = vol(klocs, kdata, dim, fov, Gm, dcf)
            
            % Assign kspace locations and data
            obj.klocs = klocs;
            obj.kdata = kdata;
            
            % Assign image space properties
            obj.dim = dim;
            obj.fov = fov;
            obj.ncoils = size(kdata,2);
            
            % If no Gmri object is passed, create it
            if isempty(Gm) || nargin < 5
                nufft_args = {dim*ones(1,3),...
                    6*ones(1,3),...
                    2*dim*ones(1,3),...
                    dim*ones(1,3)/2,...
                    'table',...
                    2^10,...
                    'minmax:kb'};
                Gm = Gmri(klocs, true(dim*ones(1,3)), ...
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
            plot3(obj.klocs(:,1),obj.klocs(:,2),obj.klocs(:,3));
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
                imc = reshape(obj.Gm'*(W*obj.kdata), [obj.dim*ones(1,3),obj.ncoils]);

                if obj.ncoils == 1 % Save single coil data
                    im = imc(:,:,:,1);
                else % RMS coil combo
                    im = sqrt(mean(imc.^2,4));
                end
            end
            
            % PCG/CP-SENSE recon
            if ~isempty(args.smap)
                
                % Make regularizer
                R = Reg1(ones(obj.dim*ones(1,3)), 'beta', 2^-12 * numel(obj.kdata)/3, ...
                    'mask', true(obj.dim*ones(1,3)));
                C = R.C;
                
                % Incorporate sensitivity encoding into system matrix
                A = Asense(obj.Gm,args.smap);
                
                % Reshape density weighting matrix
                W = Gdiag(repmat(obj.dcf(:)./obj.Gm.arg.basis.transform,1,obj.ncoils));
                
                % Recon preconditioner using conjugate-phase
                im_cp = A' * reshape(W * obj.kdata, [], 1);
                im_cp = embed(im_cp,true(obj.dim*ones(1,3)));
%                 im_cp = ir_wls_init_scale(A, obj.kdata(:), im_cp);
                
                % Recon using preconditioned conjugate gradient (iterative)
                if args.niter > 0
                    im_pcg = qpwls_pcg1(im_cp(true(obj.dim*ones(1,3))), A, 1, obj.kdata(:), 0, ...
                        'niter', nargs.iter);
                    im = embed(im_pcg,true(obj.dim*ones(1,3)));
                    
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