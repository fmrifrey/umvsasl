%% Recon using RMS coils combination
[im,imc] = rec(0,10);
lbview(im);
writenii('im',abs(im));
%% Create sensitivity map and recon using SENSE
smap = bart('ecalib -b0 -m1', fftc(imc,1:3));
im = rec(0,10,smap);
lbview(im);
writenii('im',abs(im));

function [im,imc] = rec(delay,nramp,smap,nframes)

% Load in kspace trajectory & view transformation matrices
ktraj = dir('ktraj*.txt');
ktraj = load(ktraj(1).name);
kviews = dir('kviews*.txt');
kviews = load(kviews(1).name);

% Reshape transformation matrices as an array of 3x3 matrices
T = permute(reshape(kviews(:,end-8:end)',3,3,[]),[2,1,3]);

% Load in raw data
[raw,phdr] = readpfile;
if strcmpi(phdr.image.psdname,'asl3dflex') % new sequence
	ndat = phdr.rdb.frame_size;
	nechoes = phdr.rdb.user2;
	ncoils = phdr.rdb.dab(2) - phdr.rdb.dab(1) + 1;
	ntrains = phdr.rdb.user1;
	if nargin < 4 || isempty(nframes)
		nframes = phdr.rdb.user0;
	end
	tr = phdr.image.tr*1e-3;
	dim = phdr.image.dim_X;
	fov = phdr.image.dfov/10;
else % old sequence
	ndat = phdr.rdb.frame_size;
	nechoes = phdr.rdb.nslices;
	ncoils = phdr.rdb.dab(2) - phdr.rdb.dab(1) + 1;
	ntrains = phdr.image.user0;
	nframes = phdr.rdb.nframes / ntrains;
	tr = phdr.image.tr*1e-3;
	dim = phdr.image.dim_X;
	fov = phdr.image.dfov/10;
end

% reshape: ndat x ntrains*nframes x nechoes x 1 x ncoils
%           --> ndat x ntrains x nframes x nechoes x ncoils
raw = reshape(raw,ndat,ntrains,[],nechoes,ncoils);
% permute: ndat x ntrains x nframes x nechoes x ncoils
%           --> nframes x ndat x nechoes x ntrains x ncoils
raw = permute(raw,[3,1,4,2,5]);

% Allocate space for entire trajectory
ktraj_all = zeros(ndat,3,nechoes,ntrains);

% Transform each view
for trainn = 1:ntrains
    for echon = 1:nechoes
        % Index the transformation matrix for current view
        mtxi = (trainn-1)*nechoes + echon;
        
        % Transform the trajectory
        ktraj_all(:,:,echon,trainn) = ktraj*T(:,:,mtxi)';
    end
end

% Apply delay and remove ramp points
raw = circshift(raw,[0,delay,0,0,0]);
raw = raw(1:nframes,nramp+1:end-nramp,:,:,:,:);
ktraj_all = ktraj_all(nramp+1:end-nramp,:,:,:);
    
% Create Gmri object
kspace = [reshape(ktraj_all(:,1,:,:),[],1), ...
        reshape(ktraj_all(:,2,:,:),[],1), ...
        reshape(ktraj_all(:,3,:,:),[],1)];
nufft_args = {dim*ones(1,3),...
    6*ones(1,3),...
    2*dim*ones(1,3),...
    dim*ones(1,3)/2,...
    'table',...
    2^10,...
    'minmax:kb'};
Gm = Gmri(kspace, true(dim*ones(1,3)), ...
    'fov', fov, 'basis', {'rect'}, 'nufft', nufft_args(:)');

% Calculate density weighting matrix
dcf = pipedcf(Gm.Gnufft);
W = Gdiag(dcf(:)./Gm.arg.basis.transform);

% Initialize matrices
im = zeros([dim*ones(1,3),nframes]);
imc = zeros([dim*ones(1,3),ncoils,nframes]);

if nargin < 3 || isempty(smap) % Coil-wise recon with RMS coil combo
    
    for framen = 1:nframes
        parfor coiln = 1:ncoils
            % Recon data
            data = reshape(raw(framen,:,:,:,coiln),[],1);
            imc(:,:,:,coiln,framen) = reshape(Gm' * (W*data(:)),dim,dim,dim);
        end
        im(:,:,:,framen) = sqrt(mean(imc(:,:,:,:,framen).^2,4));
    end
    
else % CG-SENSE recon
    
    % Incorporate sensitivity encoding into system matrix
    Ac = repmat({[]},ncoils,1);
    for coiln = 1:ncoils
        tmp = smap(:,:,:,coiln);
        tmp = Gdiag(tmp(true(dim*ones(1,3))),'mask',true(dim*ones(1,3)));
        Ac{coiln} = Gm * tmp;
    end
    A = block_fatrix(Ac, 'type', 'col');
    
    % Reshape density weighting matrix
    W = Gdiag(repmat(dcf(:)./Gm.arg.basis.transform,1,ncoils));
    
    % Recon data
    parfor framen = 1:nframes
        data = reshape(raw(framen,:,:,:,:),[],1);
        im(:,:,:,framen) = embed(A' * reshape(W * data, [], 1), ...
            true(dim*ones(1,3)));
    end
    
end


% Save only first frame of imc
imc = imc(:,:,:,:,1);

end
