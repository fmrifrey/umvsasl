im = rec(2,100);
    
function im = rec(delay,nramp,smap)

% Load in kspace trajectory & view transformation matrices
ktraj = load('ktraj.txt');
kviews = load('kviews.txt');

% Reshape transformation matrices as an array of 3x3 matrices
T = permute(reshape(kviews(:,end-8:end)',3,3,[]),[2,1,3]);

% Load in raw data
[raw,phdr] = readpfile;
ndat = phdr.rdb.frame_size;
nechoes = phdr.rdb.user2;
ncoils = phdr.rdb.dab(2) - phdr.rdb.dab(1) + 1;
ntrains = phdr.rdb.user1;
nframes = phdr.rdb.user0;
tr = phdr.image.tr*1e-3;
dim = 2*phdr.image.dim_X;
fov = phdr.image.dfov/10;

% reshape: ndat x ntrains*nframes x nechoes x 1 x ncoils
%           --> ndat x ntrains x nframes x nechoes x ncoils
raw = reshape(raw,ndat,ntrains,nframes + mod(nframes,2),nechoes,ncoils);
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
dcf = pipedcf(Gm.Gnufft);
W = Gdiag(dcf(:)./Gm.arg.basis.transform);

if nargin < 3 || isempty(smap)
    % Recon
    imc = zeros(dim,dim,dim,ncoils);
    for coiln = 1:ncoils
        data = reshape(raw(1,:,:,:,coiln),[],1);
        imc(:,:,:,coiln) = reshape(Gm' * (W*data(:)),dim,dim,dim);
    end
    im = sqrt(mean(imc.^2,4));
    writenii('rmsmag',abs(im));
else
    % Incorporate sensitivity encoding into system matrix
    Ac = repmat({[]},ncoils,1);
    for coiln = 1:ncoils
        tmp = smap(:,:,:,coiln);
        tmp = Gdiag(tmp(true(dim*ones(1,3))),'mask',true(dim*ones(1,3)));
        Ac{coiln} = Gm * tmp;
    end
    A = block_fatrix(Ac, 'type', 'col');
    
    W = Gdiag(repmat(dcf(:)./Gm.arg.basis.transform,1,ncoils));
    
    data = reshape(raw(1,:,:,:,:),[],1);
    im = A' * reshape(W * data, [], 1);
    im = embed(im,true(dim*ones(1,3)));
    writenii('sensemag', abs(im));
    writenii('senseang', angle(im));
end

figure, lbview(im);
title(sprintf("delay = %d", delay));
% figure, orthoview(im);

end