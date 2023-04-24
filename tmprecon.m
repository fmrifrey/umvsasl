% Load in kspace trajectory & view transformation matrices
ktraj = load('ktraj.txt');
kviews = load('kviews.txt');

% Reshape transformation matrices as an array of 3x3 matrices
T = permute(reshape(kviews(:,end-8:end)',3,3,[]),[2,1,3]);

% Load in raw data
[raw,phdr] = readpfile;
ndat = phdr.rdb.frame_size;
nechoes = phdr.rdb.nslices;
ncoils = phdr.rdb.dab(2) - phdr.rdb.dab(1) + 1;
ntrains = phdr.rdb.nechoes;
nframes = phdr.rdb.nframes - 1;
tr = phdr.image.tr*1e-3;
dim = phdr.image.dim_X;
fov = phdr.image.dfov/10;

% Correct the size
raw = raw(1:ndat,1:ntrains,1:nechoes,1:nframes,1:ncoils);

% Allocate space for entire trajectory
ktraj_all = zeros(ndat,3,ntrains,nechoes);

% Transform each view
for trainn = 1:ntrains
    for echon = 1:nechoes
        % Index the transformation matrix for current view
        mtxi = (trainn-1)*nechoes + echon;
        
        % Transform the trajectory
        ktraj_all(:,:,trainn,echon) = ktraj*T(:,:,mtxi)';
    end
end
    
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

% Recon
im = zeros(dim,dim,dim,ncoils);
for coiln = 1:ncoils
    data = reshape(raw(:,:,:,1,coiln),[],1);
    im(:,:,:,coiln) = reshape(Gm' * (W*data(:)),dim,dim,dim);
end