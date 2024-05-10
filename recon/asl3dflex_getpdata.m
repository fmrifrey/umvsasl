function [kdata,klocs,N,fov] = asl3dflex_getpdata(varargin)

% Assign defaults
defaults = struct( ...
    'pfile', './P*.7', ...
    'ndel', 0, ...
    'h5name', [], ...
    'nramp', [0,0], ...
    'frames', 'all', ...
    'echoes', 'all', ...
    'shots', 'all', ...
    'llorder', 0, ...
    'saveh5', 1 ...
    );

% Parse through arguments
args = vararginparser(defaults, varargin{:});

% Read in pfile
pfile = dir(args.pfile);
pfile = [pfile(1).folder,'/',pfile(1).name];
[raw,phdr] = readpfile(pfile); % raw = [ndat x nframes*nshots+2 x nechoes x 1 x ncoils]

% Read in kviews file
kviews = dir('./kviews*.txt');
kviews = [kviews(1).folder,'/',kviews(1).name];
kviews = load(kviews);

% Read in ktraj file
ktraj = dir('./ktraj*.txt');
ktraj = [ktraj(1).folder,'/',ktraj(1).name];
ktraj = load(ktraj);

% Save important header info
ndat = phdr.rdb.frame_size;
nframes = phdr.rdb.user1;
nshots = phdr.rdb.user2; 
nechoes = phdr.rdb.user3;
ncoils = phdr.rdb.dab(2) - phdr.rdb.dab(1) + 1;
dim = [phdr.image.dim_X,phdr.image.dim_X];
if phdr.rdb.user4 > 0 % if image is 3d
    dim = [dim, phdr.image.dim_X];
end
fov = phdr.image.dfov/10 * ones(size(dim));

% Reshape raw to make things a little easier
raw = reshape(raw(:,1:nshots*nframes,:,:,:), ...
    ndat, nshots, nframes, nechoes, ncoils); % raw = [ndat x nshots x nframes x nechoes x ncoils]
raw = permute(raw, [1,4,2,3,5]); % raw = [ndat x nechoes x nshots x nframes x ncoils]

% Loop through views and get transformed kspace locations
klocs = zeros(ndat,3,nechoes,nshots,nframes); % klocs = [ndat x 3 x nechoes x nshots x nframes]
for framen = 1:nframes
    for shotn = 1:nshots
        for echon = 1:nechoes
            % Get transformation matrix for current view
            matidx = (framen-1)*nshots*nechoes + (shotn-1)*nechoes + echon;
            T = reshape(kviews(matidx,end-8:end)',3,3)';

            % Append kspace locations and data
            klocs(:,:,echon,shotn,framen) = ktraj*T';
        end
    end
end

% Set frames to recon
if strcmpi(args.frames,'all')
    args.frames = 1:nframes;
end
nframes = length(args.frames);

% Set shots to recon
if strcmpi(args.shots,'all')
    args.shots = 1:nshots;
end
nshots = length(args.shots);

% Set echoes to recon
if strcmpi(args.echoes,'all')
    args.echoes = 1:nechoes;
end
nechoes = length(args.echoes);

% Remove unwanted data
klocs = klocs(:,:,args.echoes,args.shots,args.frames);
raw = raw(:,args.echoes,args.shots,args.frames,:);

% Apply sampling delay
raw = interp1(1:ndat,raw,(1:ndat)+args.ndel,'PCHIP',0);

% Reorder data for look-locker sequences
if args.llorder
    klocs = permute(klocs, [1,2,4,3,5]); % klocs = [ndat x 3 x nshots x nechoes x nframes]
    klocs = reshape(klocs, [ndat, 3, nshots, 1, nechoes*nframes]); % klocs = [ndat x 3 x nshots x 1 x nechoes*nframes]
    raw = permute(raw, [1,3,2,4,5]); % raw = [ndat x nshots x nechoes x nframes x ncoils]
    raw = reshape(raw, [ndat, nshots, 1, nechoes*nframes, ncoils]); % raw = [ndat x nshots x 1 x nechoes*nframes x ncoils]
    nframes = nechoes*nframes;
    nechoes = nshots;
    nshots = 1;
end

% Vectorize the data and locations
rampmask = args.nramp(1)+1:ndat-args.nramp(2);
ndat = length(rampmask);
klocs = reshape(permute(klocs(rampmask,:,:,:,:),[1,3,4,2,5]),ndat*nechoes*nshots,3,nframes);
kdata = reshape(raw(rampmask,:,:,:,:),ndat*nechoes*nshots,nframes,ncoils);

% Format the outputs
N = dim;
klocs = permute(klocs,[1,3,2]);

% save the data
if args.saveh5
   
    if isempty(args.h5name)
        fname = sprintf('%s_data.h5',pfile);
    else
        fname = [args.h5name,'.h5'];
    end
    
    if isfile(fname)
        warning('overwriting file %s\n',fname);
        system(['rm ',fname]);
    end
    
    h5create(fname,'/N',size(N));
    h5write(fname,'/N',N);
    
    h5create(fname,'/fov',size(fov));
    h5write(fname,'/fov',fov);
    
    h5create(fname,'/klocs',size(klocs));
    h5write(fname,'/klocs',klocs);
    
    h5create(fname,'/kdata_real',size(kdata));
    h5write(fname,'/kdata_real',real(kdata));
    
    h5create(fname,'/kdata_imag',size(kdata));
    h5write(fname,'/kdata_imag',imag(kdata));
    
end

end