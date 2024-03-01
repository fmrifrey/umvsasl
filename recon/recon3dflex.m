function im = recon3dflex(varargin)
% function im = recon3dflex(varargin)
%
% Part of asl3dflex project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to reconstruct images from *3dflex ASL sequence
%   using Model-based SENSE recon with Pipe & Menon Density compensation
%
%
% Notes:
%   - this function does not write images, use writenii to write the
%       outputs to file
%   - default values in help message may not be up to date - check defaults
%       structure under the function header
%
% Path dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%   - mirt (matlab version)
%       - github: JeffFessler/mirt
%       - mirt setup must have successfully ran
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'pfile'
%       - search string for pfile
%       - string containing search path for desired pfile
%       - will only be used if either 'raw' or 'info' is left blank
%       - type 'help readpfile' for more information
%       - default is 'P*.7' (uses first pfile from directory to read in
%           raw data and info)
%   - 'niter'
%       - maximum number of iterations for model-based recon
%       - integer describing number of iterations
%       - if 0 is passed, conjugate phase recon will be used
%       - if a value less than 0 is passed, NUFFT recon will be used
%       - default is 0
%   - 'frames'
%       - frames in timeseries to recon
%       - integer array containing specific frames to recon
%       - if 'all' is passed, all frames will be reconned
%       - default is 'all'
%   - 'resfac'
%       - resolution (dimension) upsampling factor
%       - float/double describing factor
%       - passing a value > 1 will result in higher output image dimension
%       - default is 1
%   - 'zoomfac'
%       - field of view zoom factor
%       - float/double describing factor
%       - passing a value > 1 will result in a smaller output image fov
%       - default is 1
%   - 'smap'
%       - coil sensitivity map for multi-coil datasets
%       - complex double/float array of dim x ncoils representing
%           sensitivity for each coil, or 'espirit' to estimate
%       - default is empty
%   - 'ccfac'
%       - coil compression factor
%       - float from 0 to 1 describing factor of reduction in # of coils
%       - default is 1
%
% Function output:
%   - im:
%       - output timeseries image (coil combined)
%       - complex array of image dimension
%

% Assign defaults
defaults = struct( ...
    'pfile', './P*.7', ...
    'ccfac', 0.25, ...
    'resfac', 1, ...
    'zoomfac', 1, ...
    'frames', 'all', ...
    'niter', 0, ...
    'smap', [] ...
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
tr = phdr.image.tr*1e-3;
te = phdr.image.te*1e-3;

% Reshape raw to make things a little easier
raw = reshape(raw(:,1:nshots*nframes,:,:,:), ...
    ndat, nshots, nframes, nechoes, ncoils); % raw = [ndat x nshots x nframes x nechoes x ncoils]
raw = permute(raw, [1,4,2,3,5]); % raw = [ndat x nechoes x nshots x nframes x ncoils]

% Loop through views and get transformed kspace locations
klocs = [];
kdata = [];
for shotn = 1:nshots
    for echon = 1:nechoes
        % Get transformation matrix for current view
        matidx = (shotn-1)*nechoes + echon;
        T = reshape(kviews(matidx,end-8:end)',3,3)';
        
        % Append kspace locations and data
        klocs = cat(1,klocs,ktraj*T');
        kdata = cat(1,kdata,squeeze(raw(:,echon,shotn,:,:)));
        
    end
end

% Set frames to recon
if strcmpi(args.frames,'all')
    args.frames = 1:nframes;
end

% Recon using CG/CP SENSE
[im,smap] = recon_cgsense(klocs, kdata(:,args.frames,:), ...
    round(args.resfac*dim), fov/args.zoomfac, ...
    'compfrac', args.ccfac, ...
    'smap', args.smap, ...
    'niter', args.niter);

% Write out image
writenii('im_mag',abs(im));
writenii('im_ang',angle(im));

% Write out sensitivity map
writenii('smap_mag',abs(smap));
writenii('smap_ang',angle(smap));

end