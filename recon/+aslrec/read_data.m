function [kdata,klocs,N,fov] = read_data(pfile)
% Function to read in the pfile and .txt file data and format it for recon

    if nargin < 1 || isempty(pfile)
        pfile = './P*.7'; % default: use first Pfile on current path
    end
    
    % find and read the pfile
    tmp = dir(pfile);
    if isempty(tmp)
        error('no pfiles found from search string: %s', pfile);
    end
    pfile = tmp(1).name;
    pdir = tmp(1).folder;
    [raw,hdr] = aslrec.ge.read_pfile([pdir,'/',pfile]);
    raw(:,all(raw == 0, [1,3:5]),:,:,:) = []; % remove empty views
    raw(:,:,all(raw == 0, [1:2,4:5]),:,:) = []; % remove empty frames
    nviews = size(raw,2);
    nframes = size(raw,3);
    ncoils = size(raw,5);
    kdata = reshape(raw,[],nviews,nframes,ncoils);
    
    % find and read the ktraj file
    tmp = dir([pdir,'/ktraj*.txt']);
    ktrajfile = tmp(1).name;
    klocs0 = load(ktrajfile);
    
    % find and read the kviews file
    tmp = dir([pdir,'/kviews*.txt']);
    kviewsfile = tmp(1).name;
    kviews = load(kviewsfile);
    
    % transform kspace locations using rotation matrices
    klocs = zeros(size(klocs0,1),3,nviews); % klocs = [N x 3 x nviews]
    for viewn = 1:nviews
        R = reshape(kviews(viewn,end-8:end)',3,3)';
        klocs(:,:,viewn) = klocs0*R';
    end
    klocs = permute(klocs,[1,3,2]); % klocs = [N x nviews x 3]
    
    % save N and fov
    N = hdr.image.dim_X * ones(1,3);
    fov = hdr.image.dfov/10 * ones(1,3);
    
end
