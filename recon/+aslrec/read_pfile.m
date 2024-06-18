function [data,hdr] = read_pfile(searchstr)
% Function to read pfile using orchestra functions
% by David Frey

    % import all ge functions
    import aslrec.ge.*

    % Set default for search string
    if nargin < 1 || isempty(searchstr)
        searchstr = 'P*.7';
    end

    % Find Pfile based on search string
    dirp = dir(searchstr);
    if size(dirp,1) > 1
        data = cell(size(dirp,1),1);
        h = cell(size(dirp,1),1);
        for n = 1:size(dirp,1) % Recurse through all pfiles
             pfile = [dirp(1).folder '/' dirp(1).name];
             [data_i,h_i] = read_pfile(pfile);
             data{n} = data_i;
             h{n} = h_i;
        end
        return;
    elseif size(dirp,1) < 1
        error('No Pfiles found for search string: %s',searchstr);
    end
    pfile = [dirp(1).folder '/' dirp(1).name];

    % Open pfile
    fid = fopen(pfile,'r','l');

    % Get revision
    if ~exist('rdbm_rev','var')
        ver = fread(fid,1,'float32');
        str = num2str(ver);
        rdbm_rev = str2double(str);
    end

    % Read rdb header
    if fseek(fid,0,'bof'), error('BOF not found'); end
    hdr.rdb = read_rdb_hdr(fid,rdbm_rev);

    % Read ps header
    if fseek(fid,hdr.rdb.off_ps,'bof')
        error('ps header offset BOF not found');
    end
    hdr.ps = read_psc_header(fid, rdbm_rev);

    % Read exam header
    if fseek(fid,hdr.rdb.off_exam,'bof')
        error('exam header offset BOF not found');
    end
    hdr.exam = read_exam_header(fid, rdbm_rev);

    % Read series header
    if fseek(fid,hdr.rdb.off_series,'bof')
        error('series header offset BOF not found');
    end
    hdr.series = read_series_header(fid,rdbm_rev);

    % Read image header
    if fseek(fid,hdr.rdb.off_image,'bof')
        error('image header offset BOF not found');
    end
    hdr.image = read_image_header(fid,rdbm_rev);

    % Read grad header (if exists)
    if isfield(hdr.rdb,'off_grad_data')
       if fseek(fid,hdr.rdb.off_grad_data,'bof')
           error('grad header offset BOF not found');
       end
       hdr.grad = read_grad_header(fid,rdbm_rev);
    end
    
    % Get data sizes from rbd header
    ndat = hdr.rdb.frame_size;
    nslices = hdr.rdb.nslices;
    nechoes = hdr.rdb.nechoes;
    nviews = hdr.rdb.nframes;
    ncoils = hdr.rdb.dab(2)-hdr.rdb.dab(1)+1;
    precision = sprintf('int%d=>int%d', ...
        hdr.rdb.point_size*8, hdr.rdb.point_size*8);
    te = hdr.image.te;
    
    % Read in data
    fseek(fid, hdr.rdb.off_data, 'bof');
    data = double(fread(fid, Inf, precision));
    data = data(1:2:end) + 1i*data(2:2:end);
    data = reshape(data,ndat,nviews+1,nslices,nechoes,ncoils);
    data = data(:,2:end,:,:,:); % Cut out the baseline

    % Close file
    fclose(fid);
    
end

