function data_rs = regrid(data,varargin)
% function im = recon3dflex(varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to quickly regrid data given zoom and
%   resolution upsampling factors
%
%
% Notes:
%   - default values in help message may not be up to date - check defaults
%       structure under the function header
%
% Path dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%
% Static Intput arguments:
%   - data:
%       - data to resample
%       - Nd float/double array with cartesian data
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'resFactor'
%       - resolution (dimension) upsampling factor
%       - float/double array describing factor along each dimension
%       - passing a value > 1 will result in higher output image dimension
%       - default is ones(ndims(data),1)
%   - 'zoomFactor'
%       - field of view zoom factor
%       - float/double array describing factor along each dimension
%       - passing a value > 1 will result in a smaller fov
%       - default is ones(ndims(data),1)
%   - 'interp'
%       - interpolation method
%       - string describing interpolation method to use in
%           griddedInterpolant()
%       - default is 'cubic'
%   - 'extrap'
%       - extrapolation method
%       - string describing extrapolation method to use in
%           griddedInterpolant()
%       - default is 'none'
%
% Function output:
%   - data_rs:
%       - resampled data
%

    % Define default arguments
    defaults = struct( ...
        'zoomFactor',   ones(1,ndims(data)), ...
        'resFactor',    ones(1,ndims(data)), ...
        'shift',        zeros(1,ndims(data)), ...
        'interp',       'cubic', ...
        'extrap',       'none');
    
    % Parse through variable inputs using matlab's built-in input parser
    args = vararg_pair(defaults,varargin);
    
    % Pad undefined dimensions with 1
    if size(args.zoomFactor(:),1) == 1
        args.zoomFactor = args.zoomFactor*ones(1,ndims(data));
    elseif size(args.zoomFactor(:),1) > ndims(data)
        error('length(zoomFactor) must be <= ndim(data)');
    elseif size(args.zoomFactor(:),1) < ndims(data)
        args.zoomFactor = [args.zoomFactor(:)', ...
            ones(1,ndims(data)-size(args.zoomFactor,2))];
    end
    if size(args.resFactor(:),1) == 1
        args.resFactor = args.resFactor*ones(1,ndims(data));
    elseif size(args.resFactor(:),1) > ndims(data)
        error('length(resFactor) must be <= ndim(data)');
    elseif size(args.resFactor(:),1) < ndims(data)
        args.resFactor = [args.resFactor(:)', ...
            ones(1,ndims(data)-size(args.resFactor,2))];
    end
    
    % Shift the data
    data = circshift(data,args.shift);

    % Make grid arrays using imgrid
    G_in = cell(1,ndims(data));
    G_out = cell(1,ndims(data));
    [G_in{:}] = imgrid(1,size(data));
    [G_out{:}] = imgrid(1./args.zoomFactor,size(data).*args.resFactor(:)');
    
    % Make interpolant and save resampled data
    F = griddedInterpolant(G_in{:},data,args.interp,args.extrap);
    data_rs = F(G_out{:});

end

