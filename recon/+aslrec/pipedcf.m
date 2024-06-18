function wi = pipedcf(G,itrmax)
% Function to generate kspace density weights for Gmri
%   reconstruction using Pipe & Menon method, as described in Pipe, J.G.,
%   Menon, P., (1999) Sampling density compensation in MRI: rationale and
%   an iterative numerical solution, Magn. Reson. Med. 41(1):179-86,
%   https://doi.org/10.1002/(sici)1522-2594(199901)41:1%3C179::aid-mrm25%3E3.0.co;2-v
%
% by David Frey
%

    % Set default for itrmax
    if nargin < 2 || isempty(itrmax)
        itrmax = 15;
    end
    
    % If G is a Gmri object, use its Gnufft object
    if isfield(G.arg,'Gnufft')
        G = G.Gnufft;
    end
    
    % Initialize weights to 1 (psf)
    wi = ones(size(G,1),1);
    
    % Loop through iterations
    for itr = 1:itrmax
        
        % Pipe algorithm: W_{i+1} = W_{i} / (G * (G' * W_{i}))
        d = real( G.arg.st.interp_table(G.arg.st, ...
            G.arg.st.interp_table_adj(G.arg.st, wi) ) );
        wi = wi ./ d;
        
    end
    
    % Normalize weights
    wi = wi / sum(abs(wi));
    
end

