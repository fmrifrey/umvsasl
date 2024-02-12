%% Construct sequence for recon
fprintf("Setting up recon...\n");
t = tic;

D = asl3df_seq('ccfac',0.25); % Compress coils down to 25%

fprintf("Done. Time elapsed = %0.3fs\n", toc(t));

%% Assign volumes for each frame
for i = 1:D.aqdims(4) % aqdims(4) is the number of frames in asl3dflex
    D = D.setup_frame(1:2, ... % frame # - average together M0 frames
        'all', ... % shot # - splice all shots together
        i ... % echo # - each echo is a new frame
        );
    D = D.setup_frame(3:2:D.aqdims(2), ... % frame # - average together label frames
        'all', ... % shot # - splice all shots together
        i ... % echo # - each echo is a new frame
        );
    D = D.setup_frame(4:2:D.aqdims(2), ... % frame # - average together control frames
        'all', ... % shot # - splice all shots together
        i ... % echo # - each echo is a new frame
        );
    fprintf("Done. Time elapsed = %0.3fs\n", toc(t));
end

%% Create coil sensitivity map
if D.ncoils > 1
    fprintf("Creating sensitivity map...\n");
    t = tic;
    
    % Create the sensitivity map
    D = D.make_sense(1);
    
    fprintf("Done. Time elapsed = %0.3fs\n", toc(t));
end

%% Recon
fprintf("Reconning...\n");
t = tic;

% Reconstruct the whole sequence
im = D.recon();

fprintf("Done. Time elapsed = %0.3fs\n", toc(t));

%% Perform pairwise subtraction
% fprintf("Subtracting...\n");
% t = tic;
% 
% % Reconstruct the whole sequence
% sub = im(:,:,:,6:2:end) - im(:,:,:,5:2:end);
% 
% fprintf("Done. Time elapsed = %0.3fs\n", toc(t));

%% Write to file
fprintf("Writing the images to file...\n");
t = tic;

% Write out
writenii('im_mag', abs(im));
writenii('im_ang', angle(im));

fprintf("Done. Time elapsed = %0.3fs\n", toc(t));
