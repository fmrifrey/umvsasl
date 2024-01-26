%% Construct sequence for recon
fprintf("Setting up recon...\n");
t = tic;

D = asl3df_seq('ccfac',1,'force3D',1); % Compress coils down to 25%

fprintf("Done. Time elapsed = %0.3fs\n", toc(t));

%% Assign volumes for each frame
for i = 1:D.aqdims(2) % aqdims(2) is the number of frames in asl3dflex
    D = D.setup_frame(i, ... % frame #
        'all', ... % shot #
        'all' ... % echo #
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
