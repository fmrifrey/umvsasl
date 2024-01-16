%% Construct sequence for recon
fprintf("Setting up recon...\n");
t = tic;

D = seq('ccfac',0.25); % Compress coils down to 25%

fprintf("Done. Time elapsed = %0.3fs\n", toc(t));

%% Assign volumes for each frame
for i = 1:D.aqdims(2) % aqdims(2) is the number of frames in asl3dflex
    D = D.setup_vol(i,'all','all');
    fprintf("Done. Time elapsed = %0.3fs\n", toc(t));
end

%% Create coil sensitivity map
if D.ncoils > 1
    fprintf("Creating sensitivity map...\n");
    t = tic;
    
    % Create the sensitivity map
    D = D.make_sense();
    
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
% writenii('sub_mag', abs(sub));
% writenii('sub_ang', angle(sub));

fprintf("Done. Time elapsed = %0.3fs\n", toc(t));
