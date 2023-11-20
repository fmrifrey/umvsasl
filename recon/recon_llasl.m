%% Construct sequence for recon
fprintf("Setting up recon...\n");
t = tic;

D = seq('ccfac',0.25); % Compress coils down to 25%

fprintf("Done. Time elapsed = %0.3fs\n", toc(t));

%% Assign volumes
for i = 1:D.aqdims(3)
    fprintf("Assigning volumes for echo %d...\n", i);
    t = tic;
    
    % Control
    D = D.setup_vol(1:2:aqdims(2),i); % All odd frame makes up the control volume for echo i
    % Label
    D = D.setup_vol(2:2:aqdims(2),i); % All even frame makes up the label volume for echo i
    
    fprintf("Done. Time elapsed = %0.3fs\n", toc(t));
end

%% Create coil sensitivity map
fprintf("Creating sensitivity map...\n");
t = tic;

% Create the sensitivity map
D = D.make_sense();

fprintf("Done. Time elapsed = %0.3fs\n", toc(t));

%% Recon
fprintf("Reconning...\n");
t = tic;

% Reconstruct the whole sequence
im = D.recon();

fprintf("Done. Time elapsed = %0.3fs\n", toc(t));

%% Perform pairwise subtraction
fprintf("Subtracting...\n");
t = tic;

% Reconstruct the whole sequence
sub = im(:,:,:,1:2:end) - im(:,:,:,2:2:end);

fprintf("Done. Time elapsed = %0.3fs\n", toc(t));

%% Write to file
fprintf("Writing the images to file...\n");
t = tic;

% Write out
writenii('im_mag', abs(im));
writenii('im_ang', angle(im));
writenii('sub_mag', abs(sub));
writenii('sub_ang', angle(sub));

fprintf("Done. Time elapsed = %0.3fs\n", toc(t));
