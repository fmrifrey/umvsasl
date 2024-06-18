function check4mirt()

    % Check for mirt
    if ~exist('Gmri','file')
        ermsg = ['mirt must be set up for this to run!'
            ' --> pull it from git@github.com:JeffFessler/mirt.git'];
        error(ermsg);
    end
    
end

