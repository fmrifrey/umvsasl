function scaninfo = read_scaninfo()

% find the scaninfo.txt file and open it
name_si = dir('scaninfo*.txt');
name_si = [name_si(1).folder,'/',name_si(1).name];
fid = fopen(name_si,'r');

% loop through lines in the file
tline = fgetl(fid);
lineCounter = 1;
while ischar(tline)
    
    % split the line at whitespaces
    cv = regexp(tline,'(\w+\.\w+)|(\w+_\w+)|\w+','match');
    if length(cv) == 2 % only lines that contain a name and value
        scaninfo.(cv{1}) = str2double(cv{2}); % save cv name and value to struct
    end
    
    % read next line
    tline = fgetl(fid);
    lineCounter = lineCounter + 1;
    
end
fclose(fid);

end