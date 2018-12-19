function [avedk fl mA binningH binningV] = readGainCal(file)

% Script to read gain cal data and provide aprox signal level for a level
% of current.
% Really it is enough to provide dark, flood and current for the last one

% Open file
fid = fopen(file);
if (fid < 0),
    error('File could not be opened');
end

% Read metadata
binningH = fread(fid,1,'int32')
binningV = fread(fid,1,'int32')
D_sd    = fread(fid,1,'float32');
%D_sD    = 68;
kV      = fread(fid,1,'int32');
n_steps = fread(fid,1,'int32');
%n_steps = 25;
current = fread(fid,n_steps,'int32');

%extTrig = fread(fid,1,'int32')
expTime = fread(fid,1,'int32')
aveDark = fread(fid,1,'float32')

% Now read images
% Set sizes
sizeX = 3072/binningH
sizeY = 1944/binningV

% Read the data
imgTmp = zeros(sizeX,sizeY,n_steps);
for iStep = 1:n_steps,
    imgTmp(:,:,iStep) = fread(fid,[sizeX sizeY],'float32');
end

%dk = squeeze(imgTmp(:,:,1));
avedk = aveDark;
fl    = squeeze(imgTmp(:,:,n_steps));
mA    = current(n_steps)/1000.0;

end
