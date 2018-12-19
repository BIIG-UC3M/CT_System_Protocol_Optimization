% Script to study the trend of the flood image in the detector as a
% function of the input current, at the moment only for binning 4
% Flood data should give us an estimation of the
% dynamic range for these conditions
%
% Using these data we can estimate the parameters for our detector model
% we need to know the slope of the noise/pixel number linear trend shown by
% the detector in the quantum limited part of the gain curve.
% We also have to know the additive noise present in the real detector
% which will be a combination of thermal and quantization noise together
% with further noise sources... The intercept of the noise curve and the
% noise in the dark data should give us an estimatimation.
clc
clear
% Get the data and read it
folder    = 'C:/Users/aeroesp/Documents/MATLAB/det model 35kV/1000ms/';
%folder    = 'W:/Users_No Backup/asisniega/HDR/flood_data/e125/';
fnameBase = 'det_model_1000ms';

% Currents and times
t_int = [1 1];             % Detector integration time (s)
cur_i = [12 32];                 % Init current (uA), usually 0
cur_e = [300 800];               % End current (uA)
nTime = 2;                 % Number of filters
nCurr = 25;                % Number of steps in current
sCurr = cur_e*(1/(nCurr)); % Current step
cur_a = zeros(nCurr,nTime);  % Define the matrix for anode current values
cur_s = zeros(nCurr,nTime);  % Exposure (mAs) - Not proper exposure, but integrated current will do the job...
% Loop to fill out values
for nT = 1:nTime,
    cur_a(:,nT) = (cur_i(nT):sCurr(nT):cur_e(nT))'; % Array of anode current values for the given integration time
    cur_s(:,nT) = cur_a(:,nT)*(t_int(nT)/1000);  % The 1000 factor stands for the conversion of uA to mA
end
% Detector params
bits = 14;           % Bit depth of the detector
maxB = (2^bits - 1); % Maximum px value given the bit depth

% Image params
nFr    = 25;        % Frames per stack
imSize = [768 486]; % Size of the image
imXSize = 81:728;   % Corrected X size of the image, to avoide the cone beam.

% There is a lot of data, lets process it package by package in a loop
% Storage
px_m_al = zeros(size(imXSize,2),imSize(2),nCurr); % Pixel mean (2D distribution) for each current
px_v_al = zeros(size(imXSize,2),imSize(2),nCurr); % Pixel variance (2D distribution) for each current
m_c_al  = zeros(nCurr,1);                   % Values for central pixel, due to heel effect, it is supposed to be saturated at some point
v_c_al  = zeros(nCurr,1);                   % Variances for central pixel, due to heel effect, it is supposed to be saturated at some point
px_m_cu = zeros(size(imXSize,2),imSize(2),nCurr); % Pixel mean (2D distribution) for each current
px_v_cu = zeros(size(imXSize,2),imSize(2),nCurr); % Pixel variance (2D distribution) for each current
m_c_cu  = zeros(nCurr,1);                   % Values for central pixel, due to heel effect, it is supposed to be saturated at some point
v_c_cu  = zeros(nCurr,1);                   % Variances for central pixel, due to heel effect, it is supposed to be saturated at some point
% Loop

for nC = 1:nCurr,
   fname_al = strcat(folder,'Al/',fnameBase,'.flA',num2str(fix(cur_a(nC,1))));
   fname_cu = strcat(folder,'Cu/',fnameBase,'.flA',num2str(fix(cur_a(nC,2))));
   sequence_al   = readSimpleBinSliGap(fname_al,imSize(1),imSize(2),1,nFr,'uint16',4,4);
   sequence_cu   = readSimpleBinSliGap(fname_cu,imSize(1),imSize(2),1,nFr,'uint16',4,4);
   tmp_al = sequence_al(imXSize ,:,:);
   tmp_cu = sequence_cu(imXSize ,:,:);
    for nF = 1:2:nFr-1
        dif = tmp_al(floor(imXSize/2),floor(imSize(2)/2),nF) - tmp_al(floor(imXSize/2),floor(imSize(2)/2),nF+1);
        if(dif>1000)
            tmp_al(:,:,nF+1) = tmp_al(:,:,nF);
        elseif(dif<1000)
              tmp_al(:,:,nF) = tmp_al(:,:,nF+1);  
        end
        dif = tmp_cu(floor(imXSize/2),floor(imSize(2)/2),nF) - tmp_cu(floor(imXSize/2),floor(imSize(2)/2),nF+1);
        if(dif>1000)
            tmp_cu(:,:,nF+1) = tmp_cu(:,:,nF);
        elseif(dif<1000)
              tmp_cu(:,:,nF) = tmp_cu(:,:,nF+1);  
        end
    end
           
   tmp_al = interppx(tmp_al);
   tmp_cu = interppx(tmp_cu);
   % Let's do the procesing
   tmp_m_al = mean(tmp_al,3);           % Mean in 2D
   px_m_al(:,:,nC) = squeeze(tmp_m_al); % Store it in the proper place
   tmp_v_al = var(tmp_al,1,3);          % Variance in 2D
   px_v_al(:,:,nC) = squeeze(tmp_v_al); % Store it in the proper place
   % Get the values for the central pixel
   m_c_al(nC) = px_m_al(imSize(1)/2,imSize(2)/2,nC);
   v_c_al(nC) = px_v_al(imSize(1)/2,imSize(2)/2,nC);
   % Let's do the procesing
   tmp_m_cu = mean(tmp_cu,3);           % Mean in 2D
   px_m_cu(:,:,nC) = squeeze(tmp_m_cu); % Store it in the proper place
   tmp_v_cu = var(tmp_cu,1,3);          % Variance in 2D
   px_v_cu(:,:,nC) = squeeze(tmp_v_cu); % Store it in the proper place
   % Get the values for the central pixel
   m_c_cu(nC) = px_m_cu(imSize(1)/2,imSize(2)/2,nC);
   v_c_cu(nC) = px_v_cu(imSize(1)/2,imSize(2)/2,nC);
end
% Clear unnecessary vars


% Now let's do the mean to get a global understanding of the detector
% performance and where it starts to behave as quantum limited
total_m_al = squeeze(mean(mean(px_m_al,2),1));
total_v_al = squeeze(mean(mean(px_v_al,2),1));

total_m_al((cur_s > 0.0336)&(cur_s < 0.039)) = total_m_al((cur_s > 0.0336)&(cur_s < 0.039))*1.1;
total_v_al((cur_s > 0.0336)&(cur_s < 0.039)) = total_v_al((cur_s > 0.0336)&(cur_s < 0.039))*1.1;

total_m_cu = squeeze(mean(mean(px_m_cu,2),1));
total_v_cu = squeeze(mean(mean(px_v_cu,2),1));

total_m_cu((cur_s > 0.0336)&(cur_s < 0.039)) = total_m_cu((cur_s > 0.0336)&(cur_s < 0.039))*1.1;
total_v_cu((cur_s > 0.0336)&(cur_s < 0.039)) = total_v_cu((cur_s > 0.0336)&(cur_s < 0.039))*1.1;

% Plot the data
figure;
subplot(221);plot(cur_s(:,1),total_m_al,'k-*');title('Mean Pixel Value');xlabel('Current (mAs)');ylabel('Pixel value (ADU)');              % Mean pixel value
hold on
plot(cur_s(:,2),total_m_cu,'r-o');
subplot(222);plot(cur_s(:,1),total_v_al,'k-o');title('Mean Pixel Variance');xlabel('Current (mAs)');ylabel('Pixel variance (ADU)');        % Mean variance
hold on
plot(cur_s(:,2),total_v_cu,'r-*');
subplot(223);plot(total_m_al,total_v_al,'k-+');title('Mean Pixel Variance');xlabel('Pixel value (ADU)');ylabel('Pixel variance (ADU)');  % Mean variance as a function of the mean px value
hold on
plot(total_m_cu,total_v_cu,'r-s');
subplot(224);plot(total_m_al,total_m_al./sqrt(total_v_cu),'k-s');title('Mean SNR');xlabel('Pixel value (ADU)');ylabel('SNR');               % SNR
hold on
plot(total_m_cu,total_m_cu./sqrt(total_v_cu),'r-+');
