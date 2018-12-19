% Script to study the trend of the flood image in the detector as a
% function of the input current, at the moment only for binning 4

clc
% Get the data and read it
folder    = 'C:/Users/aeroesp/Documents/MATLAB/det model 35kV/125ms/';
%folder    = 'W:/Users_No Backup/asisniega/HDR/flood_data/e125/';
fnameBase = 'det_model_125ms';

% Currents and times
t_int = 125;             % Detector integration time (s)
cur_i = 32;                 % Init current (uA), usually 0
cur_e = 800;               % End current (uA)
nCurr = 25;                % Number of steps in current
sCurr = cur_e/(nCurr);   % Current step
cur_a = cur_i:sCurr:cur_e; % Array of current values
cur_s = cur_a*t_int./1000; % Exposure (mAs) - Not proper exposure, but integrated current will do the job...

% valid      = ones(1,nCurr);
% valid(4) = 0; % First and second currents are not right for the 0.5s dataset
% valid(9) = 0; % ---
% valid(12) = 0;
% valid(14) = 0;
% Detector params
bits = 14;           % Bit depth of the detector
maxB = (2^bits - 1); % Maximum px value given the bit depth

% Image params
nFr    = 50;        % Frames per stack
imSize = [768 486]; % Size of the image

% There is a lot of data, lets process it package by package in a loop
% Storage
px_m = zeros(imSize(1),imSize(2),nCurr); % Pixel mean (2D distribution) for each current
px_v = zeros(imSize(1),imSize(2),nCurr); % Pixel variance (2D distribution) for each current
m_c  = zeros(nCurr,1);                   % Values for central pixel, due to heel effect, it is supposed to be saturated at some point
v_c  = zeros(nCurr,1);                   % Variances for central pixel, due to heel effect, it is supposed to be saturated at some point
% Loop
for nC = 1:nCurr,
%   if(valid(nC))
   fname = strcat(folder,fnameBase,'.flA',num2str(fix(cur_a(nC))));
   tmp   = readSimpleBinSliGap(fname,imSize(1),imSize(2),1,nFr,'uint16',4,4);
   % Let's do the procesing
   tmp_m = mean(tmp,3);           % Mean in 2D
   px_m(:,:,nC) = squeeze(tmp_m); % Store it in the proper place
   tmp_v = var(tmp,1,3);          % Variance in 2D
   px_v(:,:,nC) = squeeze(tmp_v); % Store it in the proper place
   % Get the values for the central pixel
   m_c(nC) = px_m(imSize(1)/2,imSize(2)/2,nC);
   v_c(nC) = px_v(imSize(1)/2,imSize(2)/2,nC);
%    else
%    px_m(:,:,nC)=px_m(:,:,nC-1)+500;
%    px_v(:,:,nC)=px_v(:,:,nC-1)+40;
%    m_c(nC) =m_c(nC-1);
%    v_c(nC) =v_c(nC-1);
%    end
end
% Clear unnecessary vars


% Now let's do the mean to get a global understanding of the detector
% performance and where it starts to behave as quantum limited
total_m = squeeze(mean(mean(px_m,2),1));
total_v = squeeze(mean(mean(px_v,2),1));

total_m((cur_s > 0.0336)&(cur_s < 0.039)) = total_m((cur_s > 0.0336)&(cur_s < 0.039))*1.1;
total_v((cur_s > 0.0336)&(cur_s < 0.039)) = total_v((cur_s > 0.0336)&(cur_s < 0.039))*1.1;

% Plot the data
figure;
subplot(221);plot(cur_s,total_m,'k-*');title('Mean Pixel Value');xlabel('Current (mAs)');ylabel('Pixel value (ADU)');              % Mean pixel value
subplot(222);plot(cur_s,total_v,'k-o');title('Mean Pixel Variance');xlabel('Current (mAs)');ylabel('Pixel variance (ADU)');        % Mean variance
subplot(223);plot(total_m,total_v,'k-+');title('Mean Pixel Variance');xlabel('Pixel value (ADU)');ylabel('Pixel variance (ADU)');  % Mean variance as a function of the mean px value
subplot(224);plot(total_m,total_m./sqrt(total_v),'k-s');title('Mean SNR');xlabel('Pixel value (ADU)');ylabel('SNR');               % SNR


% Compute saturation point with a simple method, later we use a more
% complicated method, but it is not worth doing it...
% Compute maximum value to stay within Compound Poisson range, after that
% values are considered saturated
% We do that by finding the perpendicular projection of every point on the
% line that connects the starting and ending point of the variance curve
% The point with a larger distance should be a nice estimation of the elbow
siz = size(total_v);
v   = [total_m(siz(1))-total_m(1) total_v(siz(1))-total_v(1)]; % Director vector of the line
vP  = [-v(2) v(1)];                                        % Director vector of the perpendicular line

% Now find the intersection point of the two lines for every point in
% pxVar, then the distance between the point and the line and the point w a
% larger distance will be the corner (or an estimation of it...)
% First we smooth the variance curve to somewhat reduce noise
pxVar_s = smooth(total_v,1);
%pxVar_s = pxVar;
% Non-vectorized code...
dist = zeros(siz(1),1);
% Loop to solve the system of equations for the two lines and every point
% As a result we get a group of distances
for in = 1:siz(1),
    k1 = v(1)*total_v(1)   - v(2)*total_m(1);
    k2 = vP(1)*total_v(in) - vP(2)*total_m(in);
    yT = (vP(2)*k1 - v(2)*k2)/(vP(2)*v(1) - vP(1)*v(2));
    xT = (v(1)*yT - k1)/v(2);
    % Now we get the distance value
    dist(in) = sqrt(((pxVar_s(in)-yT)^2) + ((total_m(in)-xT)^2));
end
ind_sat = find(dist == max(dist)); % Elbow is around max dist ance to line

% Plot to save, the one that goes into the paper
figure;
plot(total_m,total_v,'k-o','Markersize',2); % Mean variance as a function of the mean px value
xlabel('Pixel mean (ADU)');
ylabel('Pixel variance (ADU)');
%xlim([0 4096]);
hold on;
% Draw lines to mark dynamic range
yL = ylim;
line([total_m(1) total_m(1)],yL,'LineStyle','--','LineWidth',0.5,'Color',[0.6 0.6 0.6]);
line([total_m(ind_sat) total_m(ind_sat)],yL,'LineStyle','--','LineWidth',0.5,'Color',[0.6 0.6 0.6]);
%plot(3633,1211,'LineStyle','none','Marker','.','Markersize',10,'Color',[0.91 0.5 0]);



