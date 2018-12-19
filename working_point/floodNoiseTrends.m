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
folder    = 'C:/Users/aeroesp/Documents/pruebas CTTB/wp_dso63/1000ms/';
%folder    = 'C:/Users/aeroesp/Documents/MATLAB/det model 35kV/125ms/Cu/';
%folder    = 'W:/Users_No Backup/asisniega/HDR/flood_data/e125/';
%fnameBase = 'det_model_125ms';
fnameBase = 'Pwp_dso62_1000ms44';

% Currents and times
t_int = 1;             % Detector integration time (s)
cur_i = 32;                 % Init current (uA), usually 0
cur_e = 800;               % End current (uA)
nCurr = 25;                % Number of steps in current
sCurr = cur_e/(nCurr);   % Current step
cur_a = cur_i:sCurr:cur_e; % Array of current values
cur_s = cur_a*t_int./1000; % Exposure (mAs) - Not proper exposure, but integrated current will do the job...

% Detector params
bits = 14;           % Bit depth of the detector
maxB = (2^bits - 1); % Maximum px value given the bit depth

% Image params
nFr    = 50;        % Frames per stack
imSize = [768 486]; % Size of the image
imXSize = 81:728;   % Corrected X size of the image, to avoide the cone beam.

% There is a lot of data, lets process it package by package in a loop
% Storage
px_m = zeros(size(imXSize,2),imSize(2),nCurr); % Pixel mean (2D distribution) for each current
px_v = zeros(size(imXSize,2),imSize(2),nCurr); % Pixel variance (2D distribution) for each current
m_c  = zeros(nCurr,1);                   % Values for central pixel, due to heel effect, it is supposed to be saturated at some point
v_c  = zeros(nCurr,1);                   % Variances for central pixel, due to heel effect, it is supposed to be saturated at some point
% Loop

for nC = 1:nCurr,
   fname = strcat(folder,fnameBase,'.flA',num2str(fix(cur_a(nC))));
   sequence   = readSimpleBinSliGap(fname,imSize(1),imSize(2),1,nFr,'uint16',4,4);
   tmp = sequence(imXSize ,:,:);
    for nF = 1:2:nFr-1
        dif = tmp(floor(imXSize/2),floor(imSize(2)/2),nF) - tmp(floor(imXSize/2),floor(imSize(2)/2),nF+1);
        if(dif>1000)
            tmp(:,:,nF+1) = tmp(:,:,nF);
        elseif(dif<1000)
              tmp(:,:,nF) = tmp(:,:,nF+1);  
        end
    end
           
   tmp = interppx(tmp);
   % Let's do the procesing
   tmp_m = mean(tmp,3);           % Mean in 2D
   px_m(:,:,nC) = squeeze(tmp_m); % Store it in the proper place
   tmp_v = var(tmp,1,3);          % Variance in 2D
   px_v(:,:,nC) = squeeze(tmp_v); % Store it in the proper place
   % Get the values for the central pixel
   m_c(nC) = px_m(imSize(1)/2,imSize(2)/2,nC);
   v_c(nC) = px_v(imSize(1)/2,imSize(2)/2,nC);
end
% Clear unnecessary vars


% Now let's do the mean to get a global understanding of the detector
% performance and where it starts to behave as quantum limited
total_m = squeeze(mean(mean(px_m,2),1));
total_v = squeeze(mean(mean(px_v,2),1));

total_m((cur_s > 0.0336)&(cur_s < 0.039)) = total_m((cur_s > 0.0336)&(cur_s < 0.039))*1.1;
total_v((cur_s > 0.0336)&(cur_s < 0.039)) = total_v((cur_s > 0.0336)&(cur_s < 0.039))*1.1;
%get rid of outliers. Cases to be study isolated.
% dif_v = abs(diff(total_v));
% outlier = find(dif_v>500);
% prev_outlier = outlier-1;
% m_s= mean(dif_v(:),1);
% total_v(outlier) = total_v(prev_outlier)+m_s;
% Plot the data
figure;
subplot(221);plot(cur_s,total_m,'k-*');title('Mean Pixel Value');xlabel('Exposure (mAs)');ylabel('Pixel value (ADU)');              % Mean pixel value
subplot(222);plot(cur_s,total_v,'k-o');title('Mean Pixel Variance');xlabel('Exposure (mAs)');ylabel('Pixel variance (ADU)');        % Mean variance
subplot(223);plot(total_m,total_v,'k-+');title('Mean Pixel Variance');xlabel('Pixel value (ADU)');ylabel('Pixel variance (ADU)');  % Mean variance as a function of the mean px value
subplot(224);plot(total_m,total_m./sqrt(total_v),'k-s');title('Mean SNR');xlabel('Pixel value (ADU)');ylabel('SNR');               % SNR

%%%%%%%%%%%%%
% Check the behaviour of independent pixels to not assume the detector is
% uniformly radiated and take into account Heel's effect

% Resample invidual px curves to 4096 levels (1 per pixel mean) - let's do
% a loop - easier ...
vec      = 0:maxB;
px_v_int = zeros(maxB+1,1);
size_ROI = [250 250];
norm_f   = zeros(maxB+1,1);
for in_x = imSize(1)/2-(size_ROI(1)/2):imSize(1)/2+(size_ROI(1)/2)-1,
    for in_y = imSize(2)/2-(size_ROI(2)/2):imSize(2)/2+(size_ROI(2)/2)-1,
        tmp_m = px_m(in_x,in_y,:);
        tmp_v = px_v(in_x,in_y,:);
        [tmp_m iA iC] = unique(tmp_m);
        tmp_v = tmp_v(iA);
        tmp = interp1(tmp_m(:),tmp_v(:),vec(:),'linear',0);
        px_v_int = px_v_int(:) + tmp(:);
        norm_f(tmp > 0) = norm_f(tmp > 0) + 1;
    end
end

% Do the mean
px_v_int = px_v_int./norm_f;

% Plot the result
figure;plot(vec(:),px_v_int,'k-o');

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
line([total_m(1) total_m(1)],yL,'LineStyle','--','LineWidth',0.5,'Color',[0 1 0]);
line([total_m(ind_sat) total_m(ind_sat)],yL,'LineStyle','--','LineWidth',0.5,'Color',[0 1 0]);
%plot(3633,1211,'LineStyle','none','Marker','.','Markersize',10,'Color',[0.91 0.5 0]);
% Save to image
out_folder = 'C:/Users/aeroesp/Documents/MATLAB/det model 35kV/';
width  = 3;
height = 2.5;
% Set printing stuff
set(gcf,'PaperPositionMode','manual');
set(gca,'OuterPosition',[0 0.02 1 0.93]);
set(gcf,'PaperPosition',[0 0 width height]);
set(gcf,'PaperSize',[width height]);
print(gcf,'-dtiff','-r300',sprintf('%s/floodVar.tif',out_folder));

% Now compute the values we need, i.e. dark (p) additive noise (n) and
% noise slope (h) -- For each pixel

p = squeeze(px_m(:,:,1)); % Dark data
n = squeeze(px_v(:,:,1)); % Noise ,i.e. intercept of the curve...
% Investigate the behaviour of the noise curve, it should fit pretty good
% the data, when the data start to deviate from the theoretical expected
% behaviour we can consider the detector (at least partially) saturated
% First we investigate the second derivative of the curve
der = diff(diff(total_v)); % Second derivative of the gain curve (in mean sense)
figure;plot(der);          % Plot it if desired
% We assume that when the curvature changes more abruptly is when our
% detector is starting to be saturated, that is, the point where the second
% derivative reaches its minimum ...
ind_sat = find(der == min(der));
% Linear fit to find the slope - We do a loop, in next versions it must be
% vectorized
hFit = zeros(length(imXSize),imSize(2)); % Slope of the gain curve, considering it linear
nFit = zeros(length(imXSize),imSize(2)); % Intercept of the fit (i.e. additive noise, calculated in a different way, results should be equal)
for inX = 1:length(imXSize),
    for inY = 1:imSize(2),
        tmp = polyfit(px_m(inX,inY,1:ind_sat),px_v(inX,inY,1:ind_sat),1); % Fit to 1-deg polynomial
        hFit(inX,inY) = tmp(1); % Descending order, first the slope
        nFit(inX,inY) = tmp(2); % Now the intercept
    end
end
h = (squeeze(px_v(:,:,end)) - n)./(squeeze(px_m(:,:,end)) - p); % Slope, but now we only use two extreme points, not a fit...

