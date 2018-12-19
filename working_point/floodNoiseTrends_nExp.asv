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

% First of all, set graphic properties to get our current standard in
% graphics...
set(0,'DefaultFigureCreateFcn',get(0, 'FactoryFigureCreateFcn'));

% Get the data and read it
%folder    = 'D:/documents/LIM/CT/desarrollos_varios/HDR/version_2012/floodDynamicRange/';
%folder    = 'C:/Users/aeroesp/Documents/pruebas CTTB/wp_dso63/';
folder    =  'D:dr22e50/';
%folder    = 'W:/Users_No Backup/asisniega/HDR/flood_data/e125/';
%fnameBase = 'Pwp_dso62_';
fnameBase = 'DR22_fast22';
% For the conversion to exposure
mR_x_mAs  = 24.2; % mR/mAs


% Currents and times
nTime = 3;                   % Number of int times
t_int = [0.05 0.05 0.05];       % Detector integration times (s)
%cur_i = [32 8 6];             % Init current (uA), usually 0
%cur_e = [800 200 150];        % End current (uA)
cur_i = [22 22  22];             % Init current (uA), usually 0
cur_e = [550 550 550];        % End current (uA)
nCurr = 25;                  % Number of steps in current
%nCurr = 26;                  % Number of steps in current
sCurr = cur_e*(1/(nCurr)); % Current step
cur_a = zeros(nCurr,nTime);  % Define the matrix for anode current values
cur_s = zeros(nCurr,nTime);  % Exposure (mAs) - Not proper exposure, but integrated current will do the job...
% Loop to fill out values
for nT = 1:nTime,
    cur_a(:,nT) = (cur_i(nT):sCurr(nT):cur_e(nT))'; % Array of anode current values for the given integration time
    cur_s(:,nT) = cur_a(:,nT)*(t_int(nT)/1000);  % The 1000 factor stands for the conversion of uA to mA
end
% Convert to exposure
 exp_s = cur_s*mR_x_mAs;
%exp_s = cur_s; % Exposure (mAs) - Not proper exposure, but integrated current will do the job...
% Some of the values are not valid cause the current of the tube could not
% be set so low, we have to not use these values cause they are not valid,
% the source was set to full current and the detector saturated...
valid      = ones(nCurr,nT);
 valid(14,2) = 0; % First and second currents are not right for the 0.5s dataset
% valid(3,2) = 0; % ---
% valid(2,3) = 0; % First, second, third and fourth currents are not right for the 1s dataset
% valid(3,3) = 0; % ---
% valid(4,3) = 0; % ---
% valid(5,3) = 0; % ---

% get away with one out of two 
%valid(2:2:length(valid(:,1)),1) = 0;
%valid(2:2:length(valid(:,2)),2) = 0;
%valid(2:2:length(valid(:,3)),3) = 0;

% Detector params
bits = 14;           % Bit depth of the detector
maxB = (2^bits - 1); % Maximum px value given the bit depth

% Image params
nFr    = 25;        % Frames per stack
imSize = [768 486]; % Size of the image
%imXSize = 81:728;   % Corrected X size of the image, to avoide the cone beam.
% There is a lot of data, lets process it package by package in a loop
% Storage
px_m = zeros(imSize(1),imSize(2),nCurr,nTime); % Pixel mean (2D distribution) for each current
px_v = zeros(imSize(1),imSize(2),nCurr,nTime); % Pixel variance (2D distribution) for each current
m_c  = zeros(nCurr,nTime);                     % Values for central pixel, due to heel effect, it is supposed to be saturated at some point
v_c  = zeros(nCurr,nTime);                     % Variances for central pixel, due to heel effect, it is supposed to be saturated at some point
dark_aux = zeros(imSize(1),imSize(2));
dark = zeros(imSize(1),imSize(2),nFr);

% Loop
for nT = 1:nTime,
    % read dark data
    fname = strcat(folder,num2str(fix(t_int(nT)*1000)),'ms/Cu/',fnameBase,num2str(fix(t_int(nT)*1000)),'ms.dk')
    dark_aux(:,:)   = readSimpleBinSliGap(fname,imSize(1),imSize(2),1,1,'uint16',4,4);
    dark  = repmat(squeeze(dark_aux(:,:)),[1 1 nFr]);
    for nC = 1:nCurr,
       fname = strcat(folder,num2str(fix(t_int(nT)*1000)),'ms/Cu/',fnameBase,num2str(fix(t_int(nT)*1000)),'ms.flA',num2str(fix(cur_a(nC,nT))))
       tmp   = readSimpleBinSliGap(fname,imSize(1),imSize(2),1,nFr,'uint16',4,4);
       for nF = 1:2:nFr-1
        dif = tmp(floor(imSize(1)/2),floor(imSize(2)/2),nF) - tmp(floor(imSize(1)/2),floor(imSize(2)/2),nF+1);
        if(dif>1000)
            tmp(:,:,nF+1) = tmp(:,:,nF);
        elseif(dif<1000)
              tmp(:,:,nF) = tmp(:,:,nF+1);  
        end
    end
       %readSimpleBinSliGap(fname,imSize(1),imSize(2),1,nFr,'uint16',4,4)
       tmp = tmp-dark;
       tmp = interppx(tmp);
       % Let's do the procesing
       tmp_m = mean(tmp,3);           % Mean in 2D
       px_m(:,:,nC,nT) = squeeze(tmp_m); % Store it in the proper place
       tmp_v = var(tmp,1,3);          % Variance in 2D
       px_v(:,:,nC,nT) = squeeze(tmp_v); % Store it in the proper place
       % Get the values for the central pixel
       m_c(nC,nT) = px_m(imSize(1)/2,imSize(2)/2,nC,nT);
       v_c(nC,nT) = px_v(imSize(1)/2,imSize(2)/2,nC,nT);
    end
end
% Clear unnecessary vars


% Now let's do the mean to get a global understanding of the detector
% performance and where it starts to behave as quantum limited
total_m = squeeze(mean(mean(px_m,2),1));
total_v = squeeze(mean(mean(px_v,2),1));

% Compute the DR for the different integration times

% Plot the results and save the plots in the folder for images
% First plot the detector mean signal and noise as a function of dose
figure;
plot(exp_s(valid(:,1)>0.5,1),total_m(valid(:,1)>0.5,1),'-o','Color',[0 0 0],'MarkerSize',2);
hold on;
plot(exp_s(valid(:,2)>0.5,2),total_m(valid(:,2)>0.5,2),'-*','Color',[0 1 0],'MarkerSize',2);
plot(exp_s(valid(:,3)>0.5,3),total_m(valid(:,3)>0.5,3),'-s','Color',[0 0 1],'MarkerSize',2);
ax1 = gca;
hold off;
box off;
%ylim([0 4500]);
%xlim([0 0.12]);
ylabel('Mean pixel value (ADU)');
xlabel('Exposure (mAs)');
%xlabel('Exposure (mR)');
% Legends
lh = legend('T_{125}','T_{500}','T_{1000}','Location','NorthWest');
set(lh,'FontSize',6);
children = get(lh, 'Children') ;
% Resize Marker
set(children(1),'MarkerSize',4);
set(children(4),'MarkerSize',4);
set(children(7),'MarkerSize',4);
% Resize Line
XData = get(children(2), 'XData') ;
XScale = XData(2) - XData(1) ;
XData(1) = XData(1) + XScale / 4 ;
XData(2) = XData(1) + XScale / 2 ; % half location
set(children(2), 'XData', XData) ;
XData = get(children(5), 'XData') ;
XScale = XData(2) - XData(1) ;
XData(1) = XData(1) + XScale / 4 ;
XData(2) = XData(1) + XScale / 2 ; % half location
set(children(5), 'XData', XData) ;
XData = get(children(8), 'XData') ;
XScale = XData(2) - XData(1) ;
XData(1) = XData(1) + XScale / 4 ;
XData(2) = XData(1) + XScale / 2 ; % half location
set(children(8), 'XData', XData) ;
% Reduce the box
pos = get(lh,'Position');
pos(1) = pos(1) + 0.2*pos(3);
pos(3) = 0.6*pos(3);
set(lh,'Position',pos);
n='NOISE'
% Now the noise
% First create the axis
ax2 = axes('Units',get(ax1,'Units'), ...
           'Position',get(ax1,'Position'), ...
           'Parent',get(ax1,'Parent'));
% Plot the data
plot(exp_s(valid(:,1)>0.5,1),total_v(valid(:,1)>0.5,1),'--o','Color',[0 0 0],'MarkerSize',2);
hold on;
plot(exp_s(valid(:,2)>0.5,2),total_v(valid(:,2)>0.5,2),'--*','Color',[0 1 0],'MarkerSize',2);
plot(exp_s(valid(:,3)>0.5,3),total_v(valid(:,3)>0.5,3),'--s','Color',[0 0 1],'MarkerSize',2);
% Modify properties which need to be modified
set(ax2,'YAxisLocation','right','XAxisLocation','top', ...
    'Color','none', ...
    'XGrid','off','YGrid','off','Box','off', ...
    'HitTest','off');
%ylim([0 2000]);
%xlim([0 0.12]);
set(ax2,'XTickLabel',[]);
ylabel('Mean pixel variance (ADU ^2)');

% Handle resizing and this kind of stuff - same way it is made in plotyy
ax = [ax1 ax2]; % Store handles together
% From "plotyy"
% Link the "Position" and "View" properties of the axes (for plot
% manipulation purposes)
hLink = linkprop(ax,'View');
% Store a handle to the link object returned by linkprop to keep
% the listener in scope. Used for rotate3d and plot edit mode.
setappdata(ax(1),'graphicsPlotyyLinkProp',hLink);
% The position is a bit more complex. We will use a listener for that
% property:
hList(1) = handle.listener(handle(ax(1)),findprop(handle(ax(1)),'Position'),...
    'PropertyPostSet',{@globalUpdatePosition_as,ax(1),ax(2)});
hList(2) = handle.listener(handle(ax(2)),findprop(handle(ax(2)),'Position'),...
    'PropertyPostSet',{@globalUpdatePosition_as,ax(2),ax(1)});
setappdata(ax(1),'graphicsPlotyyPositionListener',hList);

% Print the figure
% out_folder = 'C:/Users/aeroesp/Documents/MATLAB/imgs';
% width  = 3;
% height = 2.5;
% % Set printing stuff
% set(gcf,'PaperPositionMode','manual');
% set(gcf,'PaperPosition',[0 0 width height]);
% set(gcf,'PaperSize',[width height]);
% %set(ax2,'OuterPosition',[0 0 1 1]);
% %set(gca,'OuterPosition',[0.15 0.1 0.7 0.9]);
% set(gca,'OuterPosition',[0.05 0.05 0.85 0.9]);
% print(gcf,'-dtiff','-r300',sprintf('%s/meanVar_vs_Exp.tif',out_folder));
% keyboard
n ='Saturation point'
% Compute saturation point with a simple method, later we use a more
% complicated method, but it is not worth doing it...
% Compute maximum value to stay within Compound Poisson range, after that
% values are considered saturated
% We do that by finding the perpendicular projection of every point on the
% line that connects the starting and ending point of the variance curve
% The point with a larger distance should be a nice estimation of the elbow
% Find the saturation index for every integration time
ind_sat = zeros(3,1);
% We do it in a loop
for nT = 1:nTime,
    % Line params
    siz = size(total_v(:,nT));
    v   = [total_m(siz(1),nT)-total_m(1,nT) total_v(siz(1),nT)-total_v(1,nT)]; % Director vector of the line
    vP  = [-v(2) v(1)];                                        % Director vector of the perpendicular line
    % Now find the intersection point of the two lines for every point in
    % pxVar, then the distance between the point and the line and the point w a
    % larger distance will be the corner (or an estimation of it...)
    % First we smooth the variance curve to somewhat reduce noise
    pxVar_s = smooth(total_v(:,nT),2);
    %pxVar_s = pxVar;
    % Non-vectorized code...
    dist = zeros(siz(1),1);
    % Loop to solve the system of equations for the two lines and every point
    % As a result we get a group of distances
    for in = 1:siz(1),
        k1 = v(1)*total_v(1,nT)   - v(2)*total_m(1,nT);
        k2 = vP(1)*total_v(in,nT) - vP(2)*total_m(in,nT);
        yT = (vP(2)*k1 - v(2)*k2)/(vP(2)*v(1) - vP(1)*v(2));
        xT = (v(1)*yT - k1)/v(2);
        % Now we get the distance value
        dist(in) = sqrt(((pxVar_s(in)-yT)^2) + ((total_m(in,nT)-xT)^2));
    end
    ind_sat(nT) = find(dist == max(dist)); % Elbow is around max distance to line

end

% Lets get the params for the linear model of pixel variance
% Additive noise
n    = zeros(3,1);
n(:) = squeeze(total_v(1,:));
% Poisson component, slope of the line - Use just the first and last points
% We just want the first one for the plot, we can adjust it later...
h    = zeros(3,1);
h(1) = (total_v(ind_sat(1),1) - total_v(1,1))/(total_m(ind_sat(1),1) - total_m(1,1));
h(2) = (total_v(ind_sat(2),2) - total_v(1,2))/(total_m(ind_sat(2),2) - total_m(1,2));
h(3) = (total_v(ind_sat(3),3) - total_v(1,3))/(total_m(ind_sat(3),3) - total_m(1,3));

% Now we have to compute the point where the SNR is equal to one
% The equation is:
%  gX = (h + sqrt(h^2 + 4*n))/2
%
% Do the calculation just for the widest range - Provide maximum DR, i.e.
% the one for the lower integration time, that have the lower P.
DR = zeros(3,1);
gX_ini = (h(1) + sqrt(h(1)^2 + 4*n(1)))/2;
gX_end = total_m(ind_sat(1),1) - total_m(1,1);
DR(1)  = gX_end - gX_ini;
gX_ini = (h(2) + sqrt(h(2)^2 + 4*n(2)))/2;
gX_end = total_m(ind_sat(2),2) - total_m(1,2);
DR(2)  = gX_end - gX_ini;
gX_ini = (h(3) + sqrt(h(3)^2 + 4*n(3)))/2;
gX_end = total_m(ind_sat(3),3) - total_m(1,3);
DR(3)  = gX_end - gX_ini;


% Plot to save, the one that goes into the paper
figure;
title('Dinamic Range');
plot(total_m(valid(:,1)>0.5,1)-total_m(1,1),total_v(valid(:,1)>0.5,1),'-o','Markersize',2,'Color',[0 0 0]); % Mean variance as a function of the mean px value
hold on;
plot(total_m(valid(:,2)>0.5,2)-total_m(1,2),total_v(valid(:,2)>0.5,2),'-*','Markersize',2,'Color',[0.33 1 0.33]); % Mean variance as a function of the mean px value
plot(total_m(valid(:,3)>0.5,3)-total_m(1,3),total_v(valid(:,3)>0.5,3),'-s','Markersize',2,'Color',[0.66 0.66 1]); % Mean variance as a function of the mean px value
xlabel('gX (ADU)');
ylabel('v(gX) (ADU ^2)');
%xlim([0 4096]);
%ylim([0 3000]);
hold on;
% Draw lines to mark dynamic range
yL = ylim;
line([gX_ini gX_ini],yL,'LineStyle','--','LineWidth',0.5,'Color',[0.5 0.5 0.5]);
line([total_m(ind_sat(1),1)-total_m(1,1) total_m(ind_sat(1),1)-total_m(1,1)],yL,'LineStyle','--','LineWidth',0.5,'Color',[0.5 0.5 0.5]);
% Legends
lh = legend('T_{125}','T_{500}','T_{1000}','Location','NorthWest');
set(lh,'FontSize',6);
children = get(lh, 'Children') ;
% Resize Marker
set(children(1),'MarkerSize',4);
set(children(4),'MarkerSize',4);
set(children(7),'MarkerSize',4);
% % Save to image
% out_folder = 'D:/documents/LIM/CT/desarrollos_varios/HDR/version_2012/paperData/imgs/fig2_mean_vs_gX';
% width  = 3;
% height = 2.5;
% % Set printing stuff
% set(gcf,'PaperPositionMode','manual');
% set(gcf,'PaperPosition',[0 0 width height]);
% set(gcf,'PaperSize',[width height]);
% set(gca,'OuterPosition',[0 0 1 1]);
% print(gcf,'-dtiff','-r300',sprintf('%s/var_vs_gX.tif',out_folder));


%%%%%%%%%%%%%
% Check the behaviour of independent pixels to not assume the detector is
% uniformly radiated and take into account Heel's effect

% Resample invidual px curves to 4096 levels (1 per pixel mean) - let's do
% a loop - easier ...
vec      = 0:maxB;
px_v_int = zeros(maxB+1,nTime);
norm_f   = zeros(maxB+1,nT);
px_v_val = ones(maxB+1,nT);
size_ROI = [128 128];
for nT = 1:nTime,
    %norm_f   = zeros(maxB+1,1);
    for in_x = imSize(1)/2-(size_ROI(1)/2):imSize(1)/2+(size_ROI(1)/2)-1,
        for in_y = imSize(2)/2-(size_ROI(2)/2):imSize(2)/2+(size_ROI(2)/2)-1,
            tmp_m = px_m(in_x,in_y,:,nT);
            tmp_v = px_v(in_x,in_y,:,nT);
            [tmp_m iA iC] = unique(tmp_m);
            tmp_v = tmp_v(iA);
            %tmp   = interp1(tmp_m(:),tmp_v(:),vec(:),'linear',0);
            %px_v_int(iA,nT)  = px_v_int(iA,nT) + tmp(:);
            siz_m = size(squeeze(tmp_m));
            for index_m = 1:siz_m,
                px_v_int(fix(tmp_m(index_m)),nT)  = px_v_int(fix(tmp_m(index_m)),nT) + tmp_v(index_m);
                norm_f(fix(tmp_m(index_m)),nT) = norm_f(fix(tmp_m(index_m)),nT) + 1;
            end
            
            %norm_f(iA) = norm_f(iA) + 1;
        end
    end
    % Do the mean
    px_v_int(:,nT) = px_v_int(:,nT)./norm_f(:,nT);
    
    % Do interpolation in the place where's there's not enough values
end



% Find useful indx and gX values for integration times
inds  = zeros(maxB+1,3,'int8');
gX_in = zeros(maxB+1,3);
gX_in(:,1) = vec(:) - total_m(1,1);
inds(:,1)  = gX_in(:,1) > 0;
gX_in(:,2) = vec(:) - total_m(1,2);
inds(:,2)  = gX_in(:,2) > 0;
gX_in(:,3) = vec(:) - total_m(1,3);
inds(:,3)  = gX_in(:,3) > 0;

% Or just find wich values have to and have not to be plot
inds(norm_f < 100) = 0;
% Now let's compute the modelled behaviour and see how it compares to the
% measured signal
sim = zeros(maxB+1,1);
for nT = 1:nTime,
    p_i    = fix(total_m(1,nT));
    maxB_i = maxB - p_i;
    for gX_i = p_i:maxB,
        sim(gX_i,nT) = varClipped(gX_i - p_i,maxB_i,h(nT),n(nT));
    end
end


% Do a simulation and see what happens
nRea = 10000;
sim_MC = zeros(129,3);
gX_sim = 0:4096/128:4096;
for nT = 1:nTime,
    p_i    = fix(total_m(1,nT));
    maxB_i = maxB - p_i;
    for in_gX = 1:129,
        mu = gX_sim(in_gX) - p_i;
        if (mu < 0),
            continue;
        end
        sigma         = sqrt(h(nT)*mu + n(nT));
        sig_s         = min(normrnd(mu,sigma,nRea,1),maxB_i);
        sim_MC(in_gX,nT) = var(sig_s,1);
    end
end
    

% Plot the result
figure;
plot(gX_in(logical(inds(:,1)),1),px_v_int(logical(inds(:,1)),1),'o','LineStyle','none','Color',[0 0 0],'MarkerSize',2); hold on;
plot(gX_in(logical(inds(:,2)),2),px_v_int(logical(inds(:,2)),2),'*','LineStyle','none','Color',[0.33 0.33 0.33],'MarkerSize',2);
plot(gX_in(logical(inds(:,3)),3),px_v_int(logical(inds(:,3)),3),'s','LineStyle','none','Color',[0.66 0.66 0.66],'MarkerSize',2);
% Simulated Data
plot(gX_in(logical(inds(:,1)),1),sim(logical(inds(:,1)),1),'r-');
plot(gX_in(logical(inds(:,2)),2),sim(logical(inds(:,2)),2),'g-');
plot(gX_in(logical(inds(:,3)),3),sim(logical(inds(:,3)),3),'y-');
% Monte Carlo
plot(gX_sim((gX_sim - total_m(1,1))>0)-total_m(1,1),sim_MC((gX_sim - total_m(1,1))>0,1),'r*');
plot(gX_sim((gX_sim - total_m(1,2))>0)-total_m(1,2),sim_MC((gX_sim - total_m(1,2))>0,2),'g*');
plot(gX_sim((gX_sim - total_m(1,3))>0)-total_m(1,3),sim_MC((gX_sim - total_m(1,3))>0,3),'y*');


% Resample the results to do an easier plot
gX_dest = cell(nTime,1);
px_v_int_resampled = cell(nTime,1);
for nT = 1:nTime,
    gX_orig     = gX_in(logical(inds(:,nT)),nT);
    gX_init     = gX_orig(1);
    gX_end      = maxB-1;
    gX_dest{nT} = gX_init:20:gX_end;
    % Resample the vector
    px_v_int_resampled{nT} = interp1(gX_orig,px_v_int(logical(inds(:,nT)),nT),gX_dest{nT}(:));
end


% Compute saturation point with a simple method, later we use a more
% complicated method, but it is not worth doing it...
% Compute maximum value to stay within Compound Poisson range, after that
% values are considered saturated
% We do that by finding the perpendicular projection of every point on the
% line that connects the starting and ending point of the variance curve
% The point with a larger distance should be a nice estimation of the elbow
% Find the saturation index for every integration time
ind_sat_indpx = zeros(3,1);
% We do it in a loop
for nT = 1:nTime,
    % Line params
    siz  = size(sim(:,nT));
    p_i  = fix(total_m(1,nT));
    gX_s = 0:maxB;
    v    = [gX_s(siz(1))-gX_s(p_i) sim(siz(1),nT)-sim(p_i,nT)]; % Director vector of the line
    vP   = [-v(2) v(1)];                                        % Director vector of the perpendicular line
    % Now find the intersection point of the two lines for every point in
    % pxVar, then the distance between the point and the line and the point w a
    % larger distance will be the corner (or an estimation of it...)
    % First we smooth the variance curve to somewhat reduce noise
    pxVar_s = smooth(sim(:,nT),2);
    %pxVar_s = pxVar;
    % Non-vectorized code...
    dist = zeros(siz(1),1);
    % Loop to solve the system of equations for the two lines and every point
    % As a result we get a group of distances
    for in = 1:siz(1),
        k1 = v(1)*sim(p_i,nT)   - v(2)*gX_s(p_i);
        k2 = vP(1)*sim(in,nT)   - vP(2)*gX_s(in);
        yT = (vP(2)*k1 - v(2)*k2)/(vP(2)*v(1) - vP(1)*v(2));
        xT = (v(1)*yT - k1)/v(2);
        % Now we get the distance value
        dist(in) = sqrt(((pxVar_s(in)-yT)^2) + ((gX_s(in)-xT)^2));
    end
    ind_sat_indpx(nT) = find(dist == max(dist)); % Elbow is around max distance to line
end

% Compute the "h" value using th individual pixel data and see if it agrees
% with previous values
% Poisson component, slope of the line - Use just the first and last points
% We just want the first one for the plot, we can adjust it later...
h_sim    = zeros(3,1);
ind1 = find(sim(:,1) > 0);
h_sim(1) = (sim(ind_sat_indpx(1),1) - sim(ind1(1),1))/(gX_in(ind_sat_indpx(1),1));% - total_m(1,1));
ind1 = find(sim(:,2) > 0);
h_sim(2) = (sim(ind_sat_indpx(2),2) - sim(ind1(1),2))/(gX_in(ind_sat_indpx(2),2));% - total_m(1,2));
ind1 = find(sim(:,3) > 0);
h_sim(3) = (sim(ind_sat_indpx(3),3) - sim(ind1(1),3))/(gX_in(ind_sat_indpx(3),3));% - total_m(1,3));

% Now the DR
DR_ind    = zeros(3,1);
DR_ind(1) = (ind_sat_indpx(1) - total_m(1,1)) - gX_ini;
DR_ind(2) = (ind_sat_indpx(2) - total_m(1,2)) - gX_ini;
DR_ind(3) = (ind_sat_indpx(3) - total_m(1,3)) - gX_ini;

% Plot the result & save it
figure;
plot(gX_dest{1},px_v_int_resampled{1},'o','LineStyle','none','Color',[0 0 0],'MarkerSize',1.5); hold on;
plot(gX_dest{2},px_v_int_resampled{2},'^','LineStyle','none','Color',[0.33 0.33 0.33],'MarkerSize',1.5);
plot(gX_dest{3},px_v_int_resampled{3},'+','LineStyle','none','Color',[0.66 0.66 0.66],'MarkerSize',1.5);
% Legends
lh = legend('T_{125}','T_{500}','T_{1000}','Location','NorthWest');
set(lh,'FontSize',6);
children = get(lh, 'Children') ;
% Resize Marker
set(children(1),'MarkerSize',4);
set(children(4),'MarkerSize',4);
set(children(7),'MarkerSize',4);
% Simulated Data
plot(gX_in(logical(inds(:,1)),1),sim(logical(inds(:,1)),1),'-','Color',[0 0 0],'LineWidth',0.7);
plot(gX_in(logical(inds(:,2)),2),sim(logical(inds(:,2)),2),'-','Color',[0.33 0.33 0.33],'LineWidth',0.7);
plot(gX_in(logical(inds(:,3)),3),sim(logical(inds(:,3)),3),'-','Color',[0.66 0.66 0.66],'LineWidth',0.7);
% Set lims
%xlim ([0 4096]);
%ylim ([0 1600]);
xlabel('gX (ADU)');
ylabel('v(gX) (ADU^2)');
% Draw lines to mark dynamic range
yL = ylim;
line([gX_ini gX_ini],yL,'LineStyle','--','LineWidth',0.5,'Color',[0.5 0.5 0.5]);
line([ind_sat_indpx(1)-total_m(1,1) ind_sat_indpx(1)-total_m(1,1)],yL,'LineStyle','--','LineWidth',0.5,'Color',[0.5 0.5 0.5]);
% Save to image
% out_folder = 'C:/Users/aeroesp/Documents/MATLAB/img/fig2_mean_vs_gX';
% width  = 4;
% height = 3;
% % Set printing stuff
% set(gcf,'PaperPositionMode','manual');
% set(gcf,'PaperPosition',[0 0 width height]);
% set(gcf,'PaperSize',[width height]);
% set(gca,'OuterPosition',[0 0 1 1]);
% print(gcf,'-dtiff','-r1200',sprintf('%s/indvar_vs_gX.tif',out_folder));
% 
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
hFit = zeros(imSize(1),imSize(2)); % Slope of the gain curve, considering it linear
nFit = zeros(imSize(1),imSize(2)); % Intercept of the fit (i.e. additive noise, calculated in a different way, results should be equal)
for inX = 1:imSize(1),
    for inY = 1:imSize(2),
        tmp = polyfit(px_m(inX,inY,1:ind_sat),px_v(inX,inY,1:ind_sat),1); % Fit to 1-deg polynomial
        hFit(inX,inY) = tmp(1); % Descending order, first the slope
        nFit(inX,inY) = tmp(2); % Now the intercept
    end
end
h = (squeeze(px_v(:,:,end)) - n)./(squeeze(px_m(:,:,end)) - p); % Slope, but now we only use two extreme points, not a fit...
