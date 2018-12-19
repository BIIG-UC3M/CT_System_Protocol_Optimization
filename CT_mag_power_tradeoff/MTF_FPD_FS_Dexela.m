% Script to compute MTF combining detector MTF and focal spot
% for the time being we get the maximum MTF 50% and 10% for a range of
% magnification factors...

% X-ray source
% The one simulated now is that requested by Ortega - 0.3 mm
%   
fs = 0.3;

% This is for the Dexela model
a_pix = 0.0748; % Pixel aperture in mm - we do not know fill factor
% Anyways we ended up copying the MTF provided by the manufacturer - but
% have only some points so far
fu       = [0 1    2    3    4    5    6    6.5];    % Frecuency points where we have valid data (lp/mm)
MTF_prev = [1 0.88 0.69 0.53 0.41 0.31 0.24 0.2];    % Values of MTF in the points provided just above

% Final frecuency array for MTF determination 
f = 0:0.01:8;

% Range of Mag values
M = 1:0.01:4;

% Resample detector MTF
MTF_fpd = interp1(fu,MTF_prev,f,'spline','extrap');

% Compute source focal spot MTF
siz = size(M);
for indM = 1:siz(2),
    MTFs_FS(indM,:) = exp((-pi)*( ((fs*((M(indM)-1)/M(indM)))^2) * (f.^2)));
end

% Now translate MTF of the detector to the object plane...
for indM = 1:siz(2),
    f_m = f/M(indM);
    MTFs_fpd(indM,:) = interp1(f,MTF_fpd,f_m,'spline');
end

% Now get the total MTF as a product of both individual MTFs
% MTF_s stands for MTF system
% Really this loop is not necessary, I just kept it to make the code
% clearer but it can and should be removed if needed
for indM = 1:siz(2),
    MTF_s(indM,:) = MTFs_FS(indM,:).*MTFs_fpd(indM,:);
end

% Now get the MTF 50%
% Same thing for the loop
for indM = 1:siz(2),
    tmp            = find(MTF_s(indM,:) <= 0.5);
    MTF_50(indM) = f(tmp(1));
end

% and MTF 10%
for indM = 1:siz(2),
    tmp          = find(MTF_s(indM,:) <= 0.2);
    MTF_20(indM) = f(tmp(1));
end

% Plot results
% 50%
figure;
plot(M,MTF_50,'k','LineWidth',2);
grid on;
xlabel('Magnification');
ylabel('Spatial Frecuency (lp/mm)');
title('MTF 50%');
% 20%
figure;
plot(M,MTF_20,'k','LineWidth',2);
grid on;
xlabel('Magnification');
ylabel('Spatial Frecuency (lp/mm)');
title('MTF 20%');