% Script to test the output provided by the detector as a function of
% spectrum.
%
% The main goal of this script is to test the consistency between
% simulated and real data and assess the validity of simulated results.

% Filtering material
% For the time being is aluminum of a given thickness
% More than one material can be used
% Al - 13; C - 6; Be - 4; W - 74 (Negative amount to compensate deviations
% from theoretical model)
% filtMat   = [13; 6; 4; 74]; % Atomic number of the materials
% filtThick = [1; 1; 0.5; -0.004];  % Thicknesses (mm)
filtMat   = [13; 6; 4]; % Atomic number of the materials
filtThick = [-0.22; 1; 0.5];  % Thicknesses (mm)

% Peak energy values to do the simu
kVps = [30 40 50 60 70];
siz  = size(kVps(:));

% Read detector response file
NvalSpec = 150;
detFile  = 'D:/documents/LIM/CT_nuevo/paper/test_protocols/spec_files/HamamatsuC7940_150um.txt';
% Load detector
fid1 = fopen(detFile,'r');
det  = fscanf(fid1,'%f',NvalSpec);
fclose(fid1);
%Test
%det  = detectorModel(56.37);

% Storage
detOut = zeros(siz(1),1);
for nEne = 1:siz(1),
    sp           = spektrSpectrum(kVps(nEne));
    spF          = spektrBeers(sp,[filtMat filtThick]);
    detOut(nEne) = sum(spF(:).*det(:));
end

% Currents to compare
A = [445 145 68 40 27];

% Done, plot it
figure;
subplot(121);
plot(kVps,detOut,'k-*');
xlabel('Energy [kVp]');
ylabel('Detector signal [AU]');
subplot(122);
plot(kVps,1./(detOut*(1/detOut(1))),'k-*');
hold on;
plot(kVps,(A*(1/A(1))),'r-o');
xlabel('Energy [kVp]');
ylabel('Detector signal ratio');