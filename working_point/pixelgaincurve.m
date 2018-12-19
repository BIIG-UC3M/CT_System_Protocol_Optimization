%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to represent the pixel gain curve 
% Multibright
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%housekeeping
clear 
clc
%Parameters
%Image size for binning H4 V4
px_x = 768; 
px_y = 486;
px_x_corr = 647;
px_y_corr = 486;
current_steps = 8;
currents = [75 85 100 120 150 200 300 600];
dark_img = zeros([px_x,px_y]);
flood_img = zeros([px_x,px_y, current_steps]);
floodavg_img = zeros(size(flood_img));
nff_img =zeros(size(flood_img));
avg_vector = zeros([1,8]);
central_px = zeros([1,6]);
%value = zeros([1,8]);
%location = zeros([1,8]);
%Process
%Loading images and computin the normalized flat-field correction images
for i = 1:current_steps
fname = sprintf('C:/Users/aeroesp/Documents/pruebas CTTB/flat_field fp0 tube10/flat_field_b44.flA%d',currents(i));
flood_img(:,:,i)=readSimpleBin(fname,px_x,px_y,1,'uint16',4);
fname = sprintf('C:/Users/aeroesp/Documents/pruebas CTTB/flat_field fp0 tube10/flat_field_b44.dk');
dark_img(:,:)=readSimpleBin(fname,px_x,px_y,1,'uint16',4);
% Average pixel value within the image
% Ib(x,y)
flood_img(:,:,i)=flood_img(:,:,i)-dark_img(:,:);
% Averaged pixel value 
avg_vector(i)=mean2(flood_img(:,:,i));
% Normalizing images
floodavg_img(:,:,i)=flood_img(:,:,i)./avg_vector(i);
%[value, location] = max(floodavg_img(:,:,i));
%[max_val(i),loc_aux]=max(value);
%max_loc(i)=location(loc_aux);
%nff_img(:,:,i)= floodavg_img(:,:,i).*max_val(i);
nff_img(:,:,i)= flood_img(:,:,i).*floodavg_img(:,:,i);
end
% Least squares fitting
currents =[75 85 100 120 150 200 ]; 
central_px(1,:)=nff_img(348,228,1:6);
ftype=fittype({'x','1'},'coeff',{'m','n'});
fresult=fit(currents()',central_px',ftype);

pendiente=fresult.m
ord_origen=fresult.n
%plotting
figure(1);
hold on;
plot(currents,central_px,'b*');
peso_fit=pendiente.*currents+ord_origen;
plot(currents,peso_fit,'r-');
%x=lsqnonneg(central_px,mean2(flood_img(:,:)));
%plot(x);
