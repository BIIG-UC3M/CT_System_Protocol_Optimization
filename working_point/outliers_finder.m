temp_ref = zeros(80,21);
dose = zeros(80,21);
dose_crop = dose;
temp_crop = temp_ref;
temp_ref = xlsread('C:/Users/aortega/Documents/MATLAB/matlab_code/Comparativa.xlsx',4);
dose = xlsread('C:/Users/aortega/Documents/MATLAB/matlab_code/Comparativa.xlsx',3);
x = reshape(dose,80*21,1);
figure
plot(x,'k-*');title('Whole sequence  ');xlabel('Time s');ylabel('Dose Gy/s');
m = mean(dose,1);
tam = size(dose);
index_first = zeros(tam(1), tam(2));
index_last = zeros(tam(1), tam(2));
dist = zeros(tam(1), tam(2));
for i=2:tam(2);
    aux = find(dose(:,i)>m(i));
    last_pos = length(aux);
    index_first(:,i) = [aux;zeros(tam(1)-last_pos,1)];
    %index_first(:,i) = [index1_aux(i,:)';zeros(tam(1)-length(index1_aux(i,:)),1)];
    if index_first(last_pos,i) < tam(1)
        aux = find(dose(index_first(last_pos,i):end,i)<2*m(i)/3);
        last_pos2 = length(aux);
        index_last(:,i) = [aux;zeros(tam(1)-last_pos2,1)];
        %index_last(:,i) = [indexL_aux(i,:);zeros(tam(1)-size(indexL_aux(i)),1)];
    end
    aux = dose(index_first(1,i):(index_first(last_pos,i)+ index_last(i,1)),i);
    last_pos = length(aux);
    dose_crop(:,i) = [aux;zeros(tam(1)-last_pos,1)];
    aux = temp_ref(index_first(1,i):(index_first(last_pos,i)+index_last(i,1)),i);
    last_pos2 = length(aux);
    temp_crop(:,i) =  [aux;zeros(tam(1)-last_pos2,1)];

% m = mean(ms150_640ma,1);
% index_first = find(ms150_640ma>m);
% index_last = find(ms150_640ma(index_first(end):end)<2*m/3);
% ms150_640ma_new = ms150_640ma(index_first(1):(index_first(end)+index_last(1)));
% %time_s = 0:length(ms500_640ma_new);
% subplot(311);plot(ms500_640ma_new,'k-*');title('500ms');xlabel('Samples');ylabel('Dose Gy/s');
% subplot(312);plot(ms150_640ma_new,'k-o');title('150ms');xlabel('Samples');ylabel('Dose Gy/s');
% subplot(313);plot(ms500_640ma_new,'g');title('Comparison');xlabel('Samples');ylabel('Dose Gy/s');
% hold on
% subplot(313);plot(ms150_640ma_new,'r');

% Compute maximum value to stay within Compound Poisson range, after that
% values are considered saturated
% We do that by finding the perpendicular projection of every point on the
% line that connects the starting and ending point of the variance curve
% The point with a larger distance should be a nice estimation of the elbow
%siz = size(dose_crop(i,:));
v   = [temp_crop(last_pos2,i)-temp_crop(1,i) dose_crop(last_pos,i)-dose_crop(1,i)]; % Director vector of the line
vP  = [-v(2) v(1)];                                        % Director vector of the perpendicular line

% Now find the intersection point of the two lines for every point in
% pxVar, then the distance between the point and the line and the point w a
% larger distance will be the corner (or an estimation of it...)
% First we smooth the dose curve to somewhat reduce noise
dose_s = smooth(dose_crop(:,i),1);
% Non-vectorized code...

% Loop to solve the system of equations for the two lines and every point
% As a result we get a group of distances
for in = 1:last_pos,
    k1 = v(1)*dose_crop(1,i)   - v(2)*temp_crop(1,i);
    k2 = vP(1)*dose_crop(in,i) - vP(2)*temp_crop(in,i);
    yT = (vP(2)*k1 - v(2)*k2)/(vP(2)*v(1) - vP(1)*v(2));
    xT = (v(1)*yT - k1)/v(2);
    % Now we get the distance value
    dist(in,i) = sqrt(((dose_s(in)-yT)^2) + ((temp_crop(in,i)-xT)^2));
end
ind_sat(:,i) = find(dist(:,i) == max(dist(:,i))); % Elbow is around max dist ance to line
ind_p = find(dose(:,i) == dose_crop(ind_sat(1,i),i)); 
%clear *crop
%clear *aux
figure 
plot(temp_ref(:,i),dose(:,i),'k-*');title(['Current Step #',num2str(i)]);xlabel('Time s');ylabel('Dose Gy/s');
hold on
plot(temp_ref(ind_p,i),dose(ind_p,i),'r-o');
keyboard
end