ms150_640ma = zeros(25000,1);
ms500_640ma = zeros(43100,1);
ms500_640ma = xlsread('C:/Users/aeroesp/Documents/MATLAB/Comparativa.xlsx',2,'a1:a43025');
ms150_640ma = xlsread('C:/Users/aeroesp/Documents/MATLAB/Comparativa.xlsx',3,'a1:a24664');
m = mean(ms500_640ma,1);
index_first = find(ms500_640ma>m);
index_last = find(ms500_640ma(index_first(end):end)<2*m/3);
ms500_640ma_new = ms500_640ma(index_first(1):(index_first(end)+index_last(1)));
m = mean(ms150_640ma,1);
index_first = find(ms150_640ma>m);
index_last = find(ms150_640ma(index_first(end):end)<2*m/3);
ms150_640ma_new = ms150_640ma(index_first(1):(index_first(end)+index_last(1)));
%time_s = 0:length(ms500_640ma_new);
subplot(311);plot(ms500_640ma_new,'k-*');title('500ms');xlabel('Samples');ylabel('Dose Gy/s');             
subplot(312);plot(ms150_640ma_new,'k-o');title('150ms');xlabel('Samples');ylabel('Dose Gy/s');
subplot(313);plot(ms500_640ma_new,'g');title('Comparison');xlabel('Samples');ylabel('Dose Gy/s');
hold on
subplot(313);plot(ms150_640ma_new,'r');