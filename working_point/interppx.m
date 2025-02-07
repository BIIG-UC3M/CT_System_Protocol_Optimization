function [ imgOut ] = interppx(imgIn,bin)
% Binning 4 image correction
% Variables
class(imgIn)
if bin==44 
        col = [32 96 160 224 288 352 385 416 449 480 481 514 544 571 577 600 608 643 672 736];
        col = [32 96 160 224 288 352 385 416 449 480 481 514 513 544 545 571 577 600 608 643 672 736];
        lines = [458];
        %px =  [33477 41531 122753 157256 195308 256509 267600 287999 307111 312598 313664 351552 351576 351616];
        pxY = [44 55 160 205 255 334 349 375 400 408 409 458 458 458];
        pxX = [453 59 641 584 236 765 336 767 679 22 320 576 600 640];
else if bin==22
        col = [63 191 319 447 575 703 831 959 1087 1198 1215 1343 1471 1284];
        lines = [915];
    else if bin==11
            col = [126 127 382 383 638 639 894 895 1150 1151 1406 1407 1662 1663 1918 1919 2174 2175 2396 2430 2431 2569 2686 2687 2942 2943];
            lines = [1830 1831];
            px=[18652    60344    62858    82547    82968    84498    86673    96540   109506   112148   120256   133462   139899   155464   156729   168813   174981   185989   240685   243996   278131   280472   293186   311920   313288   434427   448648   471106   472996   536338   581491 ...
                607139   619066   627271   649659   660059   666857   675773   679618   693297   744647   837337   840781   849930   861764   885153   939054   970787   972027  1017244  1032459  1047986  1048950  1054797  1075227  1079025  1131857  1162927  1179668  1185793  1193753  1214995 ...
                1227044  1236839  1248410  1249581  1311887  1324030  1324556  1342166  1364172  1370972  1376648  1386021  1399413  1426721  1437865  1497075  1525912  1535204  1556196  1558916  1592229  1619528  1623562  1634606  1654310  1662664  1666151  1678123  1704989  1705307  1709408 ...
                1724851  1738810  1758323  1766488  1835837  1855712  1870708  1899242  1906186  1924863  1964811  1965567  1972799  2009508  2054086  2057653  2090797  2117144  2135850  2195106  2204125  2225464  2231384  2231663  2243602  2308603  2411416  2459977  2473103  2494060  2498096 ...
                2501092  2512157  2513784  2525630  2543965  2551489  2551892  2612980  2630884  2714822  2777851  2784581  2788685  2805277  2809884  2821358  2832789  2855288  2855591  2858701  2901225  2907432  2909936  2922004  2944355  2987799  3039484  3091666  3119852  3120598  3124084 ...
                3137115  3137452  3166748  3213888  3303917  3308338  3343753  3346735  3380750  3428822  3442235  3489707  3492891  3503524  3524716  3525734  3526027  3539051  3541916  3553033  3582634  3602622  3602861  3607825  3684300  3691190  3717064  3724622  3760795  3780134  3782214 ...
                3863524  3871120  3891945  3892836  3925929  3938136  3946880  3951511  3959712  4004923  4017201  4036355  4104177  4107507  4115096  4142801  4153933  4165737  4179218  4192425  4197242  4199752  4206872  4265524  4286783  4293667  4307101  4321412  4363716  4366709  4411174 ...
                4440983  4460089  4461356  4507004  4515596  4531197  4547511  4587794  4599542  4601849  4620271  4625133  4658498  4692716  4706644  4707664  4716799  4725228  4805790  4811530  4822952  4833446  4840278  4844805  4863271  4884247  4888295  4899470  4904166  4908697  4934612 ...
                4970125  4972940  5001302  5012697  5017855  5030098  5036432  5037840  5047006  5048088  5058290  5085936  5117404  5138380  5155633  5174431  5187597  5195494  5212173  5216428  5239309  5254906  5268241  5274622  5299194  5335861  5336887  5345170  5348155  5359758  5390015 ...
                5396370  5429485  5452810  5473851  5474213  5495698  5500579  5505437  5520271  5533347  5580763  5590107  5624063  5624319  5633237  5639450  5654904  5667708  5670589  5671207  5679063  5745370  5770274  5792851  5796037  5839173  5856438  5874526  5885136  5925406  5925809 ...
                5930161  5946834  5953514  5960361];
        end
    end
end

if bin==11

imgIn(2304:2560,1831,:)=imgIn(2304:2560,1830,:);

imgIn(2304:2560,1832,:)=imgIn(2304:2560,1833,:);
end
if bin==22

imgIn(1152:1235,916,:)=(imgIn(1152:1235,915,:)+imgIn(1152:1235,915,:))./2;


end

%matriz definida como col-fil;
if bin==44
    %lines   Solo hay trozos de l�nea as� se corrigen a cap�n
%     for i=1:length(lines)
%         imgIn(:,lines(i),:)=(imgIn(:,lines(i)-1,:)+imgIn(:,lines(i)+1,:))/2;        
%     end
    imgIn(641:702,160,:)=(imgIn(641:702,159,:)+imgIn(641:702,161,:))/2;
    imgIn(577:640,458,:)=(imgIn(577:640,457,:)+imgIn(577:640,459,:))/2;
%Esto es pra corregir los pespuntes si vuelven a salir
%    imgIn(:,416,:)=(imgIn(:,415,:) +imgIn(:,417,:) )/2;

   %columns
    for i=1:length(col)
        imgIn(col(i),:,:)=(imgIn(col(i)-1,:,:)+imgIn(col(i)+1,:,:))/2;        
    end
    %pixels       
     for i=1:length(pxY)
          
        imgIn(pxX(i),pxY(i),:)= (imgIn(pxX(i)+1,pxY(i),:)+imgIn(pxX(i),pxY(i)+1,:)+imgIn(pxX(i)-1,pxY(i),:)+imgIn(pxX(i),pxY(i)-1,:))/4;    
        
     end  
   
end
imgOut=imgIn;


% 
% 
% %col = col + 1; % Since the defect map started at 0
% col_prev = col - 1;
% col_post = col + 1;
% lines = [457];
% %lines = lines + 1;
% lin_prev = lines - 1;
% lin_post = lines + 1;
% 
% %imgOut = zeros(560,586);
% %px = [ y
% %       x ]
% aux1 = 576:639;
% %aux2 = ones(size(aux1));
% 
% 
% 
% %Process
% imgOut = imgIn;
% %Auxiliar matrix with the values from the prev&post columns wrt the
% %column to correct
% aux(:,:,:,1) = imgIn(col_prev(:),:,:);
% aux(:,:,:,2) = imgIn(col_post(:),:,:);
% imgOut(col(:),:,:) = squeeze(mean(aux,4)); %Correct the columns by the mean values
% %Auxiliar matrix with the values from the prev&post lines wrt the
% %line to correct
% clear aux;
% aux(:,:,:,1) = imgIn(aux1,lin_prev(:),:);
% aux(:,:,:,2) = imgIn(aux1,lin_post(:),:);
% imgOut(aux1,lines(:),:) = squeeze(mean(aux,4));%Correct the lines by the mean values
% %Auxiliar matrix with the values from the left-right-up-down pixels wrt the
% %pixel to correct
% %clear aux;
% %aux(:,:,:,1) = imgIn(px_prev(2,:),px_prev(1,:),:);
% %aux(:,:,:,2) = imgIn(px(2,:),px_prev(1,:),:);
% %aux(:,:,:,3) = imgIn(px(2,:),px_post(1,:),:);
% %aux(:,:,:,4) = imgIn(px_post(2,:),px_post(1,:),:);
% %imgOut(px(2,:),px(1,:),:) = squeeze(mean(aux,4));%Correct the pixel by the mean value
end

