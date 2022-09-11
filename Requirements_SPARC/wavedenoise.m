function [denoised_data] = wavedenoise(data)
        lvl = 4; %set number of levels for decomposition
        [CA, LA] = wavedec(data,lvl,'sym7');
%       CA is all levels put together, LA gives length of each 
%       decomposition level 
        CA4 = CA(1:LA(1));
        CD4 = CA((LA(1)+1):(LA(1)+LA(2)));
        CD3 = CA((LA(1)+LA(2)+1):(LA(1)+LA(2)+LA(3)));
        CD2 = CA((LA(1)+LA(2)+LA(3)+1):(LA(1)+LA(2)+LA(3)+LA(4)));
        CD1 = CA((LA(1)+LA(2)+LA(3)+LA(4)+1):(LA(1)+LA(2)+LA(3)+LA(4)+LA(5)));
%       potentially rewrite this to change with # of levels
        deviation(1) = std(CA4);
        deviation(2) = std(CD4);
        deviation(3) = std(CD3);
        deviation(4) = std(CD2);
        deviation(5) = std(CD1);
%       make a threshold level for each freq interval (can be changed)
%       still need to rewrite this so that it is only using the deviation
%       from the baseline
        threshold(1) = deviation(1)*0.5;
        threshold(2) = deviation(2)*1;
        threshold(3) = deviation(3)*2;
        threshold(4) = deviation(4)*1;
        threshold(5) = deviation(5)*1;

%       threshold each interval.  hard thresholding is used to maintain signal
        TdataA4 = wthresh(CA4,'h',threshold(1));
        TdataD4 = wthresh(CD4,'h',threshold(2));
        TdataD3 = wthresh(CD3,'h',threshold(3));
        TdataD2 = wthresh(CD2,'h',threshold(4));
        TdataD1 = wthresh(CD1,'h',threshold(5));
        Tdata = [TdataA4,TdataD4,TdataD3,TdataD2,TdataD1];
        %reconstruct signal after thresholding
        denoised_data = waverec(Tdata,LA,'sym7');
%         plots both 
%         figure
%         subplot(2,1,1)
%         plot(time,data);
%         subplot(2,1,2)
%         plot(time,denoised_data)
end