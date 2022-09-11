function spikes = spikedetect(data,thresh_level,Fs)
%       finds any local maxima.
%       input starts at 30 for the window (basically to make sure you get a
%       full spike in).  outputs values of the maxima and the indices at
%       which they occur.
[pks,locs] = findpeaks(data(floor(0.0015*Fs):end),'MinPeakDistance',100); 
bckgr = std(data); %background deviation measured from baseline
% thresh_level = thresh*bckgr; %threshold level determined
% thresh_level = 3.3569; %7.458;
j = 1;
%       create a new vector that only contains the maxima that are above
%       the threshold level.
for i = 1:length(pks);
    if pks(i)>=thresh_level; 
        spikes(1,j) = pks(i);
        spikes(2,j) = locs(i);
        j = j+1;
    end
end
end
