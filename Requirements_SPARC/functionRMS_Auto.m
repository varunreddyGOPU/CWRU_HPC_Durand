function rmsMean = functionRMS_Auto(Ref_Exp, tag)

% Written by Aaron Rodrigues 11.7.21 - 9.28.22
% For outputs of read_Intan_RHD2000_file.m

data1 = regexprep(Ref_Exp, '.rhd', tag)

data10=load(data1);
data1cell=struct2cell(data10);
datamat1=cell2mat(data1cell(1:1));

%datamat1=load(data1);
%Fs = data10.samplerate (1,1);

Fs= 20000;

dataSquared = datamat1.^2;
len = length(datamat1);
T_0 = len/Fs;
t = (1:len)./Fs;

numSections = 1;
windowSeconds = 10 % Length of window in seconds
window = Fs*windowSeconds;

n = numSections;

% Overlap is the amount that windows will overlap each other; lower overlap -> faster runtime
% To add space between windows: Try -Fs/4 or -Fs
overlap = (9/10) * window;
% Ensure window > overlap, 
overlap = min([overlap, window]) 

% In each region, find the subregion of length window with the lowest RMS
for i = 1:numSections
    rmsVals(i) = realmax;
    iLowest(i) = 0;
	% Iterate over all possible windows within a section
    for j = (1+(i-1)*len/numSections):(window-overlap):((i)*len/numSections - window)
        rmsTemp = sqrt(mean(dataSquared(j:j+window)));
        if(rmsTemp < rmsVals(i))
            rmsVals(i) = rmsTemp;
            iLowest(i) = j;
        end
    end
end

tLowest = iLowest./Fs;

%%
% Print the power of each:
[rmsSorted, rmsIndex] = sort(rmsVals,'ascend');

% Pick n lowest ("best") RMS values, use for analysis
% Assumes that there will be at least n regions with a quiet subregion
for i = 1:n
    rmsMin(i) = rmsVals(rmsIndex(i));
    iMin(i) = iLowest(rmsIndex(i));
    tMin(i) = tLowest(rmsIndex(i));
end

%% Add note if output has a high variance
if (rmsSorted(1)*100 < rmsSorted(n)) 
    resultsName = regexprep(data1, '.mat', '_CHECK_results.txt')
	fid = fopen( resultsName, 'wt' );
else
	resultsName = regexprep(data1, '.mat', '_results.txt')
    fid = fopen( resultsName, 'wt' );
end

for i = 1:n
    fprintf(fid, '%s\n', "RMS: " + compose("%9.7f",rmsMin(i)));
    fprintf(fid, '%s\n', "Index: " + int2str(iMin(i)));
    fprintf(fid, '%s\n', "Time: " + compose("%9.7f",tMin(i)));
end

rmsMean = mean(rmsMin)

fprintf(fid, '%s\n', rmsMean);
fclose(fid);

threshName = regexprep(data1, '.mat', '_noiseRMS.txt')

fid2 = fopen( threshName, 'wt' );
fprintf(fid2, '%s\n', compose("%9.7f",rmsMean));
fclose(fid2);

return;