function mean_rms = functionRMS_Auto(data1)

% Written by Aaron Rodrigues 11.7.21 - 3.10.22
% For outputs of read_Intan_RHD2000_file.m

data10=load(data1);
data1cell=struct2cell(data10);
datamat1=cell2mat(data1cell(1:1));

%datamat1=load(data1);

%Fs = data10.samplerate (1,1);
Fs=20000;
%%
data_squared = datamat1.^2;
len = length(datamat1);
T_0 = len/Fs;
t = (1:len)./Fs;
%%
numSections = 10;
window = Fs/5;%Fs/4;
overlap = 0;%-Fs/4;%-Fs; % Window > overlap, lower overlap -> faster runtime

% In each region, find the subregion of length window with the lowest RMS
for i = 1:numSections
    rmsVals(i) = realmax;
    j_best(i) = 0;
    for j = (1+(i-1)*len/numSections):(window-overlap):((i)*len/numSections - window)
        rms_temp = sqrt(mean(data_squared(j:j+window)));
        if(rms_temp < rmsVals(i))
            rmsVals(i) = rms_temp;
            j_best(i) = j;
        end
    end
end

t_best = j_best./Fs;

%%
% Print the power of each:
[rms_sorted, index_rms] = sort(rmsVals,'ascend');

vals = 3;
% Pick 3 lowest ("best") RMS values, use for analysis
% Assumes that there will be at least 3 regions with a quiet subregion
for i = 2:vals+1
    rms_min(i-1) = rmsVals(index_rms(i));
    j_best_min(i-1) = j_best(index_rms(i));
    t_best_min(i-1) = t_best(index_rms(i));
end
%%
if (rms_sorted(1)*100 < rms_sorted(vals)) 
    resultsName = regexprep(data1, '.mat', '_CHECK_results.txt')
	fid = fopen( resultsName, 'wt' );
else
	resultsName = regexprep(data1, '.mat', '_results.txt')
    fid = fopen( resultsName, 'wt' );
end

for i = 1:vals
    fprintf(fid, '%s\n', rms_min(i));
    fprintf(fid, '%s\n', j_best_min(i));
    fprintf(fid, '%s\n', t_best_min(i));
end

% Calculate RMS average for all three and multiply by 10.5:
mean_rms = mean(rms_min)
noise_pkpk = mean_rms * 4.5 %3.3 * 2.5 * 3.7 * 0.95% 3 * 3.5 ; 3.3 * 2.5

fprintf(fid, '%s\n', mean_rms);
fclose(fid);

threshName = regexprep(data1, '.mat', '_noiseRMS.txt')

fid2 = fopen( threshName, 'wt' );
fprintf(fid2, '%s\n', num2str(mean_rms));
fclose(fid2);

return;