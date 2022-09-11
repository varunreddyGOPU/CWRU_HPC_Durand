function functionSpikeCounting(Ref_Exp,noiseRMS_B,noiseRMS_C)
% Written by Aaron Rodrigues 11.14.21-3.10.22
%
%takes an input of three file names: the ENG recording, the corresponding time,
%and the trigger which shows when stimulus was active (trigger input is
%optional). (in future, could input just a folder path and then it would
%load all of the files from that folder?)
% 
% Inputs in order should be:*
% 
% base file name as a string (without _data, etc), threshold level
% (as number of deviations (measured from the baseline) above zero), any 
% input to signify using saved pca, any input to signify using saved clusters,
% Name of the file that the saved pca/clusters are under,
% type of stimulus (as a string, optional),
%concentration of stimulus (if necessary/appropriate). ifdenoise is true/
% false if you want to run the wavelet denoising. last input (ifplot)
%tells if you want the function to plot the data or not.  defaults to 1
%(i.e. plot) if skipping one of these inputs, need to input ~
% 
% The type of stimulus would allow the function
% to automagically sort the output into the corresponding folder (this is
% still in progress - for now, will have separate program for different
% stimuli). the additional inputs, varargin, are stored as a cell.  

% code below checks the number of inputs passed into the function; if
% optional inputs are excluded, it just sets them to default values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants
Ref_Exp;
ifdenoise=true;

%% Set some of the parameters, loads in experimental and control data

ENG_in = regexprep(Ref_Exp, '.rhd', '_dataB_filt.mat')
ENGc_in = regexprep(Ref_Exp, '.rhd', '_dataC_filt.mat')
time_in = regexprep(Ref_Exp, '.rhd', '_time.mat')
%trigger_in = regexprep(Ref_Exp, '.rhd', '_trigger.mat')

Fs = 20000; %sampling rate
data10=load(ENG_in);
data1cell=struct2cell(data10);
dataB=cell2mat(data1cell(1:1));

data10=load(ENGc_in);
data1cell=struct2cell(data10);
dataC=cell2mat(data1cell(1:1));

time1 = load(time_in);
time = time1.time;

factor = 4.5
noise_pkpk_B = noiseRMS_B * factor
noise_pkpk_C = noiseRMS_C * factor

%% Save data
dataB_RMS = rms(dataB);
dataC_RMS = rms(dataC);
dataBPower = dataB_RMS^2;
dataCPower = dataC_RMS^2;

dataB_RMSName = regexprep(ENG_in, '.mat', '_pre_RMS.txt');
dataC_RMSName = regexprep(ENGc_in, '.mat', '_pre_RMS.txt');
dataB_PowerName = regexprep(ENG_in, '.mat', '_pre_power.txt');
dataC_PowerName = regexprep(ENGc_in, '.mat', '_pre_power.txt');

fid3 = fopen( dataB_RMSName, 'wt' );
fprintf(fid3, '%s\n', dataB_RMS);
fclose(fid3);

fid4 = fopen( dataC_RMSName, 'wt' );
fprintf(fid4, '%s\n', dataC_RMS);
fclose(fid4);

fid3 = fopen(dataB_PowerName, 'wt' );
fprintf(fid3, '%s\n', dataBPower);
fclose(fid3);

fid4 = fopen(dataC_PowerName, 'wt' );
fprintf(fid4, '%s\n', dataCPower);
fclose(fid4);
%% Figs
figure(1);
hold on;   
plot(time,dataB);
%line([time(1) time(length(time))],[noise_pkpk_B noise_pkpk_B],'Color','k','LineWidth',2);
xlabel('Time (s)')
ylabel('ENG (\muV)')
title('Spikes Identified by Cluster for B Channel')

dataA_band_figurename = regexprep(ENG_in, '.mat', '_preFilt_spikes');
%savefig(Name2); 
dataA_band_imgname = strcat(dataA_band_figurename, '.png');
saveas(gcf,dataA_band_imgname)

figure(2); 
hold on;   
plot(time,dataC);
%line([time(1) time(length(time))],[noise_pkpk_C noise_pkpk_C],'Color','k','LineWidth',2);
xlabel('Time (s)')
ylabel('ENG (\muV)')
title('Spikes Identified by Cluster for C Channel')

dataB_rate_figurename = regexprep(ENGc_in, '.mat', '_preFilt_spikes');
%savefig(Name3); 
dataB_band_imgname = strcat(dataB_rate_figurename, '.png');
saveas(gcf,dataB_band_imgname)

%% Denoise 
% Set control/experimental data to 0 for control above and below threshold
count=0;
%maxdata=max(dataB)
for i = 1:length(dataB)
    if abs(dataC(i))>=noise_pkpk_C
        dataB(i)=0;
        dataC(i)=0;
        count=count+1;
    end
end
count
%maxdataf=max(dataB)

% Wavelet denoising
if ifdenoise %if user wants data denoised, run the wavelet denoising
    dataB = wavedenoise(dataB);
    dataC = wavedenoise(dataC);
end

%% Save data
dataB_RMS = rms(dataB);
dataC_RMS = rms(dataC);
dataBPower = dataB_RMS^2;
dataCPower = dataC_RMS^2;

dataB_RMSName = regexprep(ENG_in, '.mat', '_post_RMS.txt');
dataC_RMSName = regexprep(ENGc_in, '.mat', '_post_RMS.txt');
dataB_PowerName = regexprep(ENG_in, '.mat', '_post_power.txt');
dataC_PowerName = regexprep(ENGc_in, '.mat', '_post_power.txt');

fid = fopen( dataB_RMSName, 'wt' );
fprintf(fid, '%s\n', dataB_RMS);
fclose(fid);

fid2 = fopen( dataC_RMSName, 'wt' );
fprintf(fid2, '%s\n', dataC_RMS);
fclose(fid2);

fid = fopen(dataB_PowerName, 'wt' );
fprintf(fid, '%s\n', dataBPower);
fclose(fid);

fid = fopen(dataC_PowerName, 'wt' );
fprintf(fid, '%s\n', dataCPower);
fclose(fid);

%% Plots post-filter
figure(3);
hold on;   
plot(time,dataB);
line([time(1) time(length(time))],[noise_pkpk_B noise_pkpk_B],'Color','k','LineWidth',2);

xlabel('Time (s)')
ylabel('ENG (\muV)')
title('Spikes Identified for B Channel')

dataB_band_figurename = regexprep(ENG_in, '.mat', '_postFilt_spikes');
%savefig(dataB_band_figurename); 
dataB_band_imgname = strcat(dataB_band_figurename, '.png');
saveas(gcf,dataB_band_imgname)

figure(4); 
hold on;   
plot(time,dataC);
line([time(1) time(length(time))],[noise_pkpk_C noise_pkpk_C],'Color','k','LineWidth',2);
xlabel('Time (s)')
ylabel('ENG (\muV)')
title('C Channel')

dataC_rate_figurename = regexprep(ENGc_in, '.mat', '_postFilt_spikes');
%savefig(dataC_band_figurename); 
dataC_band_imgname = strcat(dataC_rate_figurename, '.png');
saveas(gcf,dataC_band_imgname)
%% FFT
try
    L = length(dataB);    % Number of datapoints
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    f = Fs/2*linspace(0,1,NFFT/2+1);

    YB = fft(dataB,NFFT)/L;
    YC = fft(dataC,NFFT)/L;
    plot_name = 'Differential';
    ENG_min=500;
    ENG_max = 1200;

    % Plot single-sided amplitude spectrum.
    figure(6);

    subplot(2,1,1);
    plot(f,2*abs(YB(1:NFFT/2+1)),'r')
    title({plot_name 'Single-Sided Amplitude Spectrum B Channel', "Filter range: ["+ENG_min+" "+ ENG_max+"]"} );
    xlim([ENG_min-200 ENG_max+500])
    xlabel({'Frequency (Hz)'})
    ylabel('|Y(f)|')
    subplot(2,1,2);
    plot(f,2*abs(YC(1:NFFT/2+1)),'b')
    title({plot_name 'Single-Sided Amplitude Spectrum C Channel', "Filter range: ["+ENG_min+" "+ ENG_max+"]"})
    xlabel({'Frequency (Hz)';'  '; [Ref_Exp]},'Interpreter', 'none');
    ylabel('|Y(f)|')
    xlim([ENG_min-200 ENG_max+500])

    fft_figurename = regexprep(Ref_Exp, '.rhd', '_FFT_lowband');
    %savefig(dataB_rate_figurename); 
    fft_imgname = strcat(fft_figurename, '.png');
    saveas(gcf,fft_imgname)
catch
end
% %% Spike detection
% try
% 	% Runs spike detect algorithm on experimental data only
% 	spikes = spikedetect(dataB(1:(length(dataB)-0.003*Fs)),noise_pkpk_B,Fs)
%     
%     spikesName = regexprep(ENG_in, '.mat', '_spikes.mat');
%     spikesNameCSV = regexprep(ENG_in, '.mat', '_spikes.csv');
%     spikesCountName = regexprep(ENG_in, '.mat', '_spikes_count.txt');
%     SNRName = regexprep(ENG_in, '.mat', '_SNR.txt');
% 
% 	save(spikesName,'spikes')
%     
%     csvwrite(spikesNameCSV,spikes)
%     
%     fid = fopen( spikesCountName, 'wt' );
%     fprintf(fid, '%d', length(spikes));
%     fclose(fid);
%     
%     avg_spike_amplitude = mean(spikes(1,:));
%     SNR = avg_spike_amplitude / noiseRMS_B;
%     fid = fopen( SNRName, 'wt' );
%     fprintf(fid, '%d', SNR);
%     fclose(fid);
% catch
%     % No spikes detected (or some other error?)
%     spikesCountName = regexprep(ENG_in, '.mat', '_spikes_count.txt');
%     SNRName = regexprep(ENG_in, '.mat', '_SNR.txt');
% 
%     fid = fopen( spikesCountName, 'wt' );
%     fprintf(fid, '%d', 0);
%     fclose(fid);
%     
%     fid = fopen( SNRName, 'wt' );
%     fprintf(fid, '%d', 0);
%     fclose(fid);
% end
% 
% %Save
%
% figure(3);
% hold on;   
% plot(spikes(2,:)./Fs + time(1),spikes(1,:),'k*')
% 
% dataB_band_figurename = regexprep(ENG_in, '.mat', '_postFilt_spikes_stars');
% %savefig(dataB_band_figurename); 
% dataB_band_imgname = strcat(dataB_band_figurename, '.png');
% saveas(gcf,dataB_band_imgname)
% 
% %% Spike rate
% try
%     spiketimes = spikes(2,:)./Fs;
%     spikerate = zeros(1,round(length(time)/Fs));
%     spikeratetime = time(1):(time(length(time))-1);
% 
%     for i = 1:length(spiketimes)
%         index = 1 + round(spiketimes(i) - time(1)/Fs);
%         spikerate(index) = spikerate(index) + 1;
%     end
% 
%     figure(5);
%     plot(spikeratetime,spikerate)
%     xlabel('Time (s)')
%     ylabel('Impulses/s')
%     title(['Impulses/Second for B Channel'])
% 
%     dataB_rate_figurename = regexprep(ENG_in, '.mat', '_postFilt_spikerate');
%     %savefig(dataB_rate_figurename); 
%     dataB_band_imgname = strcat(dataB_rate_figurename, '.png');
%     saveas(gcf,dataB_band_imgname)
% 
%     dataB_MaxRateName = regexprep(ENG_in, '.mat', '_post_spikerateMax.txt');
%     fid = fopen( dataB_MaxRateName, 'wt' );
%     fprintf(fid, '%s\n', num2str(max(spikerate)));
%     fclose(fid);
% catch
% end

end
