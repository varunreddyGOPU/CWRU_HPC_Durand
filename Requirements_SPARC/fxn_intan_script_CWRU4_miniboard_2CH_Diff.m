function fxn_intan_script_CWRU4_miniboard_2CH_Diff(filename)
% Modified by Aaron Rodrigues 3.10.22

%For use with the CWRU4 board (4 recording channels)
% Channel A (1-2) - EEG
% Channel B (3-8) - Vagus 1
% Channel C (9-14) - Vagus 2 or control
% Channel D (15-16) - ECG

% clc
% close all
% clear


% read_Intan_RHD2000_file([directory, filename])
%if ~exist(path,'dir')
%end
%[t_amplifier, amplifier_data] = read_Intan_RHD2000_file(filename);

% [filename,t_amplifier,amplifier_data,board_adc_data] = read_Intan_RHD2000_file1(file,path);
[t_amplifier,amplifier_data,board_adc_data] = read_Intan_RHD2000_file(filename);
Fs = 20000;                    % Sampling frequency

%-----------------------------------------------------------------
% USER INPUTS
%-----------------------------------------------------------------
rec_type = 16;     % 16 = differential & 32= single-ended
include_BP = 4;    % 1 = realtime BP; 2= MAP; 3= both; 4= none
include_dataOut = 1;  % 1= dump outFiles for spike sorting; 0= no outFiles
trig_type = 0;     % 0= no_trigger, 1= automatic, 2= manual
%trig_start = 23.68;   % manual trigger start time (sec)  hyp4
%trig_stop  = 59.72;   % manual trigger stop time (sec)   hyp4
trig_start = 124.5;   % manual trigger start time (sec) hyp5
trig_stop  = 159.3;   % manual trigger stop time (sec)  hyp5
trig_scale = 10;   % Level of trigger signal plotted in ENG waveform
BP_trig_scale = .2; 
plot_FFT = 0;
plot_data = 0;


ECG_max  = 300;  % Hz
ECG_min  = 40;  % Hz
EEG_max  = 200; % Hz
ENG_min = 500; % Hz  300
ENG_max = 1200; % Hz  5000

Recording_min = 0.1;  % Hz
Recording_max = 5000; % Hz
%-----------------------------------------------------------------

% if (include_BP == 1)
%     [fileBPRAW, pathBPRAW, filterindex] = ...
%        uigetfile('*.mat', 'Select a .map Raw Blood Pressure Data File', 'MultiSelect', 'off');
%     BP_filenameRAW = [pathBPRAW,fileBPRAW];
% elseif (include_BP == 2)
%     [fileBPMAP, pathBMAP, filterindex] = ...
%        uigetfile('*.csv', 'Select a .csv MAP Blood Pressure Data TABLE', 'MultiSelect', 'off');
%     BP_filenameMAP = [pathBPMAP,fileBPMAP];
% elseif  (include_BP == 3)
%     [fileBPRAW, pathBPRAW, filterindex] = ...
%        uigetfile('*.mat', 'Select a .map Raw Blood Pressure Data File', 'MultiSelect', 'off');
%     BP_filenameRAW = [pathBPRAW,fileBPRAW];
%     [fileBPMAP, pathBPMAP, filterindex] = ...
%        uigetfile('*.csv', 'Select a .csv MAP Blood Pressure Data TABLE', 'MultiSelect', 'off');
%     BP_filenameMAP = [pathBPMAP,fileBPMAP];
% else
% end

%-----------------------------------------------------------------


if (rec_type == 16)
    plot_name = 'Differential';
    chsA = [1,2];          %A
    chsB = [3,4,5,6,7,8];  %B
    chsC = [9,10,11,12,13,14];          %C
    chsD = [15,16];  %D
else
    plot_name = 'Monopolar';
     chsAp = [1,3,5,7,9,11,13,15]; %A+
     chsAm = [2,4,6,8,10,12,14,16]; %A-
     chsBp = [17,19,21,23,25,27,29,31]; %B+  junk test dont use
     chsBm = [18,20,22,24,26,28,30,32]; %B-  junk test dont use
end

time = t_amplifier;
dataA = amplifier_data (chsA,:);
dataB = amplifier_data (chsB,:);
dataC = amplifier_data (chsC,:);
dataD = amplifier_data (chsD,:);
MA = mean(dataA);
MB = mean(dataB);
MC = mean(dataC);
MD = mean(dataD);
%% Filtering
% % data_filtB_lowband = butt_band(MB,Lowband_min,Lowband_max,Fs);
data_filtB_highband = butt_band(MB,ENG_min,ENG_max,Fs); %filtered ENG1
data_filtC_highband = butt_band(MC,ENG_min,ENG_max,Fs); %filtered ENG2 or control

d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);
data_filtA_lowband = filtfilt(d,MA);
data_filtA_lowband = lowpass(data_filtA_lowband,EEG_max,Fs); %filtered EEG
% % data_filtA_highband = butt_band(MA,Highband_min,Highband_max,Fs);
%figure;

% data_filtD_lowband = filtfilt(d,MD);
data_filtD_lowband = butt_band(MD,ECG_min,ECG_max,Fs); %filtered ECG

T = 1/Fs;                     % Sample time
% % dB_lowband = data_filtB_lowband;
dA_filt = data_filtA_lowband(1:20:end); %downsample to 1Khz 
dB_filt = data_filtB_highband;
dC_filt = data_filtC_highband;
dD_filt = data_filtD_lowband(1:20:end); %downsample to 1Khz 
% % dA_highband = data_filtA_highband;
L = length(MB);                     % Length of signal
t = (0:L-1)*T;                % Time vector

NFFT = 2^nextpow2(L); % Next power of 2 from length of y
f = Fs/2*linspace(0,1,NFFT/2+1);

% dA_lowband = downsample(dA_lowband,20); %downsampled to 1khz
MA = downsample(MA,20); %downsampled to 1khz
MD = downsample(MD,20); %downsampled to 1khz

%% Process ADC trigger input for plotting
if (trig_type == 0)
    trig_plot = zeros(1,L);
elseif (trig_type == 1)
   trigger_data = board_adc_data(1,:);
   L_trigger = length(trigger_data);
   for i=1:L_trigger
       if (trigger_data(i) > 3) 
           trig_plot(i) = 1;
       else 
           trig_plot(i) = 0;
       end
   end
else
    %L_trigger = length(board_adc_data(1,:));
    L_trigger = L;
    off_period1 = trig_start/T;
    on1  = trig_stop/T;
    on_period = on1 - (off_period1);
    off_period2 = L_trigger - (on1);
    trig_plot = [zeros(1,off_period1) ones(1,on_period) zeros(1,off_period2)];
end

trig_plot_scale = trig_plot * trig_scale;

%% FFT plotting


if plot_FFT == 1
YB = fft(dB_filt,NFFT)/L;
YC = fft(dC_filt,NFFT)/L;

% Plot single-sided amplitude spectrum.
figure(1)
subplot(2,1,1);
plot(f,2*abs(YB(1:NFFT/2+1)),'r')
title({plot_name 'Single-Sided Amplitude Spectrum B Channel', "Filter range: ["+ENG_min+" "+ ENG_max+"]"} );
xlim([ENG_min-200 ENG_max+500])
xlabel({'Frequency (Hz)'})
ylabel('|Y(f)|')
subplot(2,1,2);
plot(f,2*abs(YC(1:NFFT/2+1)),'b')
title({plot_name 'Single-Sided Amplitude Spectrum C Channel', "Filter range: ["+ENG_min+" "+ ENG_max+"]"})
xlabel({'Frequency (Hz)';'  '; [filename]},'Interpreter', 'none');
ylabel('|Y(f)|')
xlim([ENG_min-200 ENG_max+500])

fft_figurename = regexprep(filename, '.rhd', '_FFT_lowband');
%savefig(dataB_rate_figurename); 
fft_imgname = strcat(fft_figurename, '.png');
saveas(gcf,fft_imgname)

end
%%
if plot_data == 1;
figure(3) 
subplot(2,1,1); 
plot(time,MB,'g');
hold on;
plot(time,trig_plot_scale,'k');
hold off;
title({[plot_name ', 8ch Averaging, Raw, Non-Filtered Results for B Channel (min= ',num2str(Recording_min),'Hz, max= ',num2str(Recording_max),'Hz)']});
ylabel('ENG (\muV)');
subplot(2,1,2); 
plot(time,MA,'g')
hold on;
plot(time,trig_plot_scale,'k');
hold off;
title({[plot_name ', 8ch Averaging, Raw, Non-Filtered Results for A Channel (min= ',num2str(Recording_min),'Hz, max= ',num2str(Recording_max),'Hz)']});
xlabel({'Time(sec)';'  '; filename});
ylabel('ENG (\muV)');

figure(4) 
subplot(2,1,1); 
plot(time,dB_lowband,'r');
hold on;
plot(time,trig_plot_scale,'k');
hold off;
title({[plot_name ', 8ch Averaging, Lowband Filtered Results for B Channel (min= ',num2str(Lowband_min),'Hz, max= ',num2str(Lowband_max),'Hz)']});
ylabel('ENG (\muV)');
subplot(2,1,2); 
plot(time,dA_lowband,'r')
hold on;
plot(time,trig_plot_scale,'k');
hold off;
title({[plot_name ', 8ch Averaging, Lowband Filtered Results for A Channel (min= ',num2str(Lowband_min),'Hz, max= ',num2str(Lowband_max),'Hz)']});
xlabel({'Time(sec)';'  '; filename});
ylabel('ENG (\muV)');



figure(5) 
subplot(2,1,1); 
plot(time,dB_highband,'b')
hold on;
plot(time,trig_plot_scale,'k');
hold off;
title({[plot_name ', 8ch Averaging, Highband Filtered Results for B Channel (min= ',num2str(Highband_min),'Hz, max= ',num2str(Highband_max),'Hz)']}); 
ylabel('ENG (\muV)');
subplot(2,1,2); 
plot(time,dA_highband,'b')
hold on;
plot(time,trig_plot_scale,'k');
hold off;
title({[plot_name ', 8ch Averaging, Highband Filtered Results for A Channel (min= ',num2str(Highband_min),'Hz, max= ',num2str(Highband_max),'Hz)']});
xlabel({'Time(sec)';'  '; filename});
ylabel('ENG (\muV)');
end

%% Blood Pressure Plot
%-------------------------------------
if (include_BP == 1)
    load(BP_filenameRAW);
    BP_dataRAW = data(1:dataend-1000);
    BP_timeRAW = [0:1:dataend-1001] * (1/samplerate);
    figure(6)
    plot(BP_timeRAW,BP_dataRAW);
    hold on;
    plot(time,trig_plot_scale,'k');
    hold off;
    title('Abdominal Aorta Wireless Blood Pressure Measurement');
    xlabel({'Time (sec)';'  '; BP_filenameRAW})
    ylabel('Blood Pressure (mmHg)');
elseif (include_BP == 2)
    delimiterIn = ',';
    headerlinesIn = 1;
    BP_TABLE = importdata(BP_filenameMAP,delimiterIn,headerlinesIn);
    BP_timeMAP = BP_TABLE.data(:,2);
    BP_dataMAP = BP_TABLE.data(:,6);
    figure(6)
    plot(BP_timeMAP,BP_dataMAP);
    hold on;
    plot(time,((trig_plot_scale)*BP_trig_scale)+BP_dataMAP(1),'k');
    hold off;
    title('Abdominal Aorta Wireless Blood Pressure Measurement');
    xlabel({'Time (sec)';'  '; BP_filenameMAP})
    ylabel('MAP (mmHg)');
elseif (include_BP == 3)
    load(BP_filenameRAW);
    BP_dataRAW = data(1:dataend-1000);
    BP_timeRAW = [0:1:dataend-1001] * (1/samplerate);
    figure(6)
    plot(BP_timeRAW,BP_dataRAW);
    hold on;
    plot(time,trig_plot_scale,'k');
    hold off;
    title('Abdominal Aorta Wireless Blood Pressure Measurement');
    xlabel({'Time (sec)';'  '; BP_filenameRAW})
    ylabel('Blood Pressure (mmHg)');
    
    delimiterIn = ',';
    headerlinesIn = 1;
    BP_TABLE = importdata(BP_filenameMAP,delimiterIn,headerlinesIn);
    BP_timeMAP = BP_TABLE.data(:,2);
    BP_dataMAP = BP_TABLE.data(:,6);
    figure(7)
    plot(BP_timeMAP,BP_dataMAP);
    hold on;
    plot(time,((trig_plot_scale)*BP_trig_scale)+BP_dataMAP(1),'k');
    hold off;
    title('Abdominal Aorta Wireless Blood Pressure Measurement');
    xlabel({'Time (sec)';'  '; BP_filenameMAP})
    ylabel('MAP (mmHg)');
else
    
end

%% Saving Data

if (include_dataOut == 1)
%     filename = [filename(1:end-42), filename(end-23:end)];
	dataA_Lowband_filename  = regexprep(filename, '.rhd', '_dataA_filt.mat');
	dataB_Highband_filename = regexprep(filename, '.rhd', '_dataB_filt.mat');
    dataC_Highband_filename = regexprep(filename, '.rhd', '_dataC_filt.mat');
    dataD_Lowband_filename  = regexprep(filename, '.rhd', '_dataD_filt.mat');
% 
    dataA_openband_filename = regexprep(filename, '.rhd', '_dataA_openband.mat');
	dataB_openband_filename = regexprep(filename, '.rhd', '_dataB_openband.mat');
    dataC_openband_filename = regexprep(filename, '.rhd', '_dataC_openband.mat');
    dataD_openband_filename = regexprep(filename, '.rhd', '_dataD_openband.mat');
%     
	time_filename = regexprep(filename, '.rhd', '_time.mat');
% 	trig_filename = regexprep(filename, '.rhd', '_trigger.mat');
% 
	save(dataA_Lowband_filename,'dA_filt');
	save(dataB_Highband_filename,'dB_filt');
    save(dataC_Highband_filename,'dC_filt');
	save(dataD_Lowband_filename,'dD_filt');
    
    save(dataA_openband_filename,'MA');
	save(dataB_openband_filename,'MB');
    save(dataC_openband_filename,'MC');
	save(dataD_openband_filename,'MD');

    save(time_filename, 'time');
%     save(trig_filename, 'trig_plot_scale');
    
end


end