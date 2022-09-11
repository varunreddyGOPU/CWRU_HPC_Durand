function filteredRMS(Ref_Exp,thresh_level_B,thresh_level_C)
    % Written by Aaron Rodrigues 2.4.22
    %
    %% Constants
    Ref_Exp;
    saved_pca=0;
    saved_centers=0;
    saved_name=0;
    stimtype=0;
    conc=0;
    ifdenoise=1;
    ifplot=1;
    %% Set some of the parameters, loads in experimental and control data
    color = ['b','k','r','g','m','c','y','.'];

    %dataC_band_filename = strcat('_dataC_filt', '.mat');
    %dataB_band_filename = strcat('_dataB_filt', '.mat');

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
    %% save
    dataB_RMS = sqrt(mean(dataB.^2));
    dataC_RMS = sqrt(mean(dataC.^2));
    
    dataB_RMSName = regexprep(ENG_in, '.mat', '_dataB_pre_RMS.txt')
    dataC_RMSName = regexprep(ENGc_in, '.mat', '_dataC_pre_RMS.txt')

    fid3 = fopen( dataB_RMSName, 'wt' );
    fprintf(fid3, '%s\n', dataB_RMS);
    fclose(fid3);
    
    fid4 = fopen( dataC_RMSName, 'wt' );
    fprintf(fid4, '%s\n', dataC_RMS);
    fclose(fid4);
    %% Figs
    figure(5);
    hold on;   
    plot(time,dataB);
    line([time(1) time(length(time))],[thresh_level_B thresh_level_B],'Color','k','LineWidth',2);
		
    xlabel('Time (s)')
    ylabel('ENG (\muV)')
    title('Spikes Identified by Cluster for B Channel')

    dataA_band_figurename_3 = regexprep(ENG_in, '.mat', '_preFilt_spikes');
    Name2=dataA_band_figurename_3; 
    %savefig(Name2); 
    dataA_band_imgname_3 = strcat(dataA_band_figurename_3, '.png');
    saveas(gcf,dataA_band_imgname_3)
    
    figure(6); 
    hold on;   
    plot(time,dataC);
    line([time(1) time(length(time))],[thresh_level_C thresh_level_C],'Color','k','LineWidth',2);
    xlabel('Time (s)')
    ylabel('ENG (\muV)')
    title('Spikes Identified by Cluster for C Channel')

    dataB_band_figurename_3 = regexprep(ENGc_in, '.mat', '_preFilt_spikes');
    Name3=dataB_band_figurename_3; 
    %savefig(Name3); 
    dataB_band_imgname_3 = strcat(dataB_band_figurename_3, '.png');
    saveas(gcf,dataB_band_imgname_3)
    
    %% Set control/experimental data to 0 for control above and below treshhold
    count=0;
    %maxdata=max(dataB)
    for i = 1:length(dataB)
        if abs(dataC(i))>=thresh_level_C
            dataB(i)=0;
            dataC(i)=0;
            count=count+1;
        end
    end
    count
    %maxdataf=max(dataB)

    %% denoise
    if ifdenoise == 1; %if user wants data denoised, run the wavelet denoising
        dataB = wavedenoise(dataB);
        dataC = wavedenoise(dataC);
    end
    
    dataB_RMS = sqrt(mean(dataB.^2));
    dataC_RMS = sqrt(mean(dataC.^2));
    %% Save
    dataB_RMSName = regexprep(ENG_in, '.mat', '_dataB_post_RMS.txt')
    dataC_RMSName = regexprep(ENGc_in, '.mat', '_dataC_post_RMS.txt')

    fid = fopen( dataB_RMSName, 'wt' );
    fprintf(fid, '%s\n', dataB_RMS);
    fclose(fid);
    
    fid2 = fopen( dataC_RMSName, 'wt' );
    fprintf(fid2, '%s\n', dataC_RMS);
    fclose(fid2);
    %% Figs
    figure(7);
    hold on;   
    plot(time,dataB);
    line([time(1) time(length(time))],[thresh_level_B thresh_level_B],'Color','k','LineWidth',2);
		
    xlabel('Time (s)')
    ylabel('ENG (\muV)')
    title('Spikes Identified by Cluster for B Channel')

    dataA_band_figurename_3 = regexprep(ENG_in, '.mat', '_postFilt_spikes');
    Name2=dataA_band_figurename_3; 
    %savefig(Name2); 
    dataA_band_imgname_3 = strcat(dataA_band_figurename_3, '.png');
    saveas(gcf,dataA_band_imgname_3)
    
    figure(8); 
    hold on;   
    plot(time,dataC);
    line([time(1) time(length(time))],[thresh_level_C thresh_level_C],'Color','k','LineWidth',2);
    xlabel('Time (s)')
    ylabel('ENG (\muV)')
    title('Spikes Identified by Cluster for C Channel')

    dataB_band_figurename_3 = regexprep(ENGc_in, '.mat', '_postFilt_spikes');
    Name3=dataB_band_figurename_3; 
    %savefig(Name3); 
    dataB_band_imgname_3 = strcat(dataB_band_figurename_3, '.png');
    saveas(gcf,dataB_band_imgname_3)
end

