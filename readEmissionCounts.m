%% read binary file containing analog data (all data or in segments), 
% count emissions (amplitude of signal is above threshold), 
% and record segments (csv file) if emissions occur.
%% Copyright 2022, Filippe Ferrreira
% Uni Heidelberg, 08.2022
% MIT License (available at: https://opensource.org/licenses/MIT)
%_________________________________________________________________________%
%% PARAMETERS
clearvars; 
close all;
%path
folderName='/home/ffer/Desktop';
fileName='data_file_169.254.41.169_2022-08-10_08-34-03.bin';
inputFname=fullfile(folderName,fileName);
infoFname=fullfile(folderName,[fileName,'.log']);

bufferPeaksFname = fullfile(folderName,[fileName(1:end-4),'_bufferPeaks.csv']); % csv output (segments with emissions)

% Find emissions
countEmissions=1; % 1=yes, 0=no
plotPeaks = 0; % 1=yes, 0=no
% parameters to find peaks (emission counts) - > needs Signal Processing Toolbox
peakThreshold_ch1=-200; % in raw values [-16384...16384], equivalent to: [-1...1]V
peakThreshold_ch2=-360;%ch2
peakWindow_ch1= 10000; % distance until a peak is considered a new event [in number of samples]
peakWindow_ch2= 10000;%ch2
peakData_ch1=[];
peakData_ch2=[];

%write segments with peaks to a csv file
writeBufferPeaksToFile= 0; % 1=yes, 0=no

% define samples>buffer>segments
headerOffset=80; %header offset in bytes
skipOffset=92; %offset in bytes between buffers

% get sample rate and number of channels from info file 
fileId_info=fopen(infoFname);
while ~feof(fileId_info) % read the whole file, line by line, until end of file
    curline = strtrim(fgets(fileId_info)); % current line

    if startsWith(curline,'====   Data transfer')
        val =  regexp(curline,'\d+','match'); %get number from string
        timeEndAcquisition =  datetime([val{:}],'InputFormat','yyyyMMddHHmmss');
    end

    if startsWith(curline,'Current ADC speed')
        val =  regexp(curline,'\d+','match'); %get number from string
        sampleRate =  str2num(val{:});
    end

    if startsWith(curline,'The total amount')
        val =  regexp(curline,'\d+','match'); %get number from string
        nChannel =  str2num(val{:});
    end
end
fclose(fileId_info);

% get start datetime from file name
timeStartAcquisition=datetime(fileName(end-22:end-4),'InputFormat','yyyy-MM-dd_HH-mm-ss'); 
% total duration
acquisitionDuration= seconds(timeEndAcquisition-timeStartAcquisition); % duration in seconds

bufferSize = 2^14; % number of samples per buffer
nBuffers = 5; % number of buffers to read per segment (~5 looks like the sweet spot)
nSamplesPerChannel = nBuffers*bufferSize; % number of samples in a segment (ignoring offset) per channel
nSamples = nChannel*nBuffers*bufferSize; % number of samples in a segment (ignoring offset), all channels

% formatStr=[num2str(nChannel*bufferSize),'*int16=>int16']; % 16384*int8 for 8 bits? '=>int16' makes the output int16 as well
formatStr=[num2str(nChannel*bufferSize),'*int16=>single']; %find peaks needs 'single' or 'double'as output.

% Start counters
xOffset=0;

if writeBufferPeaksToFile % write headers to csv file. overwrite if file exists already
    if nChannel==1
        writematrix(['#x','CH1'],bufferPeaksFname,'WriteMode','overwrite');
    else
        writematrix(['#x','CH1','CH2'],bufferPeaksFname,'WriteMode','overwrite');
    end
end

%_________________________________________________________________________%

% GET DATA 

% get id of input file
fileID = fopen(inputFname);

% skip the headers
fseek(fileID,headerOffset,'bof');

% Single channel
if nChannel == 1
    while ~feof(fileID) % if end of file was not reached yet
    
        curSegment = fread(fileID,nSamples,formatStr,skipOffset); % 1 channel
        x=xOffset+(1:length(curSegment))';
     
        if max(curSegment)>peakWindow_ch1
            if countEmissions
                if plotPeaks
                    findpeaks(curSegment,x,'MinPeakDistance',peakWindow_ch1,'MinPeakHeight',peakThreshold_ch1);
                    ylim([-16000 16000])
                    drawnow
                    pause(0.1) % hang a bit, so we can see the data before plot is updated
                else
                    [ypos,xpos] = findpeaks(curSegment,x,'MinPeakDistance',peakWindow_ch1,'MinPeakHeight',peakThreshold_ch1);
                    peakData_ch1 = [peakData_ch1; [xpos,ypos]];
                end
            end
    
            if writeBufferPeaksToFile
                writematrix([x,curSegment],bufferPeaksFname,'WriteMode','append');
            end
        end
        xOffset=xOffset+nSamplesPerChannel;
    end
    if ~isempty(peakData_ch1)
        % plot histogram of emissions
        histogram(peakData_ch1(:,1)/sampleRate,'Normalization','cdf','NumBins',20,'DisplayName','CH 1');
        
        xlim([0,x(end)/sampleRate]) %0 to last value
        xlabel('Time [s]')
        ylabel('Emission Count')
    end

else  
    % 2 channels
    while ~feof(fileID) % if end of file was not reached yet
    
        curSegment = fread(fileID,nSamples,formatStr,skipOffset); % current segment
        x = xOffset+(1:length(curSegment)/nChannel)';
        % make every column a buffer
        curSegment_all = reshape(curSegment,[bufferSize,length(curSegment)/bufferSize]);
        % every two buffers starting from column 1 are CH1
        curSegment_ch1 = curSegment_all(:,1:2:end);
        curSegment_ch1 = curSegment_ch1(:);
        % every two buffers starting from column 2 are CH2
        curSegment_ch2 = curSegment_all(:,2:2:end);
        curSegment_ch2 = curSegment_ch2(:);
     
        if or(max(curSegment_ch1)>peakThreshold_ch1,max(curSegment_ch2)>peakThreshold_ch2) % peak in ch1 or ch2
            if countEmissions
                if plotPeaks
                    if max(curSegment_ch1)>peakThreshold_ch1
                        subplot(2,1,1)
                        findpeaks(curSegment_ch1,x,'MinPeakDistance',peakWindow_ch1,'MinPeakHeight',peakThreshold_ch1);
                        ylim([-16000 16000])
                    end
                    if max(curSegment_ch2)>peakThreshold_ch2
                        subplot(2,1,2)
                        findpeaks(curSegment_ch2,x,'MinPeakDistance',peakWindow_ch2,'MinPeakHeight',peakThreshold_ch2);
                        ylim([-16000 16000])
                    end
                    drawnow
                    pause(0.1) % hang a bit, so we can see the data before plot is updated
                else
                    if max(curSegment_ch1)>peakThreshold_ch1 % peak in ch1 or ch2
                        [ypos_ch1,xpos_ch1] = findpeaks(curSegment_ch1,x,'MinPeakDistance',peakWindow_ch1,'MinPeakHeight',peakThreshold_ch1);
                        peakData_ch1 = [peakData_ch1; [xpos_ch1,ypos_ch1]];
                    end
                    if max(curSegment_ch2)>peakThreshold_ch2
                        [ypos_ch2,xpos_ch2] = findpeaks(curSegment_ch2,x,'MinPeakDistance',peakWindow_ch2,'MinPeakHeight',peakThreshold_ch2);
                        peakData_ch2 = [peakData_ch2; [xpos_ch2,ypos_ch2]];
                    end
                end
            end
    
            if writeBufferPeaksToFile
                % write whole segment to csv file. TODO: write only in the
                % interval [peakPosition-windowSize, peakPosition+windowSize]?
                writematrix([x,curSegment_ch1,curSegment_ch2],bufferPeaksFname,'WriteMode','append');
                writematrix([nan,nan,nan],bufferPeaksFname,'WriteMode','append'); % write nan row to divide segments

            end
        end
        xOffset=xOffset+nSamplesPerChannel;
    end

    if or(~isempty(peakData_ch1), ~isempty(peakData_ch2))
        % plot histogram of emissions
        yyaxis left
        histogram(peakData_ch1(:,1)/sampleRate,'Normalization','cumcount','BinLimits',[0,x(end)/sampleRate],'NumBins',20,'DisplayName','CH 1'); % curve instead of bars: ,'DisplayStyle','stairs'
        ylabel('Cumulative Count of emissions')
        hold on
        
        yyaxis right
        histogram(peakData_ch2(:,1)/sampleRate,'Normalization','cumcount','BinLimits',[0,x(end)/sampleRate],'NumBins',20,'DisplayName','CH 2');
        hold off
        hold off

        legend('Location','northeastoutside')
        xlabel('Time [s]')
        ylabel('Cumulative Count of emissions')
    end

end

totalDurationCalculated= x(end)/sampleRate; % total duration in seconds 
if (abs(totalDurationCalculated-acquisitionDuration)/acquisitionDuration)<0.98
    warning(['The calculated Duration differs from the acquisition Duration informed in the info file by: %.2f seconds.\n', ...
        'You might need to resample your data (e.g. using resample function)\n'],(totalDurationCalculated-acquisitionDuration))
end
fclose(fileID);

%% get full data, skipping headers and offsets
% if file size is too big  (>~ max ram/2), it might raise Out of Memory error (or just crashes :))
clear all; close all;
folderName='/home/ffer/Desktop';
fileName='data_file_169.254.41.169_2022-08-08_12-45-40.bin';
inputFname=fullfile(folderName,fileName);
headerOffset=80; %header offset in bytes
skipOffset=92; %in between offset in bytes
nChannel=2; % Number of channels: 1 or 2
sampleRate= 7812500; %samples per second
bufferSize=2^14; % number of samples per buffer
formatStr=[num2str(nChannel*bufferSize),'*int16=>double']; % 16384*int8 for 8 bits?

% skip the headers
fileID = fopen(inputFname);
fseek(fileID,headerOffset,'bof');
% get data
fullBuf = fread(fileID,inf,formatStr,skipOffset);
fclose(fileID);
if nChannel == 1 
    x = (1:length(fullBuf))./sampleRate; %counts to time [in seconds] 
    plot(x,fullBuf,'DisplayName','CH1')
else
    fullBuf=reshape(fullBuf,[bufferSize,length(fullBuf)/bufferSize]);
    fullBuf_ch1=(fullBuf(:,1:2:end));
    fullBuf_ch1=fullBuf_ch1(:);
    fullBuf_ch2=(fullBuf(:,2:2:end));
    fullBuf_ch2=fullBuf_ch2(:);
    
    x = (1:length(fullBuf_ch1))./sampleRate; %counts to time [in seconds] 
    
    figure
    plot(x,fullBuf_ch1,'DisplayName','CH1')
    hold on
    plot(x,fullBuf_ch2,'DisplayName','CH2')
    hold off
end
    xlabel('Time [s]')
    ylabel('Raw Units') %counts?mV(./2^16)?
    legend ('Location','northeastoutside')
    axis tight
%% get frequencies using fft
% e.g. using fullBuf_ch1
% [s,f,t] = spectrogram(fullBuf_ch1,bufferSize,[],[],sampleRate);
nBuffers=1E2; %number of buffers
figure
spectrogram(fullBuf_ch1(1:nBuffers*bufferSize),bufferSize,[],[],sampleRate); %first 100 buffers
set(gca,'XScale','log')
%% Remove low-frequency noise
wPass=10000;% passband frequency [in Hz]. (higher will pass)
yData=fullBuf_ch1(1:(nBuffers*bufferSize));
yData=reshape(yData,[bufferSize,nBuffers]);
fullBuf_ch1_hPass = highpass(yData,wPass,sampleRate); % first 100 buffers

figure;
tiledlayout(2,1) % 2 rwos 1 column

x_hPass = (1:nBuffers*bufferSize)./sampleRate; % counts to time [in seconds] 

ax1=nexttile;
plot(x_hPass,fullBuf_ch1(1:length(x_hPass)),'DisplayName','Raw data');
ylabel('Raw Units') %counts?mV(./2^16)?
legend ('Location','northeastoutside')
axis tight
yLimAx1=get(gca,'YLim')

ax2=nexttile;
plot(x_hPass,fullBuf_ch1_hPass(:),'DisplayName',['\omega > ',num2str(wPass),'Hz']);
xlabel('Time [s]')
ylabel('Raw Units') %counts?mV(./2^16)?
legend ('Location','northeastoutside')
axis tight
set(ax2,'YLim',yLimAx1)

linkaxes([ax1,ax2],'xy')% link axes
%% get full data (including headers)
% % if file size is too big  (>~ max ram/2), it might raise Out of Memory error (or just crashes :))
% clear all; close all;
% folderName='/home/ffer/Desktop';
% fileName='data_file_169.254.41.169_2022-08-04_15-32-45.bin';
% inputFname=fullfile(folderName,fileName);
% 
% fileID = fopen(inputFname);
% fullBuf = fread(fileID,inf,'int16');
% fclose(fileID);


%% get file size in bytes, get expected number of buffer segments
% nBuffers should be smaller than totalnBuffers
% fileID = fopen(inputFname);
% fseek(fileID,0,'eof');
% filesize = ftell(fileID)
% fclose(fileID);
% totalnBuffers=round((filesize-headerOffset)/(nChannel*2^14 + skipOffset/2))-8
