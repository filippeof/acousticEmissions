%% read binary file containing analog data (all data or in segments), 
% count emissions (amplitude of signal is above threshold), 
% and record segments (csv file) if emissions occur.
%% Copyright 2022, Filippe Ferrreira
% Uni Heidelberg, 08.2022
% MIT License (available at: https://opensource.org/licenses/MIT)
%_________________________________________________________________________%
%% PARAMETERS
%path
folderName='/home/ffer/Documents/HD21/microDef/experimentalDeformation/Data/Osser';
fileName='005YP_FF_2022-03-25 151113.1_Data.txt';
inputFname=fullfile(folderName,fileName);

%format of data: 13 columns, 2nd column is string, the rest double
formatStr='%f%s%f%f%f%f%f%f%f%f%f%f%f';
griggsData=readtable(inputFname,'Delimiter','\t','Format',formatStr,'NumHeaderLines',1,'ReadVariableNames',false,'EmptyValue',nan);

% get data
time=griggsData.Var1; %time in seconds since measurement started
sig1=griggsData.Var13; %sig 1 in MPa
sig3=griggsData.Var12; %sig 3 in MPa
dispSig1=griggsData.Var4;%disp sig 1 in mm
dispSig3=griggsData.Var3;% disp sig 3 in mm
tempSample=griggsData.Var9;% temperature in deg C
%% plot whole data
figure; 

yyaxis left
plot(time,dispSig1)%,'DisplayName','Displacement')
ylabel("Displacement [mm]")

hold on
yyaxis right
plot(time,sig1)%,'DisplayName','sig1 [MPa]')
ylabel("\sigma_1 [MPa]")
xlabel("time [s]")

hold off
%% slice data based on time
t_start=3000; % starting time in seconds
t_end = 3350; % end time in seconds
id_slice=find(time>t_start & time<t_end);

figure; 

yyaxis left
plot(time(id_slice),dispSig1(id_slice))%,'DisplayName','Displacement')
ylabel("Displacement [mm]")

hold on
yyaxis right
plot(time(id_slice),sig1(id_slice))%,'DisplayName','sig1 [MPa]')

ylabel("\sigma_1 [MPa]")
xlabel("time [s]")

hold off
% legend('Location','northeastoutside')
%% filter data
nRunningMedian = 3; %number of values to get median
nRunningAverage = 5; %number of values to average

sig1_medFilter=movmedian(sig1,nRunningMedian);

sig1_medMeanFilter=movmean(sig1_medFilter,nRunningAverage);

figure; 

plot(time,sig1,'DisplayName','\sigma_1 _{raw}[MPa]')
hold on
plot(time,sig1_medMeanFilter,'DisplayName','sigma_1 _{filtered}[MPa]')

ylabel("\sigma_1 [MPa]")
xlabel("time [s]")

hold off
legend('Location','northeastoutside')

%% plot strain instead of time, + differential stress
sampleInitialLength=12; % initial sample length in mm

% slice the data starting from touch point until end of deformation stage
t_start=3000; % starting time in seconds
t_end = 3350; % end time in seconds

% slice data
id_slice=find(time>t_start & time<t_end);

% calc axial strain [%]
timeDeformation=time(id_slice)-time(id_slice(1)); %make time starts with 0
displacementDeformation=dispSig1(id_slice)-dispSig1(id_slice(1)); % make displacement start from 0 at touch point
strain=100*(displacementDeformation./sampleInitialLength); % strain in %

% calc differential stress
sig1Deformation=sig1(id_slice);
sig3Deformation=sig3(id_slice);
difStressDeformation=sig1Deformation-sig3Deformation;

figure
plot(strain,difStressDeformation)
ylabel("\sigma_1 - \sigma_3[MPa]")
xlabel("Axial Strain [%]")
%% stress curve + count of peaks
%Griggs
    % lets say you started the griggs log at:
    startGriggsRecording= datetime('10-08-2022_08:34:03','InputFormat','dd-MM-uuuu_HH:mm:ss');
%Emissions
    % lets say you started recording the emissions at:
    timeStartAcquisition = datetime('10-08-2022_09:34:03','InputFormat','dd-MM-uuuu_HH:mm:ss');
    %... and you recorded emissions during 1 hour(3600 seconds):
    emisssionRecordDuration=3600; % time in s
    %... and it thus ended recording at:
    timeEndAcquisition=timeStartAcquisition+seconds(emisssionRecordDuration);

    %this means that your emission recording time should be shifted by (+ if started after, - if started before) :
    t_shift=seconds(timeStartAcquisition- startGriggsRecording);
    
    %plot stress curve
    figure
    subplot(2,1,1)
    plot(time,sig1-sig3)
    %plot a rectangle identifying the emission count range
    pos= [t_shift, min(sig1-sig3), emisssionRecordDuration, max(sig1-sig3)];% [x y width heigth]
    hold on 
    rectangle('Position',pos)
    hold off
    % the peak data variable has two columns 1st: sample number, 2nd is intensity
    % transform number of samples to time
    sampleRate=7812500;
    peakData_ch1_s=peakData_ch1(:,1)/sampleRate;
    % add the shift
    peakData_ch1_s=peakData_ch1_s+t_shift;
    %plot
    subplot(2,1,2)
    histogram(peakData_ch1_s,'Normalization','cumcount','BinLimits',[t_shift,t_shift+emisssionRecordDuration],'NumBins',20,'DisplayName','CH 1'); % curve instead of bars: ,'DisplayStyle','stairs'
