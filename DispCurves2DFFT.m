clear all
fprintf('### STARTING THE LaseLineScanAnalysis ###\n\n')

%% LOADING PARAMETERS
folder='C:\Users\u0088749\KULeuven\PROJECTS\SARISTU\Signals\SimplePanel\DispersionCurves';
file='0deg600kHz05cycles.tdms';
% folder='C:\Users\u0088749\KULeuven\PROJECTS\Other';
% file='PlexiglassA0.tdms';
% folder='C:\Users\u0088749\KULeuven\PROJECTS\SARISTU\Signals\SimplePanel\DispersionCurves';
% file='90deg750kHz05cycles.tdms';
chNamePattern='x %f y %f/Phase0';
fileName=fullfile(folder,file);
fs=5e6;                % Sampling Frequency        
spatialStep=0.0001;     % Spatial step in [m]

%% PROCESSING PARAMETERS
normalize=1;        % Normalize the waveform from -1 to 1
envelopeOn=0;       % Enable the envelope
filterOn=1;         % Enable the predefined filter
detrendOn=1;        % Detrend the input signals
filterType='Low';   % Filter type - 'Low', 'High', 'Band'
Fstop1=160e3;       % Low stop freqeuncy
Fpass1=200e3;       % Low pass frequency
Fpass2=0.8e6;         % High pass frequency
Fstop2=1e6;       % High stop freqeuncy
TimeLims=[nan nan];  % Limits used to crop time domain signals, [nan nan]
TimeWindowOn=0;     % Time domain window
filterBackprop=0;    % Filter out backpropagating waves

%% 2D FFT PARAMETERS
zeroPad=1;          % Enable zero padding
winType='tukey';    % Common spatial and temporal window
power='lin';        % Linear or power 2 scaling of the 2D FFT
numberOfPeaks=5;    % Number of peaks to look for in abs(2DDFT)
freqUnit='kHz';     % Temporal freq. units
alpha=0.2;          % cutoff for the peak height in abs(2DFFT), alpha*max(abs(2DFFT))
figOn=1;            % Enable plotting
saveOn=1;           % Enable saving
normalize2DFFT=1;   % Normalize 2D FFT by global/column max?
FreqLims=[1e3 800e3]; % Frequency limit for displaying the results  
scalingCoef=0.3;     % scaling coefficient for waterfall plots
smoothVel=0;        % Smoothes the velocity cruves

%% VELOCITY PARAMETERS
quadFit=1;          % Enable quadratic fitting of the xcorr maxima

%% LOAD DATA
fprintf('\nLoading data ...\n');
Signals=LoadCscanDataTDMS(fileName,spatialStep*1000,spatialStep*1000,chNamePattern);
fprintf('\nSqueezing data ...\n');
Signals=squeeze(Signals)';

%% CROP THE TIME DOMAIN SIGNALS
if isnan(TimeLims(1)) == 0 && isnan(TimeLims(2)) == 0
    fprintf('\nCropping Time Domain Signals ...\n');
    Signals=Signals(TimeLims(1):TimeLims(2),1:end);
end

%% DESIGN THE FILTER
if filterOn == 1
    fprintf('\nDesigning Filter ...\n');
    switch filterType
        case 'Band'
            Hd = BandPassFIR(Fstop1,Fpass1,Fpass2,Fstop2,fs);
        case 'Low'
            Hd = LowPassFIR(Fpass2,Fstop2,fs);
    end
end

%% DETREND THE BASIC SIGNALS
if detrendOn == 1
    fprintf('\nDetrend ...\n');
    for i=1:size(Signals,2)
        Signals(:,i)=detrend(Signals(:,i),'constant');
    end
end

%% APPLY WINDOW TO THE TIME DOMAIN SIGNALS
if TimeWindowOn == 1
    fprintf('\nWindowing time domain signals ...\n');
    Signals=ApplyWindowTo2DMatrix(Signals,1,'hann',0.1);
end

%% FILTER AND NORMALIZE
if filterOn == 1
   fprintf('\nFiltering ...\n');
   Signals= filtfilt(Hd,Signals);
end
if normalize ==1
    fprintf('\nNormalize ...\n');
    for i=1:size(Signals,2)
    	Signals(:,i)=Signals(:,i)./max(abs(hilbert(Signals(:,i))));
    end
end

%% FILTER THE BACKPROPAGATING WAVES
if filterBackprop == 1
    Signals=FFT2DFilterBackwardPropWaves(Signals);
end

%% PLOT WATERFALL GRAPH
if figOn ==1
    plotResultsWaterfall(Signals(:,1:5:end),fs,detrendOn,envelopeOn,scalingCoef,'us',0,'b')
end

%% 2D FFT
if envelopeOn == 0
    warning('off','signal:findpeaks:largeMinPeakHeight');
    [Freq,X,FFTMatrix,Wvns2DFFT]=FFT2D(Signals,fs,spatialStep,zeroPad,power,winType,numberOfPeaks,...
    freqUnit,alpha,normalize2DFFT,figOn,saveOn,FreqLims,smoothVel);
end

%%
% plotSpectraWaterfall(Signals(:,1:1:end),fs,detrendOn,scalingCoef,power,0,'hann','kHz',FreqLims)
