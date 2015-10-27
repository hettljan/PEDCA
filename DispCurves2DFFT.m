% %% LOADING PARAMETERS
% folder='C:\Users\u0088749\KULeuven\PROJECTS\SARISTU\Signals\SimplePanel\DispersionCurves';
% file='0deg600kHz05cycles.tdms';

% % folder='C:\Users\u0088749\KULeuven\PROJECTS\SARISTU\Signals\SimplePanel\DispersionCurves';
% % file='90deg750kHz05cycles.tdms';
% chNamePattern='x %f y %f/Phase0';
% fileName=fullfile(folder,file);
% spatialStep=0.0001;     % Spatial step in [m]
function DispCurves2DFFT(fileName,spatialStep,chNamePattern,sigPars,disperPars,...
    varargin)
% this function preprocesses the signals for 2D FFT calculation of
% dispersion curves and further calculates them
% INPUT:
%   fileName    -   String that contains absolute path to the .tdms file with
%                   signals
%   spatialStep -   Spatial step in scanning direction in [m]
%   chNamePattern - String pattern that describes the naming of the
%                   channels
%   sigPars     -   A struct that contains the all the directives for signal 
%                   processing of the raw signals. It contains the following 
%                   fields:
%                       normalize   - Normalize the waveform from -1 to 1
%                       envelopeOn  - Enable and plot the envelope, boolean
%                                     0,1
%                       detrendOn   - Detrend the input signals, boolean 0,1
%                       filterOn    - Apply the filter to the input
%                                     signals, boolean 0,1
%                       filterType  - Type of the filter to be applied,
%                                     'Low', 'High', 'Band'
%                       Fstop1      - Stop frequency for the low cutoff [Hz]
%                       Fpass1      - Pass frequency for the low cutoff [Hz]
%                       Fpass2      - Pass frequency for the high cutoff [Hz]
%                       Fstop2      - Stop frequency for the high cutoff [Hz]
%                       TimeLims    - Vector with 2 elements. Limits used to 
%                                     crop time domain signals, default is[nan nan]
%                       TimeWindowOn- Apply hann window to time domain
%                                     signals
%   disperPars  -   A struct that contains the all the directives for 2DFFT 
%                   dispersion calculation. It contains the following fields:
%                       zeroPad     - Enable zero padding of the spatial and
%                                     temporal domain, boolean 0,1
%                       winType     - Common spatial and temporal window, string 
%                                     'tukey', 'hann' or others
%                       power       - Linear or power 2 scaling of the 2D FFT,
%                                     'lin', 'quad'
%                       nPeaks      - Number of peaks to look for in abs(2DDFT)
%                       freqUnit    - Frequency units in temporal/freq doamin
%                       alpha       - Cutoff for the peak height in abs(2DFFT),
%                                     alpha*max(abs(2DFFT))
%                       normalize2DFFT- Normalize 2D FFT by global/column max?,
%                                       boolean 0,1
% OPTIONAL:
%   filtBackprop -  Filter out backpropagating waves, boolean 0,1
%   scalingCoef  -  Scaling coefficient for waterfall plots
%   FreqLims     -  Frequency limit for displaying the results,
%                    vector of two elements [fmin fmax]
%   smoothVel    - Smoothes the velocity curves, boolean 0,1
%   figOn        -  Enable plotting of the results, boolean 0,1
%   saveOn       -  Enable saving to designated file, boolean 0,1

%% INPUT PARSING
numvarargs = length(varargin);
if numvarargs > 6
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 6 optional inputs');
end
optargs = {0, 0,[0 800e3],0,1,0};
optargs(1:numvarargs) = varargin;
[filtBackprop, scalingCoef,FreqLims,smoothVel, figOn, saveOn] = optargs{:};

%% LOAD DATA
fprintf('\nLoading data ...\n');
Signals=LoadCscanDataTDMS(fileName,spatialStep*1000,spatialStep*1000,chNamePattern);
Properties=TDMSGetChannelProperties(fileName);      % get properties from file
fs=1/Properties('wf_increment');                    % Sampling frequency [Hz] 

%%
fprintf('\nSqueezing data ...\n');
Signals=squeeze(Signals)';

%% CROP THE TIME DOMAIN SIGNALS
if isnan(sigPars.TimeLims(1)) == 0 && isnan(sigPars.TimeLims(2)) == 0
    fprintf('\nCropping Time Domain Signals ...\n');
    Signals=Signals(sigPars.TimeLims(1):sigPars.TimeLims(2),1:end);
end

%% DESIGN THE FILTER
if sigPars.filterOn == 1
    fprintf('\nDesigning Filter ...\n');
    switch sigPars.filterType
        case 'Band'
            Hd = BandPassFIR(sigPars.Fstop1,sigPars.Fpass1,sigPars.Fpass2,sigPars.Fstop2,fs);
        case 'Low'
            Hd = LowPassFIR(sigPars.Fpass2,sigPars.Fstop2,fs);
    end
end

%% DETREND THE BASIC SIGNALS
if sigPars.detrendOn == 1
    fprintf('\nDetrend ...\n');
    for i=1:size(Signals,2)
        Signals(:,i)=detrend(Signals(:,i),'constant');
    end
end

%% APPLY WINDOW TO THE TIME DOMAIN SIGNALS
if sigPars.TimeWindowOn == 1
    fprintf('\nWindowing time domain signals ...\n');
    Signals=ApplyWindowTo2DMatrix(Signals,1,'hann',0.1);
end

%% FILTER AND NORMALIZE
if sigPars.filterOn == 1
   fprintf('\nFiltering ...\n');
   Signals= filtfilt(Hd,Signals);
end
if sigPars.normalize ==1
    fprintf('\nNormalize ...\n');
    for i=1:size(Signals,2)
    	Signals(:,i)=Signals(:,i)./max(abs(hilbert(Signals(:,i))));
    end
end

%% FILTER THE BACKPROPAGATING WAVES
if filtBackprop == 1
    Signals=FFT2DFilterBackwardPropWaves(Signals);
end

%% PLOT WATERFALL GRAPH
if figOn == 1
    plotResultsWaterfall(Signals(:,1:5:end),fs,sigPars.envelopeOn,...
       scalingCoef,'us',0,'b')
end

%% 2D FFT
if sigPars.envelopeOn == 0
    warning('off','signal:findpeaks:largeMinPeakHeight');
    [~,~,~,~]=FFT2D(Signals,fs,spatialStep,disperPars.zeroPad,disperPars.power,...
        disperPars.winType,disperPars.nPeaks,disperPars.freqUnit,disperPars.alpha,...
        disperPars.normalize2DFFT,figOn,saveOn,FreqLims,smoothVel);
end

%% WATERFALL PLOT OF THE RESULTS
% if scalingCoef ~= 0
%     plotSpectraWaterfall(Signals(:,1:1:end),fs,sigPars.detrendOn,scalingCoef,...
%         disperPars.power,0,'hann','kHz',FreqLims)
% end
