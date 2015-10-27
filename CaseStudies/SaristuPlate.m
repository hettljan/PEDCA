% clear all
fprintf('### STARTING THE LaseLineScanAnalysis ###\n\n')

%% LOADING PARAMETERS
folder='C:\Users\u0088749\KULeuven\PROJECTS\SARISTU\Signals\SimplePanel\DispersionCurves';
file='0deg600kHz05cycles.tdms';
chNamePattern='x %f y %f/Phase0';
fileName=fullfile(folder,file);
spatialStep=0.0002;     % Spatial step in [m]

%% PROCESSING PARAMETERS
sigPars.normalize=1;        % Normalize the waveform from -1 to 1
sigPars.envelopeOn=0;       % Enable the envelope
sigPars.detrendOn=1;        % Detrend the input signals
sigPars.filterOn=1;         % Enable the predefined filter
sigPars.filterType='Low';   % Filter type - 'Low', 'High', 'Band'
sigPars.Fstop1=160e3;       % Low stop freqeuncy
sigPars.Fpass1=200e3;       % Low pass frequency
sigPars.Fpass2=0.8e6;       % High pass frequency
sigPars.Fstop2=1e6;         % High stop freqeuncy
sigPars.TimeLims=[nan nan]; % Limits used to crop time domain signals, [nan nan]
sigPars.TimeWindowOn=0;     % Time domain window

%% 2D FFT PARAMETERS
disperPars.zeroPad=1;               % Enable zero padding
disperPars.winType='tukey';         % Common spatial and temporal window
disperPars.power='lin';             % Linear or power 2 scaling of the 2D FFT
disperPars.nPeaks=3;                % Number of peaks to look for in abs(2DDFT)
disperPars.freqUnit='kHz';          % Temporal freq. units
disperPars.alpha=0.2;               % cutoff for the peak height in abs(2DFFT), alpha*max(abs(2DFFT))
disperPars.normalize2DFFT=1;        % Normalize 2D FFT by global/column max?
FreqLims=[1e3 500e3];    % Frequency limit for displaying the results  

%% 
filtBackprop=0;     % Filter out backpropagating waves
figOn=1;            % Enable plotting
saveOn=0;           % Enable saving
scalingCoef=0.3;    % scaling coefficient for waterfall plots
smoothVel=0;        % Smoothes the velocity curves

%%
DispCurves2DFFT(fileName,spatialStep,chNamePattern,sigPars,disperPars,...
    filtBackprop,scalingCoef,FreqLims,smoothVel,figOn,saveOn)