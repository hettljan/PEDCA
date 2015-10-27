clear all
fprintf('\nSTARTING THE LaseLineScanAnalysis ...\n')

%% LOADING PARAMETERS
folder='C:\Users\u0088749\KULeuven\PROJECTS\ALAMSA\Signals\ThinCFRPPlate\DispersionCurves\IAI_REF_A';
file='sine750kHzHann.tdms';
chNamePattern='x %f y %f/Phase0';
fileName=fullfile(folder,file);
fs=10e6;                % Sampling Frequency        
spatialStep=0.0001;     % Spatial step in [m]

%% PROCESSING PARAMETERS
normalize=1;        % Normalize the waveform from -1 to 1
envelopeOn=0;       % Enable the envelope
filterOn=1;         % Enable the predefined filter
filterType='Low';   % Filter type - 'Low', 'High', 'Band'
Fstop1=160e3;       % Low stop freqeuncy
Fpass1=200e3;       % Low pass frequency
Fpass2=0.9e6;         % High pass frequency
Fstop2=1.1e6;       % High stop freqeuncy
detrendOn=1;        % Detrend the input signals
TimeLims=[nan nan];  % Limits used to crop time domain signals, [nan nan]

%% Pencil PARAMETERS
zeroPad=1;            % Enable zero padding
N=1000;                % Number of processed samples
tolnModes=6;          % Tolerance/ Number of modes
freqUnit='kHz';       % Temporal freq. units
FreqLims=[1e3 1.1e6]; % Frequency limit for displaying the results 
VelLims=[0 10000];    % y-axis limits for phase velocity
KLims=[0 800];        % y-axis limits for wavenumber
figOn=1;              % Enable plotting
saveData=1;           % Enable data saving

%% LOAD DATA
fprintf('\nLoading data ...\n');
Signals=LoadCscanDataTDMS(fileName,spatialStep*1000,spatialStep*1000,chNamePattern);
fprintf('\nSqueezing data ...\n');
Signals=squeeze(Signals)';

%% CROP THE TIME DOMAIN SIGNALS
if isnan(TimeLims(1)) == 0 && isnan(TimeLims(2)) == 0
    fprintf('\nCropping Time Domain Signals ...\n');
    Signals=Signals(TimeLims(1):TimeLims(2),1:end-50);
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

%% SIGNAL PROCESSING - DETREND, FILTER, NORMALIZE 
if filterOn == 1
    fprintf('\nDetrend ...\n');
    Signals=detrend(Signals,'constant');
end
if filterOn == 1
    fprintf('\nFiltering ...\n');
    Signals = filtfilt(Hd,Signals);
end
if normalize == 1
    fprintf('\nNormalizing ...\n');
    for i=1:size(Signals,2)
        Signals(:,i)=Signals(:,i)./max(abs(hilbert(Signals(:,i))));
    end
end

%% GET THE WAVEFORM SIZE AND CREATE THE FREQUENCY VECTOR
[m,n]=size(Signals);
if zeroPad == 1
    m=2^nextpow2(m);
end
Freq=(-floor(m/2):floor(m/2)-1)'*fs/m;

%% CALCULATE PLAIN FFT
fprintf('\nCalculating FFT ....\n')
FFTMatrix=fftshift(fft(Signals,m),1);

%%
fprintf('\nCalculating MPM ...\n')
tic
if tolnModes<1
    PosWvns=nan(N,4);
    NegWvns=nan(N,4);
    AvgWvns=nan(N,4);
    for i=1:N
        Sig=FFTMatrix(end/2+i,:)';
        [kNeg,kPos,kAvg]=PencilMethod(Sig,spatialStep,tolnModes);
        if length(kPos)>=4
            stop=4;
        else
            stop=length(kPos);
        end
        PosWvns(i,1:stop)=kPos(1:stop);
        NegWvns(i,1:stop)=kNeg(1:stop);
        AvgWvns(i,1:stop)=kAvg(1:stop);
    end 
else
    PosWvns=nan(N,tolnModes);
    NegWvns=nan(N,tolnModes);
    AvgWvns=nan(N,tolnModes);
    parfor_progress(N); % Initialize 
    parfor i=1:N
        Sig=FFTMatrix(end/2+i,:)';
        [PosWvns(i,:),NegWvns(i,:),AvgWvns(i,:)]=PencilMethod(Sig,spatialStep,tolnModes);
        parfor_progress;
    end
    parfor_progress(0);
end
toc

%% GENERATE THE PHASE VELEOCITY MATRIX
fprintf('\n Inverting k to phase velocity\n')
VelMatrixMPM=nan(size(PosWvns));
for i=1:N
    VelMatrixMPM(i,:)=Freq(end/2+i)./AvgWvns(i,:);
end

%% GENERATE
switch freqUnit 
    case 'kHz'
        Freq=Freq./1e3;
        FreqLims=FreqLims./1e3;
    case 'MHz'
        Freq=Freq./1e6;
        FreqLims=FreqLims./1e6;
end

%% PLOTTING
fprintf('\nPlotting\n')
if figOn == 1
    figure
    subplot(2,2,1)
    plot(Freq(end/2:end/2+N-1),PosWvns,'*');
    xlim(FreqLims)
    ylim(KLims)
    xlabel(strcat('Frequency [',freqUnit,']'),'FontSize',14)
    ylabel(strcat('Wavenumber [m^{-1}]'),'FontSize',14)
    title('Positive wavenumbers')
    subplot(2,2,2)
    plot(Freq(end/2:end/2+N-1),AvgWvns,'*');
    xlim(FreqLims)
    ylim(KLims)
    xlabel(strcat('Frequency [',freqUnit,']'),'FontSize',14)
    ylabel(strcat('Wavenumber [m^{-1}]'),'FontSize',14)
    title('Back/Forw Filtered')
    
%     figure
    subplot(2,2,4)
    plot(Freq(end/2:end/2+N-1),VelMatrixMPM,'r*');
    xlabel(strcat('Frequency [',freqUnit,']'),'FontSize',14)
    ylabel(strcat('Phase velocity [ms^{-1}]'),'FontSize',14)
    xlim(FreqLims)
    ylim(VelLims)
end
   
%% 
if saveData == 1
    save('DispDataMPM', 'AvgWvns','Freq','VelMatrixMPM','PosWvns');
end