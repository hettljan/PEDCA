function [Freq,X,FFTMatrix,WvnMatrix]=FFT2D(Matrix,fs,spatialStep,varargin)
% This function calculates the dispersion curves for the ultrasonic waves
% from the 2D laser scan data. It uses the 2D FFT algorithm to do so. The 
% input is a typical 2D C-scan with time signals ordered in third dimension
% INPUT:
%   Matrix      - temporaly and spatially sampled data for 2D FFT
%                 signals are ordered column-wise
%   fs          - temporal sampling frequency [Hz]
%   spatialStep - spatial sampling step in [m]
% OPTIONAL: 
%   zeroPad     - zeropad the signals in both direction to the nearest 2^n
%   power       - linear, quadratic or log-compressed, 'lin', 'quad','log' 
%   winType     - type of the applied window - e.g. 'tukey','hann'
%   nPeaks      - deifine the number of peaks to search for in one column
%   freqUnit    - select the frequency units e.g. 'kHz'
%   alpha       - scaling coeeficient for peak finder
%   normalizeOn - normalize the final heat map by columns?
%   figOn       - turn on or of f-k map - values 1 or 0
%   saveData    - boolen, save the data
%   FreqLims    - Frequency limits [Hz]
%   smoothVel   - Smooth the velocity curves using moving avg and median
% OUTPUT:
%   Freq        - vector of the frequencies in selected units [hz, kHz or MHz]
%   X           - vector of the reciprocal values of wavelengths [1/m]
%   FFTMatrix   - results of the 2D FFT
%   WvnMatrix   - Matrix with the wavenumbers [m^{-1}]

%% PA
fprintf('\n### STARTING 2DFFT PROGRAM ###\n\n');

%% INPUT PARSING
numvarargs = length(varargin);
if numvarargs > 11
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 11 optional inputs');
end
optargs = {0,'lin','none', 1, 'Hz', 0.5, 0, 0, 0, [1 1000],0};
optargs(1:numvarargs) = varargin;
[zeroPad, power, winType, nPeaks,  freqUnit, alpha, normalizeOn, figOn,...
    saveData,FreqLims,smoothVel] = optargs{:};

%% INPUT CHECK
if size(Matrix,1) < size(Matrix,2)
    Matrix=Matrix';                 % if not columnwise then transpose
end

%% CONSTANTS
spatialFreq=1/spatialStep;
LimsY=[0 800];
MinPeakDist=10;              % min distance between the peaks in p-k space

%% PREPROCESSING
[m,n]=size(Matrix);         % get the size of the matrix
fprintf('\nOriginal size of the matrix is %d x %d\n',m,n);
fprintf('\nDetrending ...\n');
for i=1:n
    Matrix(:,i)=detrend(Matrix(:,i),'constant');  
end

%% APPLY STEEP TUKEY WINDOW TO CAN SMOOTH TRANSITION FOR ZEROPAD
% Matrix=ApplyWindowTo2DMatrix(Matrix,12,'tukey',0.1);

%% ZERO PADDING DATA
if zeroPad ==1
    [m,n]=size(Matrix);
    m=2^nextpow2(m);
    n=2^(nextpow2(n)+1);
end
fprintf('\nSize of the matrix is %d x %d\n',m,n);

%% APPLY WINDOW
if strcmp(winType,'none') ~= 1
    fprintf('\nApplying window ...\n')
    Matrix=ApplyWindowTo2DMatrix(Matrix,2,winType,0.1);
end

%% CALCULATE 2D FFT
fprintf('\nCalculating 2D FFT ...\n')
FFTMatrix=fft2(Matrix,m,n);
FFTMatrix=abs(FFTMatrix);               % magnitude
FFTMatrix=(fftshift(FFTMatrix));
FFTMatrix=fliplr(FFTMatrix');

%% X and Y AXIS VALUES
Freq=(-floor(m/2)+1:floor(m/2))'*fs/m;
X=(-floor(n/2):floor(n/2)-1)'*spatialFreq/n;

%% SCALE THE SIGNAL
fprintf('\nScaling the signal ...\n')
switch power
    case 'quad'
        FFTMatrix=FFTMatrix.^2;             % power spectrum
    case 'log'
        FFTMatrix=LogCompress(FFTMatrix,1024,0); % LOG COMPRESSION
end

%% NORMALIZE THE FFT MATRIX IN FREQUENCY SLICE
fprintf('\nNormalizing ...\n')
if normalizeOn == 1 
    for i=1:size(FFTMatrix,2)
        FFTMatrix(:,i)=FFTMatrix(:,i)./max(FFTMatrix(:,i)); % normalize by column 
    end
end

%% INVERSION FROM f-k SPACE TO PHASE VELOCITY
fprintf('\nPeak finding in f-k plane ...\n')
VelMatrix2DFFT=nan(size(FFTMatrix,2)/2,nPeaks);  % Create the matrix with velocities
WvnMatrix=nan(size(VelMatrix2DFFT));
Threshold=alpha*max(max(FFTMatrix(end/2:end,end/2:end)));
for i=1:size(FFTMatrix,2)/2
   [~,indices]=findpeaks(FFTMatrix(end/2+1:end,end/2+i),'MINPEAKHEIGHT',Threshold,...
         'MINPEAKDISTANCE',MinPeakDist,'NPEAKS',nPeaks,'SORTSTR','descend');
    TempWvn=X(end/2+indices);                   % calculate wavenumber
    WvnMatrix(i,:)=padarray(TempWvn,[nPeaks-length(TempWvn) 0],nan,'post');
end

%% CONVERSION WAVENUMBER TO PHASE VELOCITY
fprintf('\nConverting wavenumbers to phase velocity ...\n')
if smoothVel == 1
    WvnMatrix=SortDispCurves(WvnMatrix,[0 50 500],[900 900 900],35,1,...
        FreqLims./1000,[0 700]);
end
for i=1:size(FFTMatrix,2)/2
    VelMatrix2DFFT(i,:)=Freq(end/2+i)./WvnMatrix(i,:);    % calculate phase velocity
end

%% TAKE ONLY THE FORWARD PROPAGATING WAVES
fprintf('\nDeleting the backward propagating waves ...\n')
VelMatrix2DFFT(VelMatrix2DFFT<0)=nan;

%% SWITCH FREQUENCY UNITS
switch freqUnit
    case 'kHz'
        Freq=Freq./1000;
        FreqLims=FreqLims./1000;
    case 'MHz'
        Freq=Freq./1e6;
        FreqLims=FreqLims./1e6;
end

%% PLOTTING
if figOn == 1
    fprintf('\nPlotting ...\n')
    figure('units','normalized','outerposition',[0 0 1 1/2]) 
    subplot(1,2,1)
    set(gcf,'Renderer','Zbuffer')
    imagesc(Freq,X,FFTMatrix);
%     caxis([0 1])
    shading interp
    axis xy
%     axis tight
    colormap(gray)
    xlabel(strcat('Frequency [',freqUnit,']'),'FontSize',14)
    ylabel(strcat('Wavenumber [m^{-1}]'),'FontSize',14)
    xlim(FreqLims)
    ylim(LimsY)
    colorbar
    
%     figure;
    subplot(1,2,2)
    set(gcf,'Renderer','Zbuffer')
    imagesc(Freq,X,FFTMatrix);
    colormap(gray)
%     caxis([0 1])
    hold on
    plot(Freq(end/2+1:end),fliplr(WvnMatrix),'*')
    shading interp
    axis xy
    axis tight
    xlabel(strcat('Frequency [',freqUnit,']'),'FontSize',14)
    ylabel(strcat('Wavenumber [m^{-1}]'),'FontSize',14)
    xlim(FreqLims)
    ylim(LimsY)
%     colorbar
   
    figure;
%     subplot(2,2,4)
    plot(Freq(end/2+1:end),VelMatrix2DFFT,'r*');
    xlabel(strcat('Frequency [',freqUnit,']'),'FontSize',14)
    ylabel(strcat('Phase velocity [ms^{-1}]'),'FontSize',14)
    xlim(FreqLims)
    ylim([0 9000])
    tit=sprintf('FFT2D Dispersion curves\n scale: %s, winType: %s, fs: %.1f [MHz], spatialStep: %.2f [mm], alpha: %.2f, zeropad: %g',...
        power, winType,fs/1e6,spatialStep*1000,alpha,zeroPad);
    suptitle(tit);
end

%% STORE DISPMATRIX AND FREQ
if saveData == 1
    fprintf('\nSaving ...\n')
    save('DispData2DFFT','WvnMatrix','Freq','X','VelMatrix2DFFT','FFTMatrix');
end
