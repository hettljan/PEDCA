%% LOAD THE DISPERSION DATA
% load('DispData2DFFT.mat');
% Vel1=VelMatrix2DFFT;
% load('DispDataMPM.mat');
% Vel2=VelMatrixMPM;
S0=load('S0.mat');
Vel1=S0.VelMatrix2DFFT;
A0=load('DispDataMPM.mat');
Vel2=A0.VelMatrix2DFFT;
Vel1=Vel1(1:size(Vel2,1),:);              % Cut the 2D FFT data to the size of MPM 
FreqVec=Freq(end/2+1:end/2+size(Vel2,1)); % Define a corresponding freqeuncy vector

%% PARAMETERS
nCurv=size(Vel2,2);     % The total number of dispersion curves                           
nSize=100;              % Size of the neigborhood to search in for the ones
YLims=[0 8000];        % y-axis limits in ms^{-1}
freqUnit='kHz';         % Frequency units
FreqLims=[0 750];       % Frequency limits in freqUnits

%%
Mat1=zeros(10000,size(Vel1,1));
Mat2=zeros(10000,size(Vel2,1));
Mat=zeros(10000,size(Vel2,1));
for i=1:size(Vel1,1)
    for j=1:size(Vel1,2)
        if Vel1(i,j)<=10000 && round(Vel1(i,j))>0   % if Velocity is greater than 0 and smaller then 10k 
            Mat1(round(Vel1(i,j)),i)=1;
        end
    end
end
for i=1:size(Vel2,1)  
    for j=1:size(Vel2,2)
        if Vel2(i,j)<=10000 && round(Vel2(i,j))>0
            Mat2(round(Vel2(i,j)),i)=1;
        end
    end
end
[Row,Col] = find(Mat1);     % find the non-zero elements 
Row(Row<=nSize)=nan;
for i=1:size(Row)
    cont=1;
    k=-nSize;
    while isnan(Row(i))~=1 && cont == 1 && k<=nSize && Row(i)+k<=size(Mat2,1)  
        if Mat2(Row(i)+k,Col(i)) == 1       % search the -kSize to +kSize vicinity of Row(i) for ones
            newRow=round((Row(i)+(Row(i)+k))/2);     % new row index is an average of two
            Mat(newRow,Col(i))=1;
            cont=0;
        end
        k=k+1;
    end
end

%% INVERSION FROM BINARIZED IMAGE TO VELOCITY VECTOR
VelMatrix=nan(size(Vel1));
for i=1:size(Mat,2)
    k=1;
    count=1;
    while k<=size(Mat,1) && count<=nCurv
        if Mat(k,i)==1
            VelMatrix(i,count)=k;
            count=count+1;
        end
        k=k+1;
    end
end

%% PLOTTING
% figure
% subplot(1,2,1)
% % plot(repmat(FreqVec,nCurv,1),reshape(Vel1,size(Vel1,1)*nCurv,1),'bo');
% hold on
% plot(repmat(FreqVec,nCurv,1),reshape(Vel2,size(Vel2,1)*nCurv,1),'r*');
% xlim(FreqLims)
% ylim(YLims)
% xlabel(strcat('Frequency [',freqUnit,']'),'FontSize',14)
% ylabel(strcat('Phase velocity [ms^{-1}]'),'FontSize',14)
% legend('2D FFT', 'MPM')
% title('Raw DispCurves','FontSize',14)
% subplot(1,2,2)
figure
plot(FreqVec,VelMatrix,'b*')
xlabel(strcat('Frequency [',freqUnit,']'),'FontSize',14)
ylabel(strcat('Phase velocity [ms^{-1}]'),'FontSize',14)
xlim(FreqLims)
ylim(YLims)
title('Filtered disperson curve','FontSize',14)
