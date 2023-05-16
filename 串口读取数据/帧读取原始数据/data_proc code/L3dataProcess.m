function [objRange,objSpeed,objAngle]=L3dataProcess(Data_dec)
%% �״����
Tx_Number = 2;               %����������
Rx_Number = 4;               %����������
Range_Number = 128;          %�������
Doppler_Number = 128;        %������ͨ����
NChirp = Doppler_Number;     %1֡���ݵ�chirp����
NSample = Range_Number;      %ÿ��chirp ADC������
Fs = 5.209e6;                %����Ƶ�� 
c = 3.0e8;                   %����
startFreq = 77e9;            %��ʼƵ�� 
freqSlope = 60e12;           %chirp��б�� 
lambda=c/startFreq;          %�״��źŲ���
Tc = 144e-6;                 %chirp���� 
Q=180;                       %�Ƕ�FFT
d=lambda/2;                  %�������м��

%% ���ݶ�ȡ����֡����
Data_dec=hex2dec(Data_dec);
Data_zuhe=zeros(1,Tx_Number*Rx_Number*Doppler_Number*Range_Number*2); %��������洢���ݵĿվ���

for i=1:Tx_Number*Rx_Number*Doppler_Number*Range_Number*2
    
    Data_zuhe(i) = Data_dec((i-1)*2+1)+Data_dec((i-1)*2+2)*256;%�����ֽ����һ�������ڶ����ֽڳ���256�൱������8λ��
    if(Data_zuhe(i)>32767)
        Data_zuhe(i) = Data_zuhe(i) - 65536;  %���Ʒ���
    end
    
end

%% ��ӡȫ����ʵ������
Re_Data_All=zeros(1,Range_Number*Doppler_Number*Tx_Number*Rx_Number); %��������洢���ݵĿվ���
Im_Data_All=zeros(1,Range_Number*Doppler_Number*Tx_Number*Rx_Number); %��������洢���ݵĿվ���

% �鲿ʵ���ֽ�
for i=1:Tx_Number*Rx_Number*Doppler_Number*Range_Number
    Im_Data_All(i) = Data_zuhe((i-1)*2+1);
    Re_Data_All(i) = Data_zuhe((i-1)*2+2);
end
figure();
plot(Re_Data_All);
%% �鲿+ʵ����������õ����ź� 
ReIm_Data_All =complex(Re_Data_All,Im_Data_All);

%% �ַ�����
ADC_Data1=zeros(Doppler_Number,Range_Number,Rx_Number); %��������洢���ݵĿվ���
ADC_Data2=zeros(Doppler_Number,Range_Number,Rx_Number); %��������洢���ݵĿվ���

RX_num=Range_Number;
chirp_num=RX_num*Rx_Number;
TX_num=chirp_num*Doppler_Number;
for i=1:Doppler_Number
    for j=1:Range_Number
        for k=1:Rx_Number
            ADC_Data1(i,j,k) = ReIm_Data_All((i-1)*chirp_num+(k-1)*RX_num+j);
            ADC_Data2(i,j,k) = ReIm_Data_All(TX_num+(i-1)*chirp_num+(k-1)*RX_num+j);
        end
    end
end
ADC_Data = cat(3,ADC_Data1(:,:,1:Rx_Number), ADC_Data2(:,:,1:Rx_Number));
% ADC_Data=zeros(Doppler_Number,Range_Number,Tx_Number*Rx_Number); %��������洢���ݵĿվ���
% chirp_num=Range_Number;
% RX_num=chirp_num*Doppler_Number;
% for i=1:Doppler_Number
%     for j=1:Range_Number
%         for k=1:Tx_Number*Rx_Number
%             ADC_Data(i,j,k) = ReIm_Data_All((i-1)*chirp_num+(k-1)*RX_num+j);
%         end
%     end
% end

%% 1D FFT
% fft1d= zeros(Doppler_Number,Range_Number,Tx_Number*Rx_Number);
% for qq =1:Tx_Number*Rx_Number
%     for chirp_fft=1:Doppler_Number 
%         fft1d(chirp_fft,:,qq) = fft(ADC_Data(chirp_fft,:,qq));
%     end
% end
fft1d=ADC_Data;
%% ��̬�Ӳ��˳� ������ֵ����
fft1d_jingtai =zeros(Doppler_Number,Range_Number,Tx_Number*Rx_Number);
for n=1:Tx_Number*Rx_Number
    avg = sum(fft1d(:,:,n))/Doppler_Number;
    for m=1:Doppler_Number
        fft1d_jingtai(m,:,n) = fft1d(m,:,n)-avg;
    end
end
fft1d =fft1d_jingtai;

%% 2D FFT 
fft2d= zeros(Doppler_Number,Range_Number,Tx_Number*Rx_Number);
for kk=1:Tx_Number*Rx_Number
    for chirp_fft=1:Range_Number 
        fft2d(:,chirp_fft,kk) =fftshift( fft(fft1d(:,chirp_fft,kk)));  
    end
end

%% ֱ����������
fft2d(:,1,:)=0;

%% 3D FFT
fft3d = zeros(Doppler_Number,Range_Number,Q);
for qq =1:Doppler_Number
    for kk=1:Range_Number
        fft3d(qq,kk,:) = fftshift(fft(fft2d(qq,kk,:),Q));
    end
end

%% ���ƾ���FFT���
FFT1_mag=abs(fft1d(:,:,1));
figure();
mesh(FFT1_mag);
xlabel('��������');ylabel('������');zlabel('����');
xlim([0 NSample]);ylim([0 NChirp]);
title('����άFFT���');

%% ����2άFFT���
FFT2_mag=abs(fft2d(:,:,1));
[X,Y] = meshgrid(c*(0:NSample-1)*Fs/2/freqSlope/NSample,...
    (-NChirp/2:NChirp/2 - 1)*lambda/Tc/NChirp/2);
figure();
mesh(X,Y,FFT2_mag);
xlim([0 c*(NSample-1)*Fs/2/freqSlope/NSample]);
ylim([(-NChirp/2)*lambda/Tc/NChirp/2 (NChirp/2 - 1)*lambda/Tc/NChirp/2]);
xlabel('����(m)');ylabel('�ٶ�(m/s)');zlabel('����');
title('2D-FFT���');

%% ���ƽǶ�άFFT���
FFT3_mag=reshape(abs(fft3d(1,:,:)),Range_Number,Q).';
[R,Z] = meshgrid(c*(0:NSample-1)*Fs/2/freqSlope/NSample,...
    asin((-Q/2:Q/2 - 1)*lambda/d/Q)*180/pi);
figure();
mesh(R,Z,FFT3_mag);
xlim([0 c*(NSample-1)*Fs/2/freqSlope/NSample]);
ylim([asin((-Q/2)*lambda/d/Q)*180/pi asin((Q/2 - 1)*lambda/d/Q)*180/pi]);
xlabel('����ά(m)');ylabel('�Ƕ�ά(��)');zlabel('����');
title('�Ƕ�άFFT���');

%% �����ֵλ��
fft3d=abs(fft3d);
pink=max(fft3d(:));
[row,col,pag]=ind2sub(size(fft3d),find(fft3d==pink));

%% ����Ŀ����롢�ٶȡ��Ƕ�
fb = ((col-1)*Fs)/NSample;          %����Ƶ��
fd = (row-NChirp/2-1)/(NChirp*Tc);  %������Ƶ��
fw = (pag-Q/2-1)/Q;                 %�ռ�Ƶ��
R = c*fb/2/freqSlope;               %���빫ʽ
v =lambda*fd/2;                     %�ٶȹ�ʽ
theta = asin(fw*lambda/d);          %�Ƕȹ�ʽ
angle = theta*180/pi;
fprintf('Ŀ����룺 %f m\n',R);
fprintf('Ŀ���ٶȣ� %f m/s\n',v);
fprintf('Ŀ��Ƕȣ� %f��\n',angle);

%% �ٶ�ά CFAR
% %����1D CFAR ��ƽ����Ԫ���龯��
% %ȱ�㣺��Ŀ�����ڡ��Ӳ���Ե����Ƿ��
% %�ŵ㣺��ʧ����С
% %�ο���Ԫ��12
% %������Ԫ��2
% %�龯���ʣ�
% %����ֵ��

Tv=6;Pv=4;PFAv = 0.001; 

dopplerDimCfarThresholdMap = zeros(size(FFT2_mag));  %����һ����ά������dopplerάcfar��Ľ��
dopplerDimCfarResultMap    = zeros(size(FFT2_mag));

for i=1:Range_Number
   dopplerDim =  reshape(FFT2_mag(:,i),1,Doppler_Number);       %���һ������
   [cfar1D_Arr,threshold] = ac_cfar1D(Tv,Pv,PFAv,dopplerDim);   %����1D cfar
   dopplerDimCfarResultMap(:,i) = cfar1D_Arr.'; 
   dopplerDimCfarThresholdMap(:,i) = threshold.';
end

%% ����ά CFAR
%����dopplerά�ȷ���Ѱ����dopplerάcfar�о���Ϊ1�Ľ��
saveMat = zeros(size(FFT2_mag));
for range = 1:Range_Number
    indexArr = find(dopplerDimCfarResultMap(:,range));
    objDopplerArr = [indexArr;zeros(Doppler_Number - length(indexArr),1)];   %���䳤��
    saveMat(:,range) = objDopplerArr; %����doppler�±�
end
% �����������doppler����
objDopplerIndex = unique(saveMat);  % unqiue�ǲ��ظ��ķ��������е���
objDopplerIndex(objDopplerIndex==0)=[]; %�ų������е�0

% ����֮ǰdopplerά��cfar�����Ӧ���±�objDopplerIndex������Ӧ���ٶȽ���rangeά�ȵ�CFAR
Tr=6;Pr=4;PFAr = 0.003;

rangeDimCfarThresholdMap = zeros(size(FFT2_mag));  %����һ����ά������rangeάcfar��Ľ��
rangeDimCfarResultMap = zeros(size(FFT2_mag));

for i = 1:length(objDopplerIndex)
    %�����ٶ��±����range CFAR
    j = objDopplerIndex(i);     % ����������ڵ���
    rangeDim = reshape(FFT2_mag(j, :),1,Range_Number);  %���һ������
    % tip ���PFA���԰�,������õĵ�һЩ,�ڽ��з�֧�ۼ���ʱ��,���ܵĽ����û�м�⵽����
    % ��Ϊ�ڽ���rangeCFAR��ʱ��,�Ѹ��������ֵ���˵���,�������ڽ��з�ֵ�ۼ���ʱ��,�о������Ӧ�����ֵһֱ��С��
        
    [cfar1D_Avv,threshold] = ac_cfar1D(Tr,Pr,PFAr,rangeDim);  %����1D cfar
    rangeDimCfarResultMap(j,:) = cfar1D_Avv; 
    rangeDimCfarThresholdMap(j,:) = threshold;  
end

%% ����CFAR���о����
figure();
mesh(X,Y,(rangeDimCfarResultMap));
xlabel('����(m)');ylabel('�ٶ�(m/s)');zlabel('�źŷ�ֵ');
title('rangeCFAR֮���о����(��ֵ�ۼ�ǰ)');
xlim([0 c*(NSample-1)*Fs/2/freqSlope/NSample]);
ylim([(-NChirp/2)*lambda/Tc/NChirp/2 (NChirp/2 - 1)*lambda/Tc/NChirp/2]);

%% ��ֵ����
[objDprIdx,objRagIdx] = peakFocus(rangeDimCfarResultMap,FFT2_mag);% ��ֵ���� ������ٶȵ�ID��
objAglIdx=zeros(1,length(objDprIdx));%�Ƕȵ�ID��
for i=1:length(objDprIdx)
    row=objDprIdx(i);
    col=objRagIdx(i);
    [~,id]=max(fft3d(row,col,:));
    objAglIdx(i)=id;
end

%% ����Ŀ����롢�ٶȡ��Ƕȣ����
objSpeed = ( objDprIdx.' - NChirp/2 - 1)*lambda/Tc/NChirp/2;
objRange = single(c*(objRagIdx.'-1)*Fs/2/freqSlope/NSample);
objAngle = asin((objAglIdx.' - Q/2 - 1)*lambda/d/Q)*180/pi;
end

%% ���ܺ���
function [cfar1D_Arr,threshold] = ac_cfar1D(NTrain,NGuard,PFA,inputArr)
    cfar1D_Arr = zeros(size(inputArr));
    threshold = zeros(size(inputArr));

    totalNTrain = 2*(NTrain);
    a = totalNTrain*((PFA^(-1/totalNTrain))-1);
    %��ƽ��ֵ
    for i = NTrain+NGuard+1:length(inputArr)-NTrain-NGuard
        avg = mean([inputArr((i-NTrain-NGuard):(i-NGuard-1))...
            inputArr((i+NGuard+1):(i+NTrain+NGuard))]);
        threshold(1,i) = a.*avg;
        %����threshold�Ƚ�
        if(inputArr(i) < threshold(i))
            cfar1D_Arr(i) = 0;
        else
            cfar1D_Arr(i) = 1;
        end
    end
end

% ����: inputCfarResMat - ���з�ֵ�۽��Ķ�ά����,������rangeάCFAR�о���õ��Ľ������
% ���: row - �����������(��Ӧ�ٶ�)
% column - �����������(��Ӧ����)
function [row,column] = peakFocus(inputCfarResMat,FFT2_mag)
    j = 1;
    row = zeros([1 256]);
    column = zeros([1 256]);
    [d,r] = find(inputCfarResMat==1);   %Ѱ�ҽ���rangeάcfar����о�Ϊ1������
    for i = 1 : length(d)
        peakRow = d(i);
        peakColumn = r(i);
        peak = FFT2_mag(peakRow,peakColumn);  %����֤�ķ�ֵ
        % �ڸ�����3*3�����е������бȽ�,����м���������ֵ,���ж�Ϊ1  
        % ����֮ǰ���е�2��cfar,��Ϊ��TrainCell��GuardCell�����Բ���������Ե
        tempArr =[FFT2_mag(peakRow-1,peakColumn-1) , FFT2_mag(peakRow-1,peakColumn) ,  FFT2_mag(peakRow-1,peakColumn+1), ...
                  FFT2_mag(peakRow,peakColumn-1)   ,                     peak                             , FFT2_mag(peakRow,peakColumn+1), ...
                  FFT2_mag(peakRow+1,peakColumn-1) , FFT2_mag(peakRow+1,peakColumn) ,  FFT2_mag(peakRow+1,peakColumn+1)] ;    
        truePeak = max(tempArr);     % Ѱ�����ֵ
        if(truePeak == peak)         %����м�������ֵ�ͱ��浱ǰ������
            row(j) = peakRow;
            column(j) = peakColumn;
            j = j+1;
        end
    end
    row(row==0)=[]; %ȥ�������0
    column(column==0)=[];
end
