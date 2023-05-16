function [objRange,objSpeed,objAngle]=L3dataProcess(Data_dec)
%% 雷达参数
Tx_Number = 2;               %发射天线数
Rx_Number = 4;               %接收天线数
Range_Number = 128;          %距离点数
Doppler_Number = 128;        %多普勒通道数
NChirp = Doppler_Number;     %1帧数据的chirp个数
NSample = Range_Number;      %每个chirp ADC采样数
Fs = 5.209e6;                %采样频率 
c = 3.0e8;                   %光速
startFreq = 77e9;            %起始频率 
freqSlope = 60e12;           %chirp的斜率 
lambda=c/startFreq;          %雷达信号波长
Tc = 144e-6;                 %chirp周期 
Q=180;                       %角度FFT
d=lambda/2;                  %天线阵列间距

%% 数据读取、拆分、组合
Data_dec=hex2dec(Data_dec);
Data_zuhe=zeros(1,Tx_Number*Rx_Number*Doppler_Number*Range_Number*2); %建立计算存储数据的空矩阵

for i=1:Tx_Number*Rx_Number*Doppler_Number*Range_Number*2
    
    Data_zuhe(i) = Data_dec((i-1)*2+1)+Data_dec((i-1)*2+2)*256;%两个字节组成一个数，第二个字节乘以256相当于左移8位。
    if(Data_zuhe(i)>32767)
        Data_zuhe(i) = Data_zuhe(i) - 65536;  %限制幅度
    end
    
end

%% 打印全部的实虚数据
Re_Data_All=zeros(1,Range_Number*Doppler_Number*Tx_Number*Rx_Number); %建立计算存储数据的空矩阵
Im_Data_All=zeros(1,Range_Number*Doppler_Number*Tx_Number*Rx_Number); %建立计算存储数据的空矩阵

% 虚部实部分解
for i=1:Tx_Number*Rx_Number*Doppler_Number*Range_Number
    Im_Data_All(i) = Data_zuhe((i-1)*2+1);
    Re_Data_All(i) = Data_zuhe((i-1)*2+2);
end
figure();
plot(Re_Data_All);
%% 虚部+实部数据重组得到复信号 
ReIm_Data_All =complex(Re_Data_All,Im_Data_All);

%% 分放数据
ADC_Data1=zeros(Doppler_Number,Range_Number,Rx_Number); %建立计算存储数据的空矩阵
ADC_Data2=zeros(Doppler_Number,Range_Number,Rx_Number); %建立计算存储数据的空矩阵

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
% ADC_Data=zeros(Doppler_Number,Range_Number,Tx_Number*Rx_Number); %建立计算存储数据的空矩阵
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
%% 静态杂波滤除 相量均值相消
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

%% 直流分量抑制
fft2d(:,1,:)=0;

%% 3D FFT
fft3d = zeros(Doppler_Number,Range_Number,Q);
for qq =1:Doppler_Number
    for kk=1:Range_Number
        fft3d(qq,kk,:) = fftshift(fft(fft2d(qq,kk,:),Q));
    end
end

%% 绘制距离FFT结果
FFT1_mag=abs(fft1d(:,:,1));
figure();
mesh(FFT1_mag);
xlabel('采样点数');ylabel('脉冲数');zlabel('幅度');
xlim([0 NSample]);ylim([0 NChirp]);
title('距离维FFT结果');

%% 绘制2维FFT结果
FFT2_mag=abs(fft2d(:,:,1));
[X,Y] = meshgrid(c*(0:NSample-1)*Fs/2/freqSlope/NSample,...
    (-NChirp/2:NChirp/2 - 1)*lambda/Tc/NChirp/2);
figure();
mesh(X,Y,FFT2_mag);
xlim([0 c*(NSample-1)*Fs/2/freqSlope/NSample]);
ylim([(-NChirp/2)*lambda/Tc/NChirp/2 (NChirp/2 - 1)*lambda/Tc/NChirp/2]);
xlabel('距离(m)');ylabel('速度(m/s)');zlabel('幅度');
title('2D-FFT结果');

%% 绘制角度维FFT结果
FFT3_mag=reshape(abs(fft3d(1,:,:)),Range_Number,Q).';
[R,Z] = meshgrid(c*(0:NSample-1)*Fs/2/freqSlope/NSample,...
    asin((-Q/2:Q/2 - 1)*lambda/d/Q)*180/pi);
figure();
mesh(R,Z,FFT3_mag);
xlim([0 c*(NSample-1)*Fs/2/freqSlope/NSample]);
ylim([asin((-Q/2)*lambda/d/Q)*180/pi asin((Q/2 - 1)*lambda/d/Q)*180/pi]);
xlabel('距离维(m)');ylabel('角度维(°)');zlabel('幅度');
title('角度维FFT结果');

%% 计算峰值位置
fft3d=abs(fft3d);
pink=max(fft3d(:));
[row,col,pag]=ind2sub(size(fft3d),find(fft3d==pink));

%% 计算目标距离、速度、角度
fb = ((col-1)*Fs)/NSample;          %差拍频率
fd = (row-NChirp/2-1)/(NChirp*Tc);  %多普勒频率
fw = (pag-Q/2-1)/Q;                 %空间频率
R = c*fb/2/freqSlope;               %距离公式
v =lambda*fd/2;                     %速度公式
theta = asin(fw*lambda/d);          %角度公式
angle = theta*180/pi;
fprintf('目标距离： %f m\n',R);
fprintf('目标速度： %f m/s\n',v);
fprintf('目标角度： %f°\n',angle);

%% 速度维 CFAR
% %进行1D CFAR ，平均单元恒虚警率
% %缺点：多目标遮掩、杂波边缘性能欠佳
% %优点：损失率最小
% %参考单元：12
% %保护单元：2
% %虚警概率：
% %门限值：

Tv=6;Pv=4;PFAv = 0.001; 

dopplerDimCfarThresholdMap = zeros(size(FFT2_mag));  %创建一个二维矩阵存放doppler维cfar后的结果
dopplerDimCfarResultMap    = zeros(size(FFT2_mag));

for i=1:Range_Number
   dopplerDim =  reshape(FFT2_mag(:,i),1,Doppler_Number);       %变成一行数据
   [cfar1D_Arr,threshold] = ac_cfar1D(Tv,Pv,PFAv,dopplerDim);   %进行1D cfar
   dopplerDimCfarResultMap(:,i) = cfar1D_Arr.'; 
   dopplerDimCfarThresholdMap(:,i) = threshold.';
end

%% 距离维 CFAR
%沿着doppler维度方向寻找在doppler维cfar判决后为1的结果
saveMat = zeros(size(FFT2_mag));
for range = 1:Range_Number
    indexArr = find(dopplerDimCfarResultMap(:,range));
    objDopplerArr = [indexArr;zeros(Doppler_Number - length(indexArr),1)];   %补充长度
    saveMat(:,range) = objDopplerArr; %保存doppler下标
end
% 保存有物体的doppler坐标
objDopplerIndex = unique(saveMat);  % unqiue是不重复的返回数组中的数
objDopplerIndex(objDopplerIndex==0)=[]; %排除数组中的0

% 根据之前doppler维的cfar结果对应的下标objDopplerIndex，对相应的速度进行range维度的CFAR
Tr=6;Pr=4;PFAr = 0.003;

rangeDimCfarThresholdMap = zeros(size(FFT2_mag));  %创建一个二维矩阵存放range维cfar后的结果
rangeDimCfarResultMap = zeros(size(FFT2_mag));

for i = 1:length(objDopplerIndex)
    %根据速度下标进行range CFAR
    j = objDopplerIndex(i);     % 获得物体所在的行
    rangeDim = reshape(FFT2_mag(j, :),1,Range_Number);  %变成一行数据
    % tip 这个PFA很迷啊,如果设置的低一些,在进行分支聚集的时候,可能的结果是没有检测到物体
    % 因为在进行rangeCFAR的时候,把附近的最大值给滤掉了,那这样在进行峰值聚集的时候,判决结果对应的最大值一直是小的
        
    [cfar1D_Avv,threshold] = ac_cfar1D(Tr,Pr,PFAr,rangeDim);  %进行1D cfar
    rangeDimCfarResultMap(j,:) = cfar1D_Avv; 
    rangeDimCfarThresholdMap(j,:) = threshold;  
end

%% 绘制CFAR的判决结果
figure();
mesh(X,Y,(rangeDimCfarResultMap));
xlabel('距离(m)');ylabel('速度(m/s)');zlabel('信号幅值');
title('rangeCFAR之后判决结果(峰值聚集前)');
xlim([0 c*(NSample-1)*Fs/2/freqSlope/NSample]);
ylim([(-NChirp/2)*lambda/Tc/NChirp/2 (NChirp/2 - 1)*lambda/Tc/NChirp/2]);

%% 峰值检索
[objDprIdx,objRagIdx] = peakFocus(rangeDimCfarResultMap,FFT2_mag);% 峰值检索 距离和速度的ID号
objAglIdx=zeros(1,length(objDprIdx));%角度的ID号
for i=1:length(objDprIdx)
    row=objDprIdx(i);
    col=objRagIdx(i);
    [~,id]=max(fft3d(row,col,:));
    objAglIdx(i)=id;
end

%% 计算目标距离、速度、角度，输出
objSpeed = ( objDprIdx.' - NChirp/2 - 1)*lambda/Tc/NChirp/2;
objRange = single(c*(objRagIdx.'-1)*Fs/2/freqSlope/NSample);
objAngle = asin((objAglIdx.' - Q/2 - 1)*lambda/d/Q)*180/pi;
end

%% 功能函数
function [cfar1D_Arr,threshold] = ac_cfar1D(NTrain,NGuard,PFA,inputArr)
    cfar1D_Arr = zeros(size(inputArr));
    threshold = zeros(size(inputArr));

    totalNTrain = 2*(NTrain);
    a = totalNTrain*((PFA^(-1/totalNTrain))-1);
    %求平均值
    for i = NTrain+NGuard+1:length(inputArr)-NTrain-NGuard
        avg = mean([inputArr((i-NTrain-NGuard):(i-NGuard-1))...
            inputArr((i+NGuard+1):(i+NTrain+NGuard))]);
        threshold(1,i) = a.*avg;
        %根据threshold比较
        if(inputArr(i) < threshold(i))
            cfar1D_Arr(i) = 0;
        else
            cfar1D_Arr(i) = 1;
        end
    end
end

% 输入: inputCfarResMat - 进行峰值聚焦的二维矩阵,即进行range维CFAR判决后得到的结果矩阵
% 输出: row - 物体的行坐标(对应速度)
% column - 物体的列坐标(对应距离)
function [row,column] = peakFocus(inputCfarResMat,FFT2_mag)
    j = 1;
    row = zeros([1 256]);
    column = zeros([1 256]);
    [d,r] = find(inputCfarResMat==1);   %寻找进行range维cfar后的判决为1的坐标
    for i = 1 : length(d)
        peakRow = d(i);
        peakColumn = r(i);
        peak = FFT2_mag(peakRow,peakColumn);  %待验证的峰值
        % 在附近的3*3矩阵中的数进行比较,如果中间的数是最大值,就判定为1  
        % 根据之前进行的2次cfar,因为有TrainCell和GuardCell，所以不会碰到边缘
        tempArr =[FFT2_mag(peakRow-1,peakColumn-1) , FFT2_mag(peakRow-1,peakColumn) ,  FFT2_mag(peakRow-1,peakColumn+1), ...
                  FFT2_mag(peakRow,peakColumn-1)   ,                     peak                             , FFT2_mag(peakRow,peakColumn+1), ...
                  FFT2_mag(peakRow+1,peakColumn-1) , FFT2_mag(peakRow+1,peakColumn) ,  FFT2_mag(peakRow+1,peakColumn+1)] ;    
        truePeak = max(tempArr);     % 寻找最大值
        if(truePeak == peak)         %如果中间的是最大值就保存当前的坐标
            row(j) = peakRow;
            column(j) = peakColumn;
            j = j+1;
        end
    end
    row(row==0)=[]; %去掉后面的0
    column(column==0)=[];
end
