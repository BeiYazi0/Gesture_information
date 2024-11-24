function [fft1d,fft2d,fft3d]=Three_dfft(fdata)

Range_Number=512;      %距离向FFT点数
Doppler_Number=128;    %多普勒向FFT点数
Rx_Number=4;           %RX天线通道数
Tx_Number=2;           %TX天线通道数
Q = 180;               %角度FFT

%% 数据拆分、组合
% 虚部+实部数据重组得到复信号 
fileSize = size(fdata, 1);
lvds_data = zeros(1, fileSize/2);
count = 1;
for i=1:4:fileSize-5
   lvds_data(1,count) = fdata(i) + 1i*fdata(i+2); 
   lvds_data(1,count+1) = fdata(i+1)+1i*fdata(i+3); 
   count = count + 2;
end

% 将回波数据排列成3维矩阵(速度维、距离维和天线维)
ADC_Data=zeros(Doppler_Number,Range_Number,Rx_Number*Tx_Number); %建立计算存储数据的空矩阵

xSize=Range_Number;
chirpSize=xSize*Rx_Number*Tx_Number;
for i=1:Doppler_Number
    for k=1:Rx_Number*Tx_Number
        ADC_Data(i,:,k) = lvds_data((i-1)*chirpSize+(k-1)*xSize+1:(i-1)*chirpSize+k*xSize);
    end
end

%% 距离FFT
range_win = hamming(Range_Number).';   %加海明窗
range_profile = zeros(Doppler_Number,Range_Number,Rx_Number*Tx_Number);
for k=1:Rx_Number*Tx_Number
   for m=1:Doppler_Number
      temp=ADC_Data(m,:,k).*range_win;    %加窗函数
      temp_fft=fft(temp,Range_Number);    %对每个chirp做FFT
      range_profile(m,:,k)=temp_fft;
   end
end

%% 静态杂波滤除 相量均值相消
fft1d_sub =zeros(Doppler_Number,Range_Number,Tx_Number*Rx_Number);
for n=1:Tx_Number*Rx_Number
    avg = sum(range_profile(:,:,n))/Doppler_Number;
    for m=1:Doppler_Number
        fft1d_sub(m,:,n) = range_profile(m,:,n)-avg;
    end
end
range_profile =fft1d_sub;

%% 多普勒FFT
doppler_win = hamming(Doppler_Number);                %加海明窗
speed_profile = zeros(Doppler_Number,Range_Number,Rx_Number*Tx_Number);
for k=1:Rx_Number*Tx_Number
    for n=1:Range_Number
      temp=range_profile(:,n,k).*doppler_win;         %加窗函数
      temp_fft=fftshift(fft(temp,Doppler_Number));    %对rangeFFT结果进行FFT
      speed_profile(:,n,k)=temp_fft;  
    end
end

%% 直流分量抑制
range_profile(:,1,:)=0;
speed_profile(:,1,:)=0;

%% 角度FFT
angle_profile = zeros(Doppler_Number,Range_Number,Q);
for n=1:Range_Number   %range
    for m=1:Doppler_Number   %chirp
      temp=speed_profile(m,n,:);    
      temp_fft=fftshift(fft(temp,Q));    %对2D FFT结果进行Q点FFT
      angle_profile(m,n,:)=temp_fft;  
    end
end

%% 绘制距离FFT结果
% FFT1_mag=abs(range_profile(:,:,1));
% figure();
% mesh(FFT1_mag);
% xlabel('采样点数');ylabel('脉冲数');zlabel('幅度');
% xlim([0 NSample]);ylim([0 NChirp]);
% title('距离维FFT结果');

%% 绘制2D-FFT结果
% FFT2_mag=abs(speed_profile(:,:,1));
% [X,Y] = meshgrid(c*(0:NSample-1)*Fs/2/freqSlope/NSample,...
%     (-NChirp/2:NChirp/2 - 1)*lambda/Tc/NChirp/2);
% figure();
% mesh(X,Y,FFT2_mag);
% xlim([0 c*(NSample-1)*Fs/2/freqSlope/NSample]);
% ylim([(-NChirp/2)*lambda/Tc/NChirp/2 (NChirp/2 - 1)*lambda/Tc/NChirp/2]);
% xlabel('距离(m)');ylabel('速度(m/s)');zlabel('幅度');
% title('2D-FFT结果');

%% 绘制角度维FFT结果
% FFT3_mag=reshape(abs(angle_profile(1,:,:)),Range_Number,Q).';
% [R,Z] = meshgrid(c*(0:NSample-1)*Fs/2/freqSlope/NSample,...
%     asin((-Q/2:Q/2 - 1)*lambda/d/Q)*180/pi);
% figure();
% mesh(R,Z,FFT3_mag);
% xlim([0 c*(NSample-1)*Fs/2/freqSlope/NSample]);
% ylim([asin((-Q/2)*lambda/d/Q)*180/pi asin((Q/2 - 1)*lambda/d/Q)*180/pi]);
% xlabel('距离维(m)');ylabel('角度维(°)');zlabel('幅度');
% title('角度维FFT结果');

%% 结果获取
fft1d=sum(range_profile(:,:,1)).'/Doppler_Number;
[~,col]=max(fft1d);
fft2d=speed_profile(:,col,1);
[~,row]=max(fft1d);
fft3d=angle_profile(row,col,:);
end