%function frame_dataproc()
close all;clear all;clc;
%% 雷达参数（使用mmWave Studio默认参数） 
c=3.0e8;               %光速
freqSlope=30e12;       %调频斜率
Fs=5e6;                %采样率
PDF=2500;
f0=77e9;               %初始频率
lambda=c/f0;           %雷达信号波长
d=lambda/2;            %天线阵列间距
Range_Number=512;      %采样点数/脉冲
NSample=Range_Number;  %距离向FFT点数
Doppler_Number=4000;    %每帧脉冲数
NChirp=Doppler_Number; %多普勒向FFT点数
Rx_Number=4;           %RX天线通道数
Tx_Number=2;           %TX天线通道数

%% 读取回波数据
fname='C:\ti\mmwave_studio_02_01_00_00\mmWaveStudio\PostProc\adc_data.bin';
fid = fopen(fname,'rb');
bufferSize=NSample*NChirp*Rx_Number*Tx_Number*2;
%16bits，复数形式(I/Q两路)，4RX,1TX,有符号16bit，小端模式
fdata = fread(fid,bufferSize,'int16');
fclose(fid);

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
ADC_Data=zeros(Range_Number,Doppler_Number,Rx_Number*Tx_Number); %建立计算存储数据的空矩阵

xSize=Range_Number;
chirpSize=xSize*Rx_Number*Tx_Number;
for i=1:Doppler_Number
    for j=1:Rx_Number*Tx_Number
        ADC_Data(:,i,j) = lvds_data((i-1)*chirpSize+(j-1)*xSize+1:(i-1)*chirpSize+j*xSize);
    end
end

%% 距离FFT
range_win = hamming(Range_Number);   %加海明窗
range_profile = zeros(Range_Number,Doppler_Number,Rx_Number*Tx_Number);
for j=1:Rx_Number*Tx_Number
   for m=1:Doppler_Number
      temp=ADC_Data(:,m,j).*range_win;    %加窗函数
      temp_fft=fft(temp,Range_Number);    %对每个chirp做FFT
      range_profile(:,m,j)=temp_fft;
   end
end

%% 静态杂波滤除 相量均值相消
% fft1d_sub =zeros(Range_Number,Doppler_Number,Tx_Number*Rx_Number);
% for n=1:Tx_Number*Rx_Number
%     avg = sum(range_profile(:,:,n),2)/Doppler_Number;
%     for m=1:Doppler_Number
%         fft1d_sub(:,m,n) = range_profile(:,m,n)-avg;
%     end
% end
% range_profile =fft1d_sub;

%% 指数加权平均滤波
B=0;
a=0.01;
for i=1:Doppler_Number
     range_profile(:,i,:)=range_profile(:,i,:)-B;
     B=B+a*range_profile(:,i,:);
end

%% 多普勒信息提取
range_max=max(range_profile(:,:,1));
[speed_profile, F, T]=DopplerProcess(range_max,PDF);
speed_profile=fftshift(speed_profile);
F=fftshift(F);
F(F>=PDF/2)=F(F>=PDF/2)-PDF;

%% 角度信息提取
angle_profile=zeros(361,Doppler_Number);
for i=1:Doppler_Number
    X=reshape(ADC_Data(:,i,:),Range_Number,[]).';
    angle_profile(:,i)=AngleProcess(X,d);
end

%% 绘制信息提取图
t=(0:Doppler_Number-1)/PDF;
figure();
R=c*(0:NSample-1)*Fs/2/freqSlope/NSample;
mesh(t,R,abs(range_profile(:,:,1)));
% shading interp;
ylim([0 0.8]);                    %限制坐标轴
view(2);
% set(gca,'position',[0 0 1 1]);    %图像填满
% oimg =getframe(1);
% img=imstand(oimg.cdata);          %图像标准化
% imwrite(img,'range.png');
figure();
V=F*lambda/2;
mesh(T,V,abs(speed_profile));
shading interp;
ylim([-1.2 1.2]);                 %限制坐标轴
view(2);
% set(gca,'position',[0 0 1 1]);    %图像填满
% oimg = getframe(2);
% img=imstand(oimg.cdata);          %图像标准化
% imwrite(img,'doppler.png');
figure();
A=linspace(-90,90,361);
mesh(t,A,abs(angle_profile));
shading interp;
ylim([-50 50]);                   %限制坐标轴
view(2);
% set(gca,'position',[0 0 1 1]);    %图像填满
% oimg = getframe(3);
% img=imstand(oimg.cdata);          %图像标准化
% imwrite(img,'angle.png');

%% 功能函数
function [speed_profile, F, T]=DopplerProcess(range_max,fs)
win_sz = 64;
nfft = win_sz;
nooverlap = win_sz - 4;
[speed_profile, F, T] = spectrogram(range_max, win_sz, nooverlap, nfft, fs, 'yaxis');
end

function rimg=imstand(oimg)
img=imresize(oimg,[224,224]);
i=min(min(min(img)));
j=max(max(max(img)));
rimg=img-i/(j-i);
end

function Pmusic=AngleProcess(X,dd)
K=size(X,2);

Rxx=X*X'/K;                       % 计算协方差矩阵
[EV,D]=eig(Rxx);                   %特征值分解
EVA=diag(D)';                      %将特征值矩阵对角线提取并转为一行
[~,I]=sort(EVA);                 %将特征值排序 从小到大
EV=fliplr(EV(:,I));                % 对应特征矢量排序
                 
% 遍历每个角度，计算空间谱
angle = zeros(361,1);
Pmusic = zeros(361,1);
derad = pi/180;
N = 8;               % 阵元个数        
M = 1;               % 信源数目
d=0:dd:(N-1)*dd;
En=EV(:,M+1:N);                   % 取矩阵的第M+1到N列组成噪声子空间
Hm=En*En';
for iang = 1:361
    angle(iang)=(iang-181)/2;
    phim=derad*angle(iang);
    a=exp(-1i*2*pi*d*sin(phim)).';
    Pmusic(iang)=1/(a'*Hm*a);
end
Pmusic=abs(Pmusic);
Pmmax=max(Pmusic);
Pmusic=10*log10(Pmusic/Pmmax);            % 归一化处理
end
