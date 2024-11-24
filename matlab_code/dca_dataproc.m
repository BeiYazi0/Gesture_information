function dca_dataproc()
close all;clear all;clc;
%% 雷达参数（使用mmWave Studio默认参数） 
frame_num=30;          %帧数
T=0.05;                 %帧周期
c=3.0e8;               %光速
freqSlope=30e12;       %调频斜率
Tc=195.3e-6;             %chirp总周期
Fs=5e6;                %采样率
f0=77e9;               %初始频率
lambda=c/f0;           %雷达信号波长
d=lambda/2;            %天线阵列间距
Range_Number=512;      %采样点数/脉冲
NSample=Range_Number;  %距离向FFT点数
Doppler_Number=128;    %每帧脉冲数
NChirp=Doppler_Number; %多普勒向FFT点数
Rx_Number=4;           %RX天线通道数
Tx_Number=2;           %TX天线通道数
Q = 180;               %角度FFT点数

%% 读取回波数据
fname='C:\ti\mmwave_studio_02_01_00_00\mmWaveStudio\PostProc\adc_data_Raw_0.bin';
fid = fopen(fname,'rb');
bufferSize=NSample*NChirp*Rx_Number*Tx_Number*2;
%16bits，复数形式(I/Q两路)，4RX,1TX,有符号16bit，小端模式
sdata = fread(fid,frame_num*bufferSize,'int16');
fclose(fid);

%% 数据处理
Range=zeros(NSample,frame_num);
Doppler=zeros(NChirp,frame_num);
Angle=zeros(Q,frame_num);
for xx=1:frame_num
    fdata = sdata((xx-1)*bufferSize+1:xx*bufferSize);
    [fft1d,fft2d,fft3d]=Three_dfft(fdata);
    Range(:,xx)=fft1d;
    Doppler(:,xx)=fft2d;
    Angle(:,xx)=fft3d;
end
t=(0:frame_num-1)*T;
figure();
R=c*(0:NSample-1)*Fs/2/freqSlope/NSample;
pcolor(t,R,abs(Range));
shading interp;
ylim([0 0.8]);                    %限制坐标轴
set(gca,'position',[0 0 1 1]);    %图像填满
oimg =getframe(1);
img=imstand(oimg.cdata);          %图像标准化
imwrite(img,'range.png');
figure();
V=(-NChirp/2:NChirp/2 - 1)*lambda/Tc/NChirp/2;
pcolor(t,V,abs(Doppler));
shading interp;
ylim([-1.5 1.5]);                 %限制坐标轴
set(gca,'position',[0 0 1 1]);    %图像填满
oimg = getframe(2);
img=imstand(oimg.cdata);          %图像标准化
imwrite(img,'doppler.png');
figure();
A=asin((-Q/2:Q/2 - 1)*lambda/d/Q)*180/pi;
pcolor(t,A,abs(Angle));
shading interp;
ylim([-50 50]);                   %限制坐标轴
set(gca,'position',[0 0 1 1]);    %图像填满
oimg = getframe(3);
img=imstand(oimg.cdata);          %图像标准化
imwrite(img,'angle.png');
figure()
mesh(t,R,abs(Range));
figure();
mesh(t,V,abs(Doppler));
figure();
mesh(t,A,abs(Angle));
end

%% 功能函数
function rimg=imstand(oimg)
img=imresize(oimg,[224,224]);
i=min(min(min(img)));
j=max(max(max(img)));
rimg=img-i/(j-i);
end

function objres=peakfind(fft3d,lambda,freqSlope,Fs,Tc)
%% 计算峰值位置
fft3d=abs(fft3d);
pink=max(fft3d(:));
[row,col,pag]=ind2sub(size(fft3d),find(fft3d==pink));

%% 计算目标距离、速度、角度
[Doppler_Number,Range_Number,Q]=size(fft3d);
c=3.0e8;                                           %光速
fb = ((col-1)*Fs)/Range_Number;                    %差拍频率
fd = (row-Doppler_Number/2-1)/(Doppler_Number*Tc); %多普勒频率
fw = (pag-Q/2-1)/Q;                                %空间频率
d=lambda/2;                                        %天线阵列间距
R = c*(fb-fd)/2/freqSlope;                         %距离公式
v = lambda*fd/2;                                   %速度公式
theta = asin(fw*lambda/d);                         %角度公式
angle = theta*180/pi;
objres=[R;v;angle];
end