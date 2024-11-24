%function frame_dataproc()
close all;clear all;clc;
%% 文件路径
fname='C:\ti\mmwave_studio_02_01_00_00\mmWaveStudio\PostProc\adc_data_Raw_0.bin'; % 原始数据
dstpath='C:\Users\san\Desktop\data\up\b1'; % 手势信息提取图保存路径

%% 雷达参数
frame_num=20;          % 帧数
c=3.0e8;               % 光速
Rx_Number=4;           % RX天线通道数
Tx_Number=2;           % TX天线通道数
freqSlope=30e12;       % 调频斜率
Fs=5e6;                % 脉冲采样率
PDF=2500;               % 脉冲发射频率
f0=77e9;               % 初始频率
lambda=c/f0;           % 雷达信号波长
NSample=512;           % 脉冲采样点数
Range_Number=NSample;  % 距离向FFT点数
NChirp=200;            % 每帧每Tx脉冲数
Doppler_Number=NChirp*frame_num; % 总chirp数

%% 读取回波数据
fid = fopen(fname,'rb'); % 以二进制格式只读方式打开数据文件
bufferSize=Range_Number*Doppler_Number*Rx_Number*Tx_Number*2; % 读取数据大小
% 复数形式(I/Q两路)，有符号16bit，小端模式
fdata = fread(fid,bufferSize,'int16');
fclose(fid);

%% 数据拆分、组合
% 虚部+实部数据重组得到复信号 
fileSize = size(fdata, 1);
lvds_data = zeros(1, fileSize/2);
count = 1;
redata=zeros(1, fileSize/4);
for i=1:4:fileSize
   lvds_data(1,count) = fdata(i) + 1i*fdata(i+2); 
   lvds_data(1,count+1) = fdata(i+1)+1i*fdata(i+3);
   redata(count)=fdata(i);
   redata(count+1)=fdata(i+1);
   count = count + 2;
end

% 将回波数据排列成3维矩阵(距离维、速度维和天线维)
ADC_Data=zeros(Range_Number,Doppler_Number,Rx_Number*Tx_Number); % 建立计算存储数据的空矩阵
% 数据组织格式:采样点(Sample)->天线(RX)->脉冲(Chirp)
xSize=Range_Number;
chirpSize=xSize*Rx_Number*Tx_Number;
for i=1:Doppler_Number
    for j=1:Rx_Number*Tx_Number
        ADC_Data(:,i,j) = lvds_data((i-1)*chirpSize+(j-1)*xSize+1:(i-1)*chirpSize+j*xSize);
    end
end

%% 指数加权平均滤波
B=0;     % B(n):前n-1个chirp回波的加权积累
a=0.01;   % 加权系数
% s(n)为原始信号中第n个chirp回波
% 滤波:s'(n)=s(n)-B(n);加权积累:B(n+1)=B(n)+a*s'(n)
for i=1:Doppler_Number
     ADC_Data(:,i,:)=ADC_Data(:,i,:)-B;
     B=B+a*ADC_Data(:,i,:);
end

%% 距离信息提取
range_win = hamming(Range_Number);        % 加海明窗
range_profile = zeros(Range_Number,Doppler_Number,Rx_Number*Tx_Number);
for j=1:Rx_Number*Tx_Number
   for m=1:Doppler_Number
      temp=ADC_Data(:,m,j).*range_win;    % 加窗函数
      temp_fft=fft(temp,Range_Number);    % 对每个chirp做FFT
      range_profile(:,m,j)=temp_fft;
   end
end

%% 距离维滤波
rlim=0.5;                                        % 手势范围
rlimid=ceil(rlim*2*freqSlope*Range_Number/c/Fs); % 距离为rlim对应的索引
range_profile(rlimid:end,:,:)=0;                 % 将距离>rlim的杂波滤除

%% 多普勒信息提取
range_max=max(range_profile(:,:,1)); % 将距离信息中的各脉冲峰值位置数据提取出来，并重新排列成一列
win_sz=128;                           % STFT时间窗大小
[speed_profile, F, T]=DopplerProcess(range_max,PDF,win_sz);
F=fftshift(F);                       % 将零频分量移到频谱中心
F(F>=PDF/2)=F(F>=PDF/2)-PDF;         % 恢复对应频率
speed_profile=[speed_profile(win_sz/2+1:end,:);speed_profile(1:win_sz/2,:)]; % 将零频分量移到频谱中心

%% 角度信息提取
angle_profile=zeros(181,Doppler_Number);
for i=1:8:Doppler_Number
    X=reshape(ADC_Data(:,i:i+7,:),Range_Number*8,[]); % 将8个脉冲的数据拼接
    X=X.';                                            % 4(通道数)*4096(8*Range_Number)
    angle_profile(:,i:i+7)=AngleProcess(X,lambda)*ones(1,8);
end

%% 绘制信息提取图
t=(0:Doppler_Number-1)/PDF;           % 时间坐标轴    
figure();
R=c*(0:Range_Number-1)*Fs/2/freqSlope/Range_Number; % 距离坐标轴
mesh(t,R,abs(range_profile(:,:,1)));  % 绘制距离信息图
% xlabel('时间(s)');ylabel('距离(m)');
% title('距离信息提取图');
ylim([0 0.8]);                        % 限制坐标轴
view(2);                              % 转换视线
set(gca,'position',[0 0 1 1]);        % 图像填满
figure();
V=F*lambda/2;                         % 速度坐标轴
mesh(T,V,abs(speed_profile));         % 绘制速度信息图
% xlabel('时间(s)');ylabel('速度(m/s)');
% title('速度信息提取图');
ylim([-1.2 1.2]);                     % 限制坐标轴
view(2);                              % 转换视线
set(gca,'position',[0 0 1 1]);        % 图像填满
figure();
A=linspace(-90,90,181);               % 角度坐标轴
mesh(t,A,abs(angle_profile));         % 绘制角度信息图
% xlabel('时间(s)');ylabel('角度(°)');
% title('角度信息提取图');
ylim([-50 50]);                       % 限制坐标轴
view(2);                              % 转换视线
% set(gca,'position',[0 0 1 1]);        % 图像填满
figure()
[~,id]=max(angle_profile);
plot(t,id);

%% 图像保存
% oimg =getframe(1);
% img=imstand(oimg.cdata);            % 图像标准化
% imwrite(img,'range.png');
% oimg = getframe(2);
% img=imstand(oimg.cdata);            % 图像标准化
% imwrite(img,'doppler.png');
% oimg = getframe(3);
% img=imstand(oimg.cdata);            % 图像标准化
% imwrite(img,'angle.png');
% % imcopy(fname,dstpath);              % 图像拷贝到目标路径

%% 多普勒信息提取函数
function [speed_profile, F, T]=DopplerProcess(range_max,fs,win_sz)
window=hamming(win_sz);    % 使用海明窗
nfft = 2^nextpow2(win_sz); % DFT点数
nooverlap = win_sz - 4;    % 重叠单元长度
[speed_profile, F, T] = spectrogram(range_max, window, nooverlap, nfft, fs, 'yaxis');
end

%% 角度信息提取函数
function Pmusic=AngleProcess(X,lambda)
K=size(X,2);
Rxx=X*X'/K;                         % 计算协方差矩阵
[EV,D]=eig(Rxx);                    % 特征值分解
EVA=diag(D);                        % 将特征值矩阵对角线提取并转为一行
[~,I]=sort(EVA);                    % 将特征值排序 从小到大
EV=fliplr(EV(:,I));                 % 对应特征矢量排序 从大到小
                 
% 遍历每个角度，计算空间谱
Pmusic = zeros(181,1);
derad = pi/180;
N = 8;               % 阵元个数        
M = 1;               % 信源数目
dd=lambda/2;         % 天线间距
d=0:dd:(N-1)*dd;
En=EV(:,M+1:N);                     % 取矩阵的第M+1到N列组成噪声子空间
Hm=En*En';
for iang = 1:181
    angle=iang-91;
    phim=derad*angle;
    a=exp(-1i*2*pi*d/lambda*sin(phim)).';
    Pmusic(iang)=1/(a'*Hm*a);
end
end

%% 图像标准化函数
function rimg=imstand(oimg)
img=imresize(oimg,[224,224]);       % 尺寸标准化
i=min(min(min(img)));
j=max(max(max(img)));
rimg=img-i/(j-i);                   % 像素标准化
end

%% 图像拷贝函数
function imcopy(fname,dstpath)
copyfile(fname,dstpath);            % 原始数据拷贝
copyfile('Range.png',dstpath);      % 距离信息图拷贝
copyfile('Doppler.png',dstpath);    % 速度信息图拷贝
copyfile('Angle.png',dstpath);      % 角度信息图拷贝
end