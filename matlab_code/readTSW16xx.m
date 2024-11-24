function [retVal] = readTSW16xx(fileName)
%function
%输入为文件名fileName 输出为retVal
%将bin文件改为相应的格式的数据并返回做后续相应的处理
%% change based on sensor config
numADCSamples = 256; % number of ADC samples per chirp
numADCBits = 16; % number of ADC bits per sample
numRX = 4; % number of receivers
numLanes = 2; % do not change. number of lanes is always 2
isReal = 0; % set to 1 if real only data, 0 if complex data0 

%% read .bin file
fid = fopen(fileName,'r');
adcData = fread(fid, 'uint16');               

%% compensate for offset binary format
adcData = adcData - 2^15;
fclose(fid);

%% get total file
fileSize = size(adcData, 1);                %把adcdata的行数赋值给fileSize
test = adcData;

%% for complex data
adcData = reshape(adcData, numLanes, []);      %将adcData改变为二维的变量2*n
% seperate each LVDS lane into rows
LVDS = zeros(1, fileSize/2);                 %1*(1/2行数)的0矩阵
LVDS(1,:) = adcData(1, :) + sqrt(-1)*adcData(2,:);      %每一列元素由adcData的第一行元素+i*(第二行元素)
numChirps = fileSize/2/numADCSamples/numRX;        %numchirp=行数/（2*numADCSaples*numRX）     

%% organize data by receiver
adcData = zeros(numRX,numChirps*numADCSamples);      %构建一个4*（numChirps*numADCSamples）的零矩阵          
LVDS = reshape(LVDS, numADCSamples*numRX, numChirps);       %将LVDS构建成（numADCSamples*numRX）*numChirps
LVDS = LVDS.';                                              % LVDS转置为 numChirp*(numADCSamples*numRX)
for row = 1:numRX
for i = 1: numChirps
adcData(row, (i-1)*numADCSamples+1:i*numADCSamples) = LVDS(i, (row-1)*numADCSamples+1:row*numADCSamples);
end
end

%% return receiver data
retVal = adcData;

end
