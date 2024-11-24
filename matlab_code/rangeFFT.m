function [frameRange] = rangeFFT(frameSample,isAddWin,isRmClutter)
%对一帧的采样数据进行距离FFT

%frameSample ：输入数据，一个三维数组：天线高维*chirps行维*samples列维
%frameRange ： 输出

[numAnt,chirpsPerFrame,samplesPerChirp] = size(frameSample);
frameRange = zeros(numAnt,chirpsPerFrame,samplesPerChirp);

 winFun = (hann(samplesPerChirp))';
 
for ii=1:numAnt
    oneAntData = squeeze(frameSample(ii,:,:));   %删除单一维度数据，构建一个numChirp*numSample的二维数据并赋值
    if(isAddWin==1)
         oneAntData = oneAntData.*winFun; %加窗
    end      
      
 
    
    rangeData =fft(oneAntData, samplesPerChirp, 2);  %2表示对行做FFT
    rangeData=fftshift(rangeData);
    if(isRmClutter==1)
        sumFFT = sum(rangeData,1)/chirpsPerFrame; %静态杂波滤除    %sum每一行的采样点积累后求平均
        frameRange(ii,:,:) = rangeData - sumFFT;       
    else
        frameRange(ii,:,:) = rangeData;
    end
end

end

