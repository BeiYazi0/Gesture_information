function imcopy(fname,dstpath)
if nargin==0
   fname='C:\ti\mmwave_studio_02_01_00_00\mmWaveStudio\PostProc\adc_data_Raw_0.bin';
   dstpath='C:\Users\san\Desktop\data\up\b1';
end
copyfile(fname,dstpath);            % 原始数据拷贝
copyfile('Range.png',dstpath);      % 距离信息图拷贝
copyfile('Doppler.png',dstpath);    % 速度信息图拷贝
copyfile('Angle.png',dstpath);      % 角度信息图拷贝
end