function RadarRead2(loadCfg,filename,freadnums)
clc;
close all;
clear all;
if nargin==0
    loadCfg=0;
    filename='adc_data1.mat';
    freadnums=1;
end
if(loadCfg)
    delete(instrfind);           % 删除串口
    [controlSerialPort,dataSerialPort] = serial_port();     % 自动识别串口
    configurationFileName_stop ='awr1642.cfg';             
    cliCfg = readCfg(configurationFileName_stop);          % 配置文件读取
    hdataSerialPort = configureDataSport(dataSerialPort,6e6);%链接 COM串口
    mmwDemoCliPrompt = char('mmwDemo:/>');
    hControlSerialPort = configureControlPort(controlSerialPort);%链接 COM串口
 
    %Send CLI configuration to AWR16xx
    fprintf('Sending configuration from %s file to AWR16xx ...\n', configurationFileName_stop);
    for k=1:length(cliCfg)
        fprintf(hControlSerialPort, cliCfg{k});
        fprintf('%s\n', cliCfg{k});
        echo = fgetl(hControlSerialPort); % Get an echo of a command
        done = fgetl(hControlSerialPort); % Get "Done" 
        prompt = fread(hControlSerialPort, size(mmwDemoCliPrompt,2)); % Get the prompt back 
    end
    prompt_1=[];
    for k=1:freadnums
        Data_dec= fread(hdataSerialPort); % Get the prompt back
        prompt_1=[prompt_1;Data_dec];
    end
    save(filename,'prompt_1');
    fprintf(hControlSerialPort, 'sensorStop');
    fclose(hControlSerialPort);
    delete(hControlSerialPort);
    fclose(hdataSerialPort);
    delete(hdataSerialPort);
    [range,doppler,angle]=DataRead(prompt_1);
else
    adc_data=load(filename);
    Data_dec= adc_data.prompt_1;     
    [range,doppler,angle]=DataRead(Data_dec);     
end

function [range,doppler,angle]=DataRead(Data_dec)
udata=dec2hex(Data_dec,2);
[range,doppler,angle]=UartProcess(udata);
% len=size(doppler,2);
% t=(0:len-1)*100;
% index=find(range<0.5);
% n_range=interp1(index,range(index),1:len,'linear');
% n_doppler=interp1(index,doppler(index),1:len,'linear');
% n_angle=interp1(index,angle(index),1:len,'linear');
% figure();
% plot(t,range);
% xlabel('时间(ms)');ylabel('距离(m)');
% title('插值前距离图');
% exportgraphics(gcf,'r1.png','Resolution',150);
% figure();
% plot(t,n_range);
% xlabel('时间(ms)');ylabel('距离(m)');
% title('插值后距离图');
% exportgraphics(gcf,'r2.png','Resolution',150);
% figure();
% plot(t,doppler);
% xlabel('时间(ms)');ylabel('速度(m/s)');
% title('插值前速度图');
% exportgraphics(gcf,'v1.png','Resolution',150);
% figure();
% plot(t,n_doppler);
% xlabel('时间(ms)');ylabel('速度(m/s)');
% title('插值后速度图');
% exportgraphics(gcf,'v2.png','Resolution',150);
% figure();
% plot(t,angle);
% xlabel('时间(ms)');ylabel('角度(°)');
% title('插值前角度图');
% exportgraphics(gcf,'a1.png','Resolution',150);
% figure();
% plot(t,n_angle);
% xlabel('时间(ms)');ylabel('角度(°)');
% title('插值后角度图');
% exportgraphics(gcf,'a2.png','Resolution',150);

function [range,doppler,angle]=UartProcess(udata)
cur=0;
len=size(udata,1);
range=[];
doppler=[];
angle=[];
while cur<len
    framelenhex=[udata(cur+16,:),udata(cur+15,:),udata(cur+14,:),udata(cur+13,:)];
    framelen=hex2dec(framelenhex);
    next=cur+framelen;
    if next>len
        break
    end
    fdata=udata((cur+1):(cur+framelen),:);
    frame_info=FrameProcess(fdata);
    range=[range,frame_info(1)];
    doppler=[doppler,frame_info(2)];
    angle=[angle,frame_info(3)];
    cur=next;
end

function frame_info=FrameProcess(fdata)
TLV_nums=hex2dec([fdata(36,:),fdata(35,:),fdata(34,:),fdata(33,:)]);
if TLV_nums==0
    return
end
TLV_num=zeros(10,1);
cur=40;
for k=1:TLV_nums
    TLV_lenhex=[fdata(cur+8,:),fdata(cur+7,:),fdata(cur+6,:),fdata(cur+5,:)];
    TLV_len=hex2dec(TLV_lenhex);
    TLV_tag=hex2dec(fdata(cur+1,:));
    payload=fdata((cur+9):(cur+8+TLV_len),:);
    cur=cur+8+TLV_len;
    switch TLV_tag
        case 1
            [dst,velocity,angle]=DobjProcess(payload,TLV_len);
            TLV_num(1)=TLV_num(1)+1;
        case 2
            TLV_num(2)=TLV_num(2)+1;
        case 3
            TLV_num(3)=TLV_num(3)+1;
        case 4
            TLV_num(4)=TLV_num(4)+1;
        case 5
            TLV_num(5)=TLV_num(5)+1;
        case 6
            TLV_num(6)=TLV_num(6)+1;
        case 7
            TLV_num(7)=TLV_num(7)+1;
        case 8
            TLV_num(8)=TLV_num(8)+1;
        case 9
            TLV_num(9)=TLV_num(9)+1;
        case 10
            [objRange,objSpeed,objAngle]=L3dataProcess(payload)
            TLV_num(10)=TLV_num(10)+1;
        otherwise
            continue
    end
end
frame_info=[dst;velocity;angle]

function [dst,velocity,angle]=DobjProcess(payload,len)
obj_nums=len/16;
derad=180/pi;
mindst=10;
velocity=0;
angle=0;
p=typecast(uint8(hex2dec(['BF';'88';'A7';'BD'])),'single');
q=typecast(uint8(hex2dec(['85';'BE';'F8';'3C'])),'single');
r=sqrt(p^2+q^2);
for k=1:obj_nums
    cur=(k-1)*16;
    x=payload((cur+1):(cur+4),:);
    x=typecast(uint8(hex2dec(x)),'single');
    y=payload((cur+5):(cur+8),:);
    y=typecast(uint8(hex2dec(y)),'single');
    z=payload((cur+9):(cur+12),:);
    z=typecast(uint8(hex2dec(z)),'single');
    dst=sqrt(x^2+y^2+z^2);
    if dst==r
        continue
    end
    if dst<mindst
        mindst=dst;
        v=payload((cur+13):(cur+16),:);
        velocity=typecast(uint8(hex2dec(v)),'single');
        angle=atan(x/y)*derad;
    end
end
dst=mindst;

%功能函数
function [userport, dataport] = serial_port()

port =IdentifySerialComs();

user_flag = 1;
data_flag = 1;

for i = 1:length(port)
    if(user_flag)
        userportID = find(strcmp('XDS110 Class Application/User UART',  port(i,1)));        
        if(userportID)
            user_flag = 0;
            userport = port(i,2);
        end 
        
    end
    
    if(data_flag)     
        dataportID = find(strcmp('XDS110 Class Auxiliary Data Port',  port(i,1)));
        if(dataportID)
               data_flag = 0;
               dataport = port(i,2);
        end
    end  
end

userport = ['COM',num2str(cell2mat(userport))];
dataport = ['COM',num2str(cell2mat(dataport))];

function config = readCfg(filename)  %读取配置文件
    config = cell(1,100);
    fid = fopen(filename, 'r');
    if fid == -1
        fprintf('File %s not found!\n', filename);
        return;
    else
        fprintf('Opening configuration file %s ...\n', filename);
    end
    tline = fgetl(fid);
    k=1;
    while ischar(tline)
        config{k} = tline;
        tline = fgetl(fid);
        k = k + 1;
    end
    config = config(1:k-1);
    fclose(fid);

function [sphandle] = configureDataSport(comPortString, bufferSize)
  
%     comPortString = ['COM' num2str(comPortNum)];
    sphandle = serial(comPortString,'BaudRate',921600);
%     set(sphandle,'Timeout',15);
    set(sphandle,'Terminator', '');
    set(sphandle,'InputBufferSize', bufferSize);
    set(sphandle,'Timeout',10);
    set(sphandle,'ErrorFcn',@dispError);
    fopen(sphandle);
  
function [sphandle] = configureControlPort(comPortString)
    %if ~isempty(instrfind('Type','serial'))
    %    disp('Serial port(s) already open. Re-initializing...');
    %    delete(instrfind('Type','serial'));  % delete open serial ports.
    %end
    %comPortString = ['COM' num2str(comPortNum)];
    sphandle = serial(comPortString,'BaudRate',115200);
    set(sphandle,'Parity','none')    
    set(sphandle,'Terminator','LF')        
    fopen(sphandle);
   
function devices = IdentifySerialComs()

devices = [];

Skey = 'HKEY_LOCAL_MACHINE\HARDWARE\DEVICEMAP\SERIALCOMM';

[~, list] = dos(['REG QUERY ' Skey]);
if ischar(list) && strcmp('ERROR',list(1:5))  %% strcmp 两个字符串相同返回1
    disp('Error: EnumSerialComs - No SERIALCOMM registry entry')
    return;
end

list = strread(list,'%s','delimiter',' '); %#ok<FPARK> requires strread()
coms = 0;

for i = 1:numel(list)  %%numel 返回元素个数
    if strcmp(list{i}(1:3),'COM')
        if ~iscell(coms)
            coms = list(i);
        else
            coms{end+1} = list{i}; %#ok<AGROW> Loop size is always small
        end
    end
end

out = 0;
outK = 0;

for j=1:2
    switch j
        case 1
            key = 'HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Enum\USB\';
        case 2
            key = 'HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Enum\FTDIBUS\';
    end
    [~, vals] = dos(['REG QUERY ' key ' /s /f "FriendlyName" /t "REG_SZ"']);
    if ischar(vals) && strcmp('ERROR',vals(1:5))
        disp('Error: EnumSerialComs - No Enumerated USB registry entry')
        return;
    end
    vals = textscan(vals,'%s','delimiter','\t');
    vals = cat(1,vals{:});
    for i = 1:numel(vals)
        if strcmp(vals{i}(1:min(12,end)),'FriendlyName')
            if ~iscell(out)
                out = vals(i);
            else
                out{end+1} = vals{i}; %#ok<AGROW> Loop size is always small
            end
            if ~iscell(outK)
                outK = vals(i-1);
            else
                outK{end+1} = vals{i-1}; %#ok<AGROW> Loop size is always small
            end
        end
    end
end

i_dev=1;Sservices=[];
for i = 1:numel(coms)
    match = strfind(out,[coms{i},')']);
    ind = 0;
    for j = 1:numel(match)
        if ~isempty(match{j})
            ind = j;
            [~, sers] = dos(['REG QUERY "' outK{ind} '" /f "Service" /t "REG_SZ"']);
            sers = textscan(sers,'%s','delimiter','\t');
            sers = cat(1,sers{:});
            if (numel(sers)>1)
                sers=strread(sers{2},'%s','delimiter',' ');
                Sservices{i_dev} = sers{3};
                i_dev=i_dev+1;
            end
        end
    end
end
Sservices=unique(Sservices);

i_dev=1;
for ss=1:numel(Sservices)
    key = ['HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\services\' Sservices{ss} '\Enum'];
    [~, vals] = dos(['REG QUERY ' key ' /f "Count"']);
    if ischar(vals) && strcmp('ERROR',vals(1:5))
        %         disp('Error: EnumSerialComs - No Enumerated services USB registry entry')
        %         return
    end
    vals = textscan(vals,'%s','delimiter','\t');
    vals = cat(1,vals{:});

    if (numel(vals)>1)
        vals=strread(vals{2},'%s','delimiter',' ');
        Count=hex2dec(vals{3}(3:end));
        if Count>0
            [~, vals] = dos(['REG QUERY ' key]);
            vals = textscan(vals,'%s','delimiter','\t');
            vals = cat(1,vals{:});
            out=0;
            j=0;
            for i = 1:numel(vals)
                Enums=strread(vals{i},'%s','delimiter',' ');
                try nums=hex2dec(Enums{1});
                catch
                    nums=-1;
                end
                if(nums==j)
                    out=['HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Enum\' Enums{3}];
                    [~, listC] = dos(['REG QUERY "' out '" /s /f "PortName" /t "REG_SZ"']);
                    listC = textscan(listC,'%s','delimiter','\t');
                    listC = cat(1,listC{:});
                    if (numel(listC)>1)
                        listC=strread(listC{2},'%s','delimiter',' ');
                        for i = 1:numel(coms)
                            if strcmp(listC{3},coms{i})
                                [~, NameF] = dos(['REG QUERY "' out '" /s /f "FriendlyName" /t "REG_SZ"']);
                                NameF = textscan(NameF,'%s','delimiter','\t');
                                NameF = cat(1,NameF{:});
                                com = str2double(coms{i}(4:end));
                                if com > 9
                                    length = 8;
                                else
                                    length = 7;
                                end
                                devices{i_dev,1} = NameF{2}(27:end-length); %#ok<AGROW>
                                devices{i_dev,2} = com; %#ok<AGROW> Loop size is always small
                                i_dev=i_dev+1;
                            end
                        end
                    end
                    j=j+1;
                end
            end
        end
    end
end


