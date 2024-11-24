function Angle=AngleProcess(payload)
K=512;
len=length(payload)/2;
payload=typecast(uint16(hex2dec(payload)),'single');
re_data=payload(2:2:len,:);
im_data=payload(1:2:len-1,:);
data=complex(re_data,im_data);
X=reshape(data,[],K);

% 计算协方差矩阵
Rxx=X*X'/K;
% 特征值分解
[EV,D]=eig(Rxx);                   %特征值分解
EVA=diag(D)';                      %将特征值矩阵对角线提取并转为一行
[~,I]=sort(EVA);                 %将特征值排序 从小到大
EV=fliplr(EV(:,I));                % 对应特征矢量排序
                 
% 遍历每个角度，计算空间谱
angle = zeros(1,361);
Pmusic = zeros(1,361);
derad = pi/180;
N = 4;               % 阵元个数        
M = 2;               % 信源数目
dd = 0.5;            % 阵元间距 
d=0:dd:(N-1)*dd;
for iang = 1:361
    angle(iang)=(iang-181)/2;
    phim=derad*angle(iang);
    a=exp(-1i*2*pi*d*sin(phim)).';
    En=EV(:,M+1:N);                   % 取矩阵的第M+1到N列组成噪声子空间
    Pmusic(iang)=1/(a'*En*En'*a);
end
Pmusic=abs(Pmusic);
[Pmmax,index]=max(Pmusic);
Angle=angle(index);
Pmusic=10*log10(Pmusic/Pmmax);            % 归一化处理
h=plot(angle,Pmusic);
set(h,'Linewidth',2);
xlabel('入射角/(degree)');
ylabel('空间谱/(dB)');
set(gca, 'XTick',[-90:30:90]);
grid on;