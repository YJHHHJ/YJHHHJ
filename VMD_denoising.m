% VMD去噪
clc
clear all
close all

addpath('.\公式插件包\');
% 构建原始仿真冲击信号
fs = 30e3;                  % 采样频率
fn = 2e3/1;                   % 固有频率
y0 =10;                      % 位移常数
g = 0.1;                     % 阻尼系数
T = 0.005*2;                   % 重复周期
N = 4096;                  % 采样点数
NT = round(fs*T);      % 单周期采样点数
t = 0:1/fs:(N-1)/fs;      % 采样时刻
t0 = 0:1/fs:(NT-1)/fs;  % 单周期采样时刻
K = ceil(N/NT)+1;       % 重复次数
y = [];
for i = 1:K
    y = [y,y0*exp(-g*2*pi*fn*t0).*sin(2*pi*fn*sqrt(1-g^2)*t0)];
end
x = y(1:N);
figure(1)
plot(x);
ylabel('幅值A/V');%%%%纵坐标的含义
xlabel('Sample Number');%%%%横坐标的含义
   axis tight 
%  axis([0,4000,-inf,inf]);
% 信噪比和信号指定功率
amp=2;%加入随机噪声，加噪后信噪比为-2
y=noisegen(x,amp);
SNRValues1 = amp;
%%  
K=4;%分解层数
[u,residual] = vmd(y,'Nu',K);
u=u';
denoise=sum(u(3:K,:));
%%
origSignal=x;
errorSignal1=denoise-x;
signal_2 = (sum(origSignal(:).^2));
noise_2 = (sum(errorSignal1(:).^2));
SNRValues2 = 10*log10(signal_2./noise_2);
% 输出初始和最终信噪比
fprintf('初始信噪比（SNRValues1）: %.2f dB\n', SNRValues1);
fprintf('最终信噪比（SNRValues2）: %.2f dB\n', SNRValues2);


figure(5)
subplot(3,1,1)
plot(x)
title('原始无噪信号')
subplot(3,1,2)
plot(y)
title('加噪的信号')
subplot(3,1,3)
plot(denoise)
title('去噪后的信号')

sgtitle('VMD去噪')

figure(2)
[m n]=size(u);
for i=1:m
    subplot(m+1,1,i)
    plot(u(i,:))
    title(['IMF',num2str(i)])
end
subplot(m+1,1,m+1)
plot(residual')
title('残差')
sgtitle('各模态和残差')
%%
figure(5)
y0=fft(x,N);%傅里叶变换
mag0=abs(y0);
f=(0:length(y0)-1)'*fs/length(y0);
plot(f(1:N/2),mag0(1:N/2));%绘制频谱图
ylabel('幅值');xlabel('频率/HZ')
title('原始无噪频域图');
figure(6)
y1=fft(y,N);%傅里叶变换
mag1=abs(y1);
plot(f(1:N/2),mag1(1:N/2));%绘制频谱图
ylabel('幅值');xlabel('频率/HZ')
title('加噪后频域图');
figure(7)
y2=fft(denoise,N);%傅里叶变换
mag2=abs(y2);
plot(f(1:N/2),mag2(1:N/2));%绘制频谱图
ylabel('Amp');xlabel('Frequency/Hz')
title('(d)');
ylim([0 1600]);
x_min = 3000;
x_max = 5000;
y_min = 0;
y_max = 300;
% 绘制虚线框
hold on;
rectangle('Position', [x_min, y_min, x_max - x_min, y_max - y_min], 'EdgeColor', 'r', 'LineStyle', '--');
hold off;
% 绘制局部放大的图
figure(11);
plot(f(1:N/2),mag2(1:N/2));
ylabel('Amp');%%%%纵坐标的含义
xlabel('Frequency/Hz');%%%%横坐标的含义
title('(d1)');
 % 设置局部放大的区域
xlim([3000, 5000]); % 设置 X 轴范围为 [2, 4]
ylim([0, 300]); % 设置 Y 轴范围为 [0.5, 1]
figure(8)
subplot(3,1,1)
pspectrum(x,t,'spectrogram','TimeResolution', 0.1,...
    'OverlapPercent',99,'Leakage',0.85) %计算绘制时谱图,时间分辨率为0.1，相邻段之间重叠为99%，频谱泄漏量为0.85
title('原始无噪时频图');
subplot(3,1,2)
pspectrum(y,t,'spectrogram','TimeResolution', 0.1,...
    'OverlapPercent',99,'Leakage',0.85) %计算绘制时谱图,时间分辨率为0.1，相邻段之间重叠为99%，频谱泄漏量为0.85
title('加噪后时频图');
subplot(3,1,3)
pspectrum(denoise,t,'spectrogram','TimeResolution', 0.1,...
    'OverlapPercent',99,'Leakage',0.85) %计算绘制时谱图,时间分辨率为0.1，相邻段之间重叠为99%，频谱泄漏量为0.85
title('去噪后时频图');
figure(9)
plot(denoise);
ylabel('Amp');%%%%纵坐标的含义
xlabel('Sample Number');%%%%横坐标的含义
title('(d)');
   axis tight  