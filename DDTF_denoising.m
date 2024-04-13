% 清空工作区
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
xlabel('采样点个数i/个');%%%%横坐标的含义
   axis tight 
%  axis([0,4000,-inf,inf]);
% 信噪比和信号指定功率
amp=4;%加入随机噪声，加噪后信噪比为-2
y=noisegen(x,amp);
SNRValues1 = amp;
%% 参数设置
block_size=[1,64];%字典原子长度
step_size=[1,1];%每次移动步长

d_1=8;
d_2=8;
d_3=8;
d_4=8;
%%
% 字典学习
Data = im2colstep(y, block_size, step_size);
rperm = randperm(size(Data, 2));
line = round(1 * length(rperm));
Datapatch = Data(:, rperm(1:line));
learnt_dict = filter_learning(Data, 9.4);
a1=OMP(learnt_dict,Data,3);
s_n=learnt_dict * a1;
A2 = col2imstep(s_n, size(y), block_size, step_size);
Y2 = countcover(size(y), block_size, step_size);
denoise = A2 ./ Y2;
%% 
% 计算最终信噪比
origSignal = x;
errorSignal1 =x-denoise;
signal_2 = sum(origSignal(:).^2);
noise_2 = sum(errorSignal1(:).^2);
SNRValues2 = 10 * log10(signal_2 / noise_2);

% 显示字典
figure(1)
I = displayDictionaryElementsAsImage(learnt_dict, d_1, d_2, d_3, d_4);

% 输出初始和最终信噪比
fprintf('初始信噪比（SNRValues1）: %.2f dB\n', SNRValues1);
fprintf('最终信噪比（SNRValues2）: %.2f dB\n', SNRValues2);


figure(1)
I = displayDictionaryElementsAsImage(learnt_dict,d_1,d_2,d_3,d_4);


figure(2)
subplot(3,1,1)
plot(x)
title('原始无噪信号')
subplot(3,1,2)
plot(y)
title('含噪的信号')
subplot(3,1,3)
plot(denoise)
title('去噪后的信号')
sgtitle('DDTF去噪')
%%
figure(3)
y0=fft(x,N);%傅里叶变换
mag0=abs(y0);
f=(0:length(y0)-1)'*fs/length(y0);
plot(f(1:N/2),mag0(1:N/2));%绘制频谱图
ylabel('幅值');xlabel('频率/HZ')
title('原始无噪频域图');
y1=fft(y,N);%傅里叶变换
mag1=abs(y1);
figure(4)
plot(f(1:N/2),mag1(1:N/2));%绘制频谱图
ylabel('幅值');xlabel('频率/HZ')
title('加噪后频域图');
y2=fft(denoise,N);%傅里叶变换
mag2=abs(y2);
figure(5)
plot(f(1:N/2),mag2(1:N/2));%绘制频谱图
ylabel('Amp');xlabel('Frequency/Hz')
title('(e)');
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
title('(e1)');
 % 设置局部放大的区域
xlim([3000, 5000]); % 设置 X 轴范围为 [2, 4]
ylim([0, 300]); % 设置 Y 轴范围为 [0.5, 1]
figure(6)
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
figure(7)
plot(denoise);
ylabel('Amp');%%%%纵坐标的含义
xlabel('Sample Number');%%%%横坐标的含义
title('(e)');
   axis tight 