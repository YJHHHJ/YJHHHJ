clc
clear
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
ylabel('Amp');%%%%纵坐标的含义
xlabel('Sample Number');%%%%横坐标的含义
title('(a)');

   axis tight 
   ylim([-15 15]);
%  axis([0,4000,-inf,inf]);
% 信噪比和信号指定功率
amp=-10;%加入随机噪声，加噪后信噪比为-2
y=noisegen(x,amp);

SNRValues1 = amp

%%
scale=2;%小波分解层数

[c,l]=wavedec(y,scale,'db4');
[thr,sorh,keepapp]=ddencmp('den','wv',y);
%thr=29;
denoise=wdencmp('gbl',c,l,'db4',scale,thr,sorh,keepapp);
%%
%硬阈值
% o=find(abs(denoise)<6.5);
% denoise(o)=0;
%%
origSignal=x;
 errorSignal=x-denoise;
signal_2 = (sum(origSignal(:).^2));
noise_2 = (sum(errorSignal(:).^2));
SNRValues2 = 10*log10(signal_2./noise_2)

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
sgtitle('小波去噪')
%%
di=appcoef(c,l,'db4',scale);
h={};
for i=1:scale
    h1=detcoef(c,l,i);
    h{i}=h1;
end
figure(3)
for i=1:scale
    subplot(scale+1,1,i)
    plot(h{i})
    title(['高频第',num2str(i),'层小波系数'])
    xlim([1,length(h{i})])
end
subplot(scale+1,1,scale+1)
plot(di)
xlim([1,length(di)])
title('低频小波系数')
sgtitle('小波各尺度系数') 
%%
figure(5)
y0=fft(x,N);%傅里叶变换
mag0=abs(y0);
f=(0:length(y0)-1)'*fs/length(y0);
plot(f(1:N/2),mag0(1:N/2));%绘制频谱图
ylabel('Amp');xlabel('Frequency/Hz')
title('(a)');
ylim([0 2000]);
x_min = 0;
x_max = 5000;
y_min = 0;
y_max = 2000;
% 绘制虚线框
hold on;
rectangle('Position', [x_min, y_min, x_max - x_min, y_max - y_min], 'EdgeColor', 'r', 'LineStyle', '--');
hold off;
% 绘制局部放大的图
figure(6);
plot(f(1:N/2),mag0(1:N/2));
ylabel('Amp');%%%%纵坐标的含义
xlabel('Frequency/Hz');%%%%横坐标的含义
title('(a1)');
 % 设置局部放大的区域
xlim([0, 5000]); % 设置 X 轴范围为 [2, 4]
ylim([0, 2000]); % 设置 Y 轴范围为 [0.5, 1]
%%
figure(7)
y1=fft(y,N);%傅里叶变换
mag1=abs(y1);
plot(f(1:N/2),mag1(1:N/2));%绘制频谱图
ylabel('Amp');xlabel('Frequency/Hz')
title('(b)');
ylim([0 2000]);
x_min = 0;
x_max = 5000;
y_min = 0;
y_max = 2000;
% 绘制虚线框
hold on;
rectangle('Position', [x_min, y_min, x_max - x_min, y_max - y_min], 'EdgeColor', 'r', 'LineStyle', '--');
hold off;
% 绘制局部放大的图
figure(8);
plot(f(1:N/2),mag1(1:N/2));
ylabel('Amp');%%%%纵坐标的含义
xlabel('Frequency/Hz');%%%%横坐标的含义
title('(b1)');
 % 设置局部放大的区域
xlim([0, 5000]); % 设置 X 轴范围为 [2, 4]
ylim([0, 2000]); % 设置 Y 轴范围为 [0.5, 1]
figure(9)
y2=fft(denoise,N);%傅里叶变换
mag2=abs(y2);
plot(f(1:N/2),mag2(1:N/2));%绘制频谱图
ylabel('Amp');xlabel('Frequency/Hz')
title('(c)');
ylim([0 2000]);
x_min = 0;
x_max = 5000;
y_min = 0;
y_max = 2000;
% 绘制虚线框
hold on;
rectangle('Position', [x_min, y_min, x_max - x_min, y_max - y_min], 'EdgeColor', 'r', 'LineStyle', '--');
hold off;
% 绘制局部放大的图
figure(10);
plot(f(1:N/2),mag2(1:N/2));
ylabel('Amp');%%%%纵坐标的含义
xlabel('Frequency/Hz');%%%%横坐标的含义
title('(c1)');
 % 设置局部放大的区域
xlim([0, 5000]); % 设置 X 轴范围为 [2, 4]
ylim([0, 2000]); % 设置 Y 轴范围为 [0.5, 1]

figure
pspectrum(x,t,'spectrogram','TimeResolution', 0.1,...
    'OverlapPercent',99,'Leakage',0.85) %计算绘制时谱图,时间分辨率为0.1，相邻段之间重叠为99%，频谱泄漏量为0.85
title('(a)');
ylabel('Frequency/kHz');%%%%纵坐标的含义
xlabel('Time/s');%%%%横坐标的含义
figure
pspectrum(y,t,'spectrogram','TimeResolution', 0.1,...
    'OverlapPercent',99,'Leakage',0.85) %计算绘制时谱图,时间分辨率为0.1，相邻段之间重叠为99%，频谱泄漏量为0.85
title('(b)');
ylabel('Frequency/kHz');%%%%纵坐标的含义
xlabel('Time/s');%%%%横坐标的含义
figure(18)
pspectrum(denoise,t,'spectrogram','TimeResolution', 0.1,...
    'OverlapPercent',99,'Leakage',0.85) %计算绘制时谱图,时间分辨率为0.1，相邻段之间重叠为99%，频谱泄漏量为0.85
title('(c)');
ylabel('Frequency/kHz');%%%%纵坐标的含义
xlabel('Time/s');%%%%横坐标的含义
figure(12)
plot(denoise);
ylabel('Amp');%%%%纵坐标的含义
xlabel('Sample Number');%%%%横坐标的含义
title('(c)');
   axis tight 
  ylim([-15 15]);
   figure(13)
plot(y);
ylabel('Amp');%%%%纵坐标的含义
xlabel('Sample Number');%%%%横坐标的含义
title('(b)');
   axis tight 
    ylim([-15 15]);

%%
% 1. 计算数据的平均值
mean_data1 = mean(x);
mean_data2 = mean(denoise);

% 2. 计算数据的标准差
std_data1 = std(x);
std_data2 = std(denoise);

% 3. 计算协方差
covariance = cov(x, denoise);

% 4. 计算皮尔逊相关系数
pearson_corr = covariance(1,2) / (std_data1 * std_data2);

% 显示结果
disp(['皮尔逊相关系数 (Pearson correlation) = ', num2str(pearson_corr)]);
%%
% 1. 计算数据的平均值
mean_data1 = mean(x);
mean_data2 = mean(y);

% 2. 计算数据的标准差
std_data1 = std(x);
std_data2 = std(y);

% 3. 计算协方差
covariance = cov(x, y);

% 4. 计算皮尔逊相关系数
pearson_corr = covariance(1,2) / (std_data1 * std_data2);

% 显示结果
disp(['皮尔逊相关系数 (Pearson correlation)2 = ', num2str(pearson_corr)]);
