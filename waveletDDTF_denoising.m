%小波去噪
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
amp=4;%加入随机噪声，加噪后信噪比为-2
y=noisegen(x,amp);

SNRValues1 = amp;
%%
block_size=[1,64];
step_size=[1,1];
d_1=8;
d_2=8;
d_3=8;
d_4=8;
%%
scale=2;
[c1,l1]=wavedec(y,2,'db4');
x1_a2=appcoef(c1,l1,'db4',2);
x1_d2=detcoef(c1,l1,2);
x1_d1=detcoef(c1,l1,1);
%%
Datad1 = im2colstep(x1_d1, block_size, step_size);
rperm = randperm(size(Datad1, 2));
line = round(1 * length(rperm)); 
Datapatch1 = Datad1(:, rperm(1:line));
learnt_dict1 = filter_learning(Datad1, 3);%9.1
a1=OMP(learnt_dict1,Datad1,1);
s_n=learnt_dict1 * a1;
A1= col2imstep(s_n, size(x1_d1), block_size, step_size);
Y1= countcover(size(x1_d1), block_size, step_size);
yd1= A1 ./ Y1;

Datad2 = im2colstep(x1_d2, block_size, step_size);
rperm = randperm(size(Datad2, 2));
line = round(1 * length(rperm));
Datapatch2 = Datad2(:, rperm(1:line));
learnt_dict2 = filter_learning(Datad2, 3);%8
a2=OMP(learnt_dict2,Datad2,2);
s_n=learnt_dict2 * a2;
A2= col2imstep(s_n, size(x1_d2), block_size, step_size);
Y2 = countcover(size(x1_d2), block_size, step_size);
yd2 = A2 ./ Y2;

Dataa2 = im2colstep(x1_a2,block_size, step_size);
rperm = randperm(size(Dataa2, 2));
line = round(1 * length(rperm));
Datapatcha2 = Dataa2(:, rperm(1:line));
learnt_dicta2= filter_learning(Dataa2,4);%2.9
aa2=OMP(learnt_dicta2,Dataa2,4);
s_n=learnt_dicta2 * aa2;
Aa2= col2imstep(s_n, size(x1_a2), block_size, step_size);
Ya2 = countcover(size(x1_a2), block_size, step_size);
ya2 = Aa2 ./ Ya2;
%%
b2=ya2;
d2=yd2;
d1=yd1;
figure(13)
plot(d1)
  axis tight 
  ylim([-2.2,2.6])
  figure(14)
plot(d2)
  axis tight
  ylim([-4.7,4.2])
  figure(17)
plot(b2)
  axis tight 
  ylim([-9,10])
c=[b2,d2,d1];
denoise1=waverec(c,l1,'db4');
%%
[c,l]=wavedec(denoise1,scale,'db4');
[thr,sorh,keepapp]=ddencmp('den','wv',denoise1);
thr=0.35;
denoise=wdencmp('gbl',c,l1,'db4',scale,thr,sorh,keepapp);
%%
origSignal=x;
errorSignal=denoise-x;
signal_2 = (sum(origSignal(:).^2));
noise_2 = (sum(errorSignal(:).^2));
SNRValues2 = 10*log10(signal_2./noise_2);
fprintf('初始信噪比（SNRValues1）: %.2f dB\n', amp);
fprintf('最终信噪比（SNRValues2）: %.2f dB\n', SNRValues2);
     
figure(1)
I = displayDictionaryElementsAsImage(learnt_dict1,d_1, d_2,d_3,d_4);
figure(2)
I = displayDictionaryElementsAsImage(learnt_dict2,d_1, d_2,d_3,d_4);

figure(6)
I = displayDictionaryElementsAsImage(learnt_dicta2,d_1, d_2,d_3,d_4);
%%
figure(7)
subplot(3,1,1)
plot(x)
title('原始无噪信号')
subplot(3,1,2)
plot(y)
title('含噪的信号')
subplot(3,1,3)
plot(denoise)
title('去噪后的信号')
sgtitle('双稀疏去噪')
%%
di=appcoef(c,l1,'db4',scale);
h={};
for i=1:scale
    h1=detcoef(c,l1,i);
    h{i}=h1;
end
figure(8)
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
figure(9)
y0=fft(x,N);%傅里叶变换
mag0=abs(y0);
f=(0:length(y0)-1)'*fs/length(y0);
subplot(3,1,1);
plot(f(1:N/2),mag0(1:N/2));%绘制频谱图
ylabel('幅值');xlabel('频率/HZ')
title('原始无噪频域图');
y1=fft(y,N);%傅里叶变换
mag1=abs(y1);
subplot(3,1,2);
plot(f(1:N/2),mag1(1:N/2));%绘制频谱图
ylabel('幅值');xlabel('频率/HZ')
title('加噪后频域图');
y2=fft(denoise,N);%傅里叶变换
mag2=abs(y2);
subplot(3,1,3);
plot(f(1:N/2),mag2(1:N/2));%绘制频谱图
ylabel('幅值');xlabel('频率/HZ')
title('去噪后频域图');
%%
figure(10)
subplot(3,1,1)
pspectrum(x,t,'spectrogram','TimeResolution', 0.1,'OverlapPercent',99,'Leakage',0.85) %计算绘制时谱图,时间分辨率为0.1，相邻段之间重叠为99%，频谱泄漏量为0.85
title('原始无噪时频图');
subplot(3,1,2)
pspectrum(y,t,'spectrogram','TimeResolution', 0.1,'OverlapPercent',99,'Leakage',0.85) %计算绘制时谱图,时间分辨率为0.1，相邻段之间重叠为99%，频谱泄漏量为0.85
title('加噪后时频图');
subplot(3,1,3)
pspectrum(denoise,t,'spectrogram','TimeResolution', 0.1,'OverlapPercent',99,'Leakage',0.85) %计算绘制时谱图,时间分辨率为0.1，相邻段之间重叠为99%，频谱泄漏量为0.85
title('去噪后时频图');
figure(11)
plot(denoise);
ylabel('Amp');%%%%纵坐标的含义
xlabel('Sample Number');%%%%横坐标的含义
title('(f)');
   axis tight 
  figure(12) 
  y2=fft(denoise,N);%傅里叶变换
mag2=abs(y2);
   plot(f(1:N/2),mag2(1:N/2));%绘制频谱图
ylabel('Amp');xlabel('Frequency/HZ')
title('(f)');
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
title('(f1)');
 % 设置局部放大的区域
xlim([3000, 5000]); % 设置 X 轴范围为 [2, 4]
ylim([0, 300]); % 设置 Y 轴范围为 [0.5, 1]