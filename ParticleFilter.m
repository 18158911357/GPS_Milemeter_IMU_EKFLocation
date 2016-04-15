%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%程序初始化操作%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PCenter=ParticleFilter()
% clc;
% clear;
% close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%全局变量定义%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outdoor_sensor_data=361;
indoor_sensor_data=0;
sensor_data=outdoor_sensor_data+indoor_sensor_data;
d=0.1;%标准差
Theta=CreateGauss(0,d,1,sensor_data);%GPS航迹和DR航迹的夹角
ZOUT=zeros(4,outdoor_sensor_data);
ZIN=zeros(4,indoor_sensor_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%读取传感器数据%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fgps=fopen('sensor_data_041518.txt','r');%%%打开文本

for n=1:sensor_data
    gpsline=fgetl(fgps);%%%读取文本指针对应的行
    if ~ischar(gpsline) break;%%%判断是否结束
    end;
    %%%%读取室内数据
   time=sscanf(gpsline,'[Info] 2016-04-15%s(ViewController.m:%d)-[ViewController outputAccelertion:]:lat:%f;lon:%f;heading:%f;distance:%f;beacon_lat:%f;beacon_lon:%f');
   data=sscanf(gpsline,'[Info] 2016-04-15 %*s (ViewController.m:%*d)-[ViewController outputAccelertion:]:lat:%f;lon:%f;heading:%f;distance:%f;beacon_lat:%f;beacon_lon:%f');
   if(isempty(data))
       break;
   end
        result=lonLat2Mercator(data(6,1),data(5,1));
        gx(n)=result.X;%GPS经过坐标变换后的东向坐标，换算成米数
        gy(n)=result.Y;%GPS经过坐标变换后的北向坐标，换算成米数
        Phi(n)=(data(3,1)+90)*pi/180;%航向角
        dd(n)=data(4,1);%某一周期的位移
        ZIN(:,n)=[gx(n),gy(n),Phi(n),dd(n)];
end
fclose(fgps);%%%%%关闭文件指针

% 参数设置
N = 100;   %粒子总数
Q = 5;      %过程噪声
R = 5;      %测量噪声
X = zeros(2, sensor_data);    %存储系统状态
Z = zeros(2, sensor_data);    %存储系统的观测状态
P = zeros(2, N);    %建立粒子群
PCenter = zeros(2, sensor_data);  %所有粒子的中心位置
w = zeros(N, 1);         %每个粒子的权重
err = zeros(1,sensor_data);     %误差
X(:, 1) = [ZIN(1,1); ZIN(2,1)];     %初始系统状态
Z(:, 1) = [ZIN(1,1); ZIN(2,1)] ;    %初始系统的观测状态

cordinatex=round(ZIN(1,5));
cordinatey=round(ZIN(2,5));

%初始化粒子群
for i = 1 : N
    P(:, i) = [randi([cordinatex-100,cordinatex+200],1);randi([cordinatey-200,cordinatey+100],1)];
    dist = norm(P(:, i)-Z(:, 1));     %与测量位置相差的距离
    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重
end
PCenter(:, 1) = sum(P, 2) / N;      %所有粒子的几何中心位置
 
%%
err(1) = norm(X(:, 1) - PCenter(:, 1));     %粒子几何中心与系统真实状态的误差
% figure(1);
% set(gca,'FontSize',12);
% hold on
% plot(X(1, 1), X(2, 1), 'r.', 'markersize',30)   %系统状态位置
% axis([cordinatex-100 cordinatex+200 cordinatey-200 cordinatey+100]),grid on;
% plot(P(1, :), P(2, :), 'k.', 'markersize',5);   %各个粒子位置
% plot(PCenter(1, 1), PCenter(2, 1), 'b.', 'markersize',25); %所有粒子的中心位置
% legend('True State', 'Particles', 'The Center of Particles');
% title('Initial State');
% hold off

%%
%开始运动
for k = 2 : sensor_data
       
    %模拟一个弧线运动的状态
    X(:, k) = X(:, k-1) + 1 * [(sin( ZIN(3,k))); cos(ZIN(3,k))] ;     %状态方程
    Z(:, k) = ZIN(1:2, k) ;     %观测方程 
    if ZIN(1,k) == 0
        Z(:, k) = X(:, k);     %观测方程 
    end
   
    %粒子滤波
    %预测
    for i = 1 : N
        P(:, i) = P(:, i) + 1 * [sin( ZIN(3,k)); cos( ZIN(3,k))] + wgn(2, 1, 10*log10(Q));
        dist = norm(P(:, i)-Z(:, k));     %与测量位置相差的距离
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %求权重
    end
%归一化权重
    wsum = sum(w);
    for i = 1 : N
        w(i) = w(i) / wsum;
    end
   
    %重采样（更新）
    for i = 1 : N
        wmax = 2 * max(w) * rand;  %另一种重采样规则
        index = randi(N, 1);
        while(wmax > w(index))
            wmax = wmax - w(index);
            index = index + 1;
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     %得到新粒子
    end
   
    PCenter(:, k) = sum(P, 2) / N;      %所有粒子的中心位置
   
    %计算误差
    err(k) = norm(X(:, k) - PCenter(:, k));     %粒子几何中心与系统真实状态的误差
   
%     figure(2);
%     set(gca,'FontSize',12);
%     clf;
%     hold on
%     plot(X(1, k), X(2, k), 'r.', 'markersize',50);  %系统状态位置
%     axis([cordinatex-100 cordinatex+200 cordinatey-200 cordinatey+100]),grid on;
%     plot(P(1, :), P(2, :), 'k.', 'markersize',5);   %各个粒子位置
%     plot(PCenter(1, k), PCenter(2, k), 'b.', 'markersize',25); %所有粒子的中心位置
%     legend('True State', 'Particle', 'The Center of Particles');
%     hold off
%     pause(0.1);
end
%%
% figure(3);
% set(gca,'FontSize',12);
% [groundtruthx,groundtruthy]=Groud_Truth();
% plot(groundtruthx,groundtruthy,'r');hold on;
% plot( ZIN(1,1:127), ZIN(2,1:127), 'o');hold on;
% plot(PCenter(1,:), PCenter(2,:), 'g');hold off;
% axis([cordinatex-100 cordinatex+200 cordinatey-200 cordinatey+100]),grid on;
% legend('真实轨迹','观测轨迹', '粒子滤波轨迹');
% xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);
% axis equal;
