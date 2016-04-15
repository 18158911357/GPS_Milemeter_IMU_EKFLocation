%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%程序初始化操作%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%全局变量定义%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outdoor_sensor_data=370;
indoor_sensor_data=0;
sensor_data=outdoor_sensor_data+indoor_sensor_data;
d=0.1;%标准差
T=1;%采样时间
Theta=CreateGauss(0,d,1,sensor_data);%GPS航迹和DR航迹的夹角
ZOUT=zeros(4,outdoor_sensor_data);
ZIN=zeros(4,indoor_sensor_data);
Z=zeros(2,outdoor_sensor_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%读取传感器数据%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fgps=fopen('sensor_data_0415.txt','r');%%%打开文本

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
        result=lonLat2Mercator(data(2,1),data(1,1));
        gx(n)=result.X;%GPS经过坐标变换后的东向坐标，换算成米数
        gy(n)=result.Y;%GPS经过坐标变换后的北向坐标，换算成米数
        Phi(n)=(data(3,1)+90)*pi/180;%航向角
        dd(n)=data(4,1);%某一周期的位移
        dx(n)=dd(n)*sin(Phi(n))*4;%某一周期的东向位移
        dy(n)=dd(n)*cos(Phi(n))*4;%某一周期的北向位移
        Ve(n)=dd(n)*sin(Phi(n));%里程计输入的东向速度，暂时用某一周期的东向位移代替
        Vn(n)=dd(n)*cos(Phi(n));%里程计输出的北向速度，暂时用某一周期的北向位移代替
        ZOUT(:,n)=[gx(n),gy(n),dx(n),dy(n)];
        Z(:,n)=[gx(n),gy(n)];
end
fclose(fgps);%%%%%关闭文件指针

%过程向量A
A=[
    1,0,T,0,0;
    0,1,0,T,0;
    0,0,1,0,0,;
    0,0,0,1,0;
    0,0,0,0,1;
   ];
%过程噪声协方差矩阵
Q=diag([0,0,d^2,d^2,d^2]);
%观测噪声协方差矩阵
R=diag([d^2,d^2,d^2,d^2]);
P=diag([0,0,d^2,d^2,d^2]); % 滤波输出误差均方差矩阵
Xfli=[ZOUT(1,1),ZOUT(2,1),0,0,0]'; %初始条件进行估计
for k=1:sensor_data
    C=[1,0,0,0,0;
        0,1,0,0,0;
        0,0,cos(Theta(k)),-sin(Theta(k)),-Ve(k)*sin(Theta(k))-Vn(k)*cos(Theta(k));
        0,0,sin(Theta(k)),cos(Theta(k)),Ve(k)*cos(Theta(k))-Vn(k)*sin(Theta(k));
        ];
     K_location(:,k)=Xfli;
     Xest=A*Xfli; % 更新该时刻的预测值 ---kalman equation1
     %Xes=A*Xef+Gamma*W(k-1); % 预测输出误差 
     Pxe=A*P*A'+Q; % 预测误差的协方差阵 ---kalman equation
     K=Pxe*C'/(C*Pxe*C'+R); % Kalman滤波增益 ---kalman equation3
     Xfli=Xest+K*(ZOUT(:,k)-C*Xest);% k时刻Kalman滤波器的输出值 ---kalman equation4
     Px=(eye(5)-K*C)*Pxe;%滤波输出误差均方差矩阵 ---kalman equation5
end

cordinatex=ZOUT(1,5);
cordinatey=ZOUT(2,5);
%显示滤波轨迹
figure
set(gca,'FontSize',12);
[groundtruthx,groundtruthy]=Groud_Truth();
plot(groundtruthx,groundtruthy,'r');hold on;
plot(ZOUT(1,:),ZOUT(2,:),'o');hold on;
plot(K_location(1,:),K_location(2,:),'g');hold off;
axis([cordinatex-100 cordinatex+200 cordinatey-200 cordinatey+100]),grid on;
xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);
legend('真实轨迹','观测轨迹','目标滤波航迹');
axis equal;

    

