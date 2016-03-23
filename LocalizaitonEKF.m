% function XE=LocationEKF(T,offtime,d)
%LocationEKF     采用卡尔曼滤波方法，从观测数值中得到航迹的最优估计
%XE              输出x轴方向上的误差
%T              采样时间，gps观测周期
%offtime         仿真时间
%d               噪声的标准差
clc;
clear;
T=1;
offtime=10;
d=0.1;

close all
N=ceil(offtime/T); %采样点数

A=zeros(7,7);
W=zeros(7,1);
X=zeros(7,1); % 一次的滤波输出值
C=zeros(4,7);
V=zeros(4,1);
Z=zeros(4,1);%一次的观测值
XE=zeros(7,N);%所有的预测值
Ve=zeros(N,1);
Vn=zeros(N,1);
ZE=zeros(4,N);

randn('state',sum(100*clock)); % 设置随机数发生器
%%%%%%%%%%%读取文本gps值、里程计位移和航向角，并转化为直角坐标系下的坐标值%%%%%%
fgps=fopen('VelocityDetector.txt','r');%%%打开文本
n=0;
while 1
   gpsline=fgetl(fgps);%%%读取文本指针对应的行
   if ~ischar(gpsline) break;%%%判断是否结束
   end;
   n=n+1;
   time=sscanf(gpsline,'[Info] 2016-03-23%s(ViewController.m:%d)-[ViewController outputAccelertion:]:lat:%f;lon:%f;heading:%f;distance:%f');
   data=sscanf(gpsline,'[Info] 2016-03-23 %*s (ViewController.m:%*d)-[ViewController outputAccelertion:]:lat:%f;lon:%f;heading:%f;distance:%f');
   if(isempty(data))
       break;
   end
   result=lonLat2Mercator(data(2,1),data(1,1));
   gx(n)=result.X;%GPS经过坐标变换后的东向坐标，换算成米数
   gy(n)=result.Y;%GPS经过坐标变换后的北向坐标，换算成米数
   Phi(n)=deg2rad(data(3,1));%航向角
   dd(n)=data(4,1);%某一周期的位移
   dx(n)=dd(n)*sin(Phi(n));%某一周期的东向位移
   dy(n)=dd(n)*cos(Phi(n));%某一周期的北向位移
   Ve(n)=dd(n)*sin(Phi(n));%里程计输入的东向速度，暂时用某一周期的东向位移代替
   Vn(n)=dd(n)*cos(Phi(n));%里程计输出的北向速度，暂时用某一周期的北向位移代替
   ZE(:,n)=[gx(n),gy(n),dx(n),dy(n)];
end
fclose(fgps);%%%%%关闭文件指针

%过程向量A
Ae=1;
An=1;
A=[1,0,T,0,0,0,0;
    0,1,0,T,0,0,0;
    0,0,1,0,0,0,0;
    0,0,0,1,0,0,0;
    0,0,0,0,Ae,0,0;
    0,0,0,0,0,An,0;
    0,0,0,0,0,0,1;
   ];

%过程噪声
W3=CreateGauss(0,d,1,N);
W4=CreateGauss(0,d,1,N);
W5=CreateGauss(0,d,1,N);
W6=CreateGauss(0,d,1,N);
W7=CreateGauss(0,d,1,N);
W=[0,0,W3,W4,W5,W6,W7]';
%过程噪声协方差矩阵
Q=diag([0,0,d^2,d^2,d^2,d^2,d^2]);

%观测向量C

%观测噪声
nv1=CreateGauss(0,d,1,N);
nv2=CreateGauss(0,d,1,N);
nv3=CreateGauss(0,d,1,N);
nv4=CreateGauss(0,d,1,N);
V=[nv1,nv2,nv3,nv4]';

Theta=CreateGauss(0,d,1,N);%GPS航迹和DR航迹的夹角
%观测噪声协方差矩阵
R=diag([d^2,d^2,d^2,d^2]);

Xest=zeros(7,1); % 用前k-1时刻的输出值估计k时刻的预测值
Xfli=zeros(7,1); % k时刻Kalman滤波器的输出值
Xes=zeros(7,1); % 预测输出误差
Xef=zeros(7,1); % 滤波后输出的误差
Pxe=zeros(7,7); % 预测输出误差均方差矩阵 
Px=zeros(7,7);  % 滤波输出误差均方差矩阵

XE(:,1)=[gx(1),gy(1),0,0,0,0,0]';

Xfli=[gx(1),gy(1),0,0,0,0,0]'; %初始条件进行估计
%Xef=[-vx(2) Ts*W(1)/2+(vx(1)-vx(2))/Ts]'; % 滤波后输出的误差
Px=diag([0,0,d^2,d^2,d^2,d^2,d^2]); % 滤波输出误差均方差矩阵

for k=2:N
    C=[1,0,0,0,1,0,0;
        0,1,0,0,0,1,0;
        0,0,cos(Theta(k)),-sin(Theta(k)),0,0,-Ve(k)*sin(Theta(k))-Vn(k)*cos(Theta(k));
        0,0,sin(Theta(k)),cos(Theta(k)),0,0,Ve(k)*cos(Theta(k))-Vn(k)*sin(Theta(k));
        ];
    
     Xest=A*Xfli; % 更新该时刻的预测值 ---kalman equation1
     %Xes=A*Xef+Gamma*W(k-1); % 预测输出误差 
     Pxe=A*Px*A'+Q; % 预测误差的协方差阵 ---kalman equation2
     
     K=Pxe*C'/(C*Pxe*C'+R); % Kalman滤波增益 ---kalman equation3
     Z=ZE(:,k);
     Xfli=Xest+K*(Z-C*Xest);% k时刻Kalman滤波器的输出值 ---kalman equation4
     %Xef=(eye(2)-K*C)*Xes-K*vx(k);%滤波后输出的误差 
     Px=(eye(7)-K*C)*Pxe;%滤波输出误差均方差矩阵 ---kalman equation5
     
     X=Xfli;
     XE(:,k)=X;
end

randTheta=rand(1,400)';
[x,y]=TrueRoute(T,offtime);
%显示滤波轨迹
figure
plot(x,y,'r');hold on;
plot(ZE(1,:),ZE(2,:),'g');hold on;
plot(XE(1,:),XE(2,:),'b');hold off;
axis([1.351825958322876e+07-50 1.351825958322876e+07+50 3.636297153803221e+06-50 3.636297153803221e+06+50]),grid on;
legend('真实轨迹','观测轨迹','目标滤波航迹');
