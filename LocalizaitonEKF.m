% function XE=LocationEKF(T,offtime,d)
%LocationEKF     ���ÿ������˲��������ӹ۲���ֵ�еõ����������Ź���
%XE              ���x�᷽���ϵ����
%T              ����ʱ�䣬gps�۲�����
%offtime         ����ʱ��
%d               �����ı�׼��
clc;
clear;
T=1;
offtime=10;
d=0.1;

close all
N=ceil(offtime/T); %��������

A=zeros(7,7);
W=zeros(7,1);
X=zeros(7,1); % һ�ε��˲����ֵ
C=zeros(4,7);
V=zeros(4,1);
Z=zeros(4,1);%һ�εĹ۲�ֵ
XE=zeros(7,N);%���е�Ԥ��ֵ
Ve=zeros(N,1);
Vn=zeros(N,1);
ZE=zeros(4,N);

randn('state',sum(100*clock)); % ���������������
%%%%%%%%%%%��ȡ�ı�gpsֵ����̼�λ�ƺͺ���ǣ���ת��Ϊֱ������ϵ�µ�����ֵ%%%%%%
fgps=fopen('VelocityDetector.txt','r');%%%���ı�
n=0;
while 1
   gpsline=fgetl(fgps);%%%��ȡ�ı�ָ���Ӧ����
   if ~ischar(gpsline) break;%%%�ж��Ƿ����
   end;
   n=n+1;
   time=sscanf(gpsline,'[Info] 2016-03-23%s(ViewController.m:%d)-[ViewController outputAccelertion:]:lat:%f;lon:%f;heading:%f;distance:%f');
   data=sscanf(gpsline,'[Info] 2016-03-23 %*s (ViewController.m:%*d)-[ViewController outputAccelertion:]:lat:%f;lon:%f;heading:%f;distance:%f');
   if(isempty(data))
       break;
   end
   result=lonLat2Mercator(data(2,1),data(1,1));
   gx(n)=result.X;%GPS��������任��Ķ������꣬���������
   gy(n)=result.Y;%GPS��������任��ı������꣬���������
   Phi(n)=deg2rad(data(3,1));%�����
   dd(n)=data(4,1);%ĳһ���ڵ�λ��
   dx(n)=dd(n)*sin(Phi(n));%ĳһ���ڵĶ���λ��
   dy(n)=dd(n)*cos(Phi(n));%ĳһ���ڵı���λ��
   Ve(n)=dd(n)*sin(Phi(n));%��̼�����Ķ����ٶȣ���ʱ��ĳһ���ڵĶ���λ�ƴ���
   Vn(n)=dd(n)*cos(Phi(n));%��̼�����ı����ٶȣ���ʱ��ĳһ���ڵı���λ�ƴ���
   ZE(:,n)=[gx(n),gy(n),dx(n),dy(n)];
end
fclose(fgps);%%%%%�ر��ļ�ָ��

%��������A
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

%��������
W3=CreateGauss(0,d,1,N);
W4=CreateGauss(0,d,1,N);
W5=CreateGauss(0,d,1,N);
W6=CreateGauss(0,d,1,N);
W7=CreateGauss(0,d,1,N);
W=[0,0,W3,W4,W5,W6,W7]';
%��������Э�������
Q=diag([0,0,d^2,d^2,d^2,d^2,d^2]);

%�۲�����C

%�۲�����
nv1=CreateGauss(0,d,1,N);
nv2=CreateGauss(0,d,1,N);
nv3=CreateGauss(0,d,1,N);
nv4=CreateGauss(0,d,1,N);
V=[nv1,nv2,nv3,nv4]';

Theta=CreateGauss(0,d,1,N);%GPS������DR�����ļн�
%�۲�����Э�������
R=diag([d^2,d^2,d^2,d^2]);

Xest=zeros(7,1); % ��ǰk-1ʱ�̵����ֵ����kʱ�̵�Ԥ��ֵ
Xfli=zeros(7,1); % kʱ��Kalman�˲��������ֵ
Xes=zeros(7,1); % Ԥ��������
Xef=zeros(7,1); % �˲�����������
Pxe=zeros(7,7); % Ԥ���������������� 
Px=zeros(7,7);  % �˲���������������

XE(:,1)=[gx(1),gy(1),0,0,0,0,0]';

Xfli=[gx(1),gy(1),0,0,0,0,0]'; %��ʼ�������й���
%Xef=[-vx(2) Ts*W(1)/2+(vx(1)-vx(2))/Ts]'; % �˲�����������
Px=diag([0,0,d^2,d^2,d^2,d^2,d^2]); % �˲���������������

for k=2:N
    C=[1,0,0,0,1,0,0;
        0,1,0,0,0,1,0;
        0,0,cos(Theta(k)),-sin(Theta(k)),0,0,-Ve(k)*sin(Theta(k))-Vn(k)*cos(Theta(k));
        0,0,sin(Theta(k)),cos(Theta(k)),0,0,Ve(k)*cos(Theta(k))-Vn(k)*sin(Theta(k));
        ];
    
     Xest=A*Xfli; % ���¸�ʱ�̵�Ԥ��ֵ ---kalman equation1
     %Xes=A*Xef+Gamma*W(k-1); % Ԥ�������� 
     Pxe=A*Px*A'+Q; % Ԥ������Э������ ---kalman equation2
     
     K=Pxe*C'/(C*Pxe*C'+R); % Kalman�˲����� ---kalman equation3
     Z=ZE(:,k);
     Xfli=Xest+K*(Z-C*Xest);% kʱ��Kalman�˲��������ֵ ---kalman equation4
     %Xef=(eye(2)-K*C)*Xes-K*vx(k);%�˲����������� 
     Px=(eye(7)-K*C)*Pxe;%�˲��������������� ---kalman equation5
     
     X=Xfli;
     XE(:,k)=X;
end

randTheta=rand(1,400)';
[x,y]=TrueRoute(T,offtime);
%��ʾ�˲��켣
figure
plot(x,y,'r');hold on;
plot(ZE(1,:),ZE(2,:),'g');hold on;
plot(XE(1,:),XE(2,:),'b');hold off;
axis([1.351825958322876e+07-50 1.351825958322876e+07+50 3.636297153803221e+06-50 3.636297153803221e+06+50]),grid on;
legend('��ʵ�켣','�۲�켣','Ŀ���˲�����');
