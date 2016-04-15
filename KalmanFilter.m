%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�����ʼ������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ȫ�ֱ�������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outdoor_sensor_data=370;
indoor_sensor_data=0;
sensor_data=outdoor_sensor_data+indoor_sensor_data;
d=0.1;%��׼��
T=1;%����ʱ��
Theta=CreateGauss(0,d,1,sensor_data);%GPS������DR�����ļн�
ZOUT=zeros(4,outdoor_sensor_data);
ZIN=zeros(4,indoor_sensor_data);
Z=zeros(2,outdoor_sensor_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ȡ����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fgps=fopen('sensor_data_0415.txt','r');%%%���ı�

for n=1:sensor_data
    gpsline=fgetl(fgps);%%%��ȡ�ı�ָ���Ӧ����
    if ~ischar(gpsline) break;%%%�ж��Ƿ����
    end;
   %%%%��ȡ��������
   time=sscanf(gpsline,'[Info] 2016-04-15%s(ViewController.m:%d)-[ViewController outputAccelertion:]:lat:%f;lon:%f;heading:%f;distance:%f;beacon_lat:%f;beacon_lon:%f');
   data=sscanf(gpsline,'[Info] 2016-04-15 %*s (ViewController.m:%*d)-[ViewController outputAccelertion:]:lat:%f;lon:%f;heading:%f;distance:%f;beacon_lat:%f;beacon_lon:%f');
   if(isempty(data))
       break;
   end
        result=lonLat2Mercator(data(2,1),data(1,1));
        gx(n)=result.X;%GPS��������任��Ķ������꣬���������
        gy(n)=result.Y;%GPS��������任��ı������꣬���������
        Phi(n)=(data(3,1)+90)*pi/180;%�����
        dd(n)=data(4,1);%ĳһ���ڵ�λ��
        dx(n)=dd(n)*sin(Phi(n))*4;%ĳһ���ڵĶ���λ��
        dy(n)=dd(n)*cos(Phi(n))*4;%ĳһ���ڵı���λ��
        Ve(n)=dd(n)*sin(Phi(n));%��̼�����Ķ����ٶȣ���ʱ��ĳһ���ڵĶ���λ�ƴ���
        Vn(n)=dd(n)*cos(Phi(n));%��̼�����ı����ٶȣ���ʱ��ĳһ���ڵı���λ�ƴ���
        ZOUT(:,n)=[gx(n),gy(n),dx(n),dy(n)];
        Z(:,n)=[gx(n),gy(n)];
end
fclose(fgps);%%%%%�ر��ļ�ָ��

%��������A
A=[
    1,0,T,0,0;
    0,1,0,T,0;
    0,0,1,0,0,;
    0,0,0,1,0;
    0,0,0,0,1;
   ];
%��������Э�������
Q=diag([0,0,d^2,d^2,d^2]);
%�۲�����Э�������
R=diag([d^2,d^2,d^2,d^2]);
P=diag([0,0,d^2,d^2,d^2]); % �˲���������������
Xfli=[ZOUT(1,1),ZOUT(2,1),0,0,0]'; %��ʼ�������й���
for k=1:sensor_data
    C=[1,0,0,0,0;
        0,1,0,0,0;
        0,0,cos(Theta(k)),-sin(Theta(k)),-Ve(k)*sin(Theta(k))-Vn(k)*cos(Theta(k));
        0,0,sin(Theta(k)),cos(Theta(k)),Ve(k)*cos(Theta(k))-Vn(k)*sin(Theta(k));
        ];
     K_location(:,k)=Xfli;
     Xest=A*Xfli; % ���¸�ʱ�̵�Ԥ��ֵ ---kalman equation1
     %Xes=A*Xef+Gamma*W(k-1); % Ԥ�������� 
     Pxe=A*P*A'+Q; % Ԥ������Э������ ---kalman equation
     K=Pxe*C'/(C*Pxe*C'+R); % Kalman�˲����� ---kalman equation3
     Xfli=Xest+K*(ZOUT(:,k)-C*Xest);% kʱ��Kalman�˲��������ֵ ---kalman equation4
     Px=(eye(5)-K*C)*Pxe;%�˲��������������� ---kalman equation5
end

cordinatex=ZOUT(1,5);
cordinatey=ZOUT(2,5);
%��ʾ�˲��켣
figure
set(gca,'FontSize',12);
[groundtruthx,groundtruthy]=Groud_Truth();
plot(groundtruthx,groundtruthy,'r');hold on;
plot(ZOUT(1,:),ZOUT(2,:),'o');hold on;
plot(K_location(1,:),K_location(2,:),'g');hold off;
axis([cordinatex-100 cordinatex+200 cordinatey-200 cordinatey+100]),grid on;
xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);
legend('��ʵ�켣','�۲�켣','Ŀ���˲�����');
axis equal;

    

