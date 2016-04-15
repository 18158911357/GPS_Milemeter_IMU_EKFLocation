%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�����ʼ������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PCenter=ParticleFilter()
% clc;
% clear;
% close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ȫ�ֱ�������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outdoor_sensor_data=361;
indoor_sensor_data=0;
sensor_data=outdoor_sensor_data+indoor_sensor_data;
d=0.1;%��׼��
Theta=CreateGauss(0,d,1,sensor_data);%GPS������DR�����ļн�
ZOUT=zeros(4,outdoor_sensor_data);
ZIN=zeros(4,indoor_sensor_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ȡ����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fgps=fopen('sensor_data_041518.txt','r');%%%���ı�

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
        result=lonLat2Mercator(data(6,1),data(5,1));
        gx(n)=result.X;%GPS��������任��Ķ������꣬���������
        gy(n)=result.Y;%GPS��������任��ı������꣬���������
        Phi(n)=(data(3,1)+90)*pi/180;%�����
        dd(n)=data(4,1);%ĳһ���ڵ�λ��
        ZIN(:,n)=[gx(n),gy(n),Phi(n),dd(n)];
end
fclose(fgps);%%%%%�ر��ļ�ָ��

% ��������
N = 100;   %��������
Q = 5;      %��������
R = 5;      %��������
X = zeros(2, sensor_data);    %�洢ϵͳ״̬
Z = zeros(2, sensor_data);    %�洢ϵͳ�Ĺ۲�״̬
P = zeros(2, N);    %��������Ⱥ
PCenter = zeros(2, sensor_data);  %�������ӵ�����λ��
w = zeros(N, 1);         %ÿ�����ӵ�Ȩ��
err = zeros(1,sensor_data);     %���
X(:, 1) = [ZIN(1,1); ZIN(2,1)];     %��ʼϵͳ״̬
Z(:, 1) = [ZIN(1,1); ZIN(2,1)] ;    %��ʼϵͳ�Ĺ۲�״̬

cordinatex=round(ZIN(1,5));
cordinatey=round(ZIN(2,5));

%��ʼ������Ⱥ
for i = 1 : N
    P(:, i) = [randi([cordinatex-100,cordinatex+200],1);randi([cordinatey-200,cordinatey+100],1)];
    dist = norm(P(:, i)-Z(:, 1));     %�����λ�����ľ���
    w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %��Ȩ��
end
PCenter(:, 1) = sum(P, 2) / N;      %�������ӵļ�������λ��
 
%%
err(1) = norm(X(:, 1) - PCenter(:, 1));     %���Ӽ���������ϵͳ��ʵ״̬�����
% figure(1);
% set(gca,'FontSize',12);
% hold on
% plot(X(1, 1), X(2, 1), 'r.', 'markersize',30)   %ϵͳ״̬λ��
% axis([cordinatex-100 cordinatex+200 cordinatey-200 cordinatey+100]),grid on;
% plot(P(1, :), P(2, :), 'k.', 'markersize',5);   %��������λ��
% plot(PCenter(1, 1), PCenter(2, 1), 'b.', 'markersize',25); %�������ӵ�����λ��
% legend('True State', 'Particles', 'The Center of Particles');
% title('Initial State');
% hold off

%%
%��ʼ�˶�
for k = 2 : sensor_data
       
    %ģ��һ�������˶���״̬
    X(:, k) = X(:, k-1) + 1 * [(sin( ZIN(3,k))); cos(ZIN(3,k))] ;     %״̬����
    Z(:, k) = ZIN(1:2, k) ;     %�۲ⷽ�� 
    if ZIN(1,k) == 0
        Z(:, k) = X(:, k);     %�۲ⷽ�� 
    end
   
    %�����˲�
    %Ԥ��
    for i = 1 : N
        P(:, i) = P(:, i) + 1 * [sin( ZIN(3,k)); cos( ZIN(3,k))] + wgn(2, 1, 10*log10(Q));
        dist = norm(P(:, i)-Z(:, k));     %�����λ�����ľ���
        w(i) = (1 / sqrt(R) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / R);   %��Ȩ��
    end
%��һ��Ȩ��
    wsum = sum(w);
    for i = 1 : N
        w(i) = w(i) / wsum;
    end
   
    %�ز��������£�
    for i = 1 : N
        wmax = 2 * max(w) * rand;  %��һ���ز�������
        index = randi(N, 1);
        while(wmax > w(index))
            wmax = wmax - w(index);
            index = index + 1;
            if index > N
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     %�õ�������
    end
   
    PCenter(:, k) = sum(P, 2) / N;      %�������ӵ�����λ��
   
    %�������
    err(k) = norm(X(:, k) - PCenter(:, k));     %���Ӽ���������ϵͳ��ʵ״̬�����
   
%     figure(2);
%     set(gca,'FontSize',12);
%     clf;
%     hold on
%     plot(X(1, k), X(2, k), 'r.', 'markersize',50);  %ϵͳ״̬λ��
%     axis([cordinatex-100 cordinatex+200 cordinatey-200 cordinatey+100]),grid on;
%     plot(P(1, :), P(2, :), 'k.', 'markersize',5);   %��������λ��
%     plot(PCenter(1, k), PCenter(2, k), 'b.', 'markersize',25); %�������ӵ�����λ��
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
% legend('��ʵ�켣','�۲�켣', '�����˲��켣');
% xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20);
% axis equal;
