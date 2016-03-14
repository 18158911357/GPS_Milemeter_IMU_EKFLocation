function[X,Y] = TrueRoute(Ts,offtime)
%������ʵ����[X,Y]��������ʵ����ϵ�»���
%TsΪgps�۲����ڣ�ÿ��Tsʱ��õ�һ���۲�����
%offtimeΪ�����ʱ�䣬���ٴ��Ϊ0-8m/S��΢����ѧԺ��¥�Ļ��ƾ����Լ400m��������1m/s���ߣ�ʱ��Ϊ400s

if nargin>2
	error('����ı������࣬����');
end

x=zeros(offtime,1);
y=zeros(offtime,1);
X=zeros(ceil(offtime/Ts),1);
Y=zeros(ceil(offtime/Ts),1);

% t=0:100s,�ٶ�vx��vyΪ����x���y��ķ���(m/s)
x0=50;
y0=50;
vx=1;%�ٶ�Ϊ1m/s,������x��
vy=0;

for t=1:100
	x(t)=x0+vx*t;
	y(t)=y0+vy*t;
end

% t=100:200s,�ٶ�vx��vyΪ����x���y��ķ���(m/s)

vx=0;%�ٶ�Ϊ1m/s,������x��
vy=1;

for t=0:100
	x(t+101)=x(100)+vx*t;
	y(t+101)=y(100)+vy*t;
end


% t=200:300s,�ٶ�vx��vyΪ����x���y��ķ���(m/s)

vx=-1;%�ٶ�Ϊ1m/s,������x��
vy=0;

for t=0:100
	x(t+201)=x(201)+vx*t;
	y(t+201)=y(201)+vy*t;
end

% t=300:400s,�ٶ�vx��vyΪ����x���y��ķ���(m/s)

vx=0;%�ٶ�Ϊ1m/s,������x��
vy=-1;

for t=0:100
	x(t+301)=x(301)+vx*t;
	y(t+301)=y(301)+vy*t;
end

% �õ��۲�����
for n=0:Ts:offtime
    X(n/Ts+1)=x(n+1);
    Y(n/Ts+1)=y(n+1);
end

% %��ʾ��ʵ�켣
% plot(X,Y,'LineWidth',2),axis([0 200 0 200]),grid on;
% legend('Ŀ����ʵ����');
