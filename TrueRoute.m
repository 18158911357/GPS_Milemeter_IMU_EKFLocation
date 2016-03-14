function[X,Y] = TrueRoute(Ts,offtime)
%产生真实航迹[X,Y]，并在真实坐标系下绘制
%Ts为gps观测周期，每隔Ts时间得到一个观测数据
%offtime为仿真的时间，车速大概为0-8m/S，微电子学院大楼的环绕距离大约400m，按车速1m/s来走，时间为400s

if nargin>2
	error('输入的变量过多，请检查');
end

x=zeros(offtime,1);
y=zeros(offtime,1);
X=zeros(ceil(offtime/Ts),1);
Y=zeros(ceil(offtime/Ts),1);

% t=0:100s,速度vx，vy为沿着x轴和y轴的分量(m/s)
x0=50;
y0=50;
vx=1;%速度为1m/s,先沿着x轴
vy=0;

for t=1:100
	x(t)=x0+vx*t;
	y(t)=y0+vy*t;
end

% t=100:200s,速度vx，vy为沿着x轴和y轴的分量(m/s)

vx=0;%速度为1m/s,先沿着x轴
vy=1;

for t=0:100
	x(t+101)=x(100)+vx*t;
	y(t+101)=y(100)+vy*t;
end


% t=200:300s,速度vx，vy为沿着x轴和y轴的分量(m/s)

vx=-1;%速度为1m/s,先沿着x轴
vy=0;

for t=0:100
	x(t+201)=x(201)+vx*t;
	y(t+201)=y(201)+vy*t;
end

% t=300:400s,速度vx，vy为沿着x轴和y轴的分量(m/s)

vx=0;%速度为1m/s,先沿着x轴
vy=-1;

for t=0:100
	x(t+301)=x(301)+vx*t;
	y(t+301)=y(301)+vy*t;
end

% 得到观测数据
for n=0:Ts:offtime
    X(n/Ts+1)=x(n+1);
    Y(n/Ts+1)=y(n+1);
end

% %显示真实轨迹
% plot(X,Y,'LineWidth',2),axis([0 200 0 200]),grid on;
% legend('目标真实航迹');
