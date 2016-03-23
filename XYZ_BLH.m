function BLH=XYZ_BLH(X,Y,Z) 
% Transform the (X Y Z) into (B L H) under the WGS84 coordinate system 

a=6378137.00;            % 长半轴，全部参考WGS-84椭球参数表 
b=6356752.3142;          % 短半轴 
alfa=1/298.257223563;    % 扁率 
e2=0.00669437999013;     % 第一偏心率的平方 
% L=atan(Y/X); 
L=Azimuth(X,Y);          % 经度，即水平方向方位角 
N0=a; 
H0=sqrt(X^2+Y^2+Z^2)-sqrt(a*b); 
B0=atan(Z/sqrt(X^2+Y^2)/(1.0-e2*N0/(N0+H0))); 
while(true) 
    N1=a/sqrt(1.0-e2*sin(B0)^2); 
    H1=sqrt(X^2+Y^2)/cos(B0)-N0; 
    B1=atan(Z/sqrt(X^2+Y^2)/(1.0-e2*N1/(N1+H1))); 
    if (abs(H1-H0)>0.001 & abs(B1-B0)>4.848132257047973e-11)    % B的精度为0.00001" 
        N0=N1; 
        H0=H1; 
        B0=B1; 
    else 
        break; 
    end 
end 
BLH.B=B1; 
BLH.L=L; 
BLH.H=H1; 
return