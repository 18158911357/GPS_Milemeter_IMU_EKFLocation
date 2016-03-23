function G=CreateGauss(E,D,M,N)
%产生均值为E,方差为D,MxN的高斯白噪声矩阵
%
    G=randn(M,N); %生产m*n的随机数，分布为均值为0 ，方差为1
    G=G/std(G); %除以对G的列取标准差
    G=G-mean(G); %减去对列取均值
    a=E; %均值
    b=sqrt(D); %得到标准差
    G=a+b*G
    %c=var(G);
    %d=mean(G);
    %plot(G);
end
