function G=CreateGauss(E,D,M,N)
%������ֵΪE,����ΪD,MxN�ĸ�˹����������
%
    G=randn(M,N); %����m*n����������ֲ�Ϊ��ֵΪ0 ������Ϊ1
    G=G/std(G); %���Զ�G����ȡ��׼��
    G=G-mean(G); %��ȥ����ȡ��ֵ
    a=E; %��ֵ
    b=sqrt(D); %�õ���׼��
    G=a+b*G
    %c=var(G);
    %d=mean(G);
    %plot(G);
end
