function Dksvd=Dictionary11(I1,I2)
% Patchesget.m���ص�PE��PC�������������Ϣ�ḻ�Ŀ�
% ���������Ĳ�ֱͬ��Ӱ���ֵ�ѧϰʱ�䣬�����������о�������Ϊ����������ʵ���н����о�
addpath ksvdbox13;
addpath ompbox10;

k=6;    % �ظ�����Ϊ6
n0=8;     %��Ĵ�С
N=256;     
[P_E,P_C]=Patchesget(I1,I2,k,n0);    %  �п���ֱ������
P0=[P_E,P_C];

set=P0;
set=double(set/255);
params.data=(set);
[r,z]=size(set);
%params.data=X;
%params.data=X;
params.Tdata =20;
params.dictsize = N;
params.iternum =180;               % ��������Խ���ֵ�ѧϰʱ��Խ�������ս�������������һ����ϵ
params.memusage = 'high';
tic
[Dksvd,g,err] = ksvd(params,'');   %% Dictionary11.m���ص�Dksvd�������ѵ���õ����ֵ�
toc;