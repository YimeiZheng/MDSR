function Dksvd=Dictionary11(I1,I2)
% Patchesget.m返回的PE，PC就是能量块和信息丰富的块
% 迭代次数的不同直接影响字典学习时间，迭代次数的研究可以作为参数设置在实验中进行研究
addpath ksvdbox13;
addpath ompbox10;

k=6;    % 重复像素为6
n0=8;     %块的大小
N=256;     
[P_E,P_C]=Patchesget(I1,I2,k,n0);    %  切块拉直成向量
P0=[P_E,P_C];

set=P0;
set=double(set/255);
params.data=(set);
[r,z]=size(set);
%params.data=X;
%params.data=X;
params.Tdata =20;
params.dictsize = N;
params.iternum =180;               % 迭代次数越大，字典学习时间越长，最终结果与迭代次数有一定关系
params.memusage = 'high';
tic
[Dksvd,g,err] = ksvd(params,'');   %% Dictionary11.m返回的Dksvd就是最后训练得到的字典
toc;