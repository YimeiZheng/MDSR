close all;clear;clc;
Save_path='E:\MDSR\';
folderPath1 = 'E:\sourceimage\T1+T2\T1';
folderPath2 = 'E:\sourceimage\T1+T2\T2';
imageFiles1 = dir(fullfile(folderPath1)); 
imageFiles2 = dir(fullfile(folderPath2));
numImages = length(imageFiles1);
times = zeros(1, numImages);imagesA = cell(1, numImages);imagesB = cell(1, numImages);
t=0;
for c=1:numImages    
    imagePath = fullfile(folderPath1, imageFiles1(c).name);
    imagePath2 = fullfile(folderPath2, imageFiles2(c).name);
      if isequal(imageFiles1(c).name, '.') || isequal(imageFiles1(c).name, '..')
        continue;
      end
imagesA{c} = imread(imagePath2);
imagesB{c} = imread(imagePath);
%%
j=0;
I1 = imagesA{c};
I2 = imagesB{c};  

Fusion_Image=XDOGFuse(I1,I2);

figure;imshow(I1);figure;imshow(I2)
 
figure;imshow(Fusion_Image);figure;imshow(I2)
   path=strcat(Save_path,num2str(c+11),'.tif');
   imwrite(Fusion_Image,path,'TIF')
   clear  Fusion_Image   I1  I2 
end

function Fusion_image=XDOGFuse(I1,I2)
D=cell2mat(struct2cell(load('D12.mat')));  %
tic
A =im2double(I1);
B =im2double(I2);

B_YUV=ConvertRGBtoYUV(B);
B_Y=B_YUV(:,:,1); 
%% decomposition
AL = PNLS_DT(A);
BL = PNLS_DT(B);
AH=A-AL;
BH=B_Y-BL;

AL=rgb2gray(AL);
BL=rgb2gray(BL);

f1_l=AL;
f1_h=AH;
f2_l=BL;
f2_h=BH;

%% 高频
filter=pyramid_filter;
D1=downsample(D,filter);%字典采样一次
D2= downsample(D1,filter);%字典采样两次

f1_h1 = downsample(f1_h,filter);
f1_h2 = downsample(f1_h1,filter);

f2_h1 = downsample(f2_h,filter);
f2_h2 = downsample(f2_h1,filter);


%% Fh fusion
sigma=0;
if sigma>0 
    v=sigma*sigma/(255*255); 
    A=imnoise(A,'gaussian',0,v);
    B=imnoise(B,'gaussian',0,v);  
end

overlap=7;      
C=0.0035;       

if sigma==0
    epsilon=0.01;   
else
    epsilon=0.05++8*C*sigma;
end
%%
for i=1:3
     Fh1(:,:,i) =sparse_fusion3(f1_h(:,:,1),f2_h(:,:,1),D,overlap,epsilon);         %B为灰度图
end
Fh1=rgb2gray(Fh1);

 for i=1:3
     Fh2(:,:,i) =sparse_fusion3(f1_h1(:,:,1),f2_h1(:,:,1),D1,overlap,epsilon);         %B为灰度图
 end

Fh2=rgb2gray(Fh2);
Fh2=fchazhi(Fh2,[256,256]);  
 
 for i=1:3
     Fh3(:,:,i) =sparse_fusion3(f1_h2(:,:,1),f2_h2(:,:,1),D2,overlap,epsilon);         %B为灰度图
 end

Fh3=rgb2gray(Fh3);
Fh3=fchazhi(Fh3,[256,256]);  

%% 低频融合
map=abs(f1_l>f2_l);
FL=(f1_l.*map+~map.*f2_l);

%% 差图
Diff1A=abs((f1_h- Fh1));
Diff2A=abs((f1_h-Fh2));

Diff1B=abs((f2_h-Fh1));
Diff2B=abs((f2_h-Fh2));

%%

Diff1A=rgb2gray(Diff1A);
Diff1B=rgb2gray(Diff1B);
A=rgb2gray(A);

T=0.001;
[row,column]=size(A);
DiffC=zeros(row,column);
for i=1:row
    for  j=1:column
         DiffC(i,j)=abs(Diff1A(i,j)-Diff1B(i,j))+abs(Diff2A(i,j)-Diff2B(i,j));
    end
end
for i=1:row
    for  j=1:column
          if DiffC(i,j)>T
              FH(i,j)=f1_h(i,j);
          else
              FH(i,j)=f2_h(i,j); 
          end
    end
end

F=FH+FL;
Fusion_image=F;
toc;
end

function Compare_Size(image_1,image_2)
    if ( size(image_1) ~=  size(image_2)) 
        error('Input images are not of same size');
    end  
end
