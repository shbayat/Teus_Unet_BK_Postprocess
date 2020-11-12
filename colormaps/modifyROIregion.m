%% 
% % 
% points used in axial:  1:1536 averaged every 6 in xial direction
%
clear; clc;
close all;
% downsample_rate=6;
Patient=100;
CoreID=3; %CID:1091, core: 2
w = 756;
h =616;

data_dir=strcat('z:\\shared\images\ProstateVGH-2\Data\Patient',num2str(Patient));
dirlist=dir(data_dir);
data_dir = strcat(data_dir,'\',dirlist(3).name,'\');
fid = fopen([data_dir 'Bmode_A3_RMM_244.dat'],'r'); %p100
%fid = fopen([data_dir 'Bmode_A8_LML_284.dat'],'r'); %p91

 
  %fid = fopen( 'BMode_A4_RA_342.dat','r');
  % read input data. You need to specify the type that the data was saved as
   A = fread(fid,inf,'uchar');
   l =length(A);
   N = length(A)/(w*h);
   fclose(fid);
      frames = reshape(A,[w,h,N]);
   rotframes= rot90(frames, -1);
bm_backgraound=squeeze(rotframes(:,:,end));

imshow(uint8(bm_backgraound),'Colormap', gray)
img=imread([data_dir, '\BMode\ROI_Data\TestImages\BMROI_test_', int2str(CoreID),'.bmp']);
[x,y,v]=find(squeeze(img(:,:,1)==255));
min_x=min(x);
min_y=min(y);
lx=max(x)-min_x;
ly=max(y)-min_y;

hold on;
% Then, from the help:
rectangle('Position',[min_y,min_x,ly,lx],...
          'Curvature',[0.8,0.4],...
         'LineWidth',2,'LineStyle','--','EdgeColor', 'g')
RGB = insertShape(uint8(bm_backgraound), 'Rectangle', [min_y,min_x,ly,lx], 'Color',{'green'}, 'Opacity', 0.7,'LineWidth',2);
imwrite(RGB,'p100c3.bmp')