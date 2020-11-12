clear all
close all
w = 756;
h =616;

%Fill these information
Patient=100;
CoreID=8;
data_dir=strcat('z:\\shared\images\ProstateVGH-2\Data\Patient',num2str(Patient));
dirlist=dir(data_dir);
dir = strcat(data_dir,'\',dirlist(3).name,'\');
BM_file='Bmode_A8_LMM_469.dat';
save_dir=strcat(dir,'BMode\ROI_Data\WholeProstate\');


fid = fopen(strcat(dir,BM_file));
A = fread(fid,inf,'uchar');
l =length(A);
N = length(A)/(w*h);
fclose(fid);

% reshape into a matrix with shape {w,h,N}
frames = reshape(A,[w,h,N]);
rotframes= rot90(frames, -1);


h=imagesc(squeeze(rotframes(:,:,end)));colormap(gray);
title([num2str(Patient),' +++  ',num2str(CoreID)])

finished = 'NO';
i = 1;
while strcmpi(finished,'NO')
  hFH(i) = imfreehand();
  finished = questdlg('Finished?', ...
      'confirmation', ...
      'YES', 'NO', 'UNDO', 'NO');
  if strcmpi(finished, 'UNDO')
      delete(hFH(i))
      finished = 'NO';
  else
      mask{i}  = hFH(i).createMask(h);
      i = i + 1;
  end
end
xy_mask = cell2mat(mask);
indx=find(xy_mask);

BM_mask_file=strcat(save_dir,strrep(BM_file,'Bmode','BM_mask'));
BM_mask_file=strrep(BM_mask_file,'.dat','.bmp');
imwrite(xy_mask,BM_mask_file,'bmp')
