%% tune and test
% % path and file selection attn_rf20140328105051
clear; clc;
data_name = 'attn_rf20131107105335';
path_data = '/DATA/teus_philips/philips_data_flattened/';
path_pred = '/home/sedghi/';
path_mask = '/home/sedghi/projects/teus_em/data/eda/echo_prostate_mask_matlab/';
path_save = '/home/sedghi/';
if ~exist(path_save, 'dir')
   mkdir(path_save)
end

load([path_pred,data_name,'_pred.mat']);

if strcmp(data_name(1:4),'attn')
    data_name = data_name(6:end);
end

% % load data
load([path_data,data_name,'.mat']);
try
    load([path_data,data_name,'_wholetargetimage.mat']);
catch
    load(['/DATA/teus_philips/WholeTarget/',data_name,'_wholetargetimage.mat']);
end

% load([path_data,data_name,'_targetimage.mat']);
load([path_mask,data_name,'mask.mat']);

pred_benign = double(pred(:,:,1));
pred_cancer = double(pred(:,:,2));

first_frame = rf_AllFrames(:,:,1);
first_frame = (log(1+abs(first_frame)));

% downsampling
rf_x = [];
rf_y = [];
rf_first_frame = [];
for i=1:10:10*256 %ns-50+1
    temp = rf_TargetImageX(i:i+50-1,:);
    rf_x = [rf_x; mean(temp)];
    temp = rf_TargetImageY(i:i+50-1,:);
    rf_y = [rf_y; mean(temp)];
    temp = first_frame(i:i+50-1,:);
    rf_first_frame = [rf_first_frame; mean(temp)];
end

% axial tuning for rf-to-bmode transformation
n_dist = 6;
rf_first_frame = [zeros(n_dist,size(rf_first_frame,2));rf_first_frame(1:end-n_dist,:)];
rf_pred_benign = [zeros(n_dist,size(pred_benign,2));pred_benign(1:end-n_dist,:)];
rf_pred_cancer = [zeros(n_dist,size(pred_cancer,2));pred_cancer(1:end-n_dist,:)];


% echo parameters
d = echo_TargetImage;
echo_x = d(echo_TargetXY(2), :);
echo_x(1:echo_TargetXY(1)) = -echo_x(1:echo_TargetXY(1));
echo_y = d(:, echo_TargetXY(1));
echo_y(1:echo_TargetXY(2)) = -echo_y(1:echo_TargetXY(2));
[echo_x, echo_y] = meshgrid(echo_x, echo_y);
echo_d = sqrt(echo_x.^2 + echo_y.^2);

[size_y, size_x] = size(echo_d);
echo_pred_benign = zeros(size(echo_d));
echo_pred_cancer = zeros(size(echo_d));
rf_transformed = zeros(size(echo_d));

% interpolation parameter (low==sparse high==blure)
alpha = 5;
diff_x = max(abs(diff(sort(rf_x(:))))) *alpha;
diff_y = max(abs(diff(sort(rf_y(:))))) *alpha;
% diff_x = 1;
% diff_y = 1;

% rf-to-bmode transformation
parfor i_x = 1:size_x
    %disp(i_x);
    for i_y = 1:size_y
        x = echo_x(i_y,i_x);
        y = echo_y(i_y,i_x);
        rf_ind = logical((abs(rf_x-x)<diff_x).*(abs(rf_y-y)<diff_y));
        echo_pred_benign(i_y,i_x) = mean(rf_pred_benign(rf_ind));
        echo_pred_cancer(i_y,i_x) = mean(rf_pred_cancer(rf_ind));
        rf_transformed(i_y,i_x) = mean(rf_first_frame(rf_ind));
    end
end

echo_first_frame = echo_AllFrames(end:-1:1,:,1);
echo = im2double(echo_first_frame)/255;
echo_first_frame = cat(3,echo,echo,echo);

% figure; imagesc(rf_transformed); colormap(gray)
% figure; imagesc(echo_first_frame); colormap(gray)
% figure; imagesc(echo_pred_benign); colormap(jet)
% figure; imagesc(echo_pred_cancer); colormap(jet)
%%
% % prostate mask smoothing
% boundary

mask1 = bwmorph(mask,'remove');
% figure; imshow(mask1)
% image to polar
[row,col] = find(mask1);
center = [mean(row) mean(col)];
[theta,rho] = cart2pol(col-mean(col),row-mean(row));
[~,ind]=sort(theta);
theta = theta(ind);
rho = rho(ind);
% figure; plot(theta,rho)
TH = [theta-2*pi; theta; theta+2*pi];
RH = [rho; rho; rho];
ffun = fit(TH,RH,'fourier6');
RH = feval(ffun,TH);
rho1 = RH(length(rho)+1:2*length(rho),1);
[x,y] = pol2cart(theta,rho1);
x = x + mean(col);
y = y + mean(row);
% figure; plot(x,y)
ind = sub2ind(size(mask1),round(y),round(x));
mask2 = zeros(size(mask1));
mask2(ind) = 1;

% generate smooth mask
morph_str = strel('disk',10);
mask3 = imdilate(mask2, morph_str);
mask3 = imfill(mask3,'holes');
smooth_mask = imerode(mask3, morph_str);
figure; imshowpair(mask1,smooth_mask)

% % visualization and save
% map_visualize(echo_first_frame,smooth_mask,echo_pred_benign,echo_TargetXY)
% saveas(gcf,[path_save,data_name{i_data},'_benign.png'])

map_visualize(echo_first_frame,smooth_mask,echo_pred_cancer,echo_TargetXY);
% saveas(gcf,[path_save,data_name{i_data},'_cancer.png'])