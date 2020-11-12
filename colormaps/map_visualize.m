function h1 = map_visualize(echo_first_frame,smooth_mask,echo_pred_benign,echo_TargetXY)
figure;
imagesc(echo_first_frame); 
hold on;
h1 = imagesc(smooth_mask.*echo_pred_benign);
set(h1,'AlphaData',0.5*smooth_mask)
% alpha_ch = smooth_mask.*abs(echo_pred_benign-0.5)*1;
% alpha_ch = smooth_mask.*(abs(echo_pred_benign-0.5)>0)*0.5;
% set(h1,'AlphaData',alpha_ch)

% % cl0
% CMap = parula(128);
CMap = jet(128);

% % cl1
% gradient = 0:0.01:1;
% fgradient = fliplr(gradient);
% CMap = [gradient', zeros(size(gradient,2),1), fgradient'];

% % cl2
% gradient = 0.1:0.01:1;
% fgradient = fliplr(gradient);
% c1 = [zeros(size(gradient,2),1), zeros(size(gradient,2),1), fgradient'];
% c2 = [gradient', zeros(size(gradient,2),1), zeros(size(gradient,2),1)];
% % CMap = [c1; repmat([1,1,0],10,1); c2];
% CMap = [c1; c2];

colormap(CMap)
colorbar;
axis equal;
axis off;
ylim([1-0.5, size(echo_first_frame,1)+0.5])
% caxis([0 1])
a = echo_pred_benign(:); a = a(logical(smooth_mask(:)));
% caxis([min(a) max(a)]) % for min and max of each image
caxis([0 1]) % colormap for between 0 1

smooth_boundary = bwmorph(smooth_mask,'remove');
% smooth_boundary = imdilate(smooth_boundary,strel('disk',1));
h2 = imagesc(max(a)*smooth_boundary);
set(h2,'AlphaData',smooth_boundary)


viscircles(echo_TargetXY,5,'color','y','linewidth',1);
hold off;
end