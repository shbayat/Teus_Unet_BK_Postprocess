%% tune and test
% % path and file selection attn_rf20140328105051
%points used in axial:  1:1536
%downsampling rate: 1536/6= 256
clear; clc;
cut_axial=1536;
downsample_rate=6;
data_name = 'RF_DS_A0_RBL_73'; %P91, core0
data_name_mask='masked_DS_A0_RBL_73';
path_data = 'Z:\shared\images\ProstateVGH-2\Data\Patient91\2019-11-6-10-3\BMode\ROI_Data\Down_Sample\';
%% 2020_07_25_14_01_49_asd on test#105
% path_pred = 'D:\Sharareh\Prostate_Project\Preparation\DataPreparation\predict\2020_07_25_14_01_49_asd_2020_07_27_09_50_48_test\test\';
path_pred='Z:\shared\images\ProstateVGH-2\Data\Dataset\InProstate\Ds\Best_models\2020_07_25_14_01_49_asd\';
path_mask = 'Z:\shared\images\ProstateVGH-2\Data\Patient91\2019-11-6-10-3\BMode\ROI_Data\Down_Sample\';
% path_save = 'D:\Sharareh\Prostate_Project\Postprocess\Teus_Unet\colormap\Colormaps_predictions\';
% if ~exist(path_save, 'dir')
%    mkdir(path_save)
% end
%load(strcat(path_data,data_name));
RF_dims=[2344,282];
if RF_dims(2)==282 %reach 280
           left_c=2;right_c=1;
       elseif size(rf_AllFrames,2)==392 %reach 280
           left_c=57;right_c=56;
       elseif size(rf_AllFrames,2)==372 %reach 280
           left_c=47;right_c=46;
       elseif size(rf_AllFrames,2)==414 %reach 280
           left_c=68;right_c=67;
       elseif size(rf_AllFrames,2)==438 %reach 420
           left_c=10;right_c=9; DS_L=3;
       elseif size(rf_AllFrames,2)==496 %reach 420
           left_c=39;right_c=38; DS_L=3;
       elseif size(rf_AllFrames,2)==466 %reach 420
           left_c=24;right_c=23; DS_L=3;
       elseif size(rf_AllFrames,2)==530 %reach 420
           left_c=56;right_c=55; DS_L=3;
%            left_c=150;right_c=99;
       elseif size(rf_AllFrames,2)==560 % reach 420
           left_c=71;right_c=70; DS_L=3;
end

%% load predicitons
load([path_pred,'_pred.mat']);  %prediction:(105, 256, 128, 2)
pred=squeeze(prediction(8,:,:,2));
%% load Masek RF_DS
mask_RF=imread([path_mask,'masked_DS_A0_RBL_73.jpg']);  %prediction:(105, 256, 128, 2)
figure
imagesc(mask_RF)
pred = pred .* single(mask_RF(:,1:128));
pred = [pred, zeros(256,12)];

J= imresize(pred,[cut_axial RF_dims(2)-mod(RF_dims(2),140)]);
J = [zeros(cut_axial,left_c-1),J,zeros(cut_axial,right_c)];
J= [J;zeros(RF_dims(1)-cut_axial,RF_dims(2))];
figure; imagesc(J)
%Convert to Bmode
    Patient=91;
    %RF_depth should be initialized from database
    RF_depth         = 0.06;     % RF/B-mode depth (end - start) (m) ?
    %%
    RF_data_sample   = squeeze(J);                    % Get one RF sample (one 3D volume)
    %RF_data_sample=img;
    RF_data_sample   = reshape(RF_data_sample,...
                        [size(RF_data_sample,1),...
                         size(RF_data_sample,2)*size(RF_data_sample,3),...
                         size(RF_data_sample,4)]);   % reshape to [nblocks, nlines, nplanes]- I have one plane, get read of 3rd dimension
Envelope = abs(hilbert(RF_data_sample)).^.3; %Envelop detection for B-mode                     
    RF_size          = size(RF_data_sample);                               % Get the side of RF grid 
    axial_res        = 1540/2*1/7.5/1000000;%measured from the machine :0.1/1000;%Obvious RF_depth/RF_size(1);%                               % 40 MHz axial Block size
    %bug              = -0.014;                                             % Probe radius offest A-plane (m) 
    RF_start         = 841*0.086e-3;                                            % RF/B-mode start depth (m) 
    RF_ends          = 226*0.086e-3;                                               % RF/B-mode end depth (m) 
    BMode_hight      = RF_ends - RF_start;%
    Probe_width      =0.0221;%Samira: 0.02975;%Data sheet:0.0182; %measure in the machine:0.020;                                             % Transducer width (m) 
    nBlocks          = RF_dims(1);%round(RF_depth/axial_res);%              % # of samples along a axial line
    nLines           = RF_dims(2);                                         % # of axial lines of the FOV
    Rt               = 0.01185; %0.01186;%0.01105;%Data sheet% Samira: 0.01105;%reference: 0.03367202788 + RF_start + bug;                     % FOV equivalent probe radius A-plane(m) 
    %bug_2            = 0.005;                                              % Probe radius offest B-plane (m) 
    Rm               = 0; %0.00500+RF_start+bug_2;                             % Probe radius B-plane (m) 
    arc              = asin(Probe_width/2/Rt)*2* Rt;% 0.027;%measure on the machine%       0.02975;%               % Probe arc A-plane 
    dR               = RF_depth/nBlocks;                                   
    dT               = arc/nLines/Rt;   
    dP               = 0;%1.91/180*pi;  %0                                      % B-plane angle increament 
    nPhi             = 1;%RF_size(3); %1
    dX               = 0.01098594255745411*0.01;%RF_depth/RF_size(1);%axial resolution measured in machine: 0.1e-3;%Dicom: 0.009e-3;%reference 0.6e-5;                                             % 3D voxel size
    sc_CPA           = ScanConversion(Rt,...
                                      Rm,...
                                      dR,...
                                      dT,...
                                      dP,...
                                      nBlocks,...
                                      nLines,...
                                      nPhi,...
                                      616,756,1,...
                                      dX, dX, dX);
    %crop test
    RF_data_sample   = J;      % Get one RF sample (one 3D volume)  
    B_3D      = sc_CPA.ScanConvert(Envelope);
    Raw_RF_3D = sc_CPA.ScanConvert(RF_data_sample);      
    J = imrotate(B_3D , 180 ); 
% replace image3d with imagesc
     figure;imagesc(RF_data_sample);colormap gray
%     % Scan converted rf
     figure;imagesc(Raw_RF_3D);colormap gray; axis equal;
%     % B-mode Envelope
     figure;imagesc(Envelope);colormap gray
%     % Scan converted B-mode
     figure;imagesc(J);colormap gray; axis equal;
img=imread('Z:\shared\images\ProstateVGH-2\Data\Patient91\2019-11-6-10-3\BMode\ROI_Data\TestImages\BMROI_test_0.bmp');

J(J<0.6*max(J(:)))=0;
figure;imagesc(J)
img(:,:,1)=J/max(J(:))*255;
figure;imagesc(img)
