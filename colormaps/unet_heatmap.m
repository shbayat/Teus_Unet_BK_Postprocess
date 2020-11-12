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
path_save = 'D:\Sharareh\Prostate_Project\Postprocess\Teus_Unet\colormap\Colormaps_predictions\';
if ~exist(path_save, 'dir')
   mkdir(path_save)
end
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
pred=pred.*single(mask_RF(:,1:128));

pred=[pred,zeros(256,12)];

J = imresize(pred,[cut_axial RF_dims(2)-mod(RF_dims(2),140)]);

J = [zeros(cut_axial,left_c-1),J,zeros(cut_axial,right_c)];
J= [J;zeros(RF_dims(1)-cut_axial,RF_dims(2))];
figure; imagesc(J)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert RF to Bmode
    Patient=91;
    
    RF_depth         = 0.06;                                 % RF/B-mode depth (end - start) (m) 

%     data_dir=strcat('z:\\shared\images\ProstateVGH-2\Data\Patient',num2str(Patient));
%     dirlist=dir(data_dir);
%     data_dir = strcat(data_dir,'\',dirlist(3).name,'\');
%     S = dir(strcat(data_dir,'IQ*.*'));
%     CoreNo=size(S,1);
%     oldfile= strcat(data_dir,S(i).name);
%     newIQfile= strrep(oldfile,'IQ_','RF_');
% 
% 
%     %Find RF Dimension
%     fid_iq = fopen(oldfile,'rt');
%     headerstr=fgetl(fid_iq);
%     fclose(fid_iq);
% 
%     iq =   regexp(headerstr,',');
%     dims(1) = str2num(headerstr(1:iq(1)));% 496;% %hight (h)
%     iq_height=str2num(headerstr(iq(1):iq(2)));
%     if iq_height>999
%         start_ch=12;
%     else
%         start_ch=11;
%     end
%     dims(2) = 2*iq_height;%884;% %width (w)
%     
% 
%     %read IQ data
%     fid_iq = fopen(oldfile);
%     iq =   fread(fid_iq,inf,'uint8');
%     fclose(fid_iq);
% 
% 
%     iq_noheader = iq(start_ch:end,:); %take out the header - the first 11 chars are the W and H of the image and size
%     [r1,r2]= size(iq_noheader);
%     l =length(iq_noheader);
%     dims(3) = l/(dims(1)*dims(2)*2); %Number of the frames (N)
%     count=1;
%     IQ_int16=zeros(l/2,1);
%     for j = 1:2:l
%         temp = uint8(iq_noheader(j:(j+1))');
%         IQ_int16(count) = typecast(temp , 'int16');%fliplr
%         count=count+1;
%     end

%     %IQ -> RF
%     fs=30000000; %RF frequency
%     %RF= zeros(r1,r2);
%     RF= zeros(l/2,r2);
%     %RF(:,r2) = iq2rf(double(iq_noheader(:,r2)), 2, fs/4, fs); %convert IQ to RF
%     RF(:,r2) = iq2rf(double(IQ_int16), 2, fs/4, fs); %convert IQ to RF
%     % reshape into a matrix with shape {w,h,N}
%     frames = reshape(RF,[dims(2),dims(1),dims(3)]);

    %%
    RF_data_sample   = squeeze(J);                    % Get one RF sample (one 3D volume)
    %RF_data_sample=img;
    RF_data_sample   = reshape(RF_data_sample,...
                        [size(RF_data_sample,1),...
                         size(RF_data_sample,2)*size(RF_data_sample,3),...
                         size(RF_data_sample,4)]);                         % reshape to [nblocks, nlines, nplanes]- I have one plane, get read of 3rd dimension

    
    Envelope         = abs(hilbert(RF_data_sample)).^.3;                   % Envelop detection for B-mode
    RF_size          = size(RF_data_sample);                               % Get the side of RF grid 
    axial_res        = 1540/2*1/7.5/1000000;%measured from the machine :0.1/1000;%Obvious RF_depth/RF_size(1);%                               % 40 MHz axial Block size
    %bug              = -0.014;                                             % Probe radius offest A-plane (m) 
    RF_start         = 841*0.086e-3;                                            % RF/B-mode start depth (m) 
    RF_ends          = 226*0.086e-3;                                               % RF/B-mode end depth (m) 
    BMode_hight      = RF_ends - RF_start;%
    Probe_width      =0.0221;%Samira: 0.02975;%Data sheet:0.0182; %measure in the machine:0.020;                                             % Transducer width (m) 
    nBlocks          = RF_dims(1);%round(RF_depth/axial_res);%              % # of samples along a axial line
    nLines           = RF_dims(2);                                         % # of axial lines of the FOV
    Rt               = 0.01186;%0.01105;%Data sheet% Samira: 0.01105;%reference: 0.03367202788 + RF_start + bug;                     % FOV equivalent probe radius A-plane(m) 
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
    nBlocks=nBlocks-63;
    for i=1:1
        RF_data_sample   = J;%squeeze(J(1:nBlocks,:,i));                    % Get one RF sample (one 3D volume)
        
        Envelope         = abs(hilbert(RF_data_sample)).^.3;                   % Envelop detection for B-mode
        B_3D      = sc_CPA.ScanConvert(Envelope);

        Raw_RF_3D = sc_CPA.ScanConvert(RF_data_sample);

         J = flipdim(B_3D , 1 );
     end    
    %% check the 3D images
    % Raw rf
%      figure;imshow3D(RF_data_sample);colormap gray
% %     % Scan converted rf
%      figure;imshow3D(Raw_RF_3D);colormap gray; axis equal;
% %     % B-mode Envelope
%      figure;imshow3D(Envelope);colormap gray
% %     % Scan converted B-mode
%      figure;imshow3D(J);colormap gray; axis equal; 
% img=imread("\\smbhome\rcl\shared\images\ProstateVGH-2\Data\Patient91\2019-11-6-10-3\BMode\ROI_Data\TestImages\RFROI_test_8.bmp");
% img(:,:,3)=RF_data_sample*255;
% imagesc(img)

%% replace image3d with imagesc
     figure;imagesc(RF_data_sample);colormap gray
%     % Scan converted rf
     figure;imagesc(Raw_RF_3D);colormap gray; axis equal;
%     % B-mode Envelope
     figure;imagesc(Envelope);colormap gray
%     % Scan converted B-mode
     figure;imagesc(J);colormap gray; axis equal; 
img=imread('Z:\shared\images\ProstateVGH-2\Data\Patient91\2019-11-6-10-3\BMode\ROI_Data\TestImages\BMROI_test_0.bmp');
img(1:615,:,3)=J/max(J(:))*255;
figure;imagesc(img)