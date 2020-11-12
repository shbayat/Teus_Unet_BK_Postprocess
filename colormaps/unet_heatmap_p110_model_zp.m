%% 
% % 
% points used in axial:  1:1536 averaged every 6 in xial direction
%

clear; clc;
close all;

PID_test=[100, 100, 100, 100, 100, 100, 100,  91,  91,  91,  91,  91, 100,...
        92,  92,  92,  92,  92,  92,  92,  92,  93,  93,  93,  93,  93,...
        95,  95,  95,  95,  95,  95,  95,  95,  95,  96,  96,  96,  96,...
        96,  96,  96,  98,  98,  98,  98,  98,  98,  98,  98,  98,  98,...
        99,  99,  99,  99,  99,  99,  99,  99,  99, 101, 110, 110, 110,...
       106, 101, 101, 103, 103, 103, 103, 103, 103, 103, 103, 104, 104,...
       105, 105, 105, 105, 105, 105, 105, 108, 108, 108, 108, 108, 108,...
       108, 108, 108, 108, 109, 109, 109, 109, 109, 110, 110, 110, 110,...
       110];
   
%11
idcore_test=[997,  999, 1000, 1001, 1002, 1003, 1005,  907,  906,  908,  909,...
            910,  998,  911,  912,  913,  914,  916,  917,  918,  919,  922,...
            927,  928,  929,  930,  941,  942,  943,  945,  946,  947,  948,...
            949,  950,  952,  953,  954,  956,  957,  959,  961,  973,  975,...
             976,  978,  979,  980,  981,  982,  983,  984,  986,  988,  989,...
            990,  991,  992,  993,  995,  996, 1011, 1089, 1091, 1092, 1054,...
            1008, 1009, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1037,...
            1039, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1071, 1072, 1073,...
            1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084,...
            1087, 1090, 1094, 1095, 1096, 1098];
        
% ind=find(idcore_test==1096); %  p110, core7 
ind=find(idcore_test==1096); %  p110, core7         
cut_axial=1536;
% downsample_rate=6;
Patient=110;
CoreID=7; %CID:1091, core: 2
data_dir=strcat('z:\\shared\images\ProstateVGH-2\Data\Patient',num2str(Patient));
dirlist=dir(data_dir);
data_dir = strcat(data_dir,'\',dirlist(3).name,'\');

% data_name = 'RF_DS_res_A7_LML_208'; %P110, core7
% data_name_mask='masked_DS_res_A7_LML_208';
% path_data = strcat(data_dir,'\BMode\ROI_Data\Cut_axial\');
%% on test#105

path_pred='Z:\shared\images\ProstateVGH-2\Data\Dataset\InProstate\Zero-pad\Best_Models\Not_Normalized\2020_10_05_09_26_59_asd\'; %zeropad model
path_mask = strcat(data_dir,'\BMode\ROI_Data\Cut_axial\');
S = dir(strcat(path_mask,'masked1_cut_a_A*.mat'));
% path_save = 'D:\Sharareh\Prostate_Project\Postprocess\Teus_Unet\colormap\Colormaps_predictions\';
% if ~exist(path_save, 'dir')
%    mkdir(path_save)
% end
%load(strcat(path_data,data_name));
% RF_dims=[2344,282];

RF_dims=[2344,414]; %P110 core7-Benign
  

%% load predicitons
load([path_pred,'_pred.mat']);  %prediction:(105, 256, 528, 2)
img=imread([data_dir, '\BMode\ROI_Data\TestImages\BMROI_test_', int2str(CoreID),'.bmp']);

%% load mask Bmode
SS = dir([data_dir, '\BMode\ROI_Data\WholeProstate\','BM_mask_A*.bmp']);
SS(CoreID).name
if CoreID==0
    BM_mask=imread([data_dir, '\BMode\ROI_Data\WholeProstate\', SS(CoreID+1).name]);
else
    BM_mask=imread([data_dir, '\BMode\ROI_Data\WholeProstate\', SS(CoreID).name]);
end
figure(1);subplot(1,2,1); imagesc(img);title("P110-C7")
subplot(1,2,2); imagesc(BM_mask);title("mask-Bmode-P110-C7")
%%
for ii=1:2
    pred=squeeze(prediction(ind,:,:,ii));
    %% load Masek RF
%     mask_RF=imread([path_mask,S(CoreID+1).name]);  %prediction:(105, 256, 528, 2)
    mask_RF=load([path_mask,S(CoreID+1).name]);  %prediction:(105, 256, 528, 2)
    masked_RF=mask_RF.mask_RF_cut;
    figure;
    imagesc(masked_RF);colormap gray
    %% averaging every 6 samples-axial
%     k=1;
%     for i = 6:6:size(mask_RF.masked_RF,1) %1536
%     masked_RF2(k, :) = mean(masked_RF(i-5:i, :), 1); %Mean along 1st dimension
%     k = k+1;
%     end
    
    %% Zero padding
    Init=zeros(1536,528);
    Init(:,1:size(masked_RF,2))=masked_RF;
    masked_RF_zp=Init;
    
    J= imresize(pred,[cut_axial 528]);
    J = J .* single(masked_RF_zp(:,:))/max(masked_RF_zp(:));
%     pred = [pred, zeros(256,12)];

    J=J(:,1:RF_dims(2));
    
    J= [J;zeros(RF_dims(1)-cut_axial,RF_dims(2))];
    
    figure; imagesc(J);title("resized pred")
    %Convert to Bmode
     
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
         figure;imagesc(Envelope);colormap gray;title("B-mode Envelope")
         
    %     % Scan converted B-mode
         figure;imagesc(J);colormap gray; axis equal; title("Scan converted B-mode")

    

    J(J<0.5*max(J(:)))=0;
    figure;imagesc(J);title("J")
    if ii==1
        img(:,:,3)=double(img(:,:,3))+floor(J/max(J(:))*125).*BM_mask;
    else
        img(:,:,1)=double(img(:,:,1))+floor(J/max(J(:))*125).*BM_mask;
    end
    
end
    figure;imagesc(img)
    
% %% p111-140 . high inv data
% 
% %%11
% idcore_test111_140=[1103, 1121, 1129, 1101, 1122, 1123, 1125, 1126, 1127, 1128, 1130,...
%   1217, 1219, 1222, 1223, 1224, 1218, 1220, 1225, 1109, 1110, 1111,...
%   1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1131, 1132,...
%   1133, 1134, 1136, 1137, 1138, 1139, 1140, 1162, 1164, 1166, 1167,...
%   1169, 1170, 1171, 1172, 1195, 1196, 1197, 1198, 1200, 1201, 1203,...
%   1204, 1205, 1206, 1207, 1208, 1209, 1227, 1228, 1229, 1230, 1231,...
%   1232, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1241, 1242,...
%   1244, 1245, 1246, 1247, 1248, 1249, 1250, 1276, 1299, 1300, 1259,...
%   1262, 1265, 1367, 1368, 1369, 1370, 1371, 1372, 1373, 1374, 1375,...
%   1376, 1267, 1268, 1269, 1270, 1271, 1272, 1274, 1275, 1277, 1279,...
%   1284, 1285, 1286, 1287, 1288, 1289, 1290, 1291, 1292, 1293, 1294,...
%   1295, 1296, 1297, 1298, 1301, 1302, 1313, 1314, 1315, 1316, 1318,...
%   1321, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331,...
%   1332, 1334, 1335, 1336, 1337, 1338, 1339, 1340, 1341, 1342];    