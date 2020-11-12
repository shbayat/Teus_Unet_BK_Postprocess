%% p111-140 . high inv data

%%11
PID_test=[111, 113, 113, 111, 113, 113, 113, 113, 113, 113, 113, 123, 123,...
       123, 123, 123, 123, 123, 123, 112, 112, 112, 112, 112, 112, 112,...
       112, 112, 112, 112, 112, 114, 114, 114, 114, 114, 114, 114, 114,...
       114, 117, 117, 117, 117, 117, 117, 117, 117, 120, 120, 121, 121,...
       121, 121, 121, 121, 121, 121, 121, 121, 122, 124, 124, 124, 124,...
       124, 124, 124, 124, 124, 124, 124, 124, 125, 125, 125, 125, 125,...
       125, 125, 125, 125, 125, 125, 130, 132, 132, 127, 127, 127, 139,...
       139, 139, 139, 139, 139, 139, 139, 139, 139, 128, 128, 128, 128,...
       128, 128, 128, 130, 130, 130, 131, 131, 131, 131, 131, 131, 131,...
       131, 131, 131, 131, 132, 132, 132, 132, 132, 132, 134, 134, 134,...
       134, 134, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135,...
       135, 136, 136, 136, 136, 136, 136, 136, 136, 136];
   
   
idcore_test=[1103, 1121, 1129, 1101, 1122, 1123, 1125, 1126, 1127, 1128, 1130,...
  1217, 1219, 1222, 1223, 1224, 1218, 1220, 1225, 1109, 1110, 1111,...
  1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1131, 1132,...
  1133, 1134, 1136, 1137, 1138, 1139, 1140, 1162, 1164, 1166, 1167,...
  1169, 1170, 1171, 1172, 1195, 1196, 1197, 1198, 1200, 1201, 1203,...
  1204, 1205, 1206, 1207, 1208, 1209, 1227, 1228, 1229, 1230, 1231,...
  1232, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1241, 1242,...
  1244, 1245, 1246, 1247, 1248, 1249, 1250, 1276, 1299, 1300, 1259,...
  1262, 1265, 1367, 1368, 1369, 1370, 1371, 1372, 1373, 1374, 1375,...
  1376, 1267, 1268, 1269, 1270, 1271, 1272, 1274, 1275, 1277, 1279,...
  1284, 1285, 1286, 1287, 1288, 1289, 1290, 1291, 1292, 1293, 1294,...
  1295, 1296, 1297, 1298, 1301, 1302, 1313, 1314, 1315, 1316, 1318,...
  1321, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1330, 1331,...
  1332, 1334, 1335, 1336, 1337, 1338, 1339, 1340, 1341, 1342];


label=[1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,...
       1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,...
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,...
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,...
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,...
       1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0.,...
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,...
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,...
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.];
   
inv_test=[0.8 , 0.95, 0.95, 0.8 , 0.6 , 0.9 , 0.75, 0.9 , 0.95, 1.  , 1.  ,...
       0.95, 0.5 , 1.  , 0.95, 0.6 , 0.8 , 0.6 , 0.95, 0.  , 0.  , 0.  ,...
       0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,...
       0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,...
       0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,...
       0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,...
       0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,...
       0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.5 , 0.5 , 0.5 , 1.  ,...
       0.9 , 0.5 , 1.  , 1.  , 1.  , 1.  , 1.  , 1.  , 1.  , 1.  , 1.  ,...
       1.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,...
       0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,...
       0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,...
       0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,...
       0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ];   
   