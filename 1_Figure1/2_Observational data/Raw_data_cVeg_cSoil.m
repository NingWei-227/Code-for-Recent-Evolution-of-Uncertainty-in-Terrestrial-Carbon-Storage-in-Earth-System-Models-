% Benchmark analysis on cLand
% In models: cLand = cVeg+cLitter+cSoil
% Observation-derived data: 
% cVeg(aboveground and belowground): 
% Spawn et al., 2020; unit: Mg C/ha;
% NO data value: 65536;
% scaling: 0.1
% resolution: 0.0028 (~300m)           
%
% cSoil(assume including cLitter) data:
% HWSD: unit: Kg C m-2 
% resolution: 0.05
% variables: AWT_T_SOC(Area weighted topsoil carbon content) 0-30 cm
%            AWT_S_SOC(Area weighted subsoil carbon content) 30-100 cm
%
% LandGIS: unit: Kg m-2  (Kg/m2 * 10 = t/ha)
% https://zenodo.org/record/2536040#.Xy0UhSgzZjW
% resolution: 250m
% depth interval: 0每10, 10每30, 30每60, 60每100 and 100每200 cm 
% (0-100cm cSoil are used as in the benchmark analysis)
%  
% SoilGrids: unit: t/ha
% https://www.isric.org/explore/soilgrids
% resolution: 250m
% depth interval: 0每5, 5每15, 15每30, 30每60,60-100 and 100每200 cm 
% (0-100cm cSoil are used as in the benchmark analysis)
%
% NCSCD: unit: Kg/m2 (Northern circumpolar cSoil Database)
% https://bolin.su.se/data/ncscd/more.php
% resolution: 0.5 degree
% depth interval: 0每30 cm, 0每100 cm, 100每200 cm and 200每300cm depth
clear
clc
%% Observation 
% All dataset were regrided into 0.05x0.05 for consistency
% covert unit to Kg C m-2 if the unit was not 
% cVeg,
% lon -180:180; lat: 84:-61.1
% scaling: 0.1()
[cVeg_ag, cmap2]  =  imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\Biomass\ORNL\Global_Maps_C_Density_2010_1763\Global_Maps_C_Density_2010_1763\data\aboveground_biomass_carbon_2010.tif');
[cVeg_bg  cmap2] =  imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\Biomass\ORNL\Global_Maps_C_Density_2010_1763\Global_Maps_C_Density_2010_1763\data\belowground_biomass_carbon_2010.tif');
cVeg_ag(cVeg_ag==65536) =  NaN;
cVeg_bg(cVeg_bg==65536) =  NaN;
cVeg_obs = cVeg_ag + cVeg_bg;
cVeg_obs = cVeg_obs.*0.1;

cVeg_ag_un = imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\Biomass\ORNL\Global_Maps_C_Density_2010_1763\Global_Maps_C_Density_2010_1763\data\aboveground_biomass_carbon_2010_uncertainty.tif');
cVeg_bg_un = imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\Biomass\ORNL\Global_Maps_C_Density_2010_1763\Global_Maps_C_Density_2010_1763\data\belowground_biomass_carbon_2010_uncertainty.tif');
cVeg_ag_un(cVeg_ag_un==65536) =  NaN;
cVeg_bg_un(cVeg_bg_un==65536) =  NaN;
cVeg_obs_un = cVeg_ag_un + cVeg_bg_un;
cVeg_obs_un = cVeg_obs_un.*0.1;
clearvars -except cVeg_obs cVeg_obs_un

lat_interval = (84+61.1)./52201;
lon_interval = 360./129600;
nor = 6./lat_interval; 
nor_matrix(1:2159,1:129600) = NaN;
sou = (90-61.1)./lat_interval;
sou_matrix(1:10397,1:129600) = NaN;

cVeg_obs2 = [nor_matrix;cVeg_obs;sou_matrix];
cVeg_obs_un2 = [nor_matrix;cVeg_obs_un;sou_matrix];

cVeg_obs05 = imresize(cVeg_obs2,[360,720]);
cVeg_obs_un05 = imresize(cVeg_obs_un2,[360,720]);
cVeg_obs05 = double(cVeg_obs05);
cVeg_obs_un05 = double(cVeg_obs_un05);

% convert unit from Mg C/ha to Kg C/m2
cVeg_obs05kg = cVeg_obs05.*1000./10000;
cVeg_obs_un05kg = cVeg_obs_un05.*1000./10000;

cVeg_obs05kg(cVeg_obs05kg <= 0) = NaN;
cVeg_obs_un05kg(cVeg_obs_un05kg<=0) = NaN;
clearvars -except cVeg_obs05kg cVeg_obs_un05kg
save('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\Biomass\ORNL\Global_Maps_C_Density_2010_1763\cVeg05_KgCm2.mat')
%% Processing cSoil Data
clear
clc
% HWSD: unit: Kg C m-2 
% resolution: 0.05
% variables: AWT_T_SOC(Area weighted topsoil carbon content) 0-30 cm
%            AWT_S_SOC(Area weighted subsoil carbon content) 30-100 cm
HWSD_T = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\HWSD\HWSD_1247\HWSD_1247\data\AWT_T_SOC.nc4','SUM_t_c_12');
HWSD_S = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\HWSD\HWSD_1247\HWSD_1247\data\AWT_S_SOC.nc4','SUM_s_c_1');
HWSD_ST = HWSD_T + HWSD_S;
%convert to match with the global map
HWSD_ST = HWSD_ST';
HWSD_1m(1:3600,1:7200) = NaN;
for i=1:3600
    HWSD_1m(i,:) = HWSD_ST(3601-i,:); 
end

% regrid into 0.5x0.5
lat05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lat');
lon05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lon');
lat005 = -ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\HWSD\HWSD_1247\HWSD_1247\data\AWT_T_SOC.nc4','lat');
lon005 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\HWSD\HWSD_1247\HWSD_1247\data\AWT_T_SOC.nc4','lon');
[x05,y05] = meshgrid(lon05,lat05);
[x005,y005] = meshgrid(lon005,lat005);
HWSD1m_05 = interp2(x005,y005,HWSD_1m,x05,y05,'linear');
clearvars -except HWSD1m_05
save('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\HWSD\HWSD_1247\HWSDv2_05.mat')

% LandGIS: unit: Kg m-2  (Kg/m2 * 10 = t/ha)
% https://zenodo.org/record/2536040#.Xy0UhSgzZjW
% resolution: 250m
% depth interval: 0每10, 10每30, 30每60, 60每100 and 100每200 cm 
% (0-100cm cSoil are used as in the benchmark analysis)
% No data: -32768
% lon: -180:180
% lat: 87.4:-62
clear
clc
LandGIS_01 = imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\LandGIS\ArcGIS\LandGIS_05_0_10cm.tif'); 
LandGIS_03 = imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\LandGIS\ArcGIS\LandGIS_05_10_30cm.tif');
LandGIS_06 = imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\LandGIS\ArcGIS\LandGIS_05_30_60cm.tif');
LandGIS_1 = imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\LandGIS\ArcGIS\LandGIS_05_60_100cm.tif');

LandGIS_01 = double(LandGIS_01);
LandGIS_03 = double(LandGIS_03);
LandGIS_06 = double(LandGIS_06);
LandGIS_1 = double(LandGIS_1);

LandGIS_01(LandGIS_01 == -32768) = NaN;
LandGIS_03(LandGIS_03 == -32768) = NaN;
LandGIS_06(LandGIS_06 == -32768) = NaN;
LandGIS_1(LandGIS_1 == -32768) = NaN;

LandGIS_1m = LandGIS_01+LandGIS_03+LandGIS_06+LandGIS_1;

lat_interval = 0.5;
nor = (90-87.4)/lat_interval; 
nor_matrix(1:5,1:720) = NaN;
sou = (90-62)./lat_interval;
sou_matrix(1:sou,1:720) = NaN;
LandGIS_1m = [nor_matrix;LandGIS_1m;sou_matrix];
length(find(LandGIS_1m(LandGIS_1m>100)))

LandGIS_1mKg = LandGIS_1m;
clearvars -except LandGIS_1mKg
save('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\LandGIS\LandGIS05_1m.mat')

% SoilGrids: unit: t/ha (t/ha.*0.1 = Kg/m2)
% https://www.isric.org/explore/soilgrids
% resolution: 250m
% depth interval: 0每5, 5每15, 15每30, 30每60,60-100 and 100每200 cm 
% (0-100cm cSoil are used as in the benchmark analysis)
% No data: -32768
% lon: -180:180
% lat: 84:-56
clear
clc
SoilGrid_5cm = imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\SoilGrids\SoilGrid_ArcGIS\SoilGrid05_0_5cm.tif');
SoilGrid_15cm = imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\SoilGrids\SoilGrid_ArcGIS\SoilGrid05_5_15cm.tif');
SoilGrid_30cm = imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\SoilGrids\SoilGrid_ArcGIS\SoilGrid05_15_30cm.tif');
SoilGrid_60cm = imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\SoilGrids\SoilGrid_ArcGIS\SoilGrid05_30_60cm.tif');
SoilGrid_100cm = imread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\SoilGrids\SoilGrid_ArcGIS\SoilGrid05_60_100cm.tif');

SoilGrid_5cm = double(SoilGrid_5cm); SoilGrid_5cm(SoilGrid_5cm==-32768) = NaN;
SoilGrid_15cm = double(SoilGrid_15cm);SoilGrid_15cm(SoilGrid_15cm==-32768) = NaN;
SoilGrid_30cm = double(SoilGrid_30cm);SoilGrid_30cm(SoilGrid_30cm==-32768) = NaN;
SoilGrid_60cm = double(SoilGrid_60cm);SoilGrid_60cm(SoilGrid_60cm==-32768) = NaN;
SoilGrid_100cm = double(SoilGrid_100cm);SoilGrid_100cm(SoilGrid_100cm==-32768) = NaN;

SoilGrid_1m = SoilGrid_5cm + SoilGrid_15cm + SoilGrid_30cm + SoilGrid_60cm + SoilGrid_100cm;


lat_interval = 0.5;
nor = (90-84)/lat_interval; 
nor_matrix(1:nor,1:720) = NaN;
sou = (90-56)./lat_interval;
sou_matrix(1:sou,1:720) = NaN;

SoilGrid_1m = [nor_matrix;SoilGrid_1m;sou_matrix];
SoilGrid_1mkg = SoilGrid_1m.*0.1; % from t/ha to kg/m2
clearvars -except SoilGrid_1mkg
save('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\SoilGrids\SoilGrid05_1mkg.mat')

% NCSCD: unit: hg/m2.*0.1 = Kg/m2  (Northern circumpolar cSoil Database)
% https://bolin.su.se/data/ncscd/more.php
% resolution: 0.5 degree
% depth interval: 0每30 cm, 0每100 cm, 100每200 cm and 200每300cm depth
% (0-100cm cSoil are used as in the benchmark analysis)
clear
clc
NCSCD_1m = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\Circumpolar_cSoil\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_netCDF_05deg\NCSCDv2_Circumpolar_WGS84_SOCC100_05deg.nc','NCSCDv2');
NCSCD_1m = double(NCSCD_1m);
NCSCD_1m(NCSCD_1m==-32768) = NaN;

Tglobal(1:720,1:249) = NaN;
NCSCDgb_1m = cat(2,NCSCD_1m,Tglobal); 

polar_T = rot90(NCSCDgb_1m);
polar_map(1:360,1:720) = NaN;
for k=1:360
    polar_map(361-k,:) = polar_T(k,:);
end
polar_map(polar_map<=0) = NaN;

NCSCDgb_1mKg = polar_map.*0.1; % from hg/m2 to Kg/m2
clearvars -except NCSCDgb_1mKg
save('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\Circumpolar_cSoil\NCSCDv2_Circumpolar_netCDF_05deg\NCSCD05_1mkg.mat')
