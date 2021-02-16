% benchmark analysis for ecosystem residence time (TuaE)
% Data source of TuaE was obtained from Fan et al., 2020. TuaE = (Ecosystem C storage)/(GPP)
% For consistency, TuaE was calculated in the same way for CMIP56 models in the benchmark analysis
clear;clc;

%% CMIP5
% load land C storage simulated from CMIP5 models
% unit: kgC m-2
cd('G:\CMIP5\4_MatData\spatial_data\hist_rcp85')
load('cLand_cmip5_251.mat')
% calculate 2001-2005 mean
sp5_Cend5(:,:,:,1) = cLand_BCC(:,:,152:156);
sp5_Cend5(:,:,:,2) = cLand_CAN(:,:,152:156);
sp5_Cend5(:,:,:,3) = cLand_CCSM(:,:,152:156);
sp5_Cend5(:,:,:,4) = cLand_HAD(:,:,142:146);
sp5_Cend5(:,:,:,5) = cLand_IPSL(:,:,152:156);
sp5_Cend5(:,:,:,6) = cLand_MIROC(:,:,152:156);
sp5_Cend5(:,:,:,7) = cLand_MPI(:,:,152:156);
sp5_Cend5(:,:,:,8) = cLand_NOR(:,:,152:156);
sp5_Cend5(:,:,:,9) = cLand_BNU(:,:,152:156);
sp5_Cend5(:,:,:,10) = cLand_GF(:,:,141:145);
sp5_Cend5(:,:,:,11) = cLand_MRI(:,:,151:155);
sp5_Cend5_11 = nanmean(sp5_Cend5,3);
sp5_Cend5_11 = squeeze(sp5_Cend5_11);
sp5_Cend5_11(:,1:30,:) = NaN; 
% unify the terrestrial region
mask = sp5_Cend5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp5_Cend5_11 = sp5_Cend5_11.*mask;

% load GPP data and calculate 2001-2005 mean
% unit: KgC m-2 yr-1
load('GPP_cmip5_251.mat')
sp5_GPP5(:,:,:,1) = gppBCC_3D(:,:,152:156);
sp5_GPP5(:,:,:,2) = gppCAN_3D(:,:,152:156);
sp5_GPP5(:,:,:,3) = gppCCSM_3D(:,:,152:156);
sp5_GPP5(:,:,:,4) = gppHAD_3D(:,:,142:146);
sp5_GPP5(:,:,:,5) = gppIPSL_3D(:,:,152:156);
sp5_GPP5(:,:,:,6) = gppMIROC_3D(:,:,152:156);
sp5_GPP5(:,:,:,7) = gppMPI_3D(:,:,152:156);
sp5_GPP5(:,:,:,8) = gppNOR_3D(:,:,152:156);
sp5_GPP5(:,:,:,9) = gppBNU_3D(:,:,152:156);
sp5_GPP5(:,:,:,10) = gppGF_3D(:,:,141:145);
sp5_GPP5(:,:,:,11) = gppMRI_3D(:,:,151:155);
sp5_GPP5_11 = nanmean(sp5_GPP5,3);
sp5_GPP5_11 = squeeze(sp5_GPP5_11);
sp5_GPP5_11(:,1:30,:) = NaN; 
% unify the terrestrial region based on simulations from BCC
mask = sp5_GPP5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp5_GPP5_11 = sp5_GPP5_11.*mask;
Models5 = {'BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'};
% omit regions where GPP<= 0.01 KgC m-2 yr-1        
sp5_GPP5_11(sp5_GPP5_11<=0.01)=NaN;
clearvars -except sp5_Cend5_11 sp5_GPP5_11 Models5
% calculate tuaE: tuaE = Cland/GPP, unit: years
numM5 = 11;
for yr = 1:numM5
    
    tuaE_CMIP5(:,:,yr) = sp5_Cend5_11(:,:,yr)./sp5_GPP5_11(:,:,yr);
      
end
% get ride of useless variables to save memory
clearvars -except sp5_Cend5_11 sp5_GPP5_11 Models5 tuaE_CMIP5

%% CMIP6
% In CMIP6, CESM2 and NorESM2-LM simulated soil carbon storage along soil depth
% In the benchmark analysis,we only included soil carbon above 1m
% unit: KgC m-2 
cSoil_CESM1m = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6\4_1gre_cSoilAbove1m_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cSoilAbove1m');
cSoil_NOR1m = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6\10_1gre_cSoilAbove1m_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cSoilAbove1m');
cVeg6_CESM = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6\4_1gre_cVeg_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');
cVeg6_NOR = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6\10_1gre_cVeg_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');
cLit6_CESM = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6\4_1gre_cLitter_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');
cLit6_NOR = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3\CMIP6\10_1gre_cLitter_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');

cLand_CESM1m = cVeg6_CESM + cLit6_CESM + cSoil_CESM1m;
cLand_NOR1m = cVeg6_NOR + cLit6_NOR + cSoil_NOR1m;

cd('H:\CMIP56_Csink\4_MatData\spatial_data\hist_ssp585')
load('cLand_cmip6_251.mat')  
% cLand simulated from CMIP6 models
% cLand from 2001-2005
sp6_Cend6(:,:,:,1) = cLand_BCC(:,:,152:156);
sp6_Cend6(:,:,:,2) = cLand_CAN(:,:,152:156);
sp6_Cend6(:,:,:,3) = cLand_CESM1m(:,:,152:156);
sp6_Cend6(:,:,:,4) = cLand_UK(:,:,152:156);
sp6_Cend6(:,:,:,5) = cLand_IPSL(:,:,152:156);
sp6_Cend6(:,:,:,6) = cLand_MIC(:,:,152:156);
sp6_Cend6(:,:,:,7) = cLand_MPI(:,:,152:156);
sp6_Cend6(:,:,:,8) = cLand_NOR1m(:,:,152:156);
sp6_Cend6(:,:,:,9) = cLand_ASS(:,:,152:156);
sp6_Cend6(:,:,:,10) = cLand_CNRM(:,:,152:156);
sp6_Cend6(:,:,:,11) = cLand_EC(:,:,152:156);
% 2001-2005 mean cLand
sp6_Cend5_11 = nanmean(sp6_Cend6,3);
sp6_Cend5_11 = squeeze(sp6_Cend5_11);
sp6_Cend5_11(:,1:30) = NaN;
mask = sp6_Cend5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp6_Cend5_11 = sp6_Cend5_11.*mask;

% load GPP data
% unit: KgC m-2 yr-1
cd('H:\CMIP56_Csink\4_MatData\spatial_data\hist_ssp585')
load('GPP_cmip6_251.mat')             
sp6_GPP5(:,:,:,1) = gppBCC_3D(:,:,152:156);
sp6_GPP5(:,:,:,2) = gppCAN_3D(:,:,152:156);
sp6_GPP5(:,:,:,3) = gppCESM_3D(:,:,152:156);
sp6_GPP5(:,:,:,4) = gppUK_3D(:,:,152:156);
sp6_GPP5(:,:,:,5) = gppIPSL_3D(:,:,152:156);
sp6_GPP5(:,:,:,6) = gppMIC_3D(:,:,152:156);
sp6_GPP5(:,:,:,7) = gppMPI_3D(:,:,152:156);
sp6_GPP5(:,:,:,8) = gppNOR_3D(:,:,152:156);
sp6_GPP5(:,:,:,9) = gppASS_3D(:,:,152:156);
sp6_GPP5(:,:,:,10) = gppCNRM_3D(:,:,152:156);
sp6_GPP5(:,:,:,11) = gppEC_3D(:,:,152:156);
sp6_GPP5_11 = nanmean(sp6_GPP5,3);
sp6_GPP5_11 = squeeze(sp6_GPP5_11);
sp6_GPP5_11(:,1:30) = NaN;
% unify the terrestrial region based on simulations from BCC
mask = sp6_GPP5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp6_GPP5_11 = sp6_GPP5_11.*mask;

Models6 = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'};
    
clearvars -except sp5_Cend5_11 sp5_GPP5_11 Models5 tuaE_CMIP5 ...
                  sp6_Cend5_11 sp6_GPP5_11 Models6 
% omit regions where GPP<=0.01 KgC m-2 yr-1             
sp6_GPP5_11(sp6_GPP5_11<=0.01) = NaN;              
num_6M = 11;
tuaE_CMIP6(1:360,1:180,1:11) = NaN;
for yr=1:num_6M 
    tuaE_CMIP6(:,:,yr) = sp6_Cend5_11(:,:,yr)./sp6_GPP5_11(:,:,yr);
    
end
% get ride of useless variables to save memory
clearvars -except Models5 tuaE_CMIP5 ...
                  Models6 tuaE_CMIP6 

%%
% convert tuaE_CMIP5 tuaE_CMIP6 into the global map               
tuaE5_map11(1:180,1:360,1:11) = NaN;
tuaE6_map11(1:180,1:360,1:11) = NaN;
for yr =1:11
    yr
    
    tuaE5_M = tuaE_CMIP5(:,:,yr);
    tuaE6_M = tuaE_CMIP6(:,:,yr);
    
    tuaE5_t = tuaE5_M';
    tuaE6_t = tuaE6_M';
    
    tuaE5_LR(1:180,1:360) = NaN;
    tuaE6_LR(1:180,1:360) = NaN;
    
    tuaE5_LR(:,1:180) = tuaE5_t(:,181:360);
    tuaE5_LR(:,181:360) = tuaE5_t(:,1:180);
    tuaE6_LR(:,1:180) = tuaE6_t(:,181:360);
    tuaE6_LR(:,181:360) = tuaE6_t(:,1:180);
    
    map_tuaE5(1:180,1:360) = NaN;
    map_tuaE6(1:180,1:360) = NaN;
    for k=1:180
        map_tuaE5(181-k,:) = tuaE5_LR(k,:);
        map_tuaE6(181-k,:) = tuaE6_LR(k,:);
    end
    
    tuaE5_map11(:,:,yr) = map_tuaE5;
    tuaE6_map11(:,:,yr) = map_tuaE6;
      
end   
% regrid model simulations to match with data 
lat05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lat');
lon05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lon');
lat1 = 89.5:-1:-89.5; lat1 = lat1';
lon1 = -179.5:179.5;  lon1 = lon1';
[x05,y05] = meshgrid(lon05,lat05);
[x1,y1] = meshgrid(lon1,lat1);
% regrided CMIP5 and CMIP6 data into 0.5x0.5 resolution
tuaE5_sp05(1:360,1:720,1:11) = NaN;
tuaE6_sp05(1:360,1:720,1:11) = NaN;
for yr=1:11
    tuaE5_M = interp2(x1,y1,tuaE5_map11(:,:,yr),x05,y05,'linear');
    tuaE6_M = interp2(x1,y1,tuaE6_map11(:,:,yr),x05,y05,'linear');
    
    tuaE5_sp05(:,:,yr) = tuaE5_M;
    tuaE6_sp05(:,:,yr) = tuaE6_M;
end

% load observations 
tuaE_median_obs = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','tau');
tuaE_10th_obs = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','tau_q10');
tuaE_90th_obs = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','tau_q90');
% using tuaE data that only considering 0-1m soil C
tuaE_median1m_obs = tuaE_median_obs(:,:,1)';  
tuaE_10th1m_obs = tuaE_10th_obs(:,:,1)';
tuaE_90th1m_obs = tuaE_90th_obs(:,:,1)';

% match model simulations with observational dataset
tuaE_median1m_obs(~isnan(tuaE_median1m_obs)) = 1;
Mask_obs = tuaE_median1m_obs;
tuaE5_sp05 =  tuaE5_sp05.* Mask_obs;
tuaE6_sp05 =  tuaE6_sp05.* Mask_obs;
tuaE_median1m_obs = tuaE_median_obs(:,:,1)';
% get ride of useless variables to save memory
clearvars -except Models5 tuaE5_sp05...
                  Models6 tuaE6_sp05...
                  tuaE_median1m_obs tuaE_10th1m_obs tuaE_90th1m_obs
% save data used to produce Fig.2                
save E:\1_Mycase\3_CMIP56_Cland\4_Version3\1_Figure\2_Figure2_X_GPP_tauE\Model_obsData\cmip56_Obs_tuaE.mat           
              
% calculate across-model mean and standard deviation
tuaE5_sp05(tuaE5_sp05==Inf) = NaN; tuaE5_sp05(tuaE5_sp05<= 0) = NaN; tuaE5_sp05(tuaE5_sp05>10^20) = NaN;tuaE5_sp05(tuaE5_sp05==0) = NaN;
tuaE6_sp05(tuaE6_sp05==Inf) = NaN; tuaE6_sp05(tuaE6_sp05<= 0) = NaN; tuaE6_sp05(tuaE6_sp05>10^20) = NaN;tuaE6_sp05(tuaE6_sp05==0) = NaN;
tuaE_avg11_cmip5 = nanmean(tuaE5_sp05,3);
tuaE_avg11_cmip6 = nanmean(tuaE6_sp05,3);
tuaE_sd11_cmip5 = nanstd(tuaE5_sp05,0,3);
tuaE_sd11_cmip6 = nanstd(tuaE6_sp05,0,3); 
% Model-data difference in tuaE
tuaE5_bias = tuaE_avg11_cmip5 - tuaE_median1m_obs;
tuaE6_bias = tuaE_avg11_cmip6 - tuaE_median1m_obs;

%% Extended figure 4
%ax = gca;
%mycolor_tuaE_bias = colormap(ax)
%save ('E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\mycolor_tuaE_bias.mat','mycolor_tuaE_bias')

% load color code for color bar
load E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\mycolor_tuaE_bias.mat
figure
set(gcf,'position',[100 100 499.2,630])
maps = tight_subplot(2,1,[-0.1 -0.085],[0.35 -0.05],[-0.04 0.26])

% CMIP5 tuaE bias
tuaE5_Bia_M11 = tuaE5_bias;
tuaE5_Bia_M11(302:360,:) = [];
tuaE5_Bia_M11 = flipud(tuaE5_Bia_M11);
raster5_tua_bia = georasterref('RasterSize',size(tuaE5_Bia_M11),'Latlim',[-60 90],'Lonlim',[-180 180]);
cmip5 = maps(1)
axes(cmip5)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
framem('FLineWidth',1)
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.3)
framem('FLineWidth',1)
h = geoshow(tuaE5_Bia_M11,raster5_tua_bia, 'DisplayType','surface','Zdata',zeros(size(tuaE5_Bia_M11)),'CData',tuaE5_Bia_M11);
colormap(mycolor_tuaE_bias)
caxis([-100 20])
set(gca,'box','off')
setm(cmip5, 'FFaceColor', [1,1,1])
framem('FLineWidth',1,'FEdgeColor','k')
axis off
text(-2.284,2.0824, '(a) CMIP5','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);

% CMIP6 tuaE bias
tuaE6_bia_M11 = tuaE6_bias;
tuaE6_bia_M11(302:360,:) = [];
tuaE6_bia_M11 = flipud(tuaE6_bia_M11);
raster6_tuaE_bia = georasterref('RasterSize',size(tuaE6_bia_M11),'Latlim',[-60 90],'Lonlim',[-180 180]);   
cmip6 = maps(2)
axes(cmip6)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
framem('FLineWidth',1)
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.3)
framem('FLineWidth',1)
h = geoshow(tuaE6_bia_M11,raster6_tuaE_bia, 'DisplayType','surface','Zdata',zeros(size(tuaE6_bia_M11)),'CData',tuaE6_bia_M11);
colormap(mycolor_tuaE_bias)
caxis([-100 20])
set(gca,'box','off')
setm(cmip6, 'FFaceColor', [1,1,1])
framem('FLineWidth',1,'FEdgeColor','k')
axis off
text(-2.284,2.0824, '(c) CMIP6','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);    
h1 = colorbar
h1.Location = 'southoutside'
h1.Position = [0.04485,0.3864,0.609,0.0261];
h1.FontName = 'Arial'
h1.FontSize = 10;
text(-0.0033,-2.2663,'model-data difference in \tau_E (yr) ',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11)  

% Zonal mean plots for CMIP5 and CMIP6
tuaE5_sp05(tuaE5_sp05>10^5) = NaN;
tuaE6_sp05(tuaE6_sp05>10^5) = NaN;
tuaE5_zonal_11 = nanmean(tuaE5_sp05,2); tuaE5_zonal_11 = squeeze(tuaE5_zonal_11); tuaE5_zonal_11(302:360,:) = [];
tuaE6_zonal_11 = nanmean(tuaE6_sp05,2); tuaE6_zonal_11 = squeeze(tuaE6_zonal_11); tuaE6_zonal_11(302:360,:) = [];
mycolor5 = [255 0 0; 153 51 255; 237 176 33; ...%BCC CanESM2 CCSM4
            0 197 205; 0 205 0; 207 194 124;...   %Had IPSL Miroc
            255 99 71; 65 105 255;...  %mpi NorM
            0 0 0; 158 131 149; 119 136 153]./255;  %BNU, GFDL, mri 
mycolor6 = [255 0 0; 153 51 255; 237 176 33;...    %BCC-CSM2-MR, CanESM5, CESM2
           0 197 205; 0 205 0; 207 194 124;...     %UKESM1-0-LL, IPSL-CM6A-LR,MIROC-ES2L
           255 99 71; 65 105 255;...               %MPI-ESM1-2-LR,NorESM-LM
            0 0 0; 158 131 149; 119 136 153]./255; %ACCESS-ESM1-5, CNRM-ESM2-1  EC-Earth3-Veg   
tuaE_zonal_obs = nanmean(tuaE_median1m_obs,2); tuaE_zonal_obs(302:360,:) = [];  
tuaE_zonal_10th = nanmean(tuaE_10th1m_obs,2); tuaE_zonal_10th(302:360,:) = [];  
tuaE_zonal_90th = nanmean(tuaE_90th1m_obs,2); tuaE_zonal_90th(302:360,:) = [];
boudary10_90th = [tuaE_zonal_10th tuaE_zonal_90th] - tuaE_zonal_obs;
boudary10_90th(:,1) = -boudary10_90th(:,1); boudary10_90th(isnan(boudary10_90th)) = 0;

% CMIP5 zonal mean plot
panel = tight_subplot(2,1,[0.032 0.01],[0.42 0.01],[0.73 0.02])    
axes(panel(1))
hold on
X_info = gca
X_info.XScale = 'log';
X_info.XAxis.Limits = [1, 1000];
X_info.XAxis.TickValues = [10 100 1000];
lat = 90:-0.5:-60;
for yr=1:11
    CMIP5_lines(yr) = plot(tuaE5_zonal_11(:,yr),lat,'color',[0.67,0.91,1.00],'LineWidth',0.9)
end
tuaE5_zonal_avg = nanmean(tuaE5_zonal_11,2);
CMIP5_avg = plot(tuaE5_zonal_avg,lat,'LineWidth',1.8,'color',[0.02,0.64,0.91])
%tuaE5_zonal_SD = nanstd(tuaE5_zonal_11,0,2);
%CMIP5_avg = boundedline(tuaE5_zonal_avg,lat,tuaE5_zonal_SD,'alpha','transparency',0.6,...
%'orientation', 'horiz','cmap',[0.30,0.75,0.93],'nan','remove');
%set(CMIP5_avg,'LineWidth',1.8,'color',[0.02,0.64,0.91]);
Obs_line = boundedline(tuaE_zonal_obs,lat,boudary10_90th,'alpha','transparency',0.4,...
'orientation', 'horiz','cmap',[0.15,0.15,0.15],'nan','remove');
set(Obs_line,'LineWidth',1.8);
text(1.3,82, '(b)',...
        'FontName','Arial','FontSize',11)
set(gca, 'YLim',[-60 90]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xticklabels([]);
plot([35, 95],[20 20],'k-','LineWidth',1.8)
text(105, 20,'Obs','Fontname','Arial','FontSize',9)
plot([35, 95],[10 10],'LineWidth',1.8,'color',[0.02,0.64,0.91])
text(105, 10,'model','Fontname','Arial','FontSize',9)
plot([1, 1000],[0 0],'k--','LineWidth',1)


% CMIP6 zonal mean polt
axes(panel(2))
hold on
X_info = gca
X_info.XAxis.Scale = 'log';
X_info.XAxis.Limits = [1, 1000];
X_info.XAxis.TickValues = [10 100 1000];
lat = 90:-1:-59;
lat = 90:-0.5:-60;
for yr=1:11
    CMIP6_lines(yr) = plot(tuaE6_zonal_11(:,yr),lat,'color',[0.99,0.76,0.99],'LineWidth',0.9)
end
tuaE6_zonal_avg = nanmean(tuaE6_zonal_11,2);
CMIP6_avg = plot(tuaE6_zonal_avg,lat,'LineWidth',1.8,'color',[1.00,0.07,0.65]);
%tuaE6_zonal_SD = nanstd(tuaE6_zonal_11,0,2);
Obs_line = boundedline(tuaE_zonal_obs,lat,boudary10_90th,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.15,0.15,0.15],'nan','remove');
set(Obs_line,'LineWidth',1.8);
text(1.3,82, '(d)',...
        'FontName','Arial','FontSize',11)
set(gca, 'YLim',[-60 90]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
%set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xlabel('\tau_E (yr)', 'FontName','Arial','FontSize',11)  
plot([35, 95],[20 20],'k-','LineWidth',1.8)
text(105, 20,'Obs','Fontname','Arial','FontSize',9)
plot([35, 95],[10 10],'LineWidth',1.8,'color',[1.00,0.07,0.65])
text(105, 10,'model','Fontname','Arial','FontSize',9)
plot([1, 1000],[0 0],'k--','LineWidth',1)

%%  Panel (e)-(g), comparison at global, circumpolar and non-circumpolar scales
%clear
%clc
% load cLand (global, circumplar, non-circumpolar) simulated from CMIP5 and CMIP6
load E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\NOY_circumpolar\NOYpolar_cLand.mat
% load GPP (global, circumplar, non-circumpolar) simulated from CMIP5 and CMIP6
load E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\NOY_circumpolar\NOYpolar_GPP56.mat
% circumpolar residence time
tuaE5_Polar11 = cPlar5_M11./GPP5_ALLpolar_Pg
tuaE6_Polar11 = cPlar6_M11./GPP6_ALLpolar_Pg
% non-circumpolar residence time
tuaE5_NonPolar11 = cNonPlar5_M11./GPP5_NonPolarGPP
tuaE6_NonPolar11 = cNonPlar6_M11./GPP6_NonPolarGPP
% global residence time
tuaE5_end5_ag = cLand5_M11./GPP5_end5_ag'
tuaE6_end5_ag = cLand6_M11./GPP6_end5_ag'
clearvars -except tuaE5_end5_ag tuaE6_end5_ag tuaE5_Polar11 tuaE6_Polar11 tuaE5_NonPolar11 tuaE6_NonPolar11
% Models from CMIP5 and CMIP6 used in this study
leg5_str = {'BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'}
leg6_str = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'}
% define Table 
tuaE5_bar = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','PolartuaE','NoPolartuaE','GlobaltuaE','Nlimitation'});
tuaE6_bar = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','PolartuaE','NoPolartuaE','GlobaltuaE','Nlimitation'});
% Put data into the defined table
tuaE5_bar.Moldes = leg5_str';
tuaE6_bar.Moldes = leg6_str';
tuaE5_bar.PolartuaE = tuaE5_Polar11;
tuaE6_bar.PolartuaE = tuaE6_Polar11;
tuaE5_bar.NoPolartuaE = tuaE5_NonPolar11;
tuaE6_bar.NoPolartuaE = tuaE6_NonPolar11;
tuaE5_bar.GlobaltuaE = tuaE5_end5_ag;
tuaE6_bar.GlobaltuaE = tuaE6_end5_ag;
NL5 = []; NL6 = [];              % number of models considered nutrient limitation
NL5 = [0 0 1 0 0 0 0 1 0 0 0]';  % 0, model without explicit nutrient limitation 
NL6 = [0 0 1 1 1 1 1 1  1 0 1]'; % 1, model with explicit nutrient limitation
tuaE5_bar.Nlimitation = NL5;
tuaE6_bar.Nlimitation = NL6;

% Benchmarking data from Fan et al., 2020 Earth System Science Data, Table2
% cVeg (PgC)
NonPlar_Veg_obs = [358 412 399 393];
Plar_Veg_obs = [49 39 38 42];
Global_Veg_obs = [407 451 437 435];
% Csoil 0-1m (PgC)
NonPlar_soil_obs = [1215 1399 1305 764];
Plar_soil_obs = [510 796 787 568 567];
Global_soil_obs = [1725 2195 2091 1332];

% global land carbon storage was estimated based on random combination (PgC)
cLand_obs = [];
cPlar_obs = [];
cNoPlar_obs = [];
for i = 1:4
    cLand  = Global_Veg_obs(i) + Global_soil_obs ;
    C_Plar = Plar_Veg_obs(i) + Plar_soil_obs;
    C_nonPlar = NonPlar_Veg_obs(i) + NonPlar_soil_obs;
    
    cLand_obs = [cLand_obs cLand];
    cPlar_obs = [cPlar_obs C_Plar];
    cNoPlar_obs = [cNoPlar_obs C_nonPlar];
end

% GPP used in Fan et al., 2020 (PgC yr-1)
GPP_polar_mean = 7;
GPP_polar_P10 = 6;
GPP_polar_P90 = 8;

GPP_NOpolar_mean = 96;
GPP_NOpolar_P10 = 92;
GPP_NOpolar_P90 = 99;

GPP_glb_mean = 102;
GPP_glb_P10 = 99;
GPP_glb_P90 = 106;

% circumpolar residence time
obs_polar_tuaE = [cPlar_obs./GPP_polar_mean cPlar_obs./GPP_polar_P10 cPlar_obs./GPP_polar_P90]
obs_polar_tuaE_mean = nanmean(obs_polar_tuaE)
obs_polar_tuaE_SD = nanstd(obs_polar_tuaE)
% non-circumpolar residence time
obs_NOpolar_tuaE = [cNoPlar_obs./GPP_NOpolar_mean cNoPlar_obs./GPP_NOpolar_P10 cNoPlar_obs./GPP_NOpolar_P90]
obs_NOpolar_tuaE_mean = nanmean(obs_NOpolar_tuaE)
obs_NOpolar_tuaE_SD = nanstd(obs_NOpolar_tuaE)
% global terresrtial tuaE
obs_glb_tuaE = [cLand_obs./GPP_glb_mean cLand_obs./GPP_glb_P10 cLand_obs./GPP_glb_P90]
obs_glb_tuaE_mean = nanmean(obs_glb_tuaE)
obs_glb_tuaE_SD = nanstd(obs_glb_tuaE)

% save data used to produce Fig.2 
save('E:\1_Mycase\3_CMIP56_Cland\5_Version4\sp_ExtendFigure1\MataData\tauE.mat',...
    'tuaE5_bar','tuaE6_bar','obs_polar_tuaE_mean','obs_polar_tuaE_SD','obs_NOpolar_tuaE_mean',...
    'obs_NOpolar_tuaE_SD','obs_glb_tuaE_mean','obs_glb_tuaE_SD')


% Panel(e): Global terrestrial carbon residence time
panel2 = tight_subplot(1,3,[0.032 0.08],[0.07 0.70],[0.12 0.04]) 
avgtuaE5 = nanmean(tuaE5_bar.GlobaltuaE); avgtuaE6 = nanmean(tuaE6_bar.GlobaltuaE);
sdtuaE5 = nanstd(tuaE5_bar.GlobaltuaE); sdtuaE6 = nanstd(tuaE6_bar.GlobaltuaE);
axes(panel2(1))
hold on
for i = 1:11
    if tuaE5_bar.Nlimitation(i) == 0 
        leg5(i)= plot(1,tuaE5_bar.GlobaltuaE(i),'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,tuaE5_bar.GlobaltuaE(i),'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if tuaE6_bar.Nlimitation(i) == 0 
        leg6(i) = plot(2,tuaE6_bar.GlobaltuaE(i),'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,tuaE6_bar.GlobaltuaE(i),'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[avgtuaE5-sdtuaE5,avgtuaE5+sdtuaE5,avgtuaE5+sdtuaE5,avgtuaE5-sdtuaE5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgtuaE6-sdtuaE6,avgtuaE6+sdtuaE6,avgtuaE6+sdtuaE6,avgtuaE6-sdtuaE6],[1.00,0.65,0.87]);
line([0.75 1.25],[avgtuaE5 avgtuaE5],'color','k','linewidth',1.8);
line([1.75 2.25],[avgtuaE6 avgtuaE6],'color','k','linewidth',1.8);
set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.4)
x=[0.5 2.5 2.5 0.5];
y=[obs_glb_tuaE_mean-obs_glb_tuaE_SD obs_glb_tuaE_mean-obs_glb_tuaE_SD obs_glb_tuaE_mean+obs_glb_tuaE_SD obs_glb_tuaE_mean+obs_glb_tuaE_SD];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YLim',[0 30],'XLim',[0.5 2.5],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
set(gca,'YTickLabelMode','auto');
ylabel('\tau_E (yr)','Fontname','Arial','FontSize',12)
text(0.6, 27.7,'(e)','Fontname','Arial','FontSize',12)
text(1.2, -1.5,'Global','Fontname','Arial','FontSize',11)

% Panel(f): circumpolar-region terrestrial carbon residence time
avgTuaE5_polar = nanmean(tuaE5_bar.PolartuaE); avgTuaE6_polar = nanmean(tuaE6_bar.PolartuaE);
sdTuaE_polar5 = nanstd(tuaE5_bar.PolartuaE); sdTuaE6_polar = nanstd(tuaE6_bar.PolartuaE);
axes(panel2(2))
hold on
for i = 1:11
    if tuaE5_bar.Nlimitation(i) == 0 
        leg5(i)= plot(1,tuaE5_bar.PolartuaE(i),'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,tuaE5_bar.PolartuaE(i),'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if tuaE6_bar.Nlimitation(i) == 0 
        leg6(i) = plot(2,tuaE6_bar.PolartuaE(i),'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,tuaE6_bar.PolartuaE(i),'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end

H_pa5 = patch([0.75,0.75,1.25,1.25],[avgTuaE5_polar-sdTuaE_polar5,avgTuaE5_polar+sdTuaE_polar5,avgTuaE5_polar+sdTuaE_polar5,avgTuaE5_polar-sdTuaE_polar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgTuaE6_polar-sdTuaE6_polar,avgTuaE6_polar+sdTuaE6_polar,avgTuaE6_polar+sdTuaE6_polar,avgTuaE6_polar-sdTuaE6_polar],[1.00,0.65,0.87]);
line([0.75 1.25],[avgTuaE5_polar avgTuaE5_polar],'color','k','linewidth',1.8);
line([1.75 2.25],[avgTuaE6_polar avgTuaE6_polar],'color','k','linewidth',1.8);
set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.4)
set(gca,'XLim',[0.5 2.5],'YLim',[0 150],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
x=[0.5 2.5 2.5 0.5];
y=[obs_polar_tuaE_mean-obs_polar_tuaE_SD obs_polar_tuaE_mean-obs_polar_tuaE_SD obs_polar_tuaE_mean+obs_polar_tuaE_SD obs_polar_tuaE_mean+obs_polar_tuaE_SD];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
text(0.6, 140,'(f)','Fontname','Arial','FontSize',12)
text(0.9, -7,'circumpolar','Fontname','Arial','FontSize',11)

% Panel(g): non-circumpolar terrestrial carbon residence time
avgTuaE5_NOpolar = nanmean(tuaE5_bar.NoPolartuaE); avgTuaE6_NOpolar = nanmean(tuaE6_bar.NoPolartuaE);
sdTuaE_NOpolar5 = nanstd(tuaE5_bar.NoPolartuaE); sdTuaE6_NOpolar = nanstd(tuaE6_bar.NoPolartuaE);
axes(panel2(3))
hold on
for i = 1:11
    if tuaE5_bar.Nlimitation(i) == 0 
        leg5(i)= plot(1,tuaE5_bar.NoPolartuaE(i),'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,tuaE5_bar.NoPolartuaE(i),'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if tuaE6_bar.Nlimitation(i) == 0 
        leg6(i) = plot(2,tuaE6_bar.NoPolartuaE(i),'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,tuaE6_bar.NoPolartuaE(i),'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[avgTuaE5_NOpolar-sdTuaE_NOpolar5,avgTuaE5_NOpolar+sdTuaE_NOpolar5,avgTuaE5_NOpolar+sdTuaE_NOpolar5,avgTuaE5_NOpolar-sdTuaE_NOpolar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgTuaE6_NOpolar-sdTuaE6_NOpolar,avgTuaE6_NOpolar+sdTuaE6_NOpolar,avgTuaE6_NOpolar+sdTuaE6_NOpolar,avgTuaE6_NOpolar-sdTuaE6_NOpolar],[1.00,0.65,0.87]);
line([0.75 1.25],[avgTuaE5_NOpolar avgTuaE5_NOpolar],'color','k','linewidth',1.8);
line([1.75 2.25],[avgTuaE6_NOpolar avgTuaE6_NOpolar],'color','k','linewidth',1.8);
set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.4)
set(gca,'XLim',[0.5 2.5],'YLim',[0 30],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
x=[0.5 2.5 2.5 0.5];
y=[obs_NOpolar_tuaE_mean-obs_NOpolar_tuaE_SD obs_NOpolar_tuaE_mean-obs_NOpolar_tuaE_SD obs_NOpolar_tuaE_mean+obs_NOpolar_tuaE_SD obs_NOpolar_tuaE_mean+obs_NOpolar_tuaE_SD];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
text(0.6, 27.7,'(g)','Fontname','Arial','FontSize',12)
text(0.6, -1.5,'non-circumpolar','Fontname','Arial','FontSize',11)
leg_panel =legend([H_pa5 H_pa6 obs],{'CMIP5','CMIP6','Obs'},'NumColumns',3)
set(leg_panel,'FontName','Arial','FontSize',10,'Position',[0.2848,0.0026,0.5196,0.0318],...
    'color','w','EdgeColor','k')

