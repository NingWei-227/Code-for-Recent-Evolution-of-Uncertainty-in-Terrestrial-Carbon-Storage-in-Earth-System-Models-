% Figure2: data-model comparison on GPP
% observational GPP data: FLUXCOM
%                         MODIS17A2
%                         VPM
% units were unified as Kg C m-2 year-1
clear;clc;

%% load observations 
% load observations, Time: from 2001 to 2005;
% Data: FLUXCOM;  Unit: g C m-2 day-1
GPP_maps_1 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\FLUXCOM\yearlyGPP1980_2013.nc','GPP',[1 1 22],[Inf Inf 5]);
GPP_obs1(1:360,1:720,1:5) = NaN;
for i =1:5
    i
    GPPobs_M = GPP_maps_1(:,:,i);
    GPPobs_t = GPPobs_M';
    GPPobs_t = GPPobs_t.*365;  % convert unit to gC m-2 yr-1
    
    GPP_obs1(:,:,i) = GPPobs_t;   
end
% Data: MODIS17A2; Unit: gC m-2 yr-1
for i = 2001:2005
    [num,text,raw] = xlsread(['E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\MODIS17A2\MOD17A2GPPAnnual',num2str(i),'.csv']);
    for rowID=1:360
        for colID =1:720
            if raw{rowID,colID}=='NA'
                raw{rowID,colID}=nan;
            end  
        end
     end
    GPP_obs2(:,:,i-2000) = cell2mat(raw);
end
% Data: VPM; Unit: gC m-2 yr-1
for i = 2001:2005
    [num,text,raw] = xlsread(['E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\VPM\VPM_GPP_',num2str(i),'_Mean.csv']);
    for rowID=1:360
        for colID =1:720
            if raw{rowID,colID}=='NA'
                raw{rowID,colID}=nan;
            end  
        end
     end
    GPP_obs3(:,:,i-2000) = cell2mat(raw);
end
% unit: gC m-2 yr-1; 2001-2005 mean 
GPP_obs05(:,:,1) = nanmean(GPP_obs1,3);
GPP_obs05(:,:,2) = nanmean(GPP_obs2,3);
GPP_obs05(:,:,3) = nanmean(GPP_obs3,3);

%% CMIP5 and CMIP6
% load GPP from CMIP5 and CMIP6, 
% time: 2001-2005 ; 
% unit: KgC m-2 yr-1
% CMIP5
cd('G:\CMIP5\4_MatData\spatial_data\hist_rcp85')
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
% 2001-2005 mean GPP
sp5_GPP5_11 = nanmean(sp5_GPP5,3);
sp5_GPP5_11 = squeeze(sp5_GPP5_11); 
sp5_GPP5_11(:,1:30,:) = NaN;        % omit South pole
% unify terresrtial region 
mask = sp5_GPP5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp5_GPP5_11 = sp5_GPP5_11.*mask;

% CMIP6
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
% 2001-2005 mean GPP
sp6_GPP5_11 = nanmean(sp6_GPP5,3);
sp6_GPP5_11 = squeeze(sp6_GPP5_11);
sp6_GPP5_11(:,1:30) = NaN;    % omit South pole
% unify terresrtial region
mask = sp6_GPP5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp6_GPP5_11 = sp6_GPP5_11.*mask;
    
clearvars -except sp5_GPP5_11 sp6_GPP5_11
% change modeled GPP to match global map              
GPP5_map11(1:180,1:360,1:11) = NaN;
GPP6_map11(1:180,1:360,1:11) = NaN;
for i =1:11
    i
    
    GPP5_M = sp5_GPP5_11(:,:,i);
    GPP6_M = sp6_GPP5_11(:,:,i);
    
    GPP5_t = GPP5_M';
    GPP6_t = GPP6_M';
    
    GPP5_LR(1:180,1:360) = NaN;
    GPP6_LR(1:180,1:360) = NaN;
    
    GPP5_LR(:,1:180) = GPP5_t(:,181:360);
    GPP5_LR(:,181:360) = GPP5_t(:,1:180);
    GPP6_LR(:,1:180) = GPP6_t(:,181:360);
    GPP6_LR(:,181:360) = GPP6_t(:,1:180);
    
    map_GPP5(1:180,1:360) = NaN;
    map_GPP6(1:180,1:360) = NaN;
    for k=1:180
        map_GPP5(181-k,:) = GPP5_LR(k,:);
        map_GPP6(181-k,:) = GPP6_LR(k,:);
    end
    
    GPP5_map11(:,:,i) = map_GPP5;
    GPP6_map11(:,:,i) = map_GPP6;
      
end
clearvars -except GPP5_map11 GPP6_map11 GPP_obs05
save E:\1_Mycase\3_CMIP56_Cland\4_Version3\1_Figure\2_Figure2_X_GPP_tauE\Model_obsData\cmip56_GPP.mat
% regrided CMIP5 and CMIP6 data into 0.5x0.5 resolution
lat05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lat');
lon05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lon');
lat1 = 89.5:-1:-89.5; lat1 = lat1';
lon1 = -179.5:179.5;  lon1 = lon1';

[x05,y05] = meshgrid(lon05,lat05);
[x1,y1] = meshgrid(lon1,lat1);
% regrided CMIP5 and CMIP6 data into 0.5x0.5 resolution
GPP5_sp05(1:360,1:720,1:11) = NaN;
GPP6_sp05(1:360,1:720,1:11) = NaN;
for i=1:11
    GPP5_M = interp2(x1,y1,GPP5_map11(:,:,i),x05,y05,'linear');
    GPP6_M = interp2(x1,y1,GPP6_map11(:,:,i),x05,y05,'linear');
    
    GPP5_sp05(:,:,i) = GPP5_M;
    GPP6_sp05(:,:,i) = GPP6_M;
end
% For consistency, deserts and other regions where GPP<0.01 kgC m-2 yr-1 were omitted based on the criteria of Fan et al., 2019
% load Observational dataset from Fan et al., 2019 and creat mask 
NoDesert_map_obs = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','tau');
NoDesert_map_obs = NoDesert_map_obs(:,:,1)';   % considering 0-1m soil C
NoDesert_map_obs(~isnan(NoDesert_map_obs)) = 1;
mask = NoDesert_map_obs;

% convert unit to KgC m-2 yr-1
GPP_obs05KG = GPP_obs05.*10^(-3);  
% Unit of GPP(both simulated and observed): KgC m-2 yr-1
% omit region where GPP < 0.01 KgC m-2 yr-1
GPP5_sp05(GPP5_sp05<=0.01) = nan;
GPP6_sp05(GPP6_sp05<=0.01) = nan;
GPP_obs05KG(GPP_obs05KG<=0.01) = nan;

clearvars -except GPP5_map11 GPP6_map11 GPP_obs05 ...
                  GPP5_sp05 GPP6_sp05 GPP_obs05KG mask
              
GPP_obs05KG = GPP_obs05KG.*mask;
GPP5_sp05 = GPP5_sp05.*mask;
GPP6_sp05 = GPP6_sp05.*mask;            
% calculate across-model mean and standard deviation
GPP_avg11_cmip5 = nanmean(GPP5_sp05,3);
GPP_avg11_cmip6 = nanmean(GPP6_sp05,3);
GPP_sd11_cmip5 = nanstd(GPP5_sp05,0,3);
GPP_sd11_cmip6 = nanstd(GPP6_sp05,0,3);
% observational mean
GPP_obsKgC_avg = nanmean(GPP_obs05KG,3);
% Model bias in GPP 
GPP5_bias_glb = GPP_avg11_cmip5 - GPP_obsKgC_avg;
GPP6_bias_glb = GPP_avg11_cmip6 - GPP_obsKgC_avg;


%% Extended Data Fig. 3: panel(a)-(d)
%ax = gca;
%mycolor_GPPsd = colormap(ax)
%save ('E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure2_GPP\mycolor_GPPsd.mat','mycolor_GPPsd') 
%load('E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure2_GPP\mycolor_GPPsd.mat');

% load RGB code for colorbar
load E:\1_Mycase\3_CMIP56_Cland\4_Version3\1_Figure\2_Figure2_X_GPP_tauE\mycolor_cLand_bias.mat
figure
set(gcf,'position',[100 100 499.2,630])
maps = tight_subplot(2,1,[-0.1 -0.085],[0.35 -0.05],[-0.04 0.26])

% CMIP5 GPP bias
GPP5_bia_M11 = GPP5_bias_glb;
GPP5_bia_M11(302:360,:) = [];
GPP5_bia_M11 = flipud(GPP5_bia_M11);
raster5_GPP_bias = georasterref('RasterSize',size(GPP5_bia_M11),'Latlim',[-60 90],'Lonlim',[-180 180]);
cmip5 = maps(1)
axes(cmip5)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
framem('FLineWidth',1)
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.3)
framem('FLineWidth',1)
h = geoshow(GPP5_bia_M11,raster5_GPP_bias, 'DisplayType','surface','Zdata',zeros(size(GPP5_bia_M11)),'CData',GPP5_bia_M11);
colormap(mycolor_cLand_bias)
caxis([-1.2 1.2])
set(gca,'box','off')
setm(cmip5, 'FFaceColor', [1,1,1])
framem('FLineWidth',1,'FEdgeColor','k')
axis off
text(-2.284,2.0824, '(a) CMIP5','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);


% CMIP6 GPP bias
GPP6_bias_M11 = GPP6_bias_glb;
GPP6_bias_M11(302:360,:) = [];
GPP6_bias_M11 = flipud(GPP6_bias_M11);
raster6_GPP_bias = georasterref('RasterSize',size(GPP6_bias_M11),'Latlim',[-60 90],'Lonlim',[-180 180]);       
cmip6 = maps(2)
axes(cmip6)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
framem('FLineWidth',1)
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.3)
framem('FLineWidth',1)
h = geoshow(GPP6_bias_M11,raster6_GPP_bias, 'DisplayType','surface','Zdata',zeros(size(GPP6_bias_M11)),'CData',GPP6_bias_M11);
colormap(mycolor_cLand_bias)
caxis([-1.2 1.2])
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
text(-0.0033,-2.2663,'model-data difference in GPP (KgC m^-^2 yr^-^1)',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11)    
    
% Zonal mean plots for CMIP5 and CMIP6 models
GPP5_sp05(GPP5_sp05>10^3) = NaN;
GPP6_sp05(GPP6_sp05>10^3) = NaN;
GPP5_zonal_11 = nanmean(GPP5_sp05,2); GPP5_zonal_11 = squeeze(GPP5_zonal_11); GPP5_zonal_11(302:360,:) = [];
GPP6_zonal_11 = nanmean(GPP6_sp05,2); GPP6_zonal_11 = squeeze(GPP6_zonal_11); GPP6_zonal_11(302:360,:) = []
% Observational data
GPP_zonal_obs = nanmean(GPP_obs05KG,2); GPP_zonal_obs = nanmean(GPP_zonal_obs,3); GPP_zonal_obs(302:360,:) = [];  
GPP_zonal5 = nanmean(GPP_obs05KG,2); GPP_zonal5 = squeeze(GPP_zonal5); GPP_zonal5(302:360,:) = []; 
GPP_zonal_min = min(GPP_zonal5,[],2);
GPP_zonal_max = max(GPP_zonal5,[],2);
boudary_GPP = [GPP_zonal_min GPP_zonal_max] - GPP_zonal_obs;
boudary_GPP(:,1) = -boudary_GPP(:,1); boudary_GPP(isnan(boudary_GPP)) = 0;

% CMIP5
panel = tight_subplot(2,1,[0.032 0.01],[0.42 0.01],[0.73 0.02])    
axes(panel(1))
hold on
lat = 90:-0.5:-60;
for i=1:11
    CMIP5_lines(i) = plot(GPP5_zonal_11(:,i),lat,'color',[0.67,0.91,1.00],'LineWidth',0.9)
end
GPP5_zonal_avg = nanmean(GPP5_zonal_11,2);
GPP5_zonal_SD = nanstd(GPP5_zonal_11,0,2);
CMIP5_avg = plot(GPP5_zonal_avg,lat,'LineWidth',1.8,'color',[0.02,0.64,0.91])
Obs_line = boundedline(GPP_zonal_obs,lat,boudary_GPP,'alpha','transparency',0.6,...
'orientation', 'horiz','cmap',[0.15,0.15,0.15],'nan','remove');
set(Obs_line,'LineWidth',1.8);
text(0.3,80, '(b)',...
        'FontName','Arial','FontSize',11)
set(gca, 'YLim',[-60 90],'XLim',[0 6]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xticks([0 2 4 6])
xticklabels([]);
plot([0, 6],[0 0],'k--','LineWidth',1)
plot([2.8 3.6],[64 64],'k-','LineWidth',1.8)
text(3.7,64.2,'Obs','Fontname','Arial','Fontsize',10)
plot([2.8 3.6],[52 52],'LineWidth',1.8,'color',[0.02,0.64,0.91])
text(3.7,52,'model','Fontname','Arial','Fontsize',10)


% CMIP6
axes(panel(2))
hold on
lat = 90:-0.5:-60;
for i=1:11
    CMIP6_lines(i) = plot(GPP6_zonal_11(:,i),lat,'color',[0.99,0.76,0.99],'LineWidth',0.9)
end
GPP6_zonal_avg = nanmean(GPP6_zonal_11,2);
GPP6_zonal_SD = nanstd(GPP6_zonal_11,0,2);
CMIP6_avg = plot(GPP6_zonal_avg,lat,'LineWidth',1.8,'color',[1.00,0.07,0.65])
Obs_line = boundedline(GPP_zonal_obs,lat,boudary_GPP,'alpha','transparency',0.6,...
'orientation', 'horiz','cmap',[0.15,0.15,0.15],'nan','remove');
set(Obs_line,'LineWidth',1.8);
text(0.3,80, '(d)',...
        'FontName','Arial','FontSize',11)
set(gca, 'YLim',[-60 90],'XLim',[0 6]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
%set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xticks([0 2 4 6]);
plot([0, 6],[0 0],'k--','LineWidth',1)  
xlabel(['GPP',newline,'(KgC m^-^2 yr^-^1)'],'FontName','Arial','FontSize',11)
plot([2.8 3.6],[64 64],'k-','LineWidth',1.8)
text(3.7,64.2,'Obs','Fontname','Arial','Fontsize',10)
plot([2.8 3.6],[52 52],'LineWidth',1.8,'color',[1.00,0.07,0.65])
text(3.7,52,'model','Fontname','Arial','Fontsize',10)


%% Extended Data Fig. 3: panel (e)-(g), circumpolar, non-circumpolar and global terrestrial GPP
clearvars -except GPP5_sp05 GPP6_sp05 GPP_obs05
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\SoilC\Circumpolar_cSoil\NCSCDv2_Circumpolar_netCDF_05deg\NCSCD05_1mkg.mat')
% the classification of circumpolar and non-circumpolar regions is based on NCSCD data
polar_mask = NCSCDgb_1mKg;
polar_mask(~isnan(polar_mask)) = 1;     % mask for circumpolar region
Noploar_mask = polar_mask;
Noploar_mask(Noploar_mask == 1) = 0;
Noploar_mask(isnan(Noploar_mask)) = 1;
Noploar_mask(Noploar_mask == 0) = NaN;  % mask for non-circumpolar region

% global terrestrial GPP simulated from CMIP5 models
% calculate at the native resolution of models
cd('G:\CMIP5\4_MatData\temporal_data\hist_rcp85')
load('GPP_cmip5_tmp.mat')
NaNgf = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]';
GPPgf_tmp = [NaNgf; GPPgf_tmp];
NaNhad = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]';
GPPhad_tmp = [NaNhad; GPPhad_tmp];
GPPmri_tmp = [NaN; GPPmri_tmp];
GPP5_tmp(1:155,1:11) = NaN;
GPP5_tmp(:,1) = GPPbcc_tmp(2:156);
GPP5_tmp(:,2) = GPPcan_tmp(2:156);
GPP5_tmp(:,3) = GPPccsm_tmp(2:156);
GPP5_tmp(:,4) = GPPhad_tmp(2:156);
GPP5_tmp(:,5) = GPPipsl_tmp(2:156);
GPP5_tmp(:,6) = GPPmiroc_tmp(2:156);
GPP5_tmp(:,7) = GPPmpi_tmp(2:156);
GPP5_tmp(:,8) = GPPnor_tmp(2:156);
GPP5_tmp(:,9) = GPPbnu_tmp(2:156);
GPP5_tmp(:,10) = GPPgf_tmp(2:156);
GPP5_tmp(:,11) = GPPmri_tmp(2:156);
 
% global terrestrial GPP simulated from CMIP6 models
cd('H:\CMIP56_Csink\4_MatData\temporal_data\hist_ssp585')
load('2GPP_cmip6_tmp.mat')
GPP6_tmp = [GPPbcc_tmp(2:156),GPPcan_tmp(2:156),GPPcesm_tmp(2:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   GPPuk_tmp(2:156),GPPipsl_tmp(2:156),GPPmic_tmp(2:156),...          % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   GPPmpi_tmp(2:156),GPPnor_tmp(2:156),...                           % MPI-ESM1-2-LR, NorESM2
   GPPass_tmp(2:156),GPPcnrm_tmp(2:156),GPPec_tmp(2:156)];            % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg
% Models from CMIP5 and CMIP6 used in this study
Models5 = {'BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'};
Models6 = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'}; 
% global terresrtial GPP (PgC yr-1), over 2001-2005  
% models
GPP5_end5_ag = nanmean(GPP5_tmp(151:155,:),1);
GPP6_end5_ag = nanmean(GPP6_tmp(151:155,:),1); 
% observation
area05 = csvread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\CABLEGridAreaM2.csv');
GPP_obs05_Map_Pg = GPP_obs05.* area05.*10^(-15);    % Unit: PgC yr-1
GPP_obs05_Map_Pg(GPP_obs05_Map_Pg<=0) = NaN;
GPP_obs05_Pg = nansum(GPP_obs05_Map_Pg,1); GPP_obs05_Pg = nansum(GPP_obs05_Pg,2);
GPP_obs05_Pg = squeeze(GPP_obs05_Pg)

% circumpolar-region GPP
% CMIP5
GPP5_sp05_map11_Pg = GPP5_sp05.* area05.*10^(-12);  % Unit: PgC yr-1
GPP5_sp05_map11_Pg(GPP5_sp05_map11_Pg<=0) = NaN;
GPP5_sp05_polar_Pg = GPP5_sp05_map11_Pg.* polar_mask;
GPP5_ALLpolar_Pg = nansum(GPP5_sp05_polar_Pg,1); GPP5_ALLpolar_Pg = nansum(GPP5_ALLpolar_Pg,2);
GPP5_ALLpolar_Pg = squeeze(GPP5_ALLpolar_Pg)
% CMIP6
GPP6_sp05_map11_Pg = GPP6_sp05.* area05.*10^(-12);  % Unit: PgC yr-1
GPP6_sp05_map11_Pg(GPP6_sp05_map11_Pg<=0) = NaN;
GPP6_sp05_polar_Pg = GPP6_sp05_map11_Pg.* polar_mask;
GPP6_ALLpolar_Pg = nansum(GPP6_sp05_polar_Pg,1); GPP6_ALLpolar_Pg = nansum(GPP6_ALLpolar_Pg,2);
GPP6_ALLpolar_Pg = squeeze(GPP6_ALLpolar_Pg)
% Observation
GPP_obs05_polar_Pg = GPP_obs05_Map_Pg.* polar_mask;
GPP_ALLpolar_Pg = nansum(GPP_obs05_polar_Pg,1); GPP_ALLpolar_Pg = nansum(GPP_ALLpolar_Pg,2);
GPP_ALLpolar_Pg = squeeze(GPP_ALLpolar_Pg);

% non-circumpolar GPP = Global - circumpolar
% models from CMIP5 and CMIP6 
GPP5_NonPolarGPP =  GPP5_end5_ag' - GPP5_ALLpolar_Pg;
GPP6_NonPolarGPP =  GPP6_end5_ag' - GPP6_ALLpolar_Pg;
% observation
GPP_obs05_NOpolar_Pg = GPP_obs05_Map_Pg.* Noploar_mask;
GPP_ALL_NOpolar_Pg = nansum(GPP_obs05_NOpolar_Pg,1); GPP_ALL_NOpolar_Pg = nansum(GPP_ALL_NOpolar_Pg,2);
GPP_ALL_NOpolar_Pg = squeeze(GPP_ALL_NOpolar_Pg)

% save GPP5_ALLpolar_Pg GPP6_ALLpolar_Pg GPP5_NonPolarGPP GPP6_NonPolarGPP for the calculation of ecosystem carbon residence time over circumpolar and non-circumpolar regions
% used for calculation of tauE
save('E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\NOY_circumpolar\NOYpolar_GPP56.mat',...
    'GPP5_ALLpolar_Pg', 'GPP6_ALLpolar_Pg', 'GPP5_NonPolarGPP', 'GPP6_NonPolarGPP','GPP5_end5_ag','GPP6_end5_ag')

% define table
GPP5_bar = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','PolarGPP','NonPolarGPP','globalGPP','Nlimitation'});
GPP6_bar = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','PolarGPP','NonPolarGPP','globalGPP','Nlimitation'});
% put data into Table  
GPP5_bar.Moldes = Models5';
GPP6_bar.Moldes = Models6';
GPP5_bar.PolarGPP = GPP5_ALLpolar_Pg;
GPP6_bar.PolarGPP = GPP6_ALLpolar_Pg;
GPP5_bar.globalGPP = GPP5_end5_ag';
GPP6_bar.globalGPP = GPP6_end5_ag';
GPP5_bar.NonPolarGPP =  GPP5_NonPolarGPP;
GPP6_bar.NonPolarGPP =  GPP6_NonPolarGPP;
NL5 = []; NL6 = []; % number of models considered nutrient limitation
NL5 = [0 0 1 0 0 0 0 1 0 0 0]';  % 0, model without explicit nutrient limitation
NL6 = [0 0 1 1 1 1 1 1  1 0 1]'; % 1, model with explicit nutrient limitation
GPP5_bar.Nlimitation = NL5; 
GPP6_bar.Nlimitation = NL6;

clearvars -except GPP5_bar GPP6_bar GPP_ALLpolar_Pg GPP_ALL_NOpolar_Pg GPP_obs05_Pg
save('E:\1_Mycase\3_CMIP56_Cland\5_Version4\sp_ExtendFigure1\MataData\GPP.mat',...
    'GPP5_bar','GPP6_bar','GPP_ALLpolar_Pg','GPP_ALL_NOpolar_Pg','GPP_obs05_Pg')

% Panel (e): Global terrestrial GPP
panel2 = tight_subplot(1,3,[0.032 0.08],[0.07 0.70],[0.12 0.04]) 
avgGPP5 = nanmean(GPP5_bar.globalGPP); avgGPP6 = nanmean(GPP6_bar.globalGPP);
sdGPP5 = nanstd(GPP5_bar.globalGPP); sdGPP6 = nanstd(GPP6_bar.globalGPP);
axes(panel2(1))
hold on
for i = 1:11
    if GPP5_bar.Nlimitation(i) == 0 
        leg5(i)= plot(1,GPP5_bar.globalGPP(i),'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,GPP5_bar.globalGPP(i),'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if GPP6_bar.Nlimitation(i) == 0 
        leg6(i) = plot(2,GPP6_bar.globalGPP(i),'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,GPP6_bar.globalGPP(i),'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[avgGPP5-sdGPP5,avgGPP5+sdGPP5,avgGPP5+sdGPP5,avgGPP5-sdGPP5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgGPP6-sdGPP6,avgGPP6+sdGPP6,avgGPP6+sdGPP6,avgGPP6-sdGPP6],[1.00,0.65,0.87]);
line([0.75 1.25],[avgGPP5 avgGPP5],'color','k','linewidth',1.8);
line([1.75 2.25],[avgGPP6 avgGPP6],'color','k','linewidth',1.8);
set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.5)
set(gca,'YLim',[80 260],'XLim',[0.5 2.5],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
x=[0.5 2.5 2.5 0.5];
y=[min(GPP_obs05_Pg) min(GPP_obs05_Pg) max(GPP_obs05_Pg) max(GPP_obs05_Pg)];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
ylabel('GPP (PgC yr^-^1)','Fontname','Arial','FontSize',12)
text(0.6, 242,'(e)','Fontname','Arial','FontSize',11)
text(1.2, 70,'Global','Fontname','Arial','FontSize',10)

% Panel (f): circumpolar-region GPP
avgGPP5_polar = nanmean(GPP5_bar.PolarGPP); avgGPP6_polar = nanmean(GPP6_bar.PolarGPP);
sdGPP_polar5 = nanstd(GPP5_bar.PolarGPP); sdGPP6_polar = nanstd(GPP6_bar.PolarGPP);
axes(panel2(2))
hold on
for i = 1:11
    if GPP5_bar.Nlimitation(i) == 0 
        leg5(i)= plot(1,GPP5_bar.PolarGPP(i),'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,GPP5_bar.PolarGPP(i),'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if GPP6_bar.Nlimitation(i) == 0 
        leg6(i) = plot(2,GPP6_bar.PolarGPP(i),'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,GPP6_bar.PolarGPP(i),'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[avgGPP5_polar-sdGPP_polar5,avgGPP5_polar+sdGPP_polar5,avgGPP5_polar+sdGPP_polar5,avgGPP5_polar-sdGPP_polar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgGPP6_polar-sdGPP6_polar,avgGPP6_polar+sdGPP6_polar,avgGPP6_polar+sdGPP6_polar,avgGPP6_polar-sdGPP6_polar],[1.00,0.65,0.87]);
line([0.75 1.25],[avgGPP5_polar avgGPP5_polar],'color','k','linewidth',1.8);
line([1.75 2.25],[avgGPP6_polar avgGPP6_polar],'color','k','linewidth',1.8);
set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.4)
set(gca,'XLim',[0.5 2.5],'YLim',[2 23],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
x=[0.5 2.5 2.5 0.5];
y=[min(GPP_ALLpolar_Pg) min(GPP_ALLpolar_Pg) max(GPP_ALLpolar_Pg) max(GPP_ALLpolar_Pg)];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
text(0.6, 21,'(f)','Fontname','Arial','FontSize',12)
text(0.9, 0.8,'circumpolar','Fontname','Arial','FontSize',11)

% Panel (g): non-circumpolar-region GPP
avgGPP5_NONpolar = nanmean(GPP5_bar.NonPolarGPP); avgGPP6_NONpolar = nanmean(GPP6_bar.NonPolarGPP);
sdGPP_NONpolar5 = nanstd(GPP5_bar.NonPolarGPP); sdGPP6_NONpolar = nanstd(GPP6_bar.NonPolarGPP);
axes(panel2(3))
hold on
for i = 1:11
    if GPP5_bar.Nlimitation(i) == 0 
        leg5(i)= plot(1,GPP5_bar.NonPolarGPP(i),'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else 
        leg5(i)= plot(1,GPP5_bar.NonPolarGPP(i),'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
    if GPP6_bar.Nlimitation(i) == 0 
        leg6(i) = plot(2,GPP6_bar.NonPolarGPP(i),'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    else
        leg6(i) = plot(2,GPP6_bar.NonPolarGPP(i),'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',8,'LineStyle','none','LineWidth',1.5)
    end
    
end
H_pa5 = patch([0.75,0.75,1.25,1.25],[avgGPP5_NONpolar-sdGPP_NONpolar5,avgGPP5_NONpolar+sdGPP_NONpolar5,avgGPP5_NONpolar+sdGPP_NONpolar5,avgGPP5_NONpolar-sdGPP_NONpolar5],[0.00,0.77,0.80]);
H_pa6 = patch([1.75,1.75,2.25,2.25],[avgGPP6_NONpolar-sdGPP6_NONpolar,avgGPP6_NONpolar+sdGPP6_NONpolar,avgGPP6_NONpolar+sdGPP6_NONpolar,avgGPP6_NONpolar-sdGPP6_NONpolar],[1.00,0.65,0.87]);
line([0.75 1.25],[avgGPP5_NONpolar avgGPP5_NONpolar],'color','k','linewidth',1.8);
line([1.75 2.25],[avgGPP6_NONpolar avgGPP6_NONpolar],'color','k','linewidth',1.8);
set(H_pa5,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
set(H_pa6,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.6)
set(gca,'XLim',[0.5 2.5],'YLim',[60 250],'LineWidth',1,'box','on')
set(gca,'Fontname','Arial','FontSize',10)
x=[0.5 2.5 2.5 0.5];
y=[min(GPP_ALL_NOpolar_Pg) min(GPP_ALL_NOpolar_Pg) max(GPP_ALL_NOpolar_Pg) max(GPP_ALL_NOpolar_Pg)];
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.4)
set(gca,'YTickLabelMode','auto');
text(0.6, 230,'(g)','Fontname','Arial','FontSize',12)
text(0.6, 52,'non-circumpolar','Fontname','Arial','FontSize',11)
leg_panel =legend([H_pa5 H_pa6 obs],{'CMIP5','CMIP6','Obs'},'NumColumns',3)
set(leg_panel,'FontName','Arial','FontSize',10,'Position',[0.2848,0.0026,0.5196,0.0318],...
    'color','w','EdgeColor','k')



