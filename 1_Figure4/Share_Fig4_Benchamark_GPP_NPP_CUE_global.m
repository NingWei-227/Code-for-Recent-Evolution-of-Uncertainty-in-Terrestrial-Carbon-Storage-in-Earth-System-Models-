%% Fig4. Model-data comparison on global terrestrial NPP, GPP and CUE for two CMIPs
% Observational datasets
% GPP datasets:
% (1) MODIS17A2; 
% Unit: gC m-2 yr-1
% Running et al., 2015
% resolution£º 0.5 x 0.5 
%
% (2) GIMMSGPP
% reference: Smith NCC, 2016; available from: https://wkolby.org/data-code/
% unit: gC m-2 yr-1 
% resolution: 0.5 x 0.5
% standard estimation based on climate inputs from ECMWF, MERRA2 and NCEP
% 
% (3) FLUXCOM
% reference: Jung et al., 2017; https://www.bgc-jena.mpg.de/geodb/projects/Home.php
% Unit: gC m-2 day-1
% resolution: 0.5 x 0.5
%
% (4) VPM
% reference: Zhang et al., 2017;  https://doi.org/10.6084/m9.figshare.c.3789814
% Unit: gC m-2 yr-1
% resolution: 0.5 x 0.5

% NPP datasets used here are from MODIS17A2 and GIMMSGPP-NPP
% CUE (NPP/GPP) is calculated by using NPP and GPP data from the same source £¨MODIS17A2 and GIMMSGPP-NPP£©
clear;clc
%% load observational data
% Data: MODIS17A2;
GPP_obs1(1:360,1:720,1:5) = NaN;
for i = 2001:2005
    [num,text,raw] = xlsread(['E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\MODIS17A2\MOD17A2GPPAnnual',num2str(i),'.csv']);
    for rowID=1:360
        for colID =1:720
            if raw{rowID,colID}=='NA'
                raw{rowID,colID}=nan;
            end  
        end
     end
    GPP_obs1(:,:,i-2000) = cell2mat(raw);
end
NPP_obs1(1:360,1:720,1:5) = NaN;
for i=2001:2005
    nppData = imread(['E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\MOD17A3_0.5_2016fool\MOD17A3_',num2str(2000),'globalNPP.tif']);
    nppData = nppData.* 0.1; % Scale factor:0.1
    nppData(nppData<0) = NaN;
    nor_matrix(1:20,1:720) = NaN;
    sou_matrix(1:60,1:720) = NaN;
    
    nppMap = [nor_matrix; nppData; sou_matrix];
    NPP_obs1(:,:,i-2000) = nppMap;
end
CUE_obs1 =  NPP_obs1./GPP_obs1;

% Data: GIMMS GPP/NPP
% estimated based on ECMWF climate data
gpp_ECMWF = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\GIMMSNPP\GIMMSNPP\Annual\Annual\ECMWF\gpp_V4_Standard_1982_2015_Annual_GEO_30min.nc','GPP',[1 1 20],[Inf Inf 5]);
npp_ECMWF = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\GIMMSNPP\GIMMSNPP\Annual\Annual\ECMWF\npp_V4_Standard_1982_2015_Annual_GEO_30min.nc','NPP',[1 1 20],[Inf Inf 5]); 
GPP_obs2(1:360,1:720,1:5) = NaN;
for i=1:5
    gpp_data2 = gpp_ECMWF(:,:,i)';
    GPP_obs2(:,:,i) = gpp_data2;
end
NPP_obs2(1:360,1:720,1:5) = NaN;
for i=1:5
    nppData2 = npp_ECMWF(:,:,i)';
    NPP_obs2(:,:,i) = nppData2;
end
CUE_obs2 = NPP_obs2./GPP_obs2;

% estimated based on MERRA2 climate data
gpp_MERRA2 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\GIMMSNPP\GIMMSNPP\Annual\Annual\MERRA2\gpp_V4_Standard_1982_2015_Annual_GEO_30min.nc','GPP',[1 1 20],[Inf Inf 5]);
npp_MERRA2 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\GIMMSNPP\GIMMSNPP\Annual\Annual\MERRA2\npp_V4_Standard_1982_2015_Annual_GEO_30min.nc','NPP',[1 1 20],[Inf Inf 5]); 
GPP_obs3(1:360,1:720,1:5) = NaN;
for i=1:5
    gpp_data3 = gpp_MERRA2(:,:,i)';
    GPP_obs3(:,:,i) = gpp_data3;
end
NPP_obs3(1:360,1:720,1:5) = NaN;
for i=1:5
    nppData3 = npp_MERRA2(:,:,i)';
    NPP_obs3(:,:,i) = nppData3;
end
CUE_obs3 = NPP_obs3./GPP_obs3;

% estimated based on NCEP climate data
gpp_NCEP = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\GIMMSNPP\GIMMSNPP\Annual\Annual\NCEPR2\Standard\gpp_V4_Standard_1982_2015_Annual_GEO_30min.nc','GPP',[1 1 20],[Inf Inf 5]); 
npp_NCEP = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\NPP\GIMMSNPP\GIMMSNPP\Annual\Annual\NCEPR2\Standard\npp_V4_Standard_1982_2015_Annual_GEO_30min.nc','NPP',[1 1 20],[Inf Inf 5]);
GPP_obs4(1:360,1:720,1:5) = NaN;
for i=1:5
    gpp_data4 = gpp_NCEP(:,:,i)';
    GPP_obs4(:,:,i) = gpp_data4;
end
NPP_obs4(1:360,1:720,1:5) = NaN;
for i=1:5
    nppData4 = npp_NCEP(:,:,i)';
    NPP_obs4(:,:,i) = nppData4;
end
CUE_obs4 = NPP_obs4./GPP_obs4;

% other two GPP observational data
% Data: FLUXCOM;  Unit: g C m-2 day-1
GPP_maps_1 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\FLUXCOM\yearlyGPP1980_2013.nc','GPP',[1 1 22],[Inf Inf 5]);
GPP_obs5(1:360,1:720,1:5) = NaN;
for i =1:5
    i
    GPPobs_M = GPP_maps_1(:,:,i);
    GPPobs_t = GPPobs_M';
    GPPobs_t = GPPobs_t.*365;  % convert unit to gC m-2 yr-1
    
    GPP_obs5(:,:,i) = GPPobs_t;   
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
    GPP_obs6(:,:,i-2000) = cell2mat(raw);
end

% merge NPP data
% unit: gC m-2 yr-1; 
NPP_obs05_5yr(:,:,:,1) = NPP_obs1;
NPP_obs05_5yr(:,:,:,2) = NPP_obs2;
NPP_obs05_5yr(:,:,:,3) = NPP_obs3;
NPP_obs05_5yr(:,:,:,4) = NPP_obs4;

% 2001-2005 mean 
NPP_obs05(:,:,1) = nanmean(NPP_obs1,3);
NPP_obs05(:,:,2) = nanmean(NPP_obs2,3);
NPP_obs05(:,:,3) = nanmean(NPP_obs3,3);
NPP_obs05(:,:,4) = nanmean(NPP_obs4,3);

% merge GPP data
GPP_obs05_5yr(:,:,:,1) = GPP_obs1;
GPP_obs05_5yr(:,:,:,2) = GPP_obs2;
GPP_obs05_5yr(:,:,:,3) = GPP_obs3;
GPP_obs05_5yr(:,:,:,4) = GPP_obs4;
GPP_obs05_5yr(:,:,:,5) = GPP_obs5;
GPP_obs05_5yr(:,:,:,6) = GPP_obs6;

% unit: gC m-2 yr-1; 
% 2001-2005 mean
GPP_obs05(:,:,1) = nanmean(GPP_obs1,3);
GPP_obs05(:,:,2) = nanmean(GPP_obs2,3);
GPP_obs05(:,:,3) = nanmean(GPP_obs3,3);
GPP_obs05(:,:,4) = nanmean(GPP_obs4,3);
GPP_obs05(:,:,5) = nanmean(GPP_obs5,3);
GPP_obs05(:,:,6) = nanmean(GPP_obs6,3);

% merge CUE data
% 2001-2005 mean
CUE_obs05(:,:,1) = nanmean(CUE_obs1,3);
CUE_obs05(:,:,2) = nanmean(CUE_obs2,3);
CUE_obs05(:,:,3) = nanmean(CUE_obs3,3);
CUE_obs05(:,:,4) = nanmean(CUE_obs4,3);
CUE_obs05(CUE_obs05>1) = NaN; % omit unreasonable values

NPP_obs05KG = NPP_obs05.*10^(-3); % change unit from gC m-2 yr-1 into KgC m-2 yr-1
GPP_obs05KG = GPP_obs05.*10^(-3); % change unit from gC m-2 yr-1 into KgC m-2 yr-1

clearvars -except NPP_obs05KG GPP_obs05KG CUE_obs05 NPP_obs05_5yr GPP_obs05_5yr
%% model simulations from CMIP5 and CMIP6
% CMIP5 NPP 
cd('G:\CMIP5\4_MatData\spatial_data\hist_rcp85')
load('NPP_cmip5_251.mat')
% load NPP simulations from 2001 to 2005
sp5_NPP5(:,:,:,1) = nppBCC_3D(:,:,152:156);
sp5_NPP5(:,:,:,2) = nppCAN_3D(:,:,152:156);
sp5_NPP5(:,:,:,3) = nppCCSM_3D(:,:,152:156);
sp5_NPP5(:,:,:,4) = nppHAD_3D(:,:,142:146);
sp5_NPP5(:,:,:,5) = nppIPSL_3D(:,:,152:156);
sp5_NPP5(:,:,:,6) = nppMIROC_3D(:,:,152:156);
sp5_NPP5(:,:,:,7) = nppMPI_3D(:,:,152:156);
sp5_NPP5(:,:,:,8) = nppNOR_3D(:,:,152:156);
sp5_NPP5(:,:,:,9) = nppBNU_3D(:,:,152:156);
sp5_NPP5(:,:,:,10) = nppGF_3D(:,:,141:145);
sp5_NPP5(:,:,:,11) = nppMRI_3D(:,:,151:155);
% calculate 2001-2005 mean
sp5_NPP5_11 = nanmean(sp5_NPP5,3);
sp5_NPP5_11 = squeeze(sp5_NPP5_11);
sp5_NPP5_11(:,1:30,:) = NaN; 
% unify the land region based on BCC
mask = sp5_NPP5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp5_NPP5_11 = sp5_NPP5_11.*mask;
% calculate model ensemble mean and SD
sp5_NPP5_avg = nanmean(sp5_NPP5_11,3);
sp5_NPP5_std = nanstd(sp5_NPP5_11,0,3);
% put ensemble mean and SD into the matrix
sp5_NPP5_11(:,:,12) = sp5_NPP5_avg ;
sp5_NPP5_11(:,:,13) = sp5_NPP5_std ;

% GPP
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
% calculate 2001-2005 mean
sp5_GPP5_11 = nanmean(sp5_GPP5,3);
sp5_GPP5_11 = squeeze(sp5_GPP5_11);
sp5_GPP5_11(:,1:30,:) = NaN; 
% unify the land region based on BCC
mask = sp5_GPP5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp5_GPP5_11 = sp5_GPP5_11.*mask;
% calculate model ensemble mean and SD
sp5_GPP5_avg = nanmean(sp5_GPP5_11,3);
sp5_GPP5_std = nanstd(sp5_GPP5_11,0,3);
% put ensemble mean and SD into the matrix
sp5_GPP5_11(:,:,12) = sp5_GPP5_avg ;
sp5_GPP5_11(:,:,13) = sp5_GPP5_std ;

% CMIP5 CUE 
load('CUE_cmip5_251.mat')
sp5_CUE5(:,:,:,1) = cueBCC_3D(:,:,152:156);
sp5_CUE5(:,:,:,2) = cueCAN_3D(:,:,152:156);
sp5_CUE5(:,:,:,3) = cueCCSM_3D(:,:,152:156);
sp5_CUE5(:,:,:,4) = cueHAD_3D(:,:,142:146);
sp5_CUE5(:,:,:,5) = cueIPSL_3D(:,:,152:156);
sp5_CUE5(:,:,:,6) = cueMIROC_3D(:,:,152:156);
sp5_CUE5(:,:,:,7) = cueMPI_3D(:,:,152:156);
sp5_CUE5(:,:,:,8) = cueNOR_3D(:,:,152:156);
sp5_CUE5(:,:,:,9) = cueBNU_3D(:,:,152:156);
sp5_CUE5(:,:,:,10) = cueGF_3D(:,:,141:145);
sp5_CUE5(:,:,:,11) = cueMRI_3D(:,:,151:155);
% calculate 2001-2005 mean
sp5_CUE5_11 = nanmean(sp5_CUE5,3);
sp5_CUE5_11 = squeeze(sp5_CUE5_11);
sp5_CUE5_11(:,1:30,:) = NaN; 
% unify the land region based on BCC
mask = sp5_CUE5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp5_CUE5_11 = sp5_CUE5_11.*mask;
% calculate model ensemble mean and SD
sp5_CUE5_avg = nanmean(sp5_CUE5_11,3);
sp5_CUE5_std = nanstd(sp5_CUE5_11,0,3);
% put ensemble mean and SD into the matrix
sp5_CUE5_11(:,:,12) = sp5_CUE5_avg ;
sp5_CUE5_11(:,:,13) = sp5_CUE5_std ;

% CMIP6 NPP 
cd('H:\CMIP56_Csink\4_MatData\spatial_data\hist_ssp585')
load('NPP_cmip6_251.mat')
% load NPP simulations over 2001-2005              
sp6_NPP5(:,:,:,1) = nppBCC_3D(:,:,152:156);
sp6_NPP5(:,:,:,2) = nppCAN_3D(:,:,152:156);
sp6_NPP5(:,:,:,3) = nppCESM_3D(:,:,152:156);
sp6_NPP5(:,:,:,4) = nppUK_3D(:,:,152:156);
sp6_NPP5(:,:,:,5) = nppIPSL_3D(:,:,152:156);
sp6_NPP5(:,:,:,6) = nppMIC_3D(:,:,152:156);
sp6_NPP5(:,:,:,7) = nppMPI_3D(:,:,152:156);
sp6_NPP5(:,:,:,8) = nppNOR_3D(:,:,152:156);
sp6_NPP5(:,:,:,9) = nppASS_3D(:,:,152:156);
sp6_NPP5(:,:,:,10) = nppCNRM_3D(:,:,152:156);
sp6_NPP5(:,:,:,11) = nppEC_3D(:,:,152:156);
% calculate 2001-2005 mean
sp6_NPP5_11 = nanmean(sp6_NPP5,3);
sp6_NPP5_11 = squeeze(sp6_NPP5_11);
sp6_NPP5_11(:,1:30) = NaN;
% unify the land area based on BCC
mask = sp6_NPP5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp6_NPP5_11 = sp6_NPP5_11.*mask;
% calculate ensemble mean and SD
sp6_NPP5_avg = nanmean(sp6_NPP5_11,3);
sp6_NPP5_SD = nanstd(sp6_NPP5_11,0,3);
% put mean and SD into the matrix
sp6_NPP5_11(:,:,12) = sp6_NPP5_avg ;
sp6_NPP5_11(:,:,13) = sp6_NPP5_SD;

% GPP
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
% calculate 2001-2005 mean
sp6_GPP5_11 = nanmean(sp6_GPP5,3);
sp6_GPP5_11 = squeeze(sp6_GPP5_11);
sp6_GPP5_11(:,1:30) = NaN;
% unify the land area based on BCC
mask = sp6_GPP5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp6_GPP5_11 = sp6_GPP5_11.*mask;
% calculate ensemble mean and SD
sp6_GPP5_avg = nanmean(sp6_GPP5_11,3);
sp6_GPP5_SD = nanstd(sp6_GPP5_11,0,3);
% put mean and SD into the matrix
sp6_GPP5_11(:,:,12) = sp6_GPP5_avg ;
sp6_GPP5_11(:,:,13) = sp6_GPP5_SD;

% CMIP6 CUE
load('CUE_cmip6_251.mat')
sp6_CUE5(:,:,:,1) = cueBCC_3D(:,:,152:156);
sp6_CUE5(:,:,:,2) = cueCAN_3D(:,:,152:156);
sp6_CUE5(:,:,:,3) = cueCESM_3D(:,:,152:156);
sp6_CUE5(:,:,:,4) = cueUK_3D(:,:,152:156);
sp6_CUE5(:,:,:,5) = cueIPSL_3D(:,:,152:156);
sp6_CUE5(:,:,:,6) = cueMIC_3D(:,:,152:156);
sp6_CUE5(:,:,:,7) = cueMPI_3D(:,:,152:156);
sp6_CUE5(:,:,:,8) = cueNOR_3D(:,:,152:156);
sp6_CUE5(:,:,:,9) = cueASS_3D(:,:,152:156);
sp6_CUE5(:,:,:,10) = cueCNRM_3D(:,:,152:156);
sp6_CUE5(:,:,:,11) = cueEC_3D(:,:,152:156);
% calculate 2001-2005 mean
sp6_CUE5_11 = nanmean(sp6_CUE5,3);
sp6_CUE5_11 = squeeze(sp6_CUE5_11);
sp6_CUE5_11(:,1:30) = NaN;
% unify the land area based on BCC
mask = sp6_CUE5_11(:,:,1);
mask(~isnan(mask)) = 1;
sp6_CUE5_11 = sp6_CUE5_11.*mask;
% calculate ensemble mean and SD
sp6_CUE5_avg = nanmean(sp6_CUE5_11,3);
sp6_CUE5_SD = nanstd(sp6_CUE5_11,0,3);
% put mean and SD into the matrix
sp6_CUE5_11(:,:,12) = sp6_CUE5_avg ;
sp6_CUE5_11(:,:,13) = sp6_CUE5_SD;

clearvars -except NPP_obs05KG GPP_obs05KG CUE_obs05 NPP_obs05_5yr GPP_obs05_5yr...
                  sp5_NPP5_11 sp5_GPP5_11 sp5_CUE5_11 sp6_NPP5_11 sp6_GPP5_11 sp6_CUE5_11 
% model outputs to match with the global map and further regrid into 0.5x0.5 
NPP5_map11(1:180,1:360,1:13) = NaN;
NPP6_map11(1:180,1:360,1:13) = NaN;
GPP5_map11(1:180,1:360,1:13) = NaN;
GPP6_map11(1:180,1:360,1:13) = NaN;
CUE5_map11(1:180,1:360,1:13) = NaN;
CUE6_map11(1:180,1:360,1:13) = NaN;
for i = 1:13
    i
    
    NPP5_M = sp5_NPP5_11(:,:,i);
    NPP6_M = sp6_NPP5_11(:,:,i);
    GPP5_M = sp5_GPP5_11(:,:,i);
    GPP6_M = sp6_GPP5_11(:,:,i);
    CUE5_M = sp5_CUE5_11(:,:,i);
    CUE6_M = sp6_CUE5_11(:,:,i);
    
    NPP5_t = NPP5_M';
    NPP6_t = NPP6_M';
    GPP5_t = GPP5_M';
    GPP6_t = GPP6_M';
    CUE5_t = CUE5_M';
    CUE6_t = CUE6_M';
    
    NPP5_LR(1:180,1:360) = NaN;
    NPP6_LR(1:180,1:360) = NaN;
    GPP5_LR(1:180,1:360) = NaN;
    GPP6_LR(1:180,1:360) = NaN;
    CUE5_LR(1:180,1:360) = NaN; 
    CUE6_LR(1:180,1:360) = NaN; 
     
    NPP5_LR(:,1:180) = NPP5_t(:,181:360);
    NPP5_LR(:,181:360) = NPP5_t(:,1:180);
    NPP6_LR(:,1:180) = NPP6_t(:,181:360);
    NPP6_LR(:,181:360) = NPP6_t(:,1:180);
    
    GPP5_LR(:,1:180) = GPP5_t(:,181:360);
    GPP5_LR(:,181:360) = GPP5_t(:,1:180);
    GPP6_LR(:,1:180) = GPP6_t(:,181:360);
    GPP6_LR(:,181:360) = GPP6_t(:,1:180);
    
    CUE5_LR(:,1:180) = CUE5_t(:,181:360);
    CUE5_LR(:,181:360) = CUE5_t(:,1:180);
    CUE6_LR(:,1:180) = CUE6_t(:,181:360);
    CUE6_LR(:,181:360) = CUE6_t(:,1:180);
    
    map_NPP5(1:180,1:360) = NaN; 
    map_NPP6(1:180,1:360) = NaN; 
    map_GPP5(1:180,1:360) = NaN; 
    map_GPP6(1:180,1:360) = NaN;
    map_CUE5(1:180,1:360) = NaN;
    map_CUE6(1:180,1:360) = NaN;
    
    for k=1:180
        map_NPP5(181-k,:) = NPP5_LR(k,:);
        map_NPP6(181-k,:) = NPP6_LR(k,:);
        map_GPP5(181-k,:) = GPP5_LR(k,:);
        map_GPP6(181-k,:) = GPP6_LR(k,:);
        map_CUE5(181-k,:) = CUE5_LR(k,:);
        map_CUE6(181-k,:) = CUE6_LR(k,:);
        
    end
    NPP5_map11(:,:,i) = map_NPP5;
    NPP6_map11(:,:,i) = map_NPP6;
    GPP5_map11(:,:,i) = map_GPP5;
    GPP6_map11(:,:,i) = map_GPP6;
    CUE5_map11(:,:,i) = map_CUE5;
    CUE6_map11(:,:,i) = map_CUE6;
end
% regrid into 0.5x0.5 degree
lat05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lat');
lon05 = ncread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\tuaE\tau_all.nc','lon');
lat1 = 89.5:-1:-89.5; lat1 = lat1';
lon1 = -179.5:179.5;  lon1 = lon1';
[x05,y05] = meshgrid(lon05,lat05);
[x1,y1] = meshgrid(lon1,lat1);
% regrided CMIP5 and CMIP6 data into 0.5x0.5 resolution
NPP5_sp05(1:360,1:720,1:13) = NaN;
NPP6_sp05(1:360,1:720,1:13) = NaN;
GPP5_sp05(1:360,1:720,1:13) = NaN;
GPP6_sp05(1:360,1:720,1:13) = NaN;
CUE5_sp05(1:360,1:720,1:13) = NaN;
CUE6_sp05(1:360,1:720,1:13) = NaN;
for i=1:13
    NPP5_M = interp2(x1,y1,NPP5_map11(:,:,i),x05,y05,'linear');
    NPP6_M = interp2(x1,y1,NPP6_map11(:,:,i),x05,y05,'linear');
    GPP5_M = interp2(x1,y1,GPP5_map11(:,:,i),x05,y05,'linear');
    GPP6_M = interp2(x1,y1,GPP6_map11(:,:,i),x05,y05,'linear');
    CUE5_M = interp2(x1,y1,CUE5_map11(:,:,i),x05,y05,'linear');
    CUE6_M = interp2(x1,y1,CUE6_map11(:,:,i),x05,y05,'linear');
    
    NPP5_sp05(:,:,i) = NPP5_M;
    NPP6_sp05(:,:,i) = NPP6_M;
    GPP5_sp05(:,:,i) = GPP5_M;
    GPP6_sp05(:,:,i) = GPP6_M;
    CUE5_sp05(:,:,i) = CUE5_M;
    CUE6_sp05(:,:,i) = CUE6_M;
end

clearvars -except NPP_obs05KG GPP_obs05KG CUE_obs05 NPP_obs05_5yr GPP_obs05_5yr...
                  NPP5_sp05 NPP6_sp05 GPP5_sp05 GPP6_sp05 CUE5_sp05 CUE6_sp05

%% Global estimates on GPP, NPP and CUE
% CMIP5
cd('G:\CMIP5\4_MatData\temporal_data\hist_rcp85')
load('NPP_cmip5_tmp.mat')
load('GPP_cmip5_tmp.mat')
load('CUE_cmip5_tmp.mat')
% For HadGEM2-ES,GFDL-ESM2G and MRI-ESM1, their historical simulations did not fully cover the period of 1850-2005.
% We thus use NaN to supplement years without data  
NaNgf = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]';
NPPgf_tmp = [NaNgf; NPPgf_tmp]; 
GPPgf_tmp = [NaNgf; GPPgf_tmp];
CUEgf_tmp = [NaNgf; CUEgf_tmp];

NaNhad = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]';
NPPhad_tmp = [NaNhad; NPPhad_tmp];
GPPhad_tmp = [NaNhad; GPPhad_tmp];
CUEhad_tmp = [NaNhad; CUEhad_tmp];

NPPmri_tmp = [NaN; NPPmri_tmp];
GPPmri_tmp = [NaN; GPPmri_tmp];
CUEmri_tmp = [NaN; CUEmri_tmp];

NPP5_all(1:5,1:11) = NaN;
NPP5_all(:,1) = NPPbcc_tmp(152:156);
NPP5_all(:,2) = NPPcan_tmp(152:156);
NPP5_all(:,3) = NPPccsm_tmp(152:156);
NPP5_all(:,4) = NPPhad_tmp(152:156);
NPP5_all(:,5) = NPPipsl_tmp(152:156);
NPP5_all(:,6) = NPPmiroc_tmp(152:156);
NPP5_all(:,7) = NPPmpi_tmp(152:156);
NPP5_all(:,8) = NPPnor_tmp(152:156);
NPP5_all(:,9) = NPPbnu_tmp(152:156);
NPP5_all(:,10) = NPPgf_tmp(152:156);
NPP5_all(:,11) = NPPmri_tmp(152:156);

GPP5_all(1:5,1:11) = NaN;
GPP5_all(:,1) = GPPbcc_tmp(152:156);
GPP5_all(:,2) = GPPcan_tmp(152:156);
GPP5_all(:,3) = GPPccsm_tmp(152:156);
GPP5_all(:,4) = GPPhad_tmp(152:156);
GPP5_all(:,5) = GPPipsl_tmp(152:156);
GPP5_all(:,6) = GPPmiroc_tmp(152:156);
GPP5_all(:,7) = GPPmpi_tmp(152:156);
GPP5_all(:,8) = GPPnor_tmp(152:156);
GPP5_all(:,9) = GPPbnu_tmp(152:156);
GPP5_all(:,10) = GPPgf_tmp(152:156);
GPP5_all(:,11) = GPPmri_tmp(152:156);

% CMIP6
cd('H:\CMIP56_Csink\4_MatData\temporal_data\hist_ssp585')
load('2NPP_cmip6_tmp.mat')
load('2GPP_cmip6_tmp.mat')
load('2CUE_cmip6_tmp.mat')  
leg6_str = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'}    
NPP6_all = [NPPbcc_tmp(152:156),NPPcan_tmp(152:156),NPPcesm_tmp(152:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   NPPuk_tmp(152:156),NPPipsl_tmp(152:156),NPPmic_tmp(152:156),...            % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   NPPmpi_tmp(152:156),NPPnor_tmp(152:156),...                              % MPI-ESM1-2-LR, NorESM2
   NPPass_tmp(152:156),NPPcnrm_tmp(152:156),NPPec_tmp(152:156)];              % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg
GPP6_all = [GPPbcc_tmp(152:156),GPPcan_tmp(152:156),GPPcesm_tmp(152:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   GPPuk_tmp(152:156),GPPipsl_tmp(152:156),GPPmic_tmp(152:156),...          % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   GPPmpi_tmp(152:156),GPPnor_tmp(152:156),...                           % MPI-ESM1-2-LR, NorESM2
   GPPass_tmp(152:156),GPPcnrm_tmp(152:156),GPPec_tmp(152:156)];            % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

clearvars -except NPP_obs05KG GPP_obs05KG CUE_obs05 NPP_obs05_5yr GPP_obs05_5yr...
                  NPP5_sp05 NPP6_sp05 GPP5_sp05 GPP6_sp05 CUE5_sp05 CUE6_sp05 ...
                  NPP5_all GPP5_all NPP6_all GPP6_all 

area05 = csvread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\CABLEGridAreaM2.csv');
NPP_obs05_Pg = NPP_obs05_5yr.* area05.*10^(-15);    % convert Unit into PgC yr-1
NPP_obs5yr_Pg = nansum(NPP_obs05_Pg,1);  NPP_obs5yr_Pg = nansum(NPP_obs5yr_Pg,2); 
NPP_obs5yr_Pg = squeeze(NPP_obs5yr_Pg)

GPP_obs05_Pg = GPP_obs05_5yr.* area05.*10^(-15);    % Unit: PgC yr-1
GPP_obs5yr_Pg = nansum(GPP_obs05_Pg,1);  GPP_obs5yr_Pg = nansum(GPP_obs5yr_Pg,2); 
GPP_obs5yr_Pg = squeeze(GPP_obs5yr_Pg)

%% Figure
% open figure window and set position
figure
set(gcf,'position',[100 100 940,531.2])

% Panel(a): model-data comparison on GPP and NPP
Panel_a = tight_subplot(1,1,[0 0],[0.1 0.3],[0.08 0.52])
hold on
scatter(GPP5_all(:),NPP5_all(:),120,'o',...
    'MarkerFaceColor',[0.30,0.75,0.93],'MarkerEdgeColor',[0.30,0.75,0.93],...
    'MarkerFaceAlpha',0.3)
scatter(GPP6_all(:),NPP6_all(:),120,'o',...
    'MarkerFaceColor',[0.97,0.18,0.67],'MarkerEdgeColor',[0.97,0.45,0.77],...
    'MarkerFaceAlpha',0.3)
set(gca, 'YLim',[40 100],'XLim',[80 260]);
% add the observational range
% Note that, GPP and NPP data from GIMMSGPP-NPP was calculated as ensemble mean 
GPP_global_Pg = nanmean(GPP_obs5yr_Pg,1)
GIMMS_GPP = nanmean(GPP_global_Pg(2:4))
% GPP range from the min to the max
GPP_obsMax = max([GPP_global_Pg([1,5,6]),GIMMS_GPP])
GPP_obsMin = min([GPP_global_Pg([1,5,6]),GIMMS_GPP])
% NPP
NPP_global_Pg = nanmean(NPP_obs5yr_Pg,1)
GIMMS_NPP = nanmean(NPP_global_Pg(2:4))
MODIS_NPP = NPP_global_Pg(1)
ITO_NPP = [50.6 68.4] % the ranges of observation-based estimates on NPP [mean+(-)SD] estimated from Ito et al.,2011, Table2, NPP over 2000s
% NPP range from the min to the max
NPP_obsMax = max([GIMMS_NPP MODIS_NPP ITO_NPP])
NPP_obsMin = min([GIMMS_NPP MODIS_NPP ITO_NPP])

% shading area for GPP 
y=[40 100 100 40];
x=[GPP_obsMin GPP_obsMin GPP_obsMax GPP_obsMax]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.3)
% shading area for NPP
x=[80 260 260 80];
y=[NPP_obsMin NPP_obsMin NPP_obsMax NPP_obsMax]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.6,0.6,0.6],'EdgeColor','none','FaceAlpha',0.3)
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
set(gca,'linewidth',1.2,'box','on')
set(gca,'Fontname','Arial','FontSize',12);
ylabel('NPP (PgC yr^-^1)','Fontname','Arial','FontSize',13)
xlabel('GPP (PgC yr^-^1)','Fontname','Arial','FontSize',13)
set(gca,'XLim',[80 260],'YLim',[40 100]);
leg56 = legend({'CMIP5','CMIP6','GPP Obs','NPP Obs'})
text(92,95,'(a)','FontName','Arial','FontSize',12)

% GPP density plot for simulations from CMIP5 and CMIP6
Panel_Agpp = tight_subplot(1,1,[0 0],[0.705 0.06],[0.08 0.52])
hold on
[fGPP_cmip5, xiGPP_cmip5] = ksdensity(GPP5_all(:));
[fGPP_cmip6, xiGPP_cmip6] = ksdensity(GPP6_all(:));
plot(xiGPP_cmip5,fGPP_cmip5,'LineWidth',1.8,'color',[0.30,0.75,0.93]);
plot(xiGPP_cmip6,fGPP_cmip6,'LineWidth',1.8,'color',[1.00,0.07,0.65]);
set(gca,'XLim',[80 260]);
axis off
% NPP density polt for simulations from CMIP5 and CMIP6
Panel_Anpp = tight_subplot(1,1,[0 0],[0.1 0.3],[0.4807 0.38])
hold on
[fNPP_cmip5, xiNPP_cmip5] = ksdensity(NPP5_all(:));
[fNPP_cmip6, xiNPP_cmip6] = ksdensity(NPP6_all(:));
plot(fNPP_cmip5,xiNPP_cmip5,'LineWidth',1.8,'color',[0.30,0.75,0.93]);
plot(fNPP_cmip6,xiNPP_cmip6,'LineWidth',1.8,'color',[1.00,0.07,0.65]);
set(gca,'YLim',[38 102]);
axis off

%% Panel (b) and (c): Model-data comparison on CUE
% observational range of CUE
CUE_range = NPP_obs5yr_Pg./GPP_obs5yr_Pg(:,1:4)
CUE_min = min(CUE_range(:))
CUE_max = max(CUE_range(:))

% CUE in CMIP5 and CMIP6
CUE5_5yr_range = NPP5_all./GPP5_all
CUE6_5yr_range = NPP6_all./GPP6_all
% 2001-2005 mean
CUE5_5yr_avg = nanmean(CUE5_5yr_range,1)
CUE6_5yr_avg = nanmean(CUE6_5yr_range,1)
CUE5_5yr_min = min(CUE5_5yr_range)
CUE6_5yr_min = min(CUE6_5yr_range)
CUE5_5yr_max = max(CUE5_5yr_range)
CUE6_5yr_max = max(CUE6_5yr_range)

% preparing data for plotting
ID_models = 1:11;
Markers = {'s','o','^','d','s','*','o','v','o','*','d'};
Model5 = {'BCC-CSM1-1m','BNU-ESM','CanESM2','CCSM4',...
        'GFDL-ESM2G',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','MRI-ESM1','NorESM1-M'}
Model6 = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'}  
mycolor5 = [255 0 0; 0 0 255; 153 51 255; 237 176 33; ...%BCC,BNU,CanESM2,CCSM4
            139 0 139;... %GFDL
            0 197 205;... %Had
            0 205 0;...%IPSL
            207 194 124;... %Miroc2
            255 99 71; 255 20 147;...%mpi, mri
            65 105 255]./255;  %NorM    
mycolor6 = [255 0 0; 153 51 255; 237 176 33;...    %BCC-CSM2-MR, CanESM5, CESM2
           0 197 205; 0 205 0; 207 194 124;...     %UKESM1-0-LL, IPSL-CM6A-LR,MIROC-ES2L
           255 99 71; 65 105 255;...               %MPI-ESM1-2-LR,NorESM-LM
            0 0 0; 158 131 149; 119 136 153]./255; %ACCESS-ESM1-5, CNRM-ESM2-1  EC-Earth3-Veg            
CUE5_plot = [ID_models', CUE5_5yr_avg', CUE5_5yr_min', CUE5_5yr_max']
CUE6_plot = [ID_models', CUE6_5yr_avg', CUE6_5yr_min', CUE6_5yr_max']
% sort data based on CUE 
CUE5_plot_sort = sortrows(CUE5_plot,2,'descend')   
CUE6_plot_sort = sortrows(CUE6_plot,2,'descend')

% Panel(b)
Panel_bc = tight_subplot(2,1,[0.17 0],[0.18 0.08],[0.68 0.02])
axes(Panel_bc(1))
hold on
X = 1:11
for i =1:11
   h5(i) = plot(X(i),CUE5_plot_sort(i,2),'Marker',Markers{CUE5_plot_sort(i,1)},...
       'MarkerEdgeColor', mycolor5(CUE5_plot_sort(i,1),:),...
       'MarkerSize',10,'LineStyle','none','LineWidth',2)    
end
set(gca, 'YLim',[0.3 0.6],'XLim',[0 12]);
x=[0 12 12 0];
y=[CUE_min CUE_min CUE_max CUE_max]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.3)
set(gca,'linewidth',1.2,'box','on')
set(gca,'Fontname','Arial','FontSize',10);
Ylabel = gca
set(Ylabel.YAxis,'FontSize',12)
set(gca,'XTickLabelRotation',40);
set(gca,'YTickLabelMode','auto');
ylabel('CUE','Fontname','Arial','FontSize',13)
xticks(1:11);
xticklabels({Model5{CUE5_plot_sort(1:11,1)}});
text(0.5,0.57,'(b) CMIP5','FontName','Arial','FontSize',12)

% Panel(c)
axes(Panel_bc(2))
hold on
X = 1:11
for i =1:11
   h6(i) = plot(X(i),CUE6_plot_sort(i,2),'Marker',Markers{CUE6_plot_sort(i,1)},...
       'MarkerEdgeColor', mycolor6(CUE6_plot_sort(i,1),:),...
       'MarkerSize',10,'LineStyle','none','LineWidth',2);   
end
set(gca, 'YLim',[0.3 0.6],'XLim',[0 12]);
x=[0 12 12 0];
y=[CUE_min CUE_min CUE_max CUE_max]
obs = fill(x,y,'g')
set(obs,'FaceColor',[0.47,0.67,0.19],'EdgeColor','none','FaceAlpha',0.3)
set(gca,'linewidth',1.2,'box','on')
set(gca,'Fontname','Arial','FontSize',10);
Ylabel = gca
set(Ylabel.YAxis,'FontSize',12)
set(gca,'XTickLabelRotation',40);
set(gca,'YTickLabelMode','auto');
ylabel('CUE','Fontname','Arial','FontSize',13)
xticks(1:11);
xticklabels({Model6{CUE6_plot_sort(1:11,1)}});
text(0.5,0.57,'(c) CMIP6','FontName','Arial','FontSize',12)





