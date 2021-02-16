%% Extended Data Fig. 5, Model bias and inter-model spread in NPP, GPP and CUE (Evaluating at the regional scale)
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
% 
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

% merge NPP data
% unit: gC m-2 yr-1; 
% 2001-2005 mean 
NPP_obs05(:,:,1) = nanmean(NPP_obs1,3);
NPP_obs05(:,:,2) = nanmean(NPP_obs2,3);
NPP_obs05(:,:,3) = nanmean(NPP_obs3,3);
NPP_obs05(:,:,4) = nanmean(NPP_obs4,3);

% merge GPP data
% unit: gC m-2 yr-1; 
% 2001-2005 mean
GPP_obs05(:,:,1) = nanmean(GPP_obs1,3);
GPP_obs05(:,:,2) = nanmean(GPP_obs2,3);
GPP_obs05(:,:,3) = nanmean(GPP_obs3,3);
GPP_obs05(:,:,4) = nanmean(GPP_obs4,3);

% merge CUE data
% 2001-2005 mean
CUE_obs05(:,:,1) = nanmean(CUE_obs1,3);
CUE_obs05(:,:,2) = nanmean(CUE_obs2,3);
CUE_obs05(:,:,3) = nanmean(CUE_obs3,3);
CUE_obs05(:,:,4) = nanmean(CUE_obs4,3);
CUE_obs05(CUE_obs05>1) = NaN; % omit unreasonable values

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

GPP_obs05(:,:,5) = nanmean(GPP_obs5,3);
GPP_obs05(:,:,6) = nanmean(GPP_obs6,3);

 % convert unit from gC m-2 yr-1 into KgC m-2 yr-1
NPP_obs05KG = NPP_obs05.*10^(-3);
GPP_obs05KG = GPP_obs05.*10^(-3); 

clearvars -except NPP_obs05KG GPP_obs05KG CUE_obs05
save('E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure4\Obs_NPP_GPP_CUE\nppGppCUE.mat')
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

clearvars -except NPP_obs05KG GPP_obs05KG CUE_obs05...
                  sp5_NPP5_11 sp5_GPP5_11 sp5_CUE5_11 sp6_NPP5_11 sp6_GPP5_11 sp6_CUE5_11 

% convert model outputs to match with the global map and further regrid into 0.5x0.5 
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
              
%to match with observational data
CUE_mean_obs = nanmean(CUE_obs05,3);
mask_obs = nanmean(CUE_obs05,3);
mask_obs(~isnan(mask_obs)) = 1;
CUE5_avg_M11 = CUE5_sp05(:,:,12).*mask_obs; 
CUE6_avg_M11 = CUE6_sp05(:,:,12).*mask_obs; 

% estimating model bias in CUE
bis5_CUE = (CUE5_avg_M11 - CUE_mean_obs)./ CUE_mean_obs; 
bis6_CUE = (CUE6_avg_M11 - CUE_mean_obs)./ CUE_mean_obs; 

%% Figure 4: panel (a) and (e), the spatial patterns of model bias in CUE
figure
set(gcf,'position',[100 100 699.2,480])
maps = tight_subplot(2,1,[-0.13 -0.14],[0.15 -0.05],[-0.04 0.50])
load E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure4_NPP_CUE\mycolor_nppBIAS3.mat   

% CMIP5
bis5_CUE(302:360,:) = [];
bis5_CUE = flipud(bis5_CUE);
raster_bis5_CUE = georasterref('RasterSize',size(bis5_CUE),'Latlim',[-60 90],'Lonlim',[-180 180]);
cmip5 = maps(1)
axes(cmip5)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.5)
framem('FLineWidth',1)
h = geoshow(bis5_CUE,raster_bis5_CUE, 'DisplayType','surface','Zdata',zeros(size(bis5_CUE)),'CData',bis5_CUE);
colormap(mycolor_nppBIAS3)
caxis([-0.6 0.6])
set(gca,'box','off')
setm(cmip5, 'FFaceColor', [1,1,1])
framem('FLineWidth',1,'FEdgeColor','k')
axis off
text(-2.284,2.0824, '(a) CMIP5','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);

% CMIP6
bis6_CUE(302:360,:) = [];
bis6_CUE = flipud(bis6_CUE);
raster_bis6_CUE = georasterref('RasterSize',size(bis6_CUE),'Latlim',[-60 90],'Lonlim',[-180 180]);
cmip6 = maps(2)
axes(cmip6)
hold on
axesm miller
setm(gca,'MapLatLimit',[-60 90])
framem('FLineWidth',1)
framem('off')
geoshow('landareas.shp','FaceColor','none','LineWidth',0.5)
framem('FLineWidth',1)
h = geoshow(bis6_CUE,raster_bis6_CUE, 'DisplayType','surface','Zdata',zeros(size(bis6_CUE)),'CData',bis6_CUE);
colormap(mycolor_nppBIAS3)
caxis([-0.6 0.6])
set(gca,'box','off')
setm(cmip6, 'FFaceColor', [1,1,1])
framem('FLineWidth',1,'FEdgeColor','k')
axis off
text(-2.284,2.0824, '(e) CMIP6','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11);

h1 = colorbar
h1.Location = 'southoutside'
h1.Position = [0.0192,0.203,0.421,0.027];
h1.FontName = 'Arial'
h1.FontSize = 10;
h1.Ticks = [-0.6,-0.4,-0.2,0,0.2,0.4,0.6]
tt = text(0.0131,-2.4799,'$$\frac{CUE_{mod} - CUE_{obs}}{CUE_{obs}}$$',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11)  
set(tt,'Interpreter','latex');

panel = tight_subplot(2,3,[0.041 0.01],[0.235 0.028],[0.5 0.02])
%% Panel (b) and (f): the latitudinal patterns of model bias in NPP 
% omit extremely large values
NPP5_sp05(NPP5_sp05>20) = NaN;
NPP6_sp05(NPP6_sp05>20) = NaN; 
NPP_obs05KG(NPP_obs05KG>20) = NaN;
% zonal mean for individual models in CMIP5 and CMIP6
NPP5_zonal_11 = nanmean(NPP5_sp05(:,:,1:11),2); NPP5_zonal_11 = squeeze(NPP5_zonal_11); NPP5_zonal_11(302:360,:) = [];
NPP6_zonal_11 = nanmean(NPP6_sp05(:,:,1:11),2); NPP6_zonal_11 = squeeze(NPP6_zonal_11); NPP6_zonal_11(302:360,:) = [];        
% ensemble mean for observational data 
NPP_zonal_obs = nanmean(NPP_obs05KG,2); NPP_zonal_obs = nanmean(NPP_zonal_obs,3); NPP_zonal_obs(302:360,:) = [];  
% calculate model bias
npp5_zonal_bis = (NPP5_zonal_11 - NPP_zonal_obs)./ NPP_zonal_obs;
npp6_zonal_bis = (NPP6_zonal_11 - NPP_zonal_obs)./ NPP_zonal_obs;
npp5_zonal_bis(npp5_zonal_bis < -2) = NaN;
npp5_zonal_bis(npp5_zonal_bis > 2) = NaN;
npp6_zonal_bis(npp6_zonal_bis < -2) = NaN;
npp6_zonal_bis(npp6_zonal_bis > 2) = NaN;

% panel for CMIP5    
axes(panel(1))
hold on
lat = 90:-0.5:-60;
NPP5_zonal_BISavg = nanmean(npp5_zonal_bis,2); % model ensemble mean
NPP5_zonal_BISsd = nanstd(npp5_zonal_bis,0,2); % SD
% add different color for positive and negative values
NPP5_bisAVG_post(1:length(NPP5_zonal_BISavg))= NaN; 
NPP5_bisAVG_neg(1:length(NPP5_zonal_BISavg))= NaN; 
NPP5_bisSD_post(1:length(NPP5_zonal_BISavg),1:2)= NaN; 
NPP5_bisSD_neg(1:length(NPP5_zonal_BISavg),1:2)= NaN; 

NPP5_zero_over(1:length(NPP5_zonal_BISavg),1:2)= 0; 
NPP5_zero_below(1:length(NPP5_zonal_BISavg),1:2)= 0; 

for i = 1:length(NPP5_zonal_BISavg)
    
    if NPP5_zonal_BISavg(i)>0
        
       NPP5_bisAVG_post(i) = NPP5_zonal_BISavg(i);
       NPP5_bisSD_post(i,2) = NPP5_zonal_BISsd(i);
       NPP5_bisSD_post(i,1) = NPP5_zonal_BISsd(i);
       
       if NPP5_zonal_BISavg(i) - NPP5_zonal_BISsd(i) < 0
           NPP5_bisSD_post(i,1) = NPP5_zonal_BISavg(i) - 0 ;
           NPP5_zero_below(i,1) = NPP5_zonal_BISavg(i) - NPP5_zonal_BISsd(i);
           
       end
       
    end
            
    if NPP5_zonal_BISavg(i)<= 0
        
       NPP5_bisAVG_neg(i) = NPP5_zonal_BISavg(i);    
       NPP5_bisSD_neg(i,2) = NPP5_zonal_BISsd(i); 
       NPP5_bisSD_neg(i,1) = NPP5_zonal_BISsd(i);
       
       if NPP5_zonal_BISavg(i) + NPP5_zonal_BISsd(i) > 0
           NPP5_bisSD_neg(i,2) = 0 - NPP5_zonal_BISavg(i);
           NPP5_zero_over(i,2) = NPP5_zonal_BISavg(i) + NPP5_zonal_BISsd(i);
       end
       
    end
        
end
NPP5_bisSD_post(isnan(NPP5_bisSD_post)) = 0;
NPP5_bisSD_neg(isnan(NPP5_bisSD_neg)) = 0;
NPP5_bisAVG_post(isnan(NPP5_bisAVG_post)) = 0;
NPP5_bisAVG_neg(isnan(NPP5_bisAVG_neg)) = 0;
X_0(1:length(NPP5_zonal_BISavg))= 0;
CMIP5_avg = boundedline(NPP5_bisAVG_post,lat,NPP5_bisSD_post,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.88,0.15,0.02],'nan','remove');
set(CMIP5_avg,'color','none');
CMIP5_avg = boundedline(X_0,lat,-NPP5_zero_below,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.02,0.52,0.44],'nan','remove');
set(CMIP5_avg,'color','none');
CMIP5_avg2 = boundedline(NPP5_bisAVG_neg,lat,NPP5_bisSD_neg,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.02,0.52,0.44],'nan','remove');
set(CMIP5_avg2,'color','none');
CMIP5_avg2 = boundedline(X_0,lat,NPP5_zero_over,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.88,0.15,0.02],'nan','remove');
set(CMIP5_avg2,'color','none');
CMIP5_avg = plot(NPP5_zonal_BISavg,lat,'LineWidth',0.8,'color',[0.4 0.4 0.4])
plot([0 0], [-60 90], 'k-','LineWidth',0.6)
plot([-2 2], [0 0], 'k--','LineWidth',1.5)
set(gca, 'YLim',[-60 90],'XLim',[-2 2]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
%set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xticks([-2 -1 0 1 2])
xticklabels([]);
text(-1.8,80, '(b)','FontName','Arial','FontSize',11)

% panel for CMIP6          
axes(panel(4))
hold on
lat = 90:-0.5:-60;
NPP6_zonal_BISavg = nanmean(npp6_zonal_bis,2); % model ensemble mean
NPP6_zonal_BISsd = nanstd(npp6_zonal_bis,0,2); % SD
% add different color for positive and negative values
NPP6_bisAVG_post(1:length(NPP6_zonal_BISavg))= NaN; 
NPP6_bisAVG_neg(1:length(NPP6_zonal_BISavg))= NaN; 
NPP6_bisSD_post(1:length(NPP6_zonal_BISavg),1:2)= NaN; 
NPP6_bisSD_neg(1:length(NPP6_zonal_BISavg),1:2)= NaN; 

NPP6_zero_over(1:length(NPP6_zonal_BISavg),1:2)= 0; 
NPP6_zero_below(1:length(NPP6_zonal_BISavg),1:2)= 0; 

for i = 1:length(NPP6_zonal_BISavg)
    
    if NPP6_zonal_BISavg(i)>0
        
       NPP6_bisAVG_post(i) = NPP6_zonal_BISavg(i);
       NPP6_bisSD_post(i,2) = NPP6_zonal_BISsd(i);
       NPP6_bisSD_post(i,1) = NPP6_zonal_BISsd(i);
       
       if NPP6_zonal_BISavg(i) - NPP6_zonal_BISsd(i) < 0
           NPP6_bisSD_post(i,1) = NPP6_zonal_BISavg(i) - 0 ;
           NPP6_zero_below(i,1) = NPP6_zonal_BISavg(i) - NPP6_zonal_BISsd(i);
           
       end
       
    end
            
    if NPP6_zonal_BISavg(i)<= 0
        
       NPP6_bisAVG_neg(i) = NPP6_zonal_BISavg(i);    
       NPP6_bisSD_neg(i,2) = NPP6_zonal_BISsd(i); 
       NPP6_bisSD_neg(i,1) = NPP6_zonal_BISsd(i);
       
       if NPP6_zonal_BISavg(i) + NPP6_zonal_BISsd(i) > 0
           NPP6_bisSD_neg(i,2) = 0 - NPP6_zonal_BISavg(i);
           NPP6_zero_over(i,2) = NPP6_zonal_BISavg(i) + NPP6_zonal_BISsd(i);
       end
       
    end
        
end
NPP6_bisSD_post(isnan(NPP6_bisSD_post)) = 0;
NPP6_bisSD_neg(isnan(NPP6_bisSD_neg)) = 0;
NPP6_bisAVG_post(isnan(NPP6_bisAVG_post)) = 0;
NPP6_bisAVG_neg(isnan(NPP6_bisAVG_neg)) = 0;
X_0(1:length(NPP6_zonal_BISavg))= 0;
CMIP6_avg = boundedline(NPP6_bisAVG_post,lat,NPP6_bisSD_post,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.88,0.15,0.02],'nan','remove');
set(CMIP6_avg,'color','none');
CMIP6_avg = boundedline(X_0,lat,-NPP6_zero_below,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.02,0.52,0.44],'nan','remove');
set(CMIP6_avg,'color','none');
CMIP6_avg2 = boundedline(NPP6_bisAVG_neg,lat,NPP6_bisSD_neg,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.02,0.52,0.44],'nan','remove');
set(CMIP6_avg2,'color','none');
CMIP6_avg2 = boundedline(X_0,lat,NPP6_zero_over,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.88,0.15,0.02],'nan','remove');
set(CMIP6_avg2,'color','none');
CMIP6_avg = plot(NPP6_zonal_BISavg,lat,'LineWidth',0.8,'color',[0.4 0.4 0.4])
plot([0 0], [-60 90], 'k-','LineWidth',0.6)
plot([-2 2], [0 0], 'k--','LineWidth',1.5)
set(gca, 'YLim',[-60 90],'XLim',[-2 2]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
%set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels({'60S','30S','0','30N','60N','90N'});
xticks([-2 -1 0 1 2])
text(-1.8,80, '(f)','FontName','Arial','FontSize',11)  
text(-0.5,-90, 'NPP','FontName','Arial','FontSize',11) 

%%  Panel (c) and (g): the latitudinal patterns of model bias in GPP 
% omit extremely large values
GPP5_sp05(GPP5_sp05>20) = NaN; 
GPP6_sp05(GPP6_sp05>20) = NaN; 
GPP_obs05KG(GPP_obs05KG>20) = NaN;
% zonal mean for individual models in CMIP5 and CMIP6
GPP5_zonal_11 = nanmean(GPP5_sp05(:,:,1:11),2); GPP5_zonal_11 = squeeze(GPP5_zonal_11); GPP5_zonal_11(302:360,:) = [];
GPP6_zonal_11 = nanmean(GPP6_sp05(:,:,1:11),2); GPP6_zonal_11 = squeeze(GPP6_zonal_11); GPP6_zonal_11(302:360,:) = [];        
% ensemble mean for observational data 
GPP_zonal_obs = nanmean(GPP_obs05KG,2); GPP_zonal_obs = nanmean(GPP_zonal_obs,3); GPP_zonal_obs(302:360,:) = [];  
% calculate model bias
gpp5_zonal_bis = (GPP5_zonal_11 - GPP_zonal_obs)./ GPP_zonal_obs;
gpp6_zonal_bis = (GPP6_zonal_11 - GPP_zonal_obs)./ GPP_zonal_obs;
gpp5_zonal_bis(gpp5_zonal_bis < -5) = NaN;
gpp5_zonal_bis(gpp5_zonal_bis > 5) = NaN;
gpp6_zonal_bis(gpp6_zonal_bis < -5) = NaN;
gpp6_zonal_bis(gpp6_zonal_bis > 5) = NaN;


% panel for CMIP5    
axes(panel(2))
hold on
lat = 90:-0.5:-60;
GPP5_zonal_BISavg = nanmean(gpp5_zonal_bis,2); % model ensemble mean
GPP5_zonal_BISsd = nanstd(gpp5_zonal_bis,0,2); % SD
% add different color for positive and negative values
GPP5_bisAVG_post(1:length(GPP5_zonal_BISavg))= NaN; 
GPP5_bisAVG_neg(1:length(GPP5_zonal_BISavg))= NaN; 
GPP5_bisSD_post(1:length(GPP5_zonal_BISavg),1:2)= NaN; 
GPP5_bisSD_neg(1:length(GPP5_zonal_BISavg),1:2)= NaN; 

GPP5_zero_over(1:length(GPP5_zonal_BISavg),1:2)= 0; 
GPP5_zero_below(1:length(GPP5_zonal_BISavg),1:2)= 0; 
for i = 1:length(GPP5_zonal_BISavg)
    
    if GPP5_zonal_BISavg(i)>0
        
       GPP5_bisAVG_post(i) = GPP5_zonal_BISavg(i);
       GPP5_bisSD_post(i,2) = GPP5_zonal_BISsd(i);
       GPP5_bisSD_post(i,1) = GPP5_zonal_BISsd(i);
       
       if GPP5_zonal_BISavg(i) - GPP5_zonal_BISsd(i) < 0
           GPP5_bisSD_post(i,1) = GPP5_zonal_BISavg(i) - 0 ;
           GPP5_zero_below(i,1) = GPP5_zonal_BISavg(i) - GPP5_zonal_BISsd(i);
           
       end
       
    end
            
    if GPP5_zonal_BISavg(i)<= 0
        
       GPP5_bisAVG_neg(i) = GPP5_zonal_BISavg(i);    
       GPP5_bisSD_neg(i,2) = GPP5_zonal_BISsd(i); 
       GPP5_bisSD_neg(i,1) = GPP5_zonal_BISsd(i);
       
       if GPP5_zonal_BISavg(i) + GPP5_zonal_BISsd(i) > 0
           GPP5_bisSD_neg(i,2) = 0 - GPP5_zonal_BISavg(i);
           GPP5_zero_over(i,2) = GPP5_zonal_BISavg(i) + GPP5_zonal_BISsd(i);
       end
       
    end
        
end
GPP5_bisSD_post(isnan(GPP5_bisSD_post)) = 0;
GPP5_bisSD_neg(isnan(GPP5_bisSD_neg)) = 0;
GPP5_bisAVG_post(isnan(GPP5_bisAVG_post)) = 0;
GPP5_bisAVG_neg(isnan(GPP5_bisAVG_neg)) = 0;
X_0(1:length(GPP5_zonal_BISavg))= 0;
CMIP5_avg = boundedline(GPP5_bisAVG_post,lat,GPP5_bisSD_post,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.88,0.15,0.02],'nan','remove');
set(CMIP5_avg,'color','none');
CMIP5_avg = boundedline(X_0,lat,-GPP5_zero_below,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.02,0.52,0.44],'nan','remove');
set(CMIP5_avg,'color','none');
CMIP5_avg2 = boundedline(GPP5_bisAVG_neg,lat,GPP5_bisSD_neg,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.02,0.52,0.44],'nan','remove');
set(CMIP5_avg2,'color','none');
CMIP5_avg2 = boundedline(X_0,lat,GPP5_zero_over,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.88,0.15,0.02],'nan','remove');
set(CMIP5_avg2,'color','none');
CMIP5_avg = plot(GPP5_zonal_BISavg,lat,'LineWidth',0.8,'color',[0.4 0.4 0.4])
plot([0 0], [-60 90], 'k-','LineWidth',0.6)
plot([-2 2], [0 0], 'k--','LineWidth',1.5)
set(gca, 'YLim',[-60 90],'XLim',[-2 2]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
%set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels([]);
xticks([-2 -1 0 1 2])
xticklabels([]);
text(-1.8,80, '(c)','FontName','Arial','FontSize',11)

% for CMIP6        
axes(panel(5))
hold on
lat = 90:-0.5:-60;
GPP6_zonal_BISavg = nanmean(gpp6_zonal_bis,2); % model ensemble mean
GPP6_zonal_BISsd = nanstd(gpp6_zonal_bis,0,2); % SD
% add different color for positive and negative values
GPP6_bisAVG_post(1:length(GPP6_zonal_BISavg))= NaN; 
GPP6_bisAVG_neg(1:length(GPP6_zonal_BISavg))= NaN; 
GPP6_bisSD_post(1:length(GPP6_zonal_BISavg),1:2)= NaN; 
GPP6_bisSD_neg(1:length(GPP6_zonal_BISavg),1:2)= NaN; 

GPP6_zero_over(1:length(GPP6_zonal_BISavg),1:2)= 0; 
GPP6_zero_below(1:length(GPP6_zonal_BISavg),1:2)= 0; 
for i = 1:length(GPP6_zonal_BISavg)
    
    if GPP6_zonal_BISavg(i)>0
        
       GPP6_bisAVG_post(i) = GPP6_zonal_BISavg(i);
       GPP6_bisSD_post(i,2) = GPP6_zonal_BISsd(i);
       GPP6_bisSD_post(i,1) = GPP6_zonal_BISsd(i);
       
       if GPP6_zonal_BISavg(i) - GPP6_zonal_BISsd(i) < 0
           GPP6_bisSD_post(i,1) = GPP6_zonal_BISavg(i) - 0 ;
           GPP6_zero_below(i,1) = GPP6_zonal_BISavg(i) - GPP6_zonal_BISsd(i);
           
       end
       
    end
            
    if GPP6_zonal_BISavg(i)<= 0
        
       GPP6_bisAVG_neg(i) = GPP6_zonal_BISavg(i);    
       GPP6_bisSD_neg(i,2) = GPP6_zonal_BISsd(i); 
       GPP6_bisSD_neg(i,1) = GPP6_zonal_BISsd(i);
       
       if GPP6_zonal_BISavg(i) + GPP6_zonal_BISsd(i) > 0
           GPP6_bisSD_neg(i,2) = 0 - GPP6_zonal_BISavg(i);
           GPP6_zero_over(i,2) = GPP6_zonal_BISavg(i) + GPP6_zonal_BISsd(i);
       end
       
    end
        
end
GPP6_bisSD_post(isnan(GPP6_bisSD_post)) = 0;
GPP6_bisSD_neg(isnan(GPP6_bisSD_neg)) = 0;
GPP6_bisAVG_post(isnan(GPP6_bisAVG_post)) = 0;
GPP6_bisAVG_neg(isnan(GPP6_bisAVG_neg)) = 0;
X_0(1:length(GPP6_zonal_BISavg))= 0;
CMIP6_avg = boundedline(GPP6_bisAVG_post,lat,GPP6_bisSD_post,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.88,0.15,0.02],'nan','remove');
set(CMIP6_avg,'color','none');
CMIP6_avg = boundedline(X_0,lat,-GPP6_zero_below,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.02,0.52,0.44],'nan','remove');
set(CMIP6_avg,'color','none');
CMIP6_avg2 = boundedline(GPP6_bisAVG_neg,lat,GPP6_bisSD_neg,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.02,0.52,0.44],'nan','remove');
set(CMIP6_avg2,'color','none');
CMIP6_avg2 = boundedline(X_0,lat,GPP6_zero_over,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.88,0.15,0.02],'nan','remove');
set(CMIP6_avg2,'color','none');
CMIP6_avg = plot(GPP6_zonal_BISavg,lat,'LineWidth',0.8,'color',[0.4 0.4 0.4])
plot([0 0], [-60 90], 'k-','LineWidth',0.6)
plot([-2 2], [0 0], 'k--','LineWidth',1.5)
set(gca, 'YLim',[-60 90],'XLim',[-2 2]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
%set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
xticks([ -1 0 1 2])
text(-1.8,80, '(g)','FontName','Arial','FontSize',11) 
text(-0.5, -90,'GPP','FontName','Arial','FontSize',11)


%% Panel (d) and (h): the latitudinal patterns of model bias in CUE
 % omit extremely large values
CUE5_sp05(CUE5_sp05>20) = NaN;
CUE6_sp05(CUE6_sp05>20) = NaN; 
CUE_obs05(CUE_obs05>20) = NaN;
% zonal mean for individual models in CMIP5 and CMIP6
CUE5_zonal_11 = nanmean(CUE5_sp05(:,:,1:11),2); CUE5_zonal_11 = squeeze(CUE5_zonal_11); CUE5_zonal_11(302:360,:) = [];
CUE6_zonal_11 = nanmean(CUE6_sp05(:,:,1:11),2); CUE6_zonal_11 = squeeze(CUE6_zonal_11); CUE6_zonal_11(302:360,:) = [];       
% ensemble mean for observational data 
CUE_zonal_obs = nanmean(CUE_obs05,2); CUE_zonal_obs = nanmean(CUE_zonal_obs,3); CUE_zonal_obs(302:360,:) = [];  
% calculate model bias
cue5_zonal_bis = (CUE5_zonal_11 - CUE_zonal_obs)./ CUE_zonal_obs;
cue6_zonal_bis = (CUE6_zonal_11 - CUE_zonal_obs)./ CUE_zonal_obs;
cue5_zonal_bis(cue5_zonal_bis < -5) = NaN;
cue5_zonal_bis(cue5_zonal_bis > 5) = NaN;
cue6_zonal_bis(cue6_zonal_bis < -5) = NaN;
cue6_zonal_bis(cue6_zonal_bis > 5) = NaN;

% panel for CMIP5    
axes(panel(3))
hold on
lat = 90:-0.5:-60;
CUE5_zonal_BISavg = nanmean(cue5_zonal_bis,2); % model ensemble mean
CUE5_zonal_BISsd = nanstd(cue5_zonal_bis,0,2); % SD
% add different color for positive and negative values
CUE5_bisAVG_post(1:length(CUE5_zonal_BISavg))= NaN; 
CUE5_bisAVG_neg(1:length(CUE5_zonal_BISavg))= NaN; 
CUE5_bisSD_post(1:length(CUE5_zonal_BISavg),1:2)= NaN; 
CUE5_bisSD_neg(1:length(CUE5_zonal_BISavg),1:2)= NaN; 

CUE5_zero_over(1:length(CUE5_zonal_BISavg),1:2)= 0; 
CUE5_zero_below(1:length(CUE5_zonal_BISavg),1:2)= 0; 
for i = 1:length(CUE5_zonal_BISavg)
    
    if CUE5_zonal_BISavg(i)>0
        
       CUE5_bisAVG_post(i) = CUE5_zonal_BISavg(i);
       CUE5_bisSD_post(i,2) = CUE5_zonal_BISsd(i);
       CUE5_bisSD_post(i,1) = CUE5_zonal_BISsd(i);
       
       if CUE5_zonal_BISavg(i) - CUE5_zonal_BISsd(i) < 0
           CUE5_bisSD_post(i,1) = CUE5_zonal_BISavg(i) - 0 ;
           CUE5_zero_below(i,1) = CUE5_zonal_BISavg(i) - CUE5_zonal_BISsd(i);
           
       end
       
    end
            
    if CUE5_zonal_BISavg(i)<= 0
        
       CUE5_bisAVG_neg(i) = CUE5_zonal_BISavg(i);    
       CUE5_bisSD_neg(i,2) = CUE5_zonal_BISsd(i); 
       CUE5_bisSD_neg(i,1) = CUE5_zonal_BISsd(i);
       
       if CUE5_zonal_BISavg(i) + CUE5_zonal_BISsd(i) > 0
           CUE5_bisSD_neg(i,2) = 0 - CUE5_zonal_BISavg(i);
           CUE5_zero_over(i,2) = CUE5_zonal_BISavg(i) + CUE5_zonal_BISsd(i);
       end
       
    end
        
end
CUE5_bisSD_post(isnan(CUE5_bisSD_post)) = 0;
CUE5_bisSD_neg(isnan(CUE5_bisSD_neg)) = 0;
CUE5_bisAVG_post(isnan(CUE5_bisAVG_post)) = 0;
CUE5_bisAVG_neg(isnan(CUE5_bisAVG_neg)) = 0;
X_0(1:length(CUE5_zonal_BISavg))= 0;

CMIP5_avg = boundedline(CUE5_bisAVG_post,lat,CUE5_bisSD_post,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.88,0.15,0.02],'nan','remove');
set(CMIP5_avg,'color','none');
CMIP5_avg = boundedline(X_0,lat,-CUE5_zero_below,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.02,0.52,0.44],'nan','remove');
set(CMIP5_avg,'color','none');
CMIP5_avg2 = boundedline(CUE5_bisAVG_neg,lat,CUE5_bisSD_neg,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.02,0.52,0.44],'nan','remove');
set(CMIP5_avg2,'color','none');
CMIP5_avg2 = boundedline(X_0,lat,CUE5_zero_over,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.88,0.15,0.02],'nan','remove');
set(CMIP5_avg2,'color','none');
CMIP5_avg = plot(CUE5_zonal_BISavg,lat,'LineWidth',0.8,'color',[0.4 0.4 0.4])
plot([0 0], [-60 90], 'k-','LineWidth',0.6)
plot([-2 2], [0 0], 'k--','LineWidth',1.5)
set(gca, 'YLim',[-60 90],'XLim',[-1 1]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels([]);
xticks([-2 -1 0 1 2])
xticklabels([]);
text(-0.85,80, '(d)','FontName','Arial','FontSize',11)

% panel for CMIP6           
axes(panel(6))
hold on
lat = 90:-0.5:-60;
CUE6_zonal_BISavg = nanmean(cue6_zonal_bis,2); % model ensemble mean
CUE6_zonal_BISsd = nanstd(cue6_zonal_bis,0,2); % SD
% add different color for positive and negative values
CUE6_bisAVG_post(1:length(CUE6_zonal_BISavg))= NaN; 
CUE6_bisAVG_neg(1:length(CUE6_zonal_BISavg))= NaN; 
CUE6_bisSD_post(1:length(CUE6_zonal_BISavg),1:2)= NaN; 
CUE6_bisSD_neg(1:length(CUE6_zonal_BISavg),1:2)= NaN; 

CUE6_zero_over(1:length(CUE6_zonal_BISavg),1:2)= 0; 
CUE6_zero_below(1:length(CUE6_zonal_BISavg),1:2)= 0; 
for i = 1:length(CUE6_zonal_BISavg)
    
    if CUE6_zonal_BISavg(i)>0
        
       CUE6_bisAVG_post(i) = CUE6_zonal_BISavg(i);
       CUE6_bisSD_post(i,2) = CUE6_zonal_BISsd(i);
       CUE6_bisSD_post(i,1) = CUE6_zonal_BISsd(i);
       
       if CUE6_zonal_BISavg(i) - CUE6_zonal_BISsd(i) < 0
           CUE6_bisSD_post(i,1) = CUE6_zonal_BISavg(i) - 0 ;
           CUE6_zero_below(i,1) = CUE6_zonal_BISavg(i) - CUE6_zonal_BISsd(i);
           
       end
       
    end
            
    if CUE6_zonal_BISavg(i)<= 0
        
       CUE6_bisAVG_neg(i) = CUE6_zonal_BISavg(i);    
       CUE6_bisSD_neg(i,2) = CUE6_zonal_BISsd(i); 
       CUE6_bisSD_neg(i,1) = CUE6_zonal_BISsd(i);
       
       if CUE6_zonal_BISavg(i) + CUE6_zonal_BISsd(i) > 0
           CUE6_bisSD_neg(i,2) = 0 - CUE6_zonal_BISavg(i);
           CUE6_zero_over(i,2) = CUE6_zonal_BISavg(i) + CUE6_zonal_BISsd(i);
       end
       
    end
        
end
CUE6_bisSD_post(isnan(CUE6_bisSD_post)) = 0;
CUE6_bisSD_neg(isnan(CUE6_bisSD_neg)) = 0;
CUE6_bisAVG_post(isnan(CUE6_bisAVG_post)) = 0;
CUE6_bisAVG_neg(isnan(CUE6_bisAVG_neg)) = 0;
X_0(1:length(CUE6_zonal_BISavg))= 0;
CMIP6_avg = boundedline(CUE6_bisAVG_post,lat,CUE6_bisSD_post,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.88,0.15,0.02],'nan','remove');
set(CMIP6_avg,'color','none');
CMIP6_avg = boundedline(X_0,lat,-CUE6_zero_below,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.02,0.52,0.44],'nan','remove');
set(CMIP6_avg,'color','none');
CMIP6_avg2 = boundedline(CUE6_bisAVG_neg,lat,CUE6_bisSD_neg,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.02,0.52,0.44],'nan','remove');
set(CMIP6_avg2,'color','none');
CMIP6_avg2 = boundedline(X_0,lat,CUE6_zero_over,'alpha','transparency',0.5,...
'orientation', 'horiz','cmap',[0.88,0.15,0.02],'nan','remove');
set(CMIP6_avg2,'color','none');
CMIP6_avg = plot(CUE6_zonal_BISavg,lat,'LineWidth',0.8,'color',[0.4 0.4 0.4])
plot([0 0], [-60 90], 'k-','LineWidth',0.6)
plot([-2 2], [0 0], 'k--','LineWidth',1.5)
set(gca, 'YLim',[-60 90],'XLim',[-1 1]);
set(gca,'Fontname','Arial','FontSize',10);
set(gca,'lineWidth',1.2,'box', 'on')
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');
yticks([-60 -30 0 30 60 90]);
yticklabels([]);
xticks([ 0 1 ])
text(-0.85,80, '(h)','FontName','Arial','FontSize',11)   
text(-0.3,-90, 'CUE','FontName','Arial','FontSize',11)   




