% global cVeg and cSoil simulated by each model was calculated from native resolution
% CMIP5
clear;clc;
file = 'G:\CMIP5\2_DataAnalysis\3_originalYR\hist_rcp85\'
sftlf_BCC = ncread([file,'bcc-csm1-1m\sftlf_fx_bcc-csm1-1-m_r0i0p0.nc'],'sftlf');
area_BCC = ncread([file,'bcc-csm1-1m\areacella_fx_bcc-csm1-1-m_r0i0p0.nc'],'areacella');
cVeg_yr_BCC = ncread([file,'bcc-csm1-1m\cVeg_yr_bcc-csm1-1-m_historical_r1i1p1_185001-210012.nc'],'cVeg');
cLit_yr_BCC = ncread([file,'bcc-csm1-1m\clitter_yr_bcc-csm1-1-m_historical_r1i1p1_185001-210012.nc'],'cLitter');
cSOM_yr_BCC = ncread([file,'bcc-csm1-1m\cSoil_yr_bcc-csm1-1-m_historical_r1i1p1_185001-210012.nc'],'cSoil');
cSoil_yr_BCC = cLit_yr_BCC + cSOM_yr_BCC;

sftlf_BNU = ncread([file,'BNU-ESM\sftlf_fx_BNU-ESM_r0i0p0.nc'],'sftlf');
area_BNU = ncread([file,'BNU-ESM\areacella_fx_BNU-ESM_r0i0p0.nc'],'areacella');
cVeg_yr_BNU = ncread([file,'BNU-ESM\cVeg_yr_BNU-ESM_r1i1p1_185001-210012.nc'],'cVeg');
cLit_yr_BNU = ncread([file,'BNU-ESM\cLitter_yr_BNU-ESM_r1i1p1_185001-210012.nc'],'cLitter');
cSOM_yr_BNU = ncread([file,'BNU-ESM\cSoil_yr_BNU-ESM_r1i1p1_185001-210012.nc'],'cSoil');
cSoil_yr_BNU = cLit_yr_BNU + cSOM_yr_BNU;

sftlf_CAN = ncread([file,'CanESM2\sftlf_fx_CanESM2_r0i0p0.nc'],'sftlf');
area_CAN = ncread([file,'CanESM2\areacella_fx_CanESM2_r0i0p0.nc'],'areacella');
cVeg_yr_CAN = ncread([file,'CanESM2\cVeg_yr_CanESM2_r1i1p1_185001-210012.nc'],'cVeg');
cLit_yr_CAN = ncread([file,'CanESM2\cLitter_yr_CanESM2_r1i1p1_185001-210012.nc'],'cLitter');
cSOM_yr_CAN = ncread([file,'CanESM2\cSoil_yr_CanESM2_r1i1p1_185001-210012.nc'],'cSoil');
cSoil_yr_CAN = cLit_yr_CAN + cSOM_yr_CAN;

sftlf_CCSM = ncread([file,'CCSM4\sftlf_fx_CCSM4_historical_r0i0p0.nc'],'sftlf');
area_CCSM = ncread([file,'CCSM4\areacella_fx_CCSM4_historical_r0i0p0.nc'],'areacella');
cVeg_yr_CCSM = ncread([file,'CCSM4\cVeg_yr_CCSM4_r1i1p1_185001-210012.nc'],'cVeg');
cLit_yr_CCSM = ncread([file,'CCSM4\cLitter_yr_CCSM4_r1i1p1_185001-210012.nc'],'cLitter');
cSOM_yr_CCSM = ncread([file,'CCSM4\cSoil_yr_CCSM4_r1i1p1_185001-210012.nc'],'cSoil');
cSoil_yr_CCSM = cLit_yr_CCSM + cSOM_yr_CCSM;

sftlf_GF = ncread([file,'GFDL-ESM2G\sftlf_fx_GFDL-ESM2G_r0i0p0.nc'],'sftlf');
area_GF = ncread([file,'GFDL-ESM2G\areacella_fx_GFDL-ESM2G_r0i0p0.nc'],'areacella');
cVeg_yr_GF = ncread([file,'GFDL-ESM2G\cVeg_yr_GFDL-ESM2G_r1i1p1_186101-210012.nc'],'cVeg');
cSOM_yr_GF = ncread([file,'GFDL-ESM2G\cSoil_yr_GFDL-ESM2G_r1i1p1_186101-210012.nc'],'cSoil');
cSoil_yr_GF =  cSOM_yr_GF;

sftlf_HAD = ncread([file,'HadGEM2-ES\sftlf_fx_HadGEM2-ES_r0i0p0.nc'],'sftlf');
area_HAD = ncread([file,'HadGEM2-ES\areacella_fx_HadGEM2-ES_r0i0p0.nc'],'areacella');
cVeg_yr_HAD = ncread([file,'HadGEM2-ES\cVeg_yr_HadGEM2-ES_r1i1p1_186001-210012.nc'],'cVeg');
cSOM_yr_HAD = ncread([file,'HadGEM2-ES\cSoil_yr_HadGEM2-ES_r1i1p1_186001-210012.nc'],'cSoil');
cSoil_yr_HAD = cSOM_yr_HAD;

sftlf_IPSL = ncread([file,'IPSL-CM5A-MR\sftlf_fx_IPSL-CM5A-MR_historical_r0i0p0.nc'],'sftlf');
area_IPSL = ncread([file,'IPSL-CM5A-MR\areacella_fx_IPSL-CM5A-MR_historical_r0i0p0.nc'],'areacella');
cVeg_yr_IPSL = ncread([file,'IPSL-CM5A-MR\cVeg_yr_IPSL-CM5A-MR_r1i1p1_185001-210012.nc'],'cVeg');
cLit_yr_IPSL = ncread([file,'IPSL-CM5A-MR\cLitter_yr_IPSL-CM5A-MR_r1i1p1_185001-210012.nc'],'cLitter');
cSOM_yr_IPSL = ncread([file,'IPSL-CM5A-MR\cSoil_yr_IPSL-CM5A-MR_r1i1p1_185001-210012.nc'],'cSoil');
cSoil_yr_IPSL = cLit_yr_IPSL + cSOM_yr_IPSL;

sftlf_MIROC = ncread([file,'MIROC-ESM\sftlf_fx_MIROC-ESM_r0i0p0.nc'],'sftlf');
area_MIROC = ncread([file,'MIROC-ESM\areacella_fx_MIROC-ESM_r0i0p0.nc'],'areacella');
cVeg_yr_MIROC = ncread([file,'MIROC-ESM\cVeg_yr_MIROC-ESM_r1i1p1_185001-210012.nc'],'cVeg');
cLit_yr_MIROC = ncread([file,'MIROC-ESM\cLitter_yr_MIROC-ESM_r1i1p1_185001-210012.nc'],'cLitter');
cSOM_yr_MIROC = ncread([file,'MIROC-ESM\cSoil_yr_MIROC-ESM_r1i1p1_185001-210012.nc'],'cSoil');
cSoil_yr_MIROC = cLit_yr_MIROC + cSOM_yr_MIROC;

sftlf_MPI = ncread([file,'MPI-ESM-MR\sftlf_fx_MPI-ESM-MR_historical_r0i0p0.nc'],'sftlf');
area_MPI = ncread([file,'MPI-ESM-MR\areacella_fx_MPI-ESM-MR_historical_r0i0p0.nc'],'areacella');
cVeg_yr_MPI = ncread([file,'MPI-ESM-MR\cVeg_yr_MPI-ESM-MR_r1i1p1_185001-210012.nc'],'cVeg');
cLit_yr_MPI = ncread([file,'MPI-ESM-MR\cLitter_yr_MPI-ESM-MR_r1i1p1_185001-210012.nc'],'cLitter');
cSOM_yr_MPI = ncread([file,'MPI-ESM-MR\cSoil_yr_MPI-ESM-MR_r1i1p1_185001-210012.nc'],'cSoil');
cSoil_yr_MPI = cLit_yr_MPI + cSOM_yr_MPI;

sftlf_MRI = ncread([file,'MRI-ESM1\sftlf_fx_MRI-ESM1_historical_r0i0p0.nc'],'sftlf');
area_MRI = ncread([file,'MRI-ESM1\areacella_fx_MRI-ESM1_historical_r0i0p0.nc'],'areacella');
cVeg_yr_MRI = ncread([file,'MRI-ESM1\cVeg_yr_MRI-ESM1_r1i1p1_185101-210012.nc'],'cVeg');
cLit_yr_MRI = ncread([file,'MRI-ESM1\cLitter_yr_MRI-ESM1_r1i1p1_185101-210012.nc'],'cLitter');
cSOM_yr_MRI = ncread([file,'MRI-ESM1\cSoil_yr_MRI-ESM1_r1i1p1_185101-210012.nc'],'cSoil');
cSoil_yr_MRI = cLit_yr_MRI + cSOM_yr_MRI;

sftlf_NOR = ncread([file,'NorESM1-M\sftlf_fx_NorESM1-M_historical_r0i0p0.nc'],'sftlf');
area_NOR = ncread([file,'NorESM1-M\areacella_fx_NorESM1-M_historical_r0i0p0.nc'],'areacella');
cVeg_yr_NOR = ncread([file,'NorESM1-M\2_yr_cVeg_mon_NorESM1-M_r1i1p1_185001-210012.nc'],'cVeg');
cLit_yr_NOR = ncread([file,'NorESM1-M\2_yr_cLitter_mon_NorESM1-M_r1i1p1_185001-210012.nc'],'cLitter');
cSOM_yr_NOR = ncread([file,'NorESM1-M\2_yr_cSoil_mon_NorESM1-M_r1i1p1_185001-210012.nc'],'cSoil');
cSoil_yr_NOR = cLit_yr_NOR + cSOM_yr_NOR;

cVeg_yrPG_BCC = cVeg_yr_BCC.*sftlf_BCC.*area_BCC.*10^(-2).*10^(-12); cVeg_yrPG_BCC(cVeg_yrPG_BCC>10^(10)) = NaN;
cVeg_yrPG_BCC = nansum(cVeg_yrPG_BCC,1); cVeg_yrPG_BCC=nansum(cVeg_yrPG_BCC,2); cVeg_yrPG_BCC = squeeze(cVeg_yrPG_BCC);

cVeg_yrPG_BNU = cVeg_yr_BNU.*sftlf_BNU.*area_BNU.*10^(-2).*10^(-12); cVeg_yrPG_BNU(cVeg_yrPG_BNU>10^(10)) = NaN;
cVeg_yrPG_BNU = nansum(cVeg_yrPG_BNU,1); cVeg_yrPG_BNU = nansum(cVeg_yrPG_BNU,2); cVeg_yrPG_BNU = squeeze(cVeg_yrPG_BNU); 

cVeg_yrPG_CAN = cVeg_yr_CAN.*sftlf_CAN.*area_CAN.*10^(-2).*10^(-12);cVeg_yrPG_CAN(cVeg_yrPG_CAN>10^(10)) = NaN;
cVeg_yrPG_CAN = nansum(cVeg_yrPG_CAN,1); cVeg_yrPG_CAN = nansum(cVeg_yrPG_CAN,2); cVeg_yrPG_CAN = squeeze(cVeg_yrPG_CAN);

cVeg_yrPG_CCSM = cVeg_yr_CCSM.*sftlf_CCSM.*area_CCSM.*10^(-2).*10^(-12); cVeg_yrPG_CCSM(cVeg_yrPG_CCSM>10^(10)) = NaN;
cVeg_yrPG_CCSM = nansum(cVeg_yrPG_CCSM,1);cVeg_yrPG_CCSM = nansum(cVeg_yrPG_CCSM,2); cVeg_yrPG_CCSM = squeeze(cVeg_yrPG_CCSM);

cVeg_yrPG_GF = cVeg_yr_GF.*sftlf_GF.*area_GF.*10^(-2).*10^(-12); cVeg_yrPG_GF(cVeg_yrPG_GF>10^(10)) = NaN;
cVeg_yrPG_GF = nansum(cVeg_yrPG_GF,1); cVeg_yrPG_GF = nansum(cVeg_yrPG_GF,2); cVeg_yrPG_GF = squeeze(cVeg_yrPG_GF);
GF_NaN(1:11,1) = NaN;
cVeg_yrPG_GF = cat(1,GF_NaN,cVeg_yrPG_GF);

cVeg_yrPG_HAD = cVeg_yr_HAD.*sftlf_HAD.*area_HAD.*10^(-2).*10^(-12);cVeg_yrPG_HAD(cVeg_yrPG_HAD>10^(10))=NaN;
cVeg_yrPG_HAD = nansum(cVeg_yrPG_HAD,1); cVeg_yrPG_HAD = nansum(cVeg_yrPG_HAD,2); cVeg_yrPG_HAD = squeeze(cVeg_yrPG_HAD);
HAD_NaN(1:10,1) = NaN;
cVeg_yrPG_HAD = cat(1,HAD_NaN,cVeg_yrPG_HAD);

cVeg_yrPG_IPSL = cVeg_yr_IPSL.*sftlf_IPSL.*area_IPSL.*10^(-2).*10^(-12); cVeg_yrPG_IPSL(cVeg_yrPG_IPSL>10^(10))=NaN;
cVeg_yrPG_IPSL = nansum(cVeg_yrPG_IPSL,1); cVeg_yrPG_IPSL = nansum(cVeg_yrPG_IPSL,2); cVeg_yrPG_IPSL = squeeze(cVeg_yrPG_IPSL);

cVeg_yrPG_MIROC = cVeg_yr_MIROC.*sftlf_MIROC.*area_MIROC.*10^(-2).*10^(-12); cVeg_yrPG_MIROC(cVeg_yrPG_MIROC>10^(10))=NaN;
cVeg_yrPG_MIROC = nansum(cVeg_yrPG_MIROC,1); cVeg_yrPG_MIROC = nansum(cVeg_yrPG_MIROC,2); cVeg_yrPG_MIROC=squeeze(cVeg_yrPG_MIROC);

cVeg_yrPG_MPI = cVeg_yr_MPI.*sftlf_MPI.*area_MPI.*10^(-2).*10^(-12); cVeg_yrPG_MPI(cVeg_yrPG_MPI>10^(10))=NaN;
cVeg_yrPG_MPI = nansum(cVeg_yrPG_MPI,1); cVeg_yrPG_MPI = nansum(cVeg_yrPG_MPI,2); cVeg_yrPG_MPI = squeeze(cVeg_yrPG_MPI);

cVeg_yrPG_MRI = cVeg_yr_MRI.*sftlf_MRI.*area_MRI.*10^(-2).*10^(-12); cVeg_yrPG_MRI(cVeg_yrPG_MRI>10^(10))=NaN;
cVeg_yrPG_MRI = nansum(cVeg_yrPG_MRI,1); cVeg_yrPG_MRI = nansum(cVeg_yrPG_MRI,2); cVeg_yrPG_MRI = squeeze(cVeg_yrPG_MRI);
cVeg_yrPG_MRI = cat(1,NaN,cVeg_yrPG_MRI);

cVeg_yrPG_NOR = cVeg_yr_NOR.*sftlf_NOR.*area_NOR.*10^(-2).*10^(-12); cVeg_yrPG_NOR(cVeg_yrPG_NOR>10^(10))=NaN;
cVeg_yrPG_NOR = nansum(cVeg_yrPG_NOR,1); cVeg_yrPG_NOR = nansum(cVeg_yrPG_NOR,2); cVeg_yrPG_NOR = squeeze(cVeg_yrPG_NOR);

cVeg5_org_PG(1:5,1:11)=NaN;
cVeg5_org_PG(:,1) = cVeg_yrPG_BCC(152:156); 
cVeg5_org_PG(:,2) = cVeg_yrPG_CAN(152:156);  
cVeg5_org_PG(:,3) = cVeg_yrPG_CCSM(152:156); 
cVeg5_org_PG(:,4) = cVeg_yrPG_HAD(152:156);  
cVeg5_org_PG(:,5) = cVeg_yrPG_IPSL(152:156);  
cVeg5_org_PG(:,6) = cVeg_yrPG_MIROC(152:156);  
cVeg5_org_PG(:,7) = cVeg_yrPG_MPI(152:156);  
cVeg5_org_PG(:,8) = cVeg_yrPG_NOR(152:156);  
cVeg5_org_PG(:,9) = cVeg_yrPG_BNU(152:156);  
cVeg5_org_PG(:,10) = cVeg_yrPG_GF(152:156);  
cVeg5_org_PG(:,11) = cVeg_yrPG_MRI(152:156);  

cSoil_yrPG_BCC = cSoil_yr_BCC.*sftlf_BCC.*area_BCC.*10^(-2).*10^(-12);cSoil_yrPG_BCC(cSoil_yrPG_BCC>10^(10))=NaN;
cSoil_yrPG_BCC = nansum(cSoil_yrPG_BCC,1); cSoil_yrPG_BCC = nansum(cSoil_yrPG_BCC,2); cSoil_yrPG_BCC = squeeze(cSoil_yrPG_BCC);

cSoil_yrPG_BNU = cSoil_yr_BNU.*sftlf_BNU.*area_BNU.*10^(-2).*10^(-12);cSoil_yrPG_BNU(cSoil_yrPG_BNU>10^(10))=NaN;
cSoil_yrPG_BNU = nansum(cSoil_yrPG_BNU,1); cSoil_yrPG_BNU = nansum(cSoil_yrPG_BNU,2); cSoil_yrPG_BNU = squeeze(cSoil_yrPG_BNU);

cSoil_yrPG_CAN = cSoil_yr_CAN.*sftlf_CAN.*area_CAN.*10^(-2).*10^(-12);cSoil_yrPG_CAN(cSoil_yrPG_CAN>10^(10))=NaN;
cSoil_yrPG_CAN = nansum(cSoil_yrPG_CAN,1); cSoil_yrPG_CAN = nansum(cSoil_yrPG_CAN,2); cSoil_yrPG_CAN = squeeze(cSoil_yrPG_CAN);

cSoil_yrPG_CCSM = cSoil_yr_CCSM.*sftlf_CCSM.*area_CCSM.*10^(-2).*10^(-12);cSoil_yrPG_CCSM(cSoil_yrPG_CCSM>10^(10))=NaN;
cSoil_yrPG_CCSM = nansum(cSoil_yrPG_CCSM,1); cSoil_yrPG_CCSM = nansum(cSoil_yrPG_CCSM,2); cSoil_yrPG_CCSM = squeeze(cSoil_yrPG_CCSM);

cSoil_yrPG_GF = cSoil_yr_GF.*sftlf_GF.*area_GF.*10^(-2).*10^(-12);cSoil_yrPG_GF(cSoil_yrPG_GF>10^(10))=NaN;
cSoil_yrPG_GF = nansum(cSoil_yrPG_GF,1); cSoil_yrPG_GF = nansum(cSoil_yrPG_GF,2); cSoil_yrPG_GF = squeeze(cSoil_yrPG_GF);
cSoil_yrPG_GF = cat(1,GF_NaN,cSoil_yrPG_GF);

cSoil_yrPG_HAD = cSoil_yr_HAD.*sftlf_HAD.*area_HAD.*10^(-2).*10^(-12);cSoil_yrPG_HAD(cSoil_yrPG_HAD>10^(10))=NaN;
cSoil_yrPG_HAD = nansum(cSoil_yrPG_HAD,1); cSoil_yrPG_HAD = nansum(cSoil_yrPG_HAD,2); cSoil_yrPG_HAD = squeeze(cSoil_yrPG_HAD);
cSoil_yrPG_HAD = cat(1,HAD_NaN,cSoil_yrPG_HAD);

cSoil_yrPG_IPSL = cSoil_yr_IPSL.*sftlf_IPSL.*area_IPSL.*10^(-2).*10^(-12);cSoil_yrPG_IPSL(cSoil_yrPG_IPSL>10^(10))=NaN;
cSoil_yrPG_IPSL = nansum(cSoil_yrPG_IPSL,1); cSoil_yrPG_IPSL = nansum(cSoil_yrPG_IPSL,2); cSoil_yrPG_IPSL = squeeze(cSoil_yrPG_IPSL);

cSoil_yrPG_MIROC = cSoil_yr_MIROC.*sftlf_MIROC.*area_MIROC.*10^(-2).*10^(-12);cSoil_yrPG_MIROC(cSoil_yrPG_MIROC>10^(10))=NaN;
cSoil_yrPG_MIROC = nansum(cSoil_yrPG_MIROC,1); cSoil_yrPG_MIROC = nansum(cSoil_yrPG_MIROC,2); cSoil_yrPG_MIROC = squeeze(cSoil_yrPG_MIROC);

cSoil_yrPG_MPI = cSoil_yr_MPI.*sftlf_MPI.*area_MPI.*10^(-2).*10^(-12);cSoil_yrPG_MPI(cSoil_yrPG_MPI>10^(10))=NaN;
cSoil_yrPG_MPI = nansum(cSoil_yrPG_MPI,1); cSoil_yrPG_MPI = nansum(cSoil_yrPG_MPI,2); cSoil_yrPG_MPI = squeeze(cSoil_yrPG_MPI);

cSoil_yrPG_MRI = cSoil_yr_MRI.*sftlf_MRI.*area_MRI.*10^(-2).*10^(-12); cSoil_yrPG_MRI(cSoil_yrPG_MRI>10^(10))=NaN;
cSoil_yrPG_MRI = nansum(cSoil_yrPG_MRI,1); cSoil_yrPG_MRI = nansum(cSoil_yrPG_MRI,2); cSoil_yrPG_MRI = squeeze(cSoil_yrPG_MRI);
cSoil_yrPG_MRI = cat(1,NaN,cSoil_yrPG_MRI);

cSoil_yrPG_NOR = cSoil_yr_NOR.*sftlf_NOR.*area_NOR.*10^(-2).*10^(-12); cSoil_yrPG_NOR(cSoil_yrPG_NOR>10^(10))=NaN;
cSoil_yrPG_NOR = nansum(cSoil_yrPG_NOR,1); cSoil_yrPG_NOR = nansum(cSoil_yrPG_NOR,2); cSoil_yrPG_NOR = squeeze(cSoil_yrPG_NOR);

cSoil5_org_PG(1:5,1:11)=NaN;
cSoil5_org_PG(:,1) = cSoil_yrPG_BCC(152:156); 
cSoil5_org_PG(:,2) = cSoil_yrPG_CAN(152:156);  
cSoil5_org_PG(:,3) = cSoil_yrPG_CCSM(152:156); 
cSoil5_org_PG(:,4) = cSoil_yrPG_HAD(152:156);  
cSoil5_org_PG(:,5) = cSoil_yrPG_IPSL(152:156);  
cSoil5_org_PG(:,6) = cSoil_yrPG_MIROC(152:156);  
cSoil5_org_PG(:,7) = cSoil_yrPG_MPI(152:156);  
cSoil5_org_PG(:,8) = cSoil_yrPG_NOR(152:156);  
cSoil5_org_PG(:,9) = cSoil_yrPG_BNU(152:156);  
cSoil5_org_PG(:,10) = cSoil_yrPG_GF(152:156);  
cSoil5_org_PG(:,11) = cSoil_yrPG_MRI(152:156); 

% CMIP6
clearvars -except cVeg5_org_PG cSoil5_org_PG
cd('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\2_modelData\Cpool3_Native\CMIP6')
file = 'H:\CMIP56_Csink\2_DataAnalysis\3_originalYR\'
area_ASS = ncread([file,'ACCESS-ESM1-5\areacella_fx_ACCESS-ESM1-5_historical_r1i1p1f1_gn.nc'],'areacella');
sftlf_ASS = ncread([file,'ACCESS-ESM1-5\sftlf_fx_ACCESS-ESM1-5_historical_r1i1p1f1_gn.nc'],'sftlf');
cVeg_ASS = ncread([file,'ACCESS-ESM1-5\yr_cPool3_Lmon_ACCESS-ESM1-5_r1i1p1f1_gn_185001-210012.nc'],'cVeg');
cLit_ASS = ncread([file,'ACCESS-ESM1-5\yr_cPool3_Lmon_ACCESS-ESM1-5_r1i1p1f1_gn_185001-210012.nc'],'cLitter');
cSOM_ASS = ncread([file,'ACCESS-ESM1-5\yr_cPool3_Lmon_ACCESS-ESM1-5_r1i1p1f1_gn_185001-210012.nc'],'cSoil');
cSoil_ASS = cLit_ASS+cSOM_ASS;

area_BCC = ncread([file,'BCC-CSM2-MR\areacella_fx_BCC-CSM2-MR_hist-resIPO_r1i1p1f1_gn.nc'],'areacella');
sftlf_BCC = ncread([file,'BCC-CSM2-MR\sftlf_fx_BCC-CSM2-MR_hist-resIPO_r1i1p1f1_gn.nc'],'sftlf');
cVeg_BCC = ncread([file,'BCC-CSM2-MR\cPool3_BCC-CSM2-MR_r1i1p1f1_gn_185001-210012.nc'],'cVeg');
cLit_BCC = ncread([file,'BCC-CSM2-MR\cPool3_BCC-CSM2-MR_r1i1p1f1_gn_185001-210012.nc'],'cLitter');
cSOM_BCC = ncread([file,'BCC-CSM2-MR\cPool3_BCC-CSM2-MR_r1i1p1f1_gn_185001-210012.nc'],'cSoil');
cSoil_BCC = cLit_BCC+cSOM_BCC;

area_CAN = ncread([file,'CanESM5\areacella_fx_CanESM5_historical_r1i1p1f1_gn.nc'],'areacella');
sftlf_CAN = ncread([file,'CanESM5\sftlf_fx_CanESM5_historical_r1i1p1f1_gn.nc'],'sftlf');
cVeg_CAN = ncread([file,'CanESM5\cPool3_yr_CanESM5_r1i1p1f1_gn_185001-210012.nc'],'cVeg');
cLit_CAN = ncread([file,'CanESM5\cPool3_yr_CanESM5_r1i1p1f1_gn_185001-210012.nc'],'cLitter');
cSOM_CAN = ncread([file,'CanESM5\cPool3_yr_CanESM5_r1i1p1f1_gn_185001-210012.nc'],'cSoil');
cSoil_CAN = cLit_CAN+cSOM_CAN;

% CESM2 modeled soil carbon storage along soil profiles. Here, we include
% cSoilAbove1m in benchmark analysis
area_CESM = ncread([file,'CESM2\areacella_fx_CESM2_historical_r1i1p1f1_gn.nc'],'areacella');
sftlf_CESM = ncread([file,'CESM2\sftlf_fx_CESM2_historical_r1i1p1f1_gn.nc'],'sftlf');
cVeg_CESM = ncread([file,'CESM2\cPool3_yr_CESM2_r1i1p1f1_gn_185001-210012.nc'],'cVeg',[1 1 1],[Inf Inf 165]);
cLit_CESM = ncread([file,'CESM2\cPool3_yr_CESM2_r1i1p1f1_gn_185001-210012.nc'],'cLitter',[1 1 1],[Inf Inf 165]);
cSOM_CESM = ncread([file,'CESM2\2_cSoilAbove1m_Eyr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc'],'cSoilAbove1m');
cSoil_CESM = cLit_CESM+cSOM_CESM;

area_CNRM= ncread([file,'CNRM-ESM2-1\areacella_fx_CNRM-ESM2-1_amip_r1i1p1f2_gr.nc'],'areacella');
sftlf_CNRM = ncread([file,'CNRM-ESM2-1\sftlf_fx_CNRM-ESM2-1_amip_r1i1p1f2_gr.nc'],'sftlf');
cVeg_CNRM = ncread([file,'CNRM-ESM2-1\cPool3_yr_CNRM-ESM2-1_r1i1p1f2_gr_185001-210012.nc'],'cVeg');
cLit_CNRM = ncread([file,'CNRM-ESM2-1\cPool3_yr_CNRM-ESM2-1_r1i1p1f2_gr_185001-210012.nc'],'cLitter');
cSOM_CNRM = ncread([file,'CNRM-ESM2-1\cPool3_yr_CNRM-ESM2-1_r1i1p1f2_gr_185001-210012.nc'],'cSoil');
cSoil_CNRM = cLit_CNRM+cSOM_CNRM;

% EC-Earth3-Veg did not output areacella, we thus calculated based on lon_bnds and lat_bnds
sftlf_EC = ncread([file,'EC-Earth3-Veg\sftlf_fx_EC-Earth3-Veg_historical_r1i1p1f1_gr.nc'],'sftlf');
cVeg_EC= ncread([file,'EC-Earth3-Veg\cPool3_yr_EC-Earth3-Veg_r1i1p1f1_gr_185001-210012.nc'],'cVeg');
cLit_EC = ncread([file,'EC-Earth3-Veg\cPool3_yr_EC-Earth3-Veg_r1i1p1f1_gr_185001-210012.nc'],'cLitter');
cSOM_EC = ncread([file,'EC-Earth3-Veg\cPool3_yr_EC-Earth3-Veg_r1i1p1f1_gr_185001-210012.nc'],'cSoil');
cSoil_EC = cLit_EC+cSOM_EC;

file2 = 'H:\CMIP56_Csink\2_DataAnalysis\3_originalYR\EC-Earth3-Veg\'        
lon_bnds = ncread([file2,'sftlf_fx_EC-Earth3-Veg_historical_r1i1p1f1_gr.nc'],'lon_bnds');
lat_bnds = ncread([file2,'sftlf_fx_EC-Earth3-Veg_historical_r1i1p1f1_gr.nc'],'lat_bnds');
E = wgs84Ellipsoid; [~,m] = size(lat_bnds); [~,n] = size(lon_bnds);
lat_bnds = double(lat_bnds);
lon_bnds = double(lon_bnds);
for i = 1:m
    for j = 1:n
        cellarea(i,j) = areaquad(lat_bnds(1,i),lon_bnds(1,j),lat_bnds(2,i),lon_bnds(2,j),E); %m2
    end
end
area_EC = cellarea';    % unit:m2

area_IPSL= ncread([file,'IPSL-CM6A\areacella_fx_IPSL-CM6A-LR_historical_r1i1p1f1_gr.nc'],'areacella');
sftlf_IPSL = ncread([file,'IPSL-CM6A\sftlf_fx_IPSL-CM6A-LR_historical_r1i1p1f1_gr.nc'],'sftlf');
cVeg_IPSL = ncread([file,'IPSL-CM6A\2_cPool3_yr_IPSL-CM6A-LR_hist-ssp585_r1i1p1f1_gr_185001-210012.nc'],'cVeg');
cLit_IPSL = ncread([file,'IPSL-CM6A\2_cPool3_yr_IPSL-CM6A-LR_hist-ssp585_r1i1p1f1_gr_185001-210012.nc'],'cLitter');
cSOM_IPSL = ncread([file,'IPSL-CM6A\2_cPool3_yr_IPSL-CM6A-LR_hist-ssp585_r1i1p1f1_gr_185001-210012.nc'],'cSoil');
cSoil_IPSL = cLit_IPSL+cSOM_IPSL;

area_MRIOC= ncread([file,'MIROC-ES2L\areacella_fx_MIROC-ES2L_historical_r1i1p1f2_gn.nc'],'areacella');
sftlf_MRIOC = ncread([file,'MIROC-ES2L\sftlf_fx_MIROC-ES2L_historical_r1i1p1f2_gn.nc'],'sftlf');
cVeg_MRIOC = ncread([file,'MIROC-ES2L\cPool3_yr_MIROC-ES2L_r1i1p1f1_gn_185001-210012.nc'],'cVeg');
cLit_MRIOC = ncread([file,'MIROC-ES2L\cPool3_yr_MIROC-ES2L_r1i1p1f1_gn_185001-210012.nc'],'cLitter');
cSOM_MRIOC = ncread([file,'MIROC-ES2L\cPool3_yr_MIROC-ES2L_r1i1p1f1_gn_185001-210012.nc'],'cSoil');
cSoil_MRIOC = cLit_MRIOC+cSOM_MRIOC;

area_MPI= ncread([file,'MPI-ESM2-LR\areacella_fx_MPI-ESM1-2-LR_historical_r1i1p1f1_gn.nc'],'areacella');
sftlf_MPI = ncread([file,'MPI-ESM2-LR\sftlf_fx_MPI-ESM1-2-LR_historical_r1i1p1f1_gn.nc'],'sftlf');
cVeg_MPI = ncread([file,'MPI-ESM2-LR\cPool3_yr_MPI-ESM1-2-LR_r1i1p1f1_gn_185001-210012.nc'],'cVeg');
cLit_MPI = ncread([file,'MPI-ESM2-LR\cPool3_yr_MPI-ESM1-2-LR_r1i1p1f1_gn_185001-210012.nc'],'cLitter');
cSOM_MPI= ncread([file,'MPI-ESM2-LR\cPool3_yr_MPI-ESM1-2-LR_r1i1p1f1_gn_185001-210012.nc'],'cSoil');
cSoil_MPI = cLit_MPI+cSOM_MPI;

% NorESM2-LM modeled soil carbon storage along soil profiles. Here, we include
% cSoilAbove1m in benchmark analysis
area_NOR= ncread([file,'NorESM2-LM\areacella_fx_NorESM2-LM_historical_r1i1p1f1_gn.nc'],'areacella');
sftlf_NOR = ncread([file,'NorESM2-LM\sftlf_fx_NorESM2-LM_historical_r1i1p1f1_gn.nc'],'sftlf');
cVeg_NOR = ncread([file,'NorESM2-LM\cPool3_yr_NorESM2-LM_r1i1p1f1_gn_185001-210012.nc'],'cVeg',[1 1 1],[Inf Inf 165]);
cLit_NOR = ncread([file,'NorESM2-LM\cPool3_yr_NorESM2-LM_r1i1p1f1_gn_185001-210012.nc'],'cLitter',[1 1 1],[Inf Inf 165]);
cSOM_NOR = ncread([file,'NorESM2-LM\2_cSoilAbove1m_Eyr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc'],'cSoilAbove1m');
cSoil_NOR = cLit_NOR + cSOM_NOR;

area_UK= ncread([file,'UKESM1-0-LL\areacella_fx_UKESM1-0-LL_piControl_r1i1p1f2_gn.nc'],'areacella');
sftlf_UK = ncread([file,'UKESM1-0-LL\sftlf_fx_UKESM1-0-LL_piControl_r1i1p1f2_gn.nc'],'sftlf');
cVeg_UK = ncread([file,'UKESM1-0-LL\cPool2_yr__UKESM1-0-LL_r1i1p1f2_gn_185001-210012.nc'],'cVeg');
cSOM_UK = ncread([file,'UKESM1-0-LL\cPool2_yr__UKESM1-0-LL_r1i1p1f2_gn_185001-210012.nc'],'cSoil');
cSoil_UK = cSOM_UK;

% cVeg,convert unit from kgC m-2 into PgC
cVeg_ASS = cVeg_ASS.*area_ASS.*sftlf_ASS.*10^(-2).*10^(-12);cVeg_ASS(cVeg_ASS>10^(10))=NaN;
cVeg_ASS = nansum(cVeg_ASS,1); cVeg_ASS = nansum(cVeg_ASS,2); cVeg_ASS = squeeze(cVeg_ASS);

cVeg_BCC = cVeg_BCC.*area_BCC.*sftlf_BCC.*10^(-2).*10^(-12); cVeg_BCC(cVeg_BCC>10^(10))=NaN;
cVeg_BCC = nansum(cVeg_BCC,1); cVeg_BCC = nansum(cVeg_BCC,2); cVeg_BCC = squeeze(cVeg_BCC);

cVeg_CAN = cVeg_CAN.*area_CAN.*sftlf_CAN.*10^(-2).*10^(-12); cVeg_CAN(cVeg_CAN>10^(10))=NaN;
cVeg_CAN = nansum(cVeg_CAN,1); cVeg_CAN = nansum(cVeg_CAN,2); cVeg_CAN = squeeze(cVeg_CAN);

cVeg_CESM = cVeg_CESM.*area_CESM.*sftlf_CESM.*10^(-2).*10^(-12);  cVeg_CESM(cVeg_CESM>10^(10))=NaN;
cVeg_CESM = nansum(cVeg_CESM,1); cVeg_CESM = nansum(cVeg_CESM,2); cVeg_CESM = squeeze(cVeg_CESM);

cVeg_CNRM = cVeg_CNRM.*area_CNRM.*sftlf_CNRM.*10^(-2).*10^(-12); cVeg_CNRM(cVeg_CNRM>10^(10))=NaN; 
cVeg_CNRM = nansum(cVeg_CNRM,1); cVeg_CNRM = nansum(cVeg_CNRM,2); cVeg_CNRM = squeeze(cVeg_CNRM);

cVeg_EC = cVeg_EC.*area_EC.*sftlf_EC.*10^(-2).*10^(-12); cVeg_E(cVeg_EC>10^(10))=NaN;
cVeg_EC = nansum(cVeg_EC,1); cVeg_EC = nansum(cVeg_EC,2); cVeg_EC = squeeze(cVeg_EC);

cVeg_IPSL = cVeg_IPSL.*area_IPSL.*sftlf_IPSL.*10^(-2).*10^(-12); cVeg_IPSL(cVeg_IPSL>10^(10))=NaN;
cVeg_IPSL = nansum(cVeg_IPSL,1); cVeg_IPSL = nansum(cVeg_IPSL,2); cVeg_IPSL = squeeze(cVeg_IPSL);

cVeg_MRIOC = cVeg_MRIOC.*area_MRIOC.*sftlf_MRIOC.*10^(-2).*10^(-12);  cVeg_MRIOC(cVeg_MRIOC>10^(10))=NaN; 
cVeg_MRIOC = nansum(cVeg_MRIOC,1); cVeg_MRIOC = nansum(cVeg_MRIOC,2); cVeg_MRIOC = squeeze(cVeg_MRIOC);

cVeg_MPI = cVeg_MPI.*area_MPI.*sftlf_MPI.*10^(-2).*10^(-12); cVeg_MPI(cVeg_MPI>10^(10))=NaN; 
cVeg_MPI = nansum(cVeg_MPI,1); cVeg_MPI = nansum(cVeg_MPI,2); cVeg_MPI = squeeze(cVeg_MPI);

cVeg_NOR = cVeg_NOR.* area_NOR.* sftlf_NOR.*10^(-2).*10^(-12); cVeg_NOR(cVeg_NOR>10^(10))=NaN;
cVeg_NOR = nansum(cVeg_NOR,1); cVeg_NOR = nansum(cVeg_NOR,2); cVeg_NOR = squeeze(cVeg_NOR);

cVeg_UK = cVeg_UK.*area_UK.*sftlf_UK.*10^(-2).*10^(-12);   cVeg_UK(cVeg_UK>10^(10))=NaN; 
cVeg_UK = nansum(cVeg_UK,1); cVeg_UK = nansum(cVeg_UK,2); cVeg_UK = squeeze(cVeg_UK);


% cSoil, convert unit from kgC m-2 into PgC
cSoil_ASS = cSoil_ASS.*area_ASS.*sftlf_ASS.*10^(-2).*10^(-12); cSoil_ASS(cSoil_ASS>10^(10))=NaN;
cSoil_ASS = nansum(cSoil_ASS,1); cSoil_ASS = nansum(cSoil_ASS,2); cSoil_ASS = squeeze(cSoil_ASS);

cSoil_BCC = cSoil_BCC.*area_BCC.*sftlf_BCC.*10^(-2).*10^(-12); cSoil_BCC(cSoil_BCC>10^(10))=NaN;
cSoil_BCC = nansum(cSoil_BCC,1); cSoil_BCC = nansum(cSoil_BCC,2); cSoil_BCC = squeeze(cSoil_BCC);

cSoil_CAN = cSoil_CAN.*area_CAN.*sftlf_CAN.*10^(-2).*10^(-12); cSoil_CAN(cSoil_CAN>10^(10))=NaN;
cSoil_CAN = nansum(cSoil_CAN,1); cSoil_CAN = nansum(cSoil_CAN,2); cSoil_CAN = squeeze(cSoil_CAN);

cSoil_CESM  = cSoil_CESM.*area_CESM.*sftlf_CESM.*10^(-2).*10^(-12); cSoil_CESM(cSoil_CESM>10^(10))=NaN; 
cSoil_CESM = nansum(cSoil_CESM,1); cSoil_CESM = nansum(cSoil_CESM,2); cSoil_CESM = squeeze(cSoil_CESM);

cSoil_CNRM = cSoil_CNRM.*area_CNRM.*sftlf_CNRM.*10^(-2).*10^(-12); cSoil_CNRM(cSoil_CNRM>10^(10))=NaN; 
cSoil_CNRM = nansum(cSoil_CNRM,1); cSoil_CNRM = nansum(cSoil_CNRM,2); cSoil_CNRM = squeeze(cSoil_CNRM);

cSoil_EC = cSoil_EC.*area_EC.*sftlf_EC.*10^(-2).*10^(-12); cSoil_EC(cSoil_EC>10^(10))=NaN;
cSoil_EC = nansum(cSoil_EC,1); cSoil_EC = nansum(cSoil_EC,2); cSoil_EC = squeeze(cSoil_EC);

cSoil_IPSL = cSoil_IPSL.*area_IPSL.*sftlf_IPSL.*10^(-2).*10^(-12); cSoil_IPSL(cSoil_IPSL>10^(10))=NaN;
cSoil_IPSL = nansum(cSoil_IPSL,1); cSoil_IPSL = nansum(cSoil_IPSL,2); cSoil_IPSL = squeeze(cSoil_IPSL);

cSoil_MRIOC = cSoil_MRIOC.*area_MRIOC.*sftlf_MRIOC.*10^(-2).*10^(-12); cSoil_MRIOC(cSoil_MRIOC>10^(10))=NaN; 
cSoil_MRIOC  = nansum(cSoil_MRIOC ,1); cSoil_MRIOC = nansum(cSoil_MRIOC,2); cSoil_MRIOC  = squeeze(cSoil_MRIOC);

cSoil_MPI = cSoil_MPI.*area_MPI.*sftlf_MPI.*10^(-2).*10^(-12); cSoil_MPI(cSoil_MPI>10^(10))=NaN; 
cSoil_MPI = nansum(cSoil_MPI,1); cSoil_MPI = nansum(cSoil_MPI,2); cSoil_MPI = squeeze(cSoil_MPI);

cSoil_NOR = cSoil_NOR.* area_NOR.* sftlf_NOR.*10^(-2).*10^(-12); cSoil_NOR(cSoil_NOR>10^(10))=NaN;
cSoil_NOR = nansum(cSoil_NOR,1); cSoil_NOR = nansum(cSoil_NOR,2); cSoil_NOR = squeeze(cSoil_NOR);

cSoil_UK = cSoil_UK.*area_UK.*sftlf_UK.*10^(-2).*10^(-12); cSoil_UK(cSoil_UK>10^(10))=NaN; 
cSoil_UK = nansum(cSoil_UK,1); cSoil_UK = nansum(cSoil_UK,2); cSoil_UK = squeeze(cSoil_UK);

cVeg6_org_PG(1:5,1:11)=NaN;
cVeg6_org_PG(:,1) = cVeg_BCC(152:156); 
cVeg6_org_PG(:,2) = cVeg_CAN(152:156);  
cVeg6_org_PG(:,3) = cVeg_CESM(152:156); 
cVeg6_org_PG(:,4) = cVeg_UK(152:156);  
cVeg6_org_PG(:,5) = cVeg_IPSL(152:156);  
cVeg6_org_PG(:,6) = cVeg_MRIOC(152:156);  
cVeg6_org_PG(:,7) = cVeg_MPI(152:156);  
cVeg6_org_PG(:,8) = cVeg_NOR(152:156);  
cVeg6_org_PG(:,9) = cVeg_ASS(152:156);  
cVeg6_org_PG(:,10) = cVeg_CNRM(152:156);  
cVeg6_org_PG(:,11) = cVeg_EC(152:156); 

cSoil6_org_PG(1:5,1:11)=NaN;
cSoil6_org_PG(:,1) = cSoil_BCC(152:156); 
cSoil6_org_PG(:,2) = cSoil_CAN(152:156);  
cSoil6_org_PG(:,3) = cSoil_CESM(152:156); 
cSoil6_org_PG(:,4) = cSoil_UK(152:156);  
cSoil6_org_PG(:,5) = cSoil_IPSL(152:156);  
cSoil6_org_PG(:,6) = cSoil_MRIOC(152:156);  
cSoil6_org_PG(:,7) = cSoil_MPI(152:156);  
cSoil6_org_PG(:,8) = cSoil_NOR(152:156);  
cSoil6_org_PG(:,9) = cSoil_ASS(152:156);  
cSoil6_org_PG(:,10) = cSoil_CNRM(152:156);  
cSoil6_org_PG(:,11) = cSoil_EC(152:156); 

% CMIP6
clearvars -except cVeg5_org_PG cSoil5_org_PG cVeg6_org_PG cSoil6_org_PG
save E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\4_MatData\cVeg_cSoil_CMIP56_native.mat











