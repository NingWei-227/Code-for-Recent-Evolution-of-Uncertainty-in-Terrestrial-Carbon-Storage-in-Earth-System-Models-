% write .mat data into csv for further her.part analysis in R
clear;clc
cd('file path of matData')

load('cLand_cmip5_tmp.mat') % C storage for 11 models
load('Xc_cmip5_tmp.mat')    % C storage capacity for 11 models
load('Xp_cmip5_tmp.mat')    % C storage potential for 11 models
load('NPP_cmip5_tmp.mat')
load('tuaE_cmip5_tmp.mat')

load('GPP_cmip5_tmp.mat')
load('CUE_cmip5_tmp.mat')
load('TuaEbase_CMIP5.mat')
load('tuaE_scaler_CMIP5.mat')
load('tuaEop_CMIP5.mat')

leg5_str = {'Year','BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'}
Years = 1851:2005;

NaNgf = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]';
CLgf_tmp = [NaNgf; CLgf_tmp];
XcGf_tmp = [NaNgf; XcGf_tmp];
XpGf_tmp = [NaNgf; XpGf_tmp];
NPPgf_tmp = [NaNgf; NPPgf_tmp];
tuaEgf_tmp = [NaNgf; tuaEgf_tmp];
GPPgf_tmp = [NaNgf; GPPgf_tmp];
CUEgf_tmp = [NaNgf; CUEgf_tmp];
tuaE_scaler_gf2G = [[NaNgf NaNgf NaNgf];tuaE_scaler_gf2G'];
tuaE_scaler_gf2G = tuaE_scaler_gf2G';
tuaE_opGf2G = [NaNgf' tuaE_opGf2G];

NaNhad = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]';
CLhad_tmp = [NaNhad; CLhad_tmp];
XcHad_tmp = [NaNhad; XcHad_tmp];
XpHad_tmp = [NaNhad; XpHad_tmp];
NPPhad_tmp = [NaNhad; NPPhad_tmp];
tuaEhad_tmp = [NaNhad; tuaEhad_tmp];
GPPhad_tmp = [NaNhad; GPPhad_tmp];
CUEhad_tmp = [NaNhad; CUEhad_tmp];
tuaE_scaler_hadES = [[NaNhad NaNhad NaNhad];tuaE_scaler_hadES'];
tuaE_scaler_hadES = tuaE_scaler_hadES';
tuaE_opHadES = [NaNhad' tuaE_opHadES];

CLmri_tmp = [NaN; CLmri_tmp];
XcMri_tmp = [NaN; XcMri_tmp;];
XpMri_tmp = [NaN; XpMri_tmp;];
NPPmri_tmp = [NaN; NPPmri_tmp];
tuaEmri_tmp = [NaN; tuaEmri_tmp];
GPPmri_tmp = [NaN; GPPmri_tmp];
CUEmri_tmp = [NaN; CUEmri_tmp];
tuaE_scaler_mri= [[NaN NaN NaN];tuaE_scaler_mri'];
tuaE_scaler_mri = tuaE_scaler_mri';
tuaE_opMri = [NaN' tuaE_opMri];

X5_all(1:155,1:11) = NaN;
X5_all(:,1) = CLbcc_tmp(2:156);
X5_all(:,2) = CLcan_tmp(2:156);
X5_all(:,3) = CLccsm_tmp(2:156);
X5_all(:,4) = CLhad_tmp(2:156);
X5_all(:,5) = CLipsl_tmp(2:156);
X5_all(:,6) = CLmiroc_tmp(2:156);
X5_all(:,7) = CLmpi_tmp(2:156);
X5_all(:,8) = CLnor_tmp(2:156);
X5_all(:,9) = CLbnu_tmp(2:156);
X5_all(:,10) = CLgf_tmp(2:156);
X5_all(:,11) = CLmri_tmp(2:156);

Xc5_all(1:155,1:11) = NaN;
Xc5_all(:,1) = XcBcc_tmp(2:156);
Xc5_all(:,2) = XcCan_tmp(2:156);
Xc5_all(:,3) = XcCcsm_tmp(2:156);
Xc5_all(:,4) = XcHad_tmp(2:156);
Xc5_all(:,5) = XcIpsl_tmp(2:156);
Xc5_all(:,6) = XcMiroc_tmp(2:156);
Xc5_all(:,7) = XcMpi_tmp(2:156);
Xc5_all(:,8) = XcNor_tmp(2:156);
Xc5_all(:,9) = XcBnu_tmp(2:156);
Xc5_all(:,10) = XcGf_tmp(2:156);
Xc5_all(:,11) = XcMri_tmp(2:156);

Xp5_all(1:155,1:11) = NaN;
Xp5_all(:,1) = XpBcc_tmp(2:156);
Xp5_all(:,2) = XpCan_tmp(2:156);
Xp5_all(:,3) = XpCcsm_tmp(2:156);
Xp5_all(:,4) = XpHad_tmp(2:156);
Xp5_all(:,5) = XpIpsl_tmp(2:156);
Xp5_all(:,6) = XpMiroc_tmp(2:156);
Xp5_all(:,7) = XpMpi_tmp(2:156);
Xp5_all(:,8) = XpNor_tmp(2:156);
Xp5_all(:,9) = XpBnu_tmp(2:156);
Xp5_all(:,10) = XpGf_tmp(2:156);
Xp5_all(:,11) = XpMri_tmp(2:156);

NPP5_all(1:155,1:11) = NaN;
NPP5_all(:,1) = NPPbcc_tmp(2:156);
NPP5_all(:,2) = NPPcan_tmp(2:156);
NPP5_all(:,3) = NPPccsm_tmp(2:156);
NPP5_all(:,4) = NPPhad_tmp(2:156);
NPP5_all(:,5) = NPPipsl_tmp(2:156);
NPP5_all(:,6) = NPPmiroc_tmp(2:156);
NPP5_all(:,7) = NPPmpi_tmp(2:156);
NPP5_all(:,8) = NPPnor_tmp(2:156);
NPP5_all(:,9) = NPPbnu_tmp(2:156);
NPP5_all(:,10) = NPPgf_tmp(2:156);
NPP5_all(:,11) = NPPmri_tmp(2:156);

tuaE5_all(1:155,1:11) = NaN;
tuaE5_all(:,1) = tuaEbcc_tmp(2:156);
tuaE5_all(:,2) = tuaEcan_tmp(2:156);
tuaE5_all(:,3) = tuaEccsm_tmp(2:156);
tuaE5_all(:,4) = tuaEhad_tmp(2:156);
tuaE5_all(:,5) = tuaEipsl_tmp(2:156);
tuaE5_all(:,6) = tuaEmiroc_tmp(2:156);
tuaE5_all(:,7) = tuaEmpi_tmp(2:156);
tuaE5_all(:,8) = tuaEnor_tmp(2:156);
tuaE5_all(:,9) = tuaEbnu_tmp(2:156);
tuaE5_all(:,10) = tuaEgf_tmp(2:156);
tuaE5_all(:,11) = tuaEmri_tmp(2:156);

GPP5_all(1:155,1:11) = NaN;
GPP5_all(:,1) = GPPbcc_tmp(2:156);
GPP5_all(:,2) = GPPcan_tmp(2:156);
GPP5_all(:,3) = GPPccsm_tmp(2:156);
GPP5_all(:,4) = GPPhad_tmp(2:156);
GPP5_all(:,5) = GPPipsl_tmp(2:156);
GPP5_all(:,6) = GPPmiroc_tmp(2:156);
GPP5_all(:,7) = GPPmpi_tmp(2:156);
GPP5_all(:,8) = GPPnor_tmp(2:156);
GPP5_all(:,9) = GPPbnu_tmp(2:156);
GPP5_all(:,10) = GPPgf_tmp(2:156);
GPP5_all(:,11) = GPPmri_tmp(2:156);

CUE5_all(1:155,1:11) = NaN;
CUE5_all(:,1) = CUEbcc_tmp(2:156);   %BCC-CSM1-1m
CUE5_all(:,2) = CUEcan_tmp(2:156);   %CanESM2
CUE5_all(:,3) = CUEccsm_tmp(2:156);  %CCSM4
CUE5_all(:,4) = CUEhad_tmp(2:156);   %HadGEM2-ES (the modeled historical period: 1861-2005)
CUE5_all(:,5) = CUEipsl_tmp(2:156);  %IPSL-CM5A-MR
CUE5_all(:,6) = CUEmiroc_tmp(2:156); %MIROC-ESM
CUE5_all(:,7) = CUEmpi_tmp(2:156);   %MPI-ESM-MR 
CUE5_all(:,8) = CUEnor_tmp(2:156);   %NorESM1-M
CUE5_all(:,9) = CUEbnu_tmp(2:156);   %BNU-ESM
CUE5_all(:,10) = CUEgf_tmp(2:156);   %GFDL-ESM2G (the modeled historical period: 1862-2005)
CUE5_all(:,11) = CUEmri_tmp(2:156);  %MRI_ESM1

TuaEbase5_all(1:155,1:11) = NaN;
TuaEbase5_all(:,1) = baseTuaE_bcc;
TuaEbase5_all(:,2) = baseTuaE_can;
TuaEbase5_all(:,3) = baseTuaE_ccsm;
TuaEbase5_all(11:155,4) = baseTuaE_hadES;
TuaEbase5_all(:,5) = baseTuaE_ipsl5A;
TuaEbase5_all(:,6) = baseTuaE_miroc;
TuaEbase5_all(:,7) = baseTuaE_mpiM;
TuaEbase5_all(:,8) = baseTuaE_norM;
TuaEbase5_all(:,9) = baseTuaE_bnu;
TuaEbase5_all(12:155,10) = baseTuaE_gf2G;
TuaEbase5_all(2:155,11) = baseTuaE_mri;  

TuaE5_scaler_tas(1:155,1:11) = NaN;
TuaE5_scaler_tas(:,1) = tuaE_scaler_bcc(1,:);
TuaE5_scaler_tas(:,2) = tuaE_scaler_can(1,:);
TuaE5_scaler_tas(:,3) = tuaE_scaler_ccsm(1,:);
TuaE5_scaler_tas(:,4) = tuaE_scaler_hadES(1,:);
TuaE5_scaler_tas(:,5) = tuaE_scaler_ipsl5A(1,:);
TuaE5_scaler_tas(:,6) = tuaE_scaler_miroc(1,:);
TuaE5_scaler_tas(:,7) = tuaE_scaler_mpiM(1,:);
TuaE5_scaler_tas(:,8) = tuaE_scaler_norM(1,:);
TuaE5_scaler_tas(:,9) = tuaE_scaler_bnu(1,:);
TuaE5_scaler_tas(:,10) = tuaE_scaler_gf2G(1,:);
TuaE5_scaler_tas(:,11) = tuaE_scaler_mri(1,:); 

TuaE5_scaler_pr(1:155,1:11) = NaN;
TuaE5_scaler_pr(:,1) = tuaE_scaler_bcc(2,:);
TuaE5_scaler_pr(:,2) = tuaE_scaler_can(2,:);
TuaE5_scaler_pr(:,3) = tuaE_scaler_ccsm(2,:);
TuaE5_scaler_pr(:,4) = tuaE_scaler_hadES(2,:);
TuaE5_scaler_pr(:,5) = tuaE_scaler_ipsl5A(2,:);
TuaE5_scaler_pr(:,6) = tuaE_scaler_miroc(2,:);
TuaE5_scaler_pr(:,7) = tuaE_scaler_mpiM(2,:);
TuaE5_scaler_pr(:,8) = tuaE_scaler_norM(2,:);
TuaE5_scaler_pr(:,9) = tuaE_scaler_bnu(2,:);
TuaE5_scaler_pr(:,10) = tuaE_scaler_gf2G(2,:);
TuaE5_scaler_pr(:,11) = tuaE_scaler_mri(2,:); 

tuaE5_opt(1:155,1:11) = NaN;
tuaE5_opt(:,1) = tuaE_opBcc;
tuaE5_opt(:,2) = tuaE_opCan;
tuaE5_opt(:,3) = tuaE_opCcsm;
tuaE5_opt(:,4) = tuaE_opHadES;
tuaE5_opt(:,5) = tuaE_opIpsl5A;
tuaE5_opt(:,6) = tuaE_opMiroc;
tuaE5_opt(:,7) = tuaE_opMpiM;
tuaE5_opt(:,8) = tuaE_opNorM;
tuaE5_opt(:,9) = tuaE_opBnu;
tuaE5_opt(:,10) = tuaE_opGf2G;
tuaE5_opt(:,11) = tuaE_opMri;


X5_all = [Years', X5_all];
X5_all = [leg5_str; num2cell(X5_all)];

Xc5_all = [Years', Xc5_all];
Xc5_all = [leg5_str; num2cell(Xc5_all)];

Xp5_all = [Years', Xp5_all];
Xp5_all = [leg5_str; num2cell(Xp5_all)];

NPP5_all = [Years', NPP5_all];
NPP5_all = [leg5_str; num2cell(NPP5_all)];
tuaE5_all = [Years', tuaE5_all];
tuaE5_all = [leg5_str; num2cell(tuaE5_all)];

GPP5_all = [Years', GPP5_all];
GPP5_all = [leg5_str; num2cell(GPP5_all)];
CUE5_all = [Years', CUE5_all];
CUE5_all = [leg5_str; num2cell(CUE5_all)];

TuaEbase5_all = [Years', TuaEbase5_all];
TuaEbase5_all = [leg5_str; num2cell(TuaEbase5_all)];

TuaE5_scaler_tas = [Years', TuaE5_scaler_tas];
TuaE5_scaler_tas = [leg5_str; num2cell(TuaE5_scaler_tas)];
TuaE5_scaler_pr = [Years', TuaE5_scaler_pr];
TuaE5_scaler_pr = [leg5_str; num2cell(TuaE5_scaler_pr)];
tuaE5_opt = [Years', tuaE5_opt];
tuaE5_opt = [leg5_str; num2cell(tuaE5_opt)];

xlswrite('CMIP5_temp11.xlsx',X5_all,'X')
xlswrite('CMIP5_temp11.xlsx',Xc5_all,'Xc')
xlswrite('CMIP5_temp11.xlsx',Xp5_all,'Xp')
xlswrite('CMIP5_temp11.xlsx',NPP5_all,'NPP')
xlswrite('CMIP5_temp11.xlsx',tuaE5_all,'tuaE')
xlswrite('CMIP5_temp11.xlsx',GPP5_all,'GPP')
xlswrite('CMIP5_temp11.xlsx',CUE5_all,'CUE')
   
xlswrite('CMIP5_temp11.xlsx',TuaE5_scaler_tas ,'scaler_tas')    
xlswrite('CMIP5_temp11.xlsx',TuaE5_scaler_pr,'scaler_pr') 
xlswrite('CMIP5_temp11.xlsx',TuaEbase5_all,'tuaE_Base') 
xlswrite('CMIP5_temp11.xlsx',tuaE5_opt,'tuaE5_opt')    









