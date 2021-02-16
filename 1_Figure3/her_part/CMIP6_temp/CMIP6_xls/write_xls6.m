% write .mat data into csv for further her.part analysis in R
% 0-1m soil C are used for models with soil multiple profiles
clear;clc
cd('file_path\CMIP6_temp')

load('2cLand_cmip6_tmp.mat') % C storage for 11 models
load('2Xc_cmip6_tmp.mat')    % C storage capacity for 11 models
load('2Xp_cmip6_tmp.mat')    % C storage potential for 11 models
load('2NPP_cmip6_tmp.mat')
load('2tuaE_cmip6_tmp.mat')

load('2GPP_cmip6_tmp.mat')
load('2CUE_cmip6_tmp.mat')
load('TuaEbase_CMIP6.mat')
load('tuaE_scaler_CMIP6.mat')
load('tuaEop_CMIP6.mat')

leg6_str = {'Years','BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'}
Years = 1851:2005;


X6_all = [CLbcc_tmp(2:156),CLcan_tmp(2:156),CLcesm1m_tmp(2:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   CLuk_tmp(2:156),CLipsl_tmp(2:156),CLmic_tmp(2:156),...          % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   CLmpi_tmp(2:156),CLnor1m_tmp(2:156),...                           % MPI-ESM1-2-LR, NorESM2
   CLass_tmp(2:156),CLcnrm_tmp(2:156),CLec_tmp(2:156)];            % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

Xc6_all = [XcBcc_tmp(2:156),XcCan_tmp(2:156),XcCesm_tmp1m(2:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   XcUk_tmp(2:156),XcIpsl_tmp(2:156),XcMic_tmp(2:156),...          % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   XcMpi_tmp(2:156),XcNor_tmp1m(2:156),...                           % MPI-ESM1-2-LR, NorESM2
   XcAss_tmp(2:156),XcCnrm_tmp(2:156),XcEc_tmp(2:156)];            % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

Xp6_all = [XpBcc_tmp(2:156),XpCan_tmp(2:156),XpCesm_tmp1m(2:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   XpUk_tmp(2:156),XpIpsl_tmp(2:156),XpMic_tmp(2:156),...          % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   XpMpi_tmp(2:156),XpNor_tmp1m(2:156),...                           % MPI-ESM1-2-LR, NorESM2
   XpAss_tmp(2:156),XpCnrm_tmp(2:156),XpEc_tmp(2:156)];            % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

NPP6_all = [NPPbcc_tmp(2:156),NPPcan_tmp(2:156),NPPcesm_tmp(2:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   NPPuk_tmp(2:156),NPPipsl_tmp(2:156),NPPmic_tmp(2:156),...            % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   NPPmpi_tmp(2:156),NPPnor_tmp(2:156),...                              % MPI-ESM1-2-LR, NorESM2
   NPPass_tmp(2:156),NPPcnrm_tmp(2:156),NPPec_tmp(2:156)];              % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

tuaE6_all = [tuaEbcc_tmp(2:156),tuaEcan_tmp(2:156),tuaEcesm_tmp1m(2:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   tuaEuk_tmp(2:156),tuaEipsl_tmp(2:156),tuaEmic_tmp(2:156),...             % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   tuaEmpi_tmp(2:156),tuaEnor_tmp1m(2:156),...                                % MPI-ESM1-2-LR, NorESM2
   tuaEass_tmp(2:156),tuaEcnrm_tmp(2:156),tuaEec_tmp(2:156)];               % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

GPP6_all = [GPPbcc_tmp(2:156),GPPcan_tmp(2:156),GPPcesm_tmp(2:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   GPPuk_tmp(2:156),GPPipsl_tmp(2:156),GPPmic_tmp(2:156),...          % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   GPPmpi_tmp(2:156),GPPnor_tmp(2:156),...                           % MPI-ESM1-2-LR, NorESM2
   GPPass_tmp(2:156),GPPcnrm_tmp(2:156),GPPec_tmp(2:156)];            % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

CUE6_all = [CUEbcc_tmp(2:156),CUEcan_tmp(2:156),CUEcesm_tmp(2:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   CUEuk_tmp(2:156),CUEipsl_tmp(2:156),CUEmic_tmp(2:156),...          % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   CUEmpi_tmp(2:156),CUEnor_tmp(2:156),...                           % MPI-ESM1-2-LR, NorESM2
   CUEass_tmp(2:156),CUEcnrm_tmp(2:156),CUEec_tmp(2:156)];            % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

TuaEbase6_all(1:155,1:11) = NaN;
TuaEbase6_all(:,1) = baseTuaE_bccC;
TuaEbase6_all(:,2) = baseTuaE_can;
TuaEbase6_all(:,3) = baseTuaE_cesm;
TuaEbase6_all(:,4) = baseTuaE_uk;
TuaEbase6_all(:,5) = baseTuaE_ipsl;
TuaEbase6_all(:,6) = baseTuaE_miroc;
TuaEbase6_all(:,7) = baseTuaE_mpi;
TuaEbase6_all(:,8) = baseTuaE_nor;
TuaEbase6_all(:,9) = baseTuaE_ASS;
TuaEbase6_all(:,10) = baseTuaE_cnrm;
TuaEbase6_all(:,11) = baseTuaE_ec;

TuaE6_scaler_tas = [tuaE_scaler_bccC(1,:)',tuaE_scaler_can(1,:)',tuaE_scaler_cesm(1,:)',...% BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   tuaE_scaler_uk(1,:)',tuaE_scaler_ipsl(1,:)',tuaE_scaler_miroc(1,:)',...             % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   tuaE_scaler_mpi(1,:)',tuaE_scaler_nor(1,:)',...                                    % MPI-ESM1-2-LR, NorESM2
   tuaE_scaler_ASS(1,:)',tuaE_scaler_cnrm(1,:)',tuaE_scaler_ec(1,:)'];                 % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

TuaE6_scaler_pr = [tuaE_scaler_bccC(2,:)',tuaE_scaler_can(2,:)',tuaE_scaler_cesm(2,:)',...% BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   tuaE_scaler_uk(2,:)',tuaE_scaler_ipsl(2,:)',tuaE_scaler_miroc(2,:)',...             % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   tuaE_scaler_mpi(2,:)',tuaE_scaler_nor(2,:)',...                                    % MPI-ESM1-2-LR, NorESM2
   tuaE_scaler_ASS(2,:)',tuaE_scaler_cnrm(2,:)',tuaE_scaler_ec(2,:)'];                 % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

% TuaE_op = basedTuaE * pr_scaler * tas_scaler
tuaEop6_all = [tuaE_opBccC',tuaE_opCan', tuaE_opCesm',tuaE_opUk',...
               tuaE_opIpsl', tuaE_opMiroc', tuaE_opMpi', tuaE_opNor',...
               tuaE_opASS', tuaE_opCnrm', tuaE_opEc'];

X6_all = [Years', X6_all];
X6_all = [leg6_str; num2cell(X6_all)];

Xc6_all = [Years', Xc6_all];
Xc6_all = [leg6_str; num2cell(Xc6_all)];

Xp6_all = [Years', Xp6_all];
Xp6_all = [leg6_str; num2cell(Xp6_all)];

NPP6_all = [Years', NPP6_all];
NPP6_all = [leg6_str; num2cell(NPP6_all)];
tuaE6_all = [Years', tuaE6_all];
tuaE6_all = [leg6_str; num2cell(tuaE6_all)];

GPP6_all = [Years', GPP6_all];
GPP6_all = [leg6_str; num2cell(GPP6_all)];
CUE6_all = [Years', CUE6_all];
CUE6_all = [leg6_str; num2cell(CUE6_all)];

TuaEbase6_all = [Years', TuaEbase6_all];
TuaEbase6_all = [leg6_str; num2cell(TuaEbase6_all)];

TuaE6_scaler_tas = [Years', TuaE6_scaler_tas];
TuaE6_scaler_tas = [leg6_str; num2cell(TuaE6_scaler_tas)];
TuaE6_scaler_pr = [Years', TuaE6_scaler_pr];
TuaE6_scaler_pr = [leg6_str; num2cell(TuaE6_scaler_pr)];
tuaEop6_all = [Years', tuaEop6_all];
tuaEop6_all = [leg6_str; num2cell(tuaEop6_all)];

xlswrite('CMIP6_temp11.xlsx',X6_all,'X')
xlswrite('CMIP6_temp11.xlsx',Xc6_all,'Xc')
xlswrite('CMIP6_temp11.xlsx',Xp6_all,'Xp')
xlswrite('CMIP6_temp11.xlsx',NPP6_all,'NPP')
xlswrite('CMIP6_temp11.xlsx',tuaE6_all,'tuaE')
xlswrite('CMIP6_temp11.xlsx',GPP6_all,'GPP')
xlswrite('CMIP6_temp11.xlsx',CUE6_all,'CUE')
   
xlswrite('CMIP6_temp11.xlsx',TuaE6_scaler_tas ,'scaler_tas')    
xlswrite('CMIP6_temp11.xlsx',TuaE6_scaler_pr,'scaler_pr') 
xlswrite('CMIP6_temp11.xlsx',TuaEbase6_all,'tuaE_Base')  
xlswrite('CMIP6_temp11.xlsx',tuaEop6_all,'tuaE_opt')    















