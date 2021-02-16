% In CMIP6, CESM2 and NorESM2-LM simulated soil carbon storage along soil depth
% In the traceability analysis, we used soil carbon above 1m
% unit: KgC m-2 
clear
clc
cd('file_path')
load('2NPP_cmip6_tmp.mat')
clearvars -except NPPcesm_tmp NPPnor_tmp
cSoil_CESM1m_month = ncread('cSoilAbove1m_Emon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cSoilAbove1m');
cSoil_CESM1m(1:288,1:192,1:165) = NaN;
for yr = 1:165
    yr
    interval_yr = (yr-1).*12+1: yr.*12;
    cSoil_cesm1m_12 = cSoil_CESM1m_month(:,:,interval_yr);
    cSoil_CESM1m(:,:,yr) =  nanmean(cSoil_cesm1m_12,3);    
end
cSoil_NOR1m = ncread('2_cSoilAbove1m_Eyr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cSoilAbove1m');

cVeg6_CESM = ncread('2_cVeg_yr_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');
cVeg6_NOR = ncread('2_cVeg_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cVeg');

cLit6_CESM = ncread('2_cLitter_Lmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');
cLit6_NOR = ncread('2_cLitter_yr_NorESM2-LM_historical_r1i1p1f1_gn_185001-201412.nc','cLitter');

cLand_CESM1m = cVeg6_CESM + cLit6_CESM + cSoil_CESM1m; % unit: KgC m-2
cLand_NOR1m = cVeg6_NOR + cLit6_NOR + cSoil_NOR1m;     % unit: KgC m-2

areacell_cesm = ncread('areacella_fx_CESM2_historical_r1i1p1f1_gn.nc','areacella');    %unit: m2
sftlf_cesm = ncread('sftlf_fx_CESM2_historical_r1i1p1f1_gn.nc','sftlf');               %unit: %
cellarea_CESM = areacell_cesm.* sftlf_cesm*10^(-2);                                    %unit: m2
CLcesm1m_pg = cLand_CESM1m.* cellarea_CESM.* 10^(-12);  % unit: PgC
CLcesm1m_tmp = nansum(CLcesm1m_pg,1); CLcesm1m_tmp = nansum(CLcesm1m_tmp,2);
CLcesm1m_tmp = squeeze(CLcesm1m_tmp);

areacell_nor = ncread('areacella_fx_NorESM2-LM_historical_r1i1p1f1_gn.nc','areacella');  %unit: m2
sftlf_nor = ncread('sftlf_fx_NorESM2-LM_historical_r1i1p1f1_gn.nc','sftlf');             %unit: %
cellarea_nor = areacell_nor.* sftlf_nor*10^(-2);                                                                            %unit: m2
CLnor1m_pg = cLand_NOR1m.* cellarea_nor.* 10^(-12);  % unit: PgC
CLnor1m_tmp = nansum(CLnor1m_pg,1); CLnor1m_tmp = nansum(CLnor1m_tmp,2);
CLnor1m_tmp = squeeze(CLnor1m_tmp);

nyear = 165;
for yr = 1: nyear-1
    NetCcesm_tmp1m(yr+1,1) = CLcesm1m_tmp(yr+1) - CLcesm1m_tmp(yr);
    NetCnor_tmp1m(yr+1,1) = CLnor1m_tmp(yr+1) - CLnor1m_tmp(yr);
    
    tuaEcesm_tmp1m(yr+1,1) = CLcesm1m_tmp(yr+1)./(NPPcesm_tmp(yr+1) - NetCcesm_tmp1m(yr+1));
    tuaEnor_tmp1m(yr+1,1) = CLnor1m_tmp(yr+1)./(NPPnor_tmp(yr+1) - NetCnor_tmp1m(yr+1));
    
end

XcCesm_tmp1m = tuaEcesm_tmp1m.*NPPcesm_tmp(1:165);
XpCesm_tmp1m = XcCesm_tmp1m - CLcesm1m_tmp;

XcNor_tmp1m = tuaEnor_tmp1m.*NPPnor_tmp(1:165);
XpNor_tmp1m = XcNor_tmp1m - CLnor1m_tmp;

cd('file path to write results into .mat datasets')
save('2cLand_cmip6_tmp.mat','CLcesm1m_tmp','-append')
save('2tuaE_cmip6_tmp.mat','tuaEcesm_tmp1m','-append')
save('2Xc_cmip6_tmp.mat','XcCesm_tmp1m','-append')
save('2Xp_cmip6_tmp.mat','XpCesm_tmp1m','-append')

save('2cLand_cmip6_tmp.mat','CLnor1m_tmp','-append')
save('2tuaE_cmip6_tmp.mat','tuaEnor_tmp1m','-append')
save('2Xc_cmip6_tmp.mat','XcNor_tmp1m','-append')
save('2Xp_cmip6_tmp.mat','XpNor_tmp1m','-append')

%% optimize analysis 
clear
clc
cd('file path of matData')
load('2pr_cmip6_tmp.mat')
load('2tas_cmip6_tmp.mat')
load('2tuaE_cmip6_tmp.mat')
cd('set the file path to save results')
%% ACCESS-ESM1-5
% put pre tas tuaE data into inputs
inputs_ASS(1,:) = tasASS_tmp(2:156,:);
inputs_ASS(2,:) =  prASS_tmp(2:156,:);
inputs_ASS(3,:) =  tuaEass_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_ASS,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_ASS),[1,20],[],[],[],[],...
    [0,0],[5,max(inputs_ASS(3,:))],@(x)myfunc_tuaE(x,inputs_ASS));

Q10_ASS = Q10base_ASS (1); % x(1): Q10
baseTuaE_ASS = Q10base_ASS(2); %x(2): baseTuaE
tuaE_opASS = cal_tuaE_op (Q10base_ASS, inputs_ASS); % calculate the estimated optimal tuaE
tuaE_scaler_ASS = cal_scaler(Q10base_ASS, inputs_ASS); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_ASS = cal_R2rmse(Q10base_ASS, inputs_ASS) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP6.mat','Q10_ASS');
save('TuaEbase_CMIP6.mat','baseTuaE_ASS')
save('tuaEop_CMIP6.mat','tuaE_opASS')
save('tuaE_scaler_CMIP6.mat','tuaE_scaler_ASS')
save('tuaE_R2rmse_CMIP6.mat','tuaE_R2rmse_ASS')

%% BCC-CSM2-MR
% put pre tas tuaE data into inputs
inputs_bccC(1,:) = tasBCC_tmp(2:156,:);
inputs_bccC(2,:) =  prBCC_tmp(2:156,:);
inputs_bccC(3,:) =  tuaEbcc_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_bccC,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_bccC),[1,20],[],[],[],[],...
    [0,0],[5,max(inputs_bccC(3,:))],@(x)myfunc_tuaE(x,inputs_bccC));

Q10_bccC = Q10base_bccC (1); % x(1): Q10
baseTuaE_bccC = Q10base_bccC(2); %x(2): baseTuaE
tuaE_opBccC = cal_tuaE_op (Q10base_bccC, inputs_bccC); % calculate the estimated optimal tuaE
tuaE_scaler_bccC = cal_scaler(Q10base_bccC, inputs_bccC); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_bccC = cal_R2rmse(Q10base_bccC, inputs_bccC) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP6.mat','Q10_bccC','-append');
save('TuaEbase_CMIP6.mat','baseTuaE_bccC','-append')
save('tuaEop_CMIP6.mat','tuaE_opBccC','-append')
save('tuaE_scaler_CMIP6.mat','tuaE_scaler_bccC','-append')
save('tuaE_R2rmse_CMIP6.mat','tuaE_R2rmse_bccC','-append')

%% CanESM5
inputs_can(1,:) =  tasCAN_tmp(2:156,:);
inputs_can(2,:) =  prCAN_tmp(2:156,:);
inputs_can(3,:) =  tuaEcan_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_can,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_can),[1,20],[],[],[],[],...
    [0,0],[15,max(inputs_can(3,:))],@(x)myfunc_tuaE(x,inputs_can));

Q10_can = Q10base_can (1); % x(1): Q10
baseTuaE_can = Q10base_can(2); %x(2): baseTuaE
tuaE_opCan = cal_tuaE_op (Q10base_can, inputs_can); % calculate the estimated optimal tuaE
tuaE_scaler_can = cal_scaler(Q10base_can, inputs_can); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_can = cal_R2rmse(Q10base_can, inputs_can) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP6.mat','Q10_can','-append');
save('TuaEbase_CMIP6.mat','baseTuaE_can','-append')
save('tuaEop_CMIP6.mat','tuaE_opCan','-append')
save('tuaE_scaler_CMIP6.mat','tuaE_scaler_can','-append')
save('tuaE_R2rmse_CMIP6.mat','tuaE_R2rmse_can','-append')

%% CESM2
inputs_cesm(1,:) =  tasCESM_tmp(2:156,:);
inputs_cesm(2,:) =  prCESM_tmp(2:156,:);
inputs_cesm(3,:) =  tuaEcesm_tmp1m(2:156,:);

[Q10base_cesm,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_cesm),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_cesm(3,:))],@(x)myfunc_tuaE(x,inputs_cesm));

Q10_cesm = Q10base_cesm (1); % x(1): Q10
baseTuaE_cesm = Q10base_cesm(2); %x(2): baseTuaE
tuaE_opCesm = cal_tuaE_op (Q10base_cesm, inputs_cesm); % calculate the estimated optimal tuaE
tuaE_scaler_cesm = cal_scaler(Q10base_cesm, inputs_cesm); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_cesm = cal_R2rmse(Q10base_cesm, inputs_cesm) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP6.mat','Q10_cesm','-append');
save('TuaEbase_CMIP6.mat','baseTuaE_cesm','-append')
save('tuaEop_CMIP6.mat','tuaE_opCesm','-append')
save('tuaE_scaler_CMIP6.mat','tuaE_scaler_cesm','-append')
save('tuaE_R2rmse_CMIP6.mat','tuaE_R2rmse_cesm','-append')


%% CNRM-ESM2-1
inputs_cnrm(1,:) =  tasCNRM_tmp(2:156,:);
inputs_cnrm(2,:) =  prCNRM_tmp(2:156,:);
inputs_cnrm(3,:) =  tuaEcnrm_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_cnrm,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_cnrm),[1,20],[],[],[],[],...
    [0,0],[5,max(inputs_cnrm(3,:))],@(x)myfunc_tuaE(x,inputs_cnrm));

Q10_cnrm = Q10base_cnrm (1); % x(1): Q10
baseTuaE_cnrm = Q10base_cnrm(2); %x(2): baseTuaE
tuaE_opCnrm = cal_tuaE_op (Q10base_cnrm, inputs_cnrm); % calculate the estimated optimal tuaE
tuaE_scaler_cnrm = cal_scaler(Q10base_cnrm, inputs_cnrm); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_cnrm = cal_R2rmse(Q10base_cnrm, inputs_cnrm) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP6.mat','Q10_cnrm','-append');
save('TuaEbase_CMIP6.mat','baseTuaE_cnrm','-append')
save('tuaEop_CMIP6.mat','tuaE_opCnrm','-append')
save('tuaE_scaler_CMIP6.mat','tuaE_scaler_cnrm','-append')
save('tuaE_R2rmse_CMIP6.mat','tuaE_R2rmse_cnrm','-append')

%% EC-Earth3-Veg
inputs_ec(1,:) =  tasEC_tmp(2:156,:);
inputs_ec(2,:) =  prEC_tmp(2:156,:);
inputs_ec(3,:) =  tuaEec_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_ec,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_ec),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_ec(3,:))],@(x)myfunc_tuaE(x,inputs_ec));

Q10_ec = Q10base_ec (1); % x(1): Q10
baseTuaE_ec = Q10base_ec(2); %x(2): baseTuaE
tuaE_opEc = cal_tuaE_op (Q10base_ec, inputs_ec); % calculate the estimated optimal tuaE
tuaE_scaler_ec = cal_scaler(Q10base_ec, inputs_ec); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_ec = cal_R2rmse(Q10base_ec, inputs_ec) % calculate the R2 and RMSE


save('TuaE_Q10_CMIP6.mat','Q10_ec','-append');
save('TuaEbase_CMIP6.mat','baseTuaE_ec','-append')
save('tuaEop_CMIP6.mat','tuaE_opEc','-append')
save('tuaE_scaler_CMIP6.mat','tuaE_scaler_ec','-append')
save('tuaE_R2rmse_CMIP6.mat','tuaE_R2rmse_ec','-append')

%% IPSL-CM6A-LR
inputs_ipsl(1,:) =  tasIPSL_tmp(2:156,:);
inputs_ipsl(2,:) =  prIPSL_tmp(2:156,:);
inputs_ipsl(3,:) =  tuaEipsl_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_ipsl,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_ipsl),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_ipsl(3,:))],@(x)myfunc_tuaE(x,inputs_ipsl));

Q10_ipsl = Q10base_ipsl (1); % x(1): Q10
baseTuaE_ipsl = Q10base_ipsl(2); %x(2): baseTuaE
tuaE_opIpsl = cal_tuaE_op (Q10base_ipsl, inputs_ipsl); % calculate the estimated optimal tuaE
tuaE_scaler_ipsl = cal_scaler(Q10base_ipsl, inputs_ipsl); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_ipsl = cal_R2rmse(Q10base_ipsl, inputs_ipsl) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP6.mat','Q10_ipsl','-append');
save('TuaEbase_CMIP6.mat','baseTuaE_ipsl','-append')
save('tuaEop_CMIP6.mat','tuaE_opIpsl','-append')
save('tuaE_scaler_CMIP6.mat','tuaE_scaler_ipsl','-append')
save('tuaE_R2rmse_CMIP6.mat','tuaE_R2rmse_ipsl','-append')

%% MIROC-ES2L
inputs_miroc(1,:) =  tasMIC_tmp(2:156,:);
inputs_miroc(2,:) =  prMIC_tmp(2:156,:);
inputs_miroc(3,:) =  tuaEmic_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_miroc,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_miroc),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_miroc(3,:))],@(x)myfunc_tuaE(x,inputs_miroc));

Q10_miroc = Q10base_miroc (1); % x(1): Q10
baseTuaE_miroc = Q10base_miroc(2); %x(2): baseTuaE
tuaE_opMiroc = cal_tuaE_op (Q10base_miroc, inputs_miroc); % calculate the estimated optimal tuaE
tuaE_scaler_miroc = cal_scaler(Q10base_miroc, inputs_miroc); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_miroc = cal_R2rmse(Q10base_miroc, inputs_miroc) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP6.mat','Q10_miroc','-append');
save('TuaEbase_CMIP6.mat','baseTuaE_miroc','-append')
save('tuaEop_CMIP6.mat','tuaE_opMiroc','-append')
save('tuaE_scaler_CMIP6.mat','tuaE_scaler_miroc','-append')
save('tuaE_R2rmse_CMIP6.mat','tuaE_R2rmse_miroc','-append')

%% MPI-ESM1-2-LR
inputs_mpi(1,:) =  tasMPI_tmp(2:156,:);
inputs_mpi(2,:) =  prMPI_tmp(2:156,:);
inputs_mpi(3,:) =  tuaEmpi_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_mpi,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_mpi),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_mpi(3,:))],@(x)myfunc_tuaE(x,inputs_mpi));

Q10_mpi = Q10base_mpi (1); % x(1): Q10
baseTuaE_mpi = Q10base_mpi(2); %x(2): baseTuaE
tuaE_opMpi = cal_tuaE_op (Q10base_mpi, inputs_mpi); % calculate the estimated optimal tuaE
tuaE_scaler_mpi = cal_scaler(Q10base_mpi, inputs_mpi); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_mpi = cal_R2rmse(Q10base_mpi, inputs_mpi) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP6.mat','Q10_mpi','-append');
save('TuaEbase_CMIP6.mat','baseTuaE_mpi','-append')
save('tuaEop_CMIP6.mat','tuaE_opMpi','-append')
save('tuaE_scaler_CMIP6.mat','tuaE_scaler_mpi','-append')
save('tuaE_R2rmse_CMIP6.mat','tuaE_R2rmse_mpi','-append')

%% NorESM2-LM
inputs_nor(1,:) =  tasNOR_tmp(2:156,:);
inputs_nor(2,:) =  prNOR_tmp(2:156,:);
inputs_nor(3,:) =  tuaEnor_tmp1m(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_nor,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_nor),[1,25],[],[],[],[],...
    [0,0],[5,max(inputs_nor(3,:))],@(x)myfunc_tuaE(x,inputs_nor));

Q10_nor = Q10base_nor (1); % x(1): Q10
baseTuaE_nor = Q10base_nor(2); %x(2): baseTuaE
tuaE_opNor = cal_tuaE_op (Q10base_nor, inputs_nor); % calculate the estimated optimal tuaE
tuaE_scaler_nor = cal_scaler(Q10base_nor, inputs_nor); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_nor = cal_R2rmse(Q10base_nor, inputs_nor) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP6.mat','Q10_nor','-append');
save('TuaEbase_CMIP6.mat','baseTuaE_nor','-append')
save('tuaEop_CMIP6.mat','tuaE_opNor','-append')
save('tuaE_scaler_CMIP6.mat','tuaE_scaler_nor','-append')
save('tuaE_R2rmse_CMIP6.mat','tuaE_R2rmse_nor','-append')

%% UKESM1-0-LL
inputs_uk(1,:) =  tasUK_tmp(2:156,:);
inputs_uk(2,:) =  prUK_tmp(2:156,:);
inputs_uk(3,:) =  tuaEuk_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_uk,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_uk),[1,25],[],[],[],[],...
    [0,0],[5,max(inputs_uk(3,:))],@(x)myfunc_tuaE(x,inputs_uk));

Q10_uk = Q10base_uk (1); % x(1): Q10
baseTuaE_uk = Q10base_uk(2); %x(2): baseTuaE
tuaE_opUk = cal_tuaE_op (Q10base_uk, inputs_uk); % calculate the estimated optimal tuaE
tuaE_scaler_uk = cal_scaler(Q10base_uk, inputs_uk); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_uk = cal_R2rmse(Q10base_uk, inputs_uk) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP6.mat','Q10_uk','-append');
save('TuaEbase_CMIP6.mat','baseTuaE_uk','-append')
save('tuaEop_CMIP6.mat','tuaE_opUk','-append')
save('tuaE_scaler_CMIP6.mat','tuaE_scaler_uk','-append')
save('tuaE_R2rmse_CMIP6.mat','tuaE_R2rmse_uk','-append')

%% Figure:
% 1:1 line: the modeled tuaE to that estimated from optimizing process
clear
clc
load('2tuaE_cmip6_tmp.mat')
load('tuaEop_CMIP6.mat')

tuaE6_all = [tuaEbcc_tmp(2:156),tuaEcan_tmp(2:156),tuaEcesm_tmp(2:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   tuaEuk_tmp(2:156),tuaEipsl_tmp(2:156),tuaEmic_tmp(2:156),...             % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   tuaEmpi_tmp(2:156),tuaEnor_tmp(2:156),...                                % MPI-ESM1-2-LR, NorESM2
   tuaEass_tmp(2:156),tuaEcnrm_tmp(2:156),tuaEec_tmp(2:156)];               % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg

% TuaE_op = basedTuaE * pr_scaler * tas_scaler
tuaEop6_all = [tuaE_opBccC; tuaE_opCan; tuaE_opCesm; tuaE_opUk;...
               tuaE_opIpsl; tuaE_opMiroc; tuaE_opMpi; tuaE_opNor;...
               tuaE_opASS; tuaE_opCnrm; tuaE_opEc;];

tuaEop6_all = tuaEop6_all';

mycolor6 = [255 0 0; 153 51 255; 237 176 33;...    %BCC-CSM2-MR, CanESM5, CESM2
           0 197 205; 0 205 0; 207 194 124;...     %UKESM1-0-LL, IPSL-CM6A-LR,MIROC-ES2L
           255 99 71; 65 105 255;...               %MPI-ESM1-2-LR,NorESM-LM
            0 0 0; 158 131 149; 119 136 153]./255; %ACCESS-ESM1-5, CNRM-ESM2-1  EC-Earth3-Veg  

clearvars -except tuaE6_all  tuaEop6_all  mycolor6

Figure
hold on              
for i=1:11
    plot(tuaE6_all(:,i),tuaEop6_all(:,i),'Marker','o','MarkerEdgeColor', mycolor6(i,:),...
        'MarkerSize',5,'LineStyle','none','LineWidth',1.2)
end
set(gca,'linewidth',1.2,'box','on');
set(gca,'XLim',[0 80]);
set(gca,'YLim',[0 80]);
set(gca,'Fontname','Arial','FontSize',12)
set(gca,'YTickLabelMode','auto');
set(gca,'XTickLabelMode','auto');

xlabel('Residence time simulated from models (year)','Fontname','Arial','FontSize',14)
ylabel('Residemce time reproduced (year)','Fontname','Arial','FontSize',14)

X=linspace(0,80,100);
Y = X;
plot(X,Y,'k:') 

leg6_str = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'}    
leg6 = legend(leg6_str)    
set(leg6,'color','none','EdgeColor','none','Fontname','Arial','Fontsize',10,'location','northwest')   
text(55, 10,'(b) CMIP6','Fontname','Arial','Fontsize',14)








