clear
clc

cd('file path of CMIP5 matData')
load('tuaE_cmip5_tmp.mat')
load('pr_cmip5_tmp.mat')
load('tas_cmip5_tmp.mat')

cd('file path to save results as .mat')
%% BCC-CSM1-1m
% put pre tas tuaE data into inputs
inputs_bcc(1,:) =  tasBCC_tmp(2:156,:);  % input1: near surface temperature
inputs_bcc(2,:) =  prBCC_tmp(2:156,:);   % input2: precipitation
inputs_bcc(3,:) =  tuaEbcc_tmp(2:156,:); % input3: residence time
% x1: Q10, x2: baseTuaE
[Q10base_bcc,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_bcc),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_bcc(3,:))],@(x)myfunc_tuaE(x,inputs_bcc));

Q10_bcc = Q10base_bcc (1);     % x(1): Q10
baseTuaE_bcc = Q10base_bcc(2); %x(2): baseTuaE
tuaE_opBcc = cal_tuaE_op (Q10base_bcc, inputs_bcc);    % calculate the estimated optimal tuaE
tuaE_scaler_bcc = cal_scaler(Q10base_bcc, inputs_bcc); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_bcc = cal_R2rmse(Q10base_bcc, inputs_bcc)  % calculate the R2 and RMSE

save('TuaE_Q10_CMIP5.mat','Q10_bcc');
save('TuaEbase_CMIP5.mat','baseTuaE_bcc')
save('tuaEop_CMIP5.mat','tuaE_opBcc')
save('tuaE_scaler_CMIP5.mat','tuaE_scaler_bcc')
save('tuaE_R2rmse_CMIP5.mat','tuaE_R2rmse_bcc')

%% BNU-ESM
inputs_bnu(1,:) = tasBNU_tmp(2:156,:);
inputs_bnu(2,:) =  prBNU_tmp(2:156,:);
inputs_bnu(3,:) =  tuaEbnu_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_bnu,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_bnu),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_bnu(3,:))],@(x)myfunc_tuaE(x,inputs_bnu));

Q10_bnu = Q10base_bnu (1); % x(1): Q10
baseTuaE_bnu = Q10base_bnu(2); %x(2): baseTuaE
tuaE_opBnu = cal_tuaE_op (Q10base_bnu, inputs_bnu); % calculate the estimated optimal tuaE
tuaE_scaler_bnu = cal_scaler(Q10base_bnu, inputs_bnu); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_bnu = cal_R2rmse(Q10base_bnu, inputs_bnu) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP5.mat','Q10_bnu','-append');
save('TuaEbase_CMIP5.mat','baseTuaE_bnu','-append')
save('tuaEop_CMIP5.mat','tuaE_opBnu','-append')
save('tuaE_scaler_CMIP5.mat','tuaE_scaler_bnu','-append')
save('tuaE_R2rmse_CMIP5.mat','tuaE_R2rmse_bnu','-append')

%% CanESM2
inputs_can(1,:) =  tasCAN_tmp(2:156,:);
inputs_can(2,:) =  prCAN_tmp(2:156,:);
inputs_can(3,:) =  tuaEcan_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_can,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_can),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_can(3,:))],@(x)myfunc_tuaE(x,inputs_can));

Q10_can = Q10base_can (1); % x(1): Q10
baseTuaE_can = Q10base_can(2); %x(2): baseTuaE
tuaE_opCan = cal_tuaE_op (Q10base_can, inputs_can); % calculate the estimated optimal tuaE
tuaE_scaler_can = cal_scaler(Q10base_can, inputs_can); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_can = cal_R2rmse(Q10base_can, inputs_can) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP5.mat','Q10_can','-append');
save('TuaEbase_CMIP5.mat','baseTuaE_can','-append')
save('tuaEop_CMIP5.mat','tuaE_opCan','-append')
save('tuaE_scaler_CMIP5.mat','tuaE_scaler_can','-append')
save('tuaE_R2rmse_CMIP5.mat','tuaE_R2rmse_can','-append')

%% CCSM4
inputs_ccsm(1,:) =  tasCCSM_tmp(2:156,:);
inputs_ccsm(2,:) =  prCCSM_tmp(2:156,:);
inputs_ccsm(3,:) =  tuaEccsm_tmp(2:156,:);

[Q10base_ccsm,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_ccsm),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_ccsm(3,:))],@(x)myfunc_tuaE(x,inputs_ccsm));

Q10_ccsm = Q10base_ccsm (1); % x(1): Q10
baseTuaE_ccsm = Q10base_ccsm(2); %x(2): baseTuaE
tuaE_opCcsm = cal_tuaE_op (Q10base_ccsm, inputs_ccsm); % calculate the estimated optimal tuaE
tuaE_scaler_ccsm = cal_scaler(Q10base_ccsm, inputs_ccsm); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_ccsm = cal_R2rmse(Q10base_ccsm, inputs_ccsm) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP5.mat','Q10_ccsm','-append');
save('TuaEbase_CMIP5.mat','baseTuaE_ccsm','-append')
save('tuaEop_CMIP5.mat','tuaE_opCcsm','-append')
save('tuaE_scaler_CMIP5.mat','tuaE_scaler_ccsm','-append')
save('tuaE_R2rmse_CMIP5.mat','tuaE_R2rmse_ccsm','-append')


%% GFDL-ESM2G (the modeled historical period: 1861-2005)
inputs_gf2G(1,:) =  tasGF_tmp(2:145,:);
inputs_gf2G(2,:) =  prGF_tmp(2:145,:);
inputs_gf2G(3,:) =  tuaEgf_tmp(2:145,:);

% x1: Q10, x2: baseTuaE
[Q10base_gf2G,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_gf2G),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_gf2G(3,:))],@(x)myfunc_tuaE(x,inputs_gf2G));

Q10_gf2G = Q10base_gf2G (1); % x(1): Q10
baseTuaE_gf2G = Q10base_gf2G(2); %x(2): baseTuaE
tuaE_opGf2G = cal_tuaE_op (Q10base_gf2G, inputs_gf2G); % calculate the estimated optimal tuaE
tuaE_scaler_gf2G = cal_scaler(Q10base_gf2G, inputs_gf2G); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_gf2G = cal_R2rmse(Q10base_gf2G, inputs_gf2G) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP5.mat','Q10_gf2G','-append');
save('TuaEbase_CMIP5.mat','baseTuaE_gf2G','-append')
save('tuaEop_CMIP5.mat','tuaE_opGf2G','-append')
save('tuaE_scaler_CMIP5.mat','tuaE_scaler_gf2G','-append')
save('tuaE_R2rmse_CMIP5.mat','tuaE_R2rmse_gf2G','-append')

%% HadGEM2-ES(the modeled historical period: 1860-2005)
inputs_hadES(1,:) =  tasHAD_tmp(2:146,:);
inputs_hadES(2,:) =  prHAD_tmp(2:146,:);
inputs_hadES(3,:) =  tuaEhad_tmp(2:146,:);

% x1: Q10, x2: baseTuaE
[Q10base_hadES,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_hadES),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_hadES(3,:))],@(x)myfunc_tuaE(x,inputs_hadES));

Q10_hadES = Q10base_hadES (1); % x(1): Q10
baseTuaE_hadES = Q10base_hadES(2); %x(2): baseTuaE
tuaE_opHadES = cal_tuaE_op (Q10base_hadES, inputs_hadES); % calculate the estimated optimal tuaE
tuaE_scaler_hadES = cal_scaler(Q10base_hadES, inputs_hadES); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_hadES = cal_R2rmse(Q10base_hadES, inputs_hadES) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP5.mat','Q10_hadES','-append');
save('TuaEbase_CMIP5.mat','baseTuaE_hadES','-append')
save('tuaEop_CMIP5.mat','tuaE_opHadES','-append')
save('tuaE_scaler_CMIP5.mat','tuaE_scaler_hadES','-append')
save('tuaE_R2rmse_CMIP5.mat','tuaE_R2rmse_hadES','-append')

%% IPSL-CM5A-MR
inputs_ipsl5A(1,:) =  tasIPSL_tmp(2:156,:);
inputs_ipsl5A(2,:) =  prIPSL_tmp(2:156,:);
inputs_ipsl5A(3,:) =  tuaEipsl_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_ipsl5A,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_ipsl5A),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_ipsl5A(3,:))],@(x)myfunc_tuaE(x,inputs_ipsl5A));

Q10_ipsl5A = Q10base_ipsl5A (1); % x(1): Q10
baseTuaE_ipsl5A = Q10base_ipsl5A(2); %x(2): baseTuaE
tuaE_opIpsl5A = cal_tuaE_op (Q10base_ipsl5A, inputs_ipsl5A); % calculate the estimated optimal tuaE
tuaE_scaler_ipsl5A = cal_scaler(Q10base_ipsl5A, inputs_ipsl5A); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_ipsl5A = cal_R2rmse(Q10base_ipsl5A, inputs_ipsl5A) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP5.mat','Q10_ipsl5A','-append');
save('TuaEbase_CMIP5.mat','baseTuaE_ipsl5A','-append')
save('tuaEop_CMIP5.mat','tuaE_opIpsl5A','-append')
save('tuaE_scaler_CMIP5.mat','tuaE_scaler_ipsl5A','-append')
save('tuaE_R2rmse_CMIP5.mat','tuaE_R2rmse_ipsl5A','-append')

%%  MIROC-ESM
inputs_miroc(1,:) =  tasMIROC_tmp(2:156,:);
inputs_miroc(2,:) =  prMIROC_tmp(2:156,:);
inputs_miroc(3,:) =  tuaEmiroc_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_miroc,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_miroc),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_miroc(3,:))],@(x)myfunc_tuaE(x,inputs_miroc));

Q10_miroc = Q10base_miroc (1); % x(1): Q10
baseTuaE_miroc = Q10base_miroc(2); %x(2): baseTuaE
tuaE_opMiroc = cal_tuaE_op (Q10base_miroc, inputs_miroc); % calculate the estimated optimal tuaE
tuaE_scaler_miroc = cal_scaler(Q10base_miroc, inputs_miroc); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_miroc = cal_R2rmse(Q10base_miroc, inputs_miroc) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP5.mat','Q10_miroc','-append');
save('TuaEbase_CMIP5.mat','baseTuaE_miroc','-append')
save('tuaEop_CMIP5.mat','tuaE_opMiroc','-append')
save('tuaE_scaler_CMIP5.mat','tuaE_scaler_miroc','-append')
save('tuaE_R2rmse_CMIP5.mat','tuaE_R2rmse_miroc','-append')


%%  MPI-ESM-MR
inputs_mpiM(1,:) =  tasMPI_tmp(2:156,:);
inputs_mpiM(2,:) =  prMPI_tmp(2:156,:);
inputs_mpiM(3,:) =  tuaEmpi_tmp(2:156,:);

% x1: Q10, x2: baseTuaE
[Q10base_mpiM,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_mpiM),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_mpiM(3,:))],@(x)myfunc_tuaE(x,inputs_mpiM));

Q10_mpiM = Q10base_mpiM (1); % x(1): Q10
baseTuaE_mpiM = Q10base_mpiM(2); %x(2): baseTuaE
tuaE_opMpiM = cal_tuaE_op (Q10base_mpiM, inputs_mpiM); % calculate the estimated optimal tuaE
tuaE_scaler_mpiM = cal_scaler(Q10base_mpiM, inputs_mpiM); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_mpiM = cal_R2rmse(Q10base_mpiM, inputs_mpiM) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP5.mat','Q10_mpiM','-append');
save('TuaEbase_CMIP5.mat','baseTuaE_mpiM','-append')
save('tuaEop_CMIP5.mat','tuaE_opMpiM','-append')
save('tuaE_scaler_CMIP5.mat','tuaE_scaler_mpiM','-append')
save('tuaE_R2rmse_CMIP5.mat','tuaE_R2rmse_mpiM','-append')

%%  MRI-ESM1(the modeled historical period: 1851-2005)
inputs_mri(1,:) =  tasMRI_tmp(2:155,:);
inputs_mri(2,:) =  prMRI_tmp(2:155,:);
inputs_mri(3,:) =  tuaEmri_tmp(2:155,:);

% x1: Q10, x2: baseTuaE
[Q10base_mri,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_mri),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_mri(3,:))],@(x)myfunc_tuaE(x,inputs_mri));

Q10_mri = Q10base_mri (1); % x(1): Q10
baseTuaE_mri = Q10base_mri(2); %x(2): baseTuaE
tuaE_opMri = cal_tuaE_op (Q10base_mri, inputs_mri); % calculate the estimated optimal tuaE
tuaE_scaler_mri = cal_scaler(Q10base_mri, inputs_mri); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_mri = cal_R2rmse(Q10base_mri, inputs_mri) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP5.mat','Q10_mri','-append');
save('TuaEbase_CMIP5.mat','baseTuaE_mri','-append')
save('tuaEop_CMIP5.mat','tuaE_opMri','-append')
save('tuaE_scaler_CMIP5.mat','tuaE_scaler_mri','-append')
save('tuaE_R2rmse_CMIP5.mat','tuaE_R2rmse_mri','-append')

%%  NorESM1-M
inputs_norM(1,:) =  tasNOR_tmp(2:156,:);
inputs_norM(2,:) =  prNOR_tmp(2:156,:);
inputs_norM(3,:) =  tuaEnor_tmp(2:156,:);
inputs_norM = double(inputs_norM);

% x1: Q10, x2: baseTuaE
[Q10base_norM,b,c]=fmincon(@(x)rmseR2_tuaE(x,inputs_norM),[1,20],[],[],[],[],...
    [0,0],[10,max(inputs_norM(3,:))],@(x)myfunc_tuaE(x,inputs_norM));

Q10_norM = Q10base_norM (1); % x(1): Q10
baseTuaE_norM = Q10base_norM(2); %x(2): baseTuaE
tuaE_opNorM = cal_tuaE_op (Q10base_norM, inputs_norM); % calculate the estimated optimal tuaE
tuaE_scaler_norM = cal_scaler(Q10base_norM, inputs_norM); % calculate the environmental scalers: scaler_tem, scaler_pre, T_scaler;
tuaE_R2rmse_norM = cal_R2rmse(Q10base_norM, inputs_norM) % calculate the R2 and RMSE

save('TuaE_Q10_CMIP5.mat','Q10_norM','-append');
save('TuaEbase_CMIP5.mat','baseTuaE_norM','-append')
save('tuaEop_CMIP5.mat','tuaE_opNorM','-append')
save('tuaE_scaler_CMIP5.mat','tuaE_scaler_norM','-append')
save('tuaE_R2rmse_CMIP5.mat','tuaE_R2rmse_norM','-append')


%% Figure:
% 1:1 line: the modeled tuaE to that estimated from optimizing process
clear
clc
load('tuaE_cmip5_tmp.mat')
load('tuaEop_CMIP5.mat')

NaNgf = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]';
tuaEgf_tmp = [NaNgf; tuaEgf_tmp];
tuaE_opGf2G = [NaNgf' tuaE_opGf2G];

NaNhad = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]';
tuaEhad_tmp = [NaNhad; tuaEhad_tmp];
tuaE_opHadES = [NaNhad' tuaE_opHadES];

tuaEmri_tmp = [NaN; tuaEmri_tmp];
tuaE_opMri = [NaN' tuaE_opMri];

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

mycolor5 = [255 0 0; 153 51 255; 237 176 33; ...%BCC CanESM2 CCSM4
            0 197 205; 0 205 0; 207 194 124;...   %Had IPSL Miroc
            255 99 71; 65 105 255;...  %mpi NorM
            0 0 0; 158 131 149; 119 136 153]./255;  %BNU, GFDL, mri         
clearvars -except tuaE5_all tuaE5_opt mycolor5   

figure 
hold on
for i=1:11
    plot(tuaE5_all(:,i),tuaE5_opt(:,i),'Marker','o','MarkerEdgeColor', mycolor5(i,:),...
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

leg5_str = {'BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'};    
leg5 = legend(leg5_str)    
set(leg5,'color','none','EdgeColor','none','Fontname','Arial','Fontsize',10,'location','northwest')
text(55,10,'(a) CMIP5','Fontname','Arial','Fontsize',14)


