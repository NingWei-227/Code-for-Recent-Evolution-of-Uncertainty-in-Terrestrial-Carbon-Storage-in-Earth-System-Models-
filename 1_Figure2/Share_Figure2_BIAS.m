% Fig.2 Evaluating model bias and inter-model spread in:
% Global terrestrial carbon storage, GPP and ¦ÓE for CMIP5 and CMIP6
% Observation-based estimates (at the global scale) include: 
%         cLand (Pg C):   cVeg: from Spawn et al., 2020;
%                         cSoil: estimated based on  Fan et al., 2020 Earth System Science Data, Table2
%         GPP(Pg C yr-1): estimated based on, FLUXCOM, MODIS17A2 and VPM
%         ¦ÓE (yr):       estimated from Fan et al., 2020 Earth System Science Data.
% the reference value used to estimate model bias is calculate as ensemble mean 
clear;clc

%% global land C storage estimated from Data, CMIP5 and CMIP6
% simulations from CMIP56
% load global terrestrial C stocks in soil and vegetation
load E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\4_MatData\cVeg_cSoil_CMIP56_native.mat
cVeg5_org_avg11 = nanmean(cVeg5_org_PG);
cSoil5_org_avg11 = nanmean(cSoil5_org_PG);
cVeg6_org_avg11 = nanmean(cVeg6_org_PG);
cSoil6_org_avg11 = nanmean(cSoil6_org_PG);
% calculate global terrestrial carbon storage for CMIP5 and CMIP6
cLand5_Pg_11M = cVeg5_org_avg11 + cSoil5_org_avg11;
cLand6_Pg_11M = cVeg6_org_avg11 + cSoil6_org_avg11; 

% estimates from data
% load cVeg data and calculate the global scale cVeg
load('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\cLand\Biomass\ORNL\Global_Maps_C_Density_2010_1763\cVeg05_KgCm2.mat')
area05 = csvread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\CABLEGridAreaM2.csv');              
Global_Veg_obs = cVeg_obs05kg.*area05.*10^(-12);   % convert unit into Pg C
Global_Veg_obs = nansum(Global_Veg_obs(:))         
Global_VegUN_obs = cVeg_obs_un05kg.*area05.*10^(-12); % uncertianty 
Global_VegUN_obs = nansum(Global_VegUN_obs(:))
% Benchmarking data from Fan et al., 2020 Earth System Science Data, Table2
% Csoil 0-1m
Global_soil_obs = [2195 2091 1332];  % unit: Pg C

% Global cLand estimated from data
cGlobal_obs = Global_soil_obs + Global_Veg_obs;
cGlobal_obs_sd1 = cGlobal_obs + Global_VegUN_obs;
cGlobal_obs_sd2 = cGlobal_obs - Global_VegUN_obs;
cGlobal_obs = [cGlobal_obs_sd1 cGlobal_obs cGlobal_obs_sd2];
cLand_obs05_Pg = nanmean(cGlobal_obs)

% Model-data difference
cLand_bias5 = (cLand5_Pg_11M - cLand_obs05_Pg)./cLand_obs05_Pg
cLand_bias6 = (cLand6_Pg_11M - cLand_obs05_Pg)./cLand_obs05_Pg

%nanmean(cLand_bias5)
%nanmean(cLand_bias6)

% get ride of useless variables to save memory 
clearvars -except cLand_bias5 cLand_bias6
%% global GPP estimated from Data, CMIP5 and CMIP6
% load CMIP5 data
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
% Load CMIP6 data
cd('H:\CMIP56_Csink\4_MatData\temporal_data\hist_ssp585')
load('2GPP_cmip6_tmp.mat')
GPP6_tmp = [GPPbcc_tmp(2:156),GPPcan_tmp(2:156),GPPcesm_tmp(2:156),...  % BCC_CSM2_MR, CanESM5 CESM2, CESM2 
   GPPuk_tmp(2:156),GPPipsl_tmp(2:156),GPPmic_tmp(2:156),...          % UKESM1_0_LL, IPSL_CM6A_LR, MIROC_ES2L  
   GPPmpi_tmp(2:156),GPPnor_tmp(2:156),...                           % MPI-ESM1-2-LR, NorESM2
   GPPass_tmp(2:156),GPPcnrm_tmp(2:156),GPPec_tmp(2:156)];            % ACCESS-ESM1-5, CNRM_ESM2_1, EC_Earth3_Veg
% 2001-2005 mean   
GPP5_end5_ag = nanmean(GPP5_tmp(151:155,:),1);
GPP6_end5_ag = nanmean(GPP6_tmp(151:155,:),1); 

% load Observational data
load E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure4\Obs_NPP_GPP_CUE\nppGppCUE.mat
area05 = csvread('E:\1_Mycase\3_CMIP56_Cland\2_Benchmark\1_obsData\GPP\FiveAnnualGPP2WeiN\FiveAnnualGPP2WeiN\CABLEGridAreaM2.csv');
GPP_obs05_Pg = GPP_obs05KG.* area05.*10^(-12);    % convert unit to PgC yr-1
GPP_obs5yr_Pg = nansum(GPP_obs05_Pg,1);  GPP_obs5yr_Pg = nansum(GPP_obs5yr_Pg,2); 
GPP_obs5yr_Pg = squeeze(GPP_obs5yr_Pg)           % global terrestrial GPP estimated from three data products
% ensemble mean of data products
GPP_obsPg = nanmean(GPP_obs5yr_Pg([1,5,6]));
% Model-data difference
GPP_bias5 = (GPP5_end5_ag - GPP_obsPg)./GPP_obsPg;
GPP_bias6 = (GPP6_end5_ag - GPP_obsPg)./GPP_obsPg;

%nanmean(GPP_bias5)
%nanmean(GPP_bias6)

clearvars -except cLand_bias5 cLand_bias6 ...
    GPP_bias5 GPP_bias6 GPP_obs5yr_Pg GPP5_end5_ag GPP6_end5_ag...
    CUE_obs05 NPP_obs05KG 

%% Global tauE estimated from Data, CMIP5 and CMIP6
% load processed data
load E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\NOY_circumpolar\NOYpolar_cLand.mat
load E:\1_Mycase\3_CMIP56_Cland\3_Cstorage_New\1_Figure\1_Main\Figure3_tauE\NOY_circumpolar\NOYpolar_GPP56.mat
% global tauE simulated from models
tuaE5_end5_ag = cLand5_M11./GPP5_end5_ag'
tuaE6_end5_ag = cLand6_M11./GPP6_end5_ag'

% Benchmarking data from Fan et al., 2020 Earth System Science Data, Table2
% cVeg (PgC)
Global_Veg_obs = [407 451 437 435];
% Csoil 0-1m (PgC)
Global_soil_obs = [1725 2195 2091 1332];
% global land carbon storage was estimated based on random combination (PgC)
cLand_obs = [];
for i = 1:4
    cLand  = Global_Veg_obs(i) + Global_soil_obs ;
    cLand_obs = [cLand_obs cLand];
    
end
% GPP used in Fan et al., 2020 (PgC yr-1)
GPP_glb_mean = 102;
% global terresrtial tuaE
obs_glb_tuaE = cLand_obs./GPP_glb_mean
obs_tuaE_mean = nanmean(obs_glb_tuaE)

% Model-data differnce
tauE_bias5 = (tuaE5_end5_ag - obs_tuaE_mean)./obs_tuaE_mean
tauE_bias6 = (tuaE6_end5_ag - obs_tuaE_mean)./obs_tuaE_mean

%nanmean(tauE_bias5)
%nanmean(tauE_bias6)
clearvars -except cLand_bias5 cLand_bias6 ...
    GPP_bias5 GPP_bias6 GPP_obs5yr_Pg GPP5_end5_ag GPP6_end5_ag...
    CUE_obs05 NPP_obs05KG ...
    tauE_bias5 tauE_bias6

%% Figures
% define table for storing data from CMIP5 and CMIP6
cmip5_3vars = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','cLand','GPP','tauE','Nlimitation'});
cmip6_3vars = table('Size',[11 5],'VariableType',{'string','double','double','double','double'},...
    'VariableName',{'Moldes','cLand','GPP','tauE','Nlimitation'});
% CMIP5 and CMIP6 models used in this study
Models5 = {'BCC-CSM1-1m','CanESM2','CCSM4',...
        'HadGEM2-ES','IPSL-CM5A-MR',...
        'MIROC-ESM','MPI-ESM-MR','NorESM1-M',...
        'BNU-ESM','GFDL-ESM2G','MRI-ESM1'};
Models6 = {'BCC-CSM2-MR', 'CanESM5', 'CESM2', ...
        'UKESM1-0-LL', 'IPSL-CM6A-LR', 'MIROC-ES2L',...
        'MPI-ESM1-2-LR', 'NorESM2-LM',...
        'ACCESS-ESM1-5', 'CNRM-ESM2-1','EC-Earth3-Veg'};
% put data into Table 
cmip5_3vars.Moldes = Models5';
cmip6_3vars.Moldes = Models6';
cmip5_3vars.cLand = cLand_bias5';
cmip6_3vars.cLand = cLand_bias6';
cmip5_3vars.GPP = GPP_bias5';
cmip6_3vars.GPP = GPP_bias6';
cmip5_3vars.tauE =  tauE_bias5;
cmip6_3vars.tauE =  tauE_bias6;
NL5 = [0 0 1 0 0 0 0 1 0 0 0]';  % 0, model without explicit nutrient limitation
NL6 = [0 0 1 1 1 1 1 1  1 0 1]'; % 1, model with explicit nutrient limitation
cmip5_3vars.Nlimitation = NL5;
cmip6_3vars.Nlimitation = NL6;

% RGB color and creat the figure window 
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
figure 
set(gcf,'position',[500,200,885,416])

%% Panel(a) for CMIP5
cmip56(1) = subtightplot(1,2,1,[0.03 0.02],[0.2 0.015],[0.08 0.01])
hold on
set(gca,'Fontname','Arial','FontSize',11)
set(gca,'linewidth',1,'box','off')
h5_X = raincloud_plot(cLand_bias5,'box_on', 0,...
                     'color', [0.97,0.63,0.59], 'cloud_edge_col', [1.00,0.47,0.41],...
                     'bandwidth', 12,'density_type', 'ks','alpha', 0.5);
set(h5_X{2},'Marker','none')                                  
h5_GPP = raincloud_plot(GPP_bias5,'box_on', 0, ...
    'color', [0.80,1.00,0.47], 'cloud_edge_col', [0.56,0.80,0.16],...
    'bandwidth', 12,'density_type', 'ks','alpha', 0.5);
set(h5_GPP{2},'Marker','none')
h5_tauE = raincloud_plot(tauE_bias5,'box_on', 0,...
    'color', [0.36,0.92,0.82], 'cloud_edge_col', [0.36,0.92,0.82],...
    'bandwidth', 12,'density_type', 'ks','alpha', 0.2);
set(h5_tauE{2},'Marker','none')
set(gca, 'YLim',[-8 4],'XLim',[-1.2 1.2]);
plot([0 0], [-8 4],'k--','LineWidth',1.5)
plot( [-1.5 1.5],[0 0],'-','LineWidth',1.8,'color',[0.65 0.65 0.65])

% model bias for each model
avgcLand5 = nanmean(cmip5_3vars.cLand); 
avgGPP5 = nanmean(cmip5_3vars.GPP); 
avgtauE5 = nanmean(cmip5_3vars.tauE); 
for i = 1:11
    
    if cmip5_3vars.Nlimitation(i) == 0 
        leg5_X(i)= plot(cmip5_3vars.cLand(i),-2,'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg5_GPP(i)= plot(cmip5_3vars.GPP(i),-4.5,'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg5_tauE(i)= plot(cmip5_3vars.tauE(i),-6.5,'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
    
    else 
        leg5_X(i)= plot(cmip5_3vars.cLand(i),-2,'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg5_GPP(i)= plot(cmip5_3vars.GPP(i),-4.5,'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg5_tauE(i)= plot(cmip5_3vars.tauE(i),-6.5,'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
    end
    
end

H_pa5_cLand = patch([min(cmip5_3vars.cLand),max(cmip5_3vars.cLand),max(cmip5_3vars.cLand),min(cmip5_3vars.cLand)],[-2.5,-2.5,-1.5,-1.5],cb(4,:));
set(H_pa5_cLand,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
line([avgcLand5 avgcLand5],[-2.5 -1.5],'color','k','linewidth',2);

H_pa5_GPP = patch([min(cmip5_3vars.GPP),max(cmip5_3vars.GPP),max(cmip5_3vars.GPP),min(cmip5_3vars.GPP)],[-5,-5,-4,-4],[0.56,0.80,0.16]);
set(H_pa5_GPP,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
line([avgGPP5 avgGPP5],[-5 -4],'color','k','linewidth',2);

H_pa5_tauE = patch([min(cmip5_3vars.tauE),max(cmip5_3vars.tauE),max(cmip5_3vars.tauE),min(cmip5_3vars.tauE)],[-7,-7,-6,-6],[0.00,0.77,0.80]);
set(H_pa5_tauE,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.35)
line([avgtauE5 avgtauE5],[-7 -6],'color','k','linewidth',2);

yticks([-6 -4 -2]);
yticklabels({'\tau_E','GPP','cLand'});
Ylabel = gca
set(Ylabel.YAxis,'FontSize',14,'FontWeight','bold')
set(gca,'Fontname','Arial','FontSize',12)
plot([-1.2, 1.2],[4,4],'k-','LineWidth',1)
text(-1.1,3.3,'(a) CMIP5','Fontname','Arial','FontSize',13)

annotation('arrow',[0.304 0.5289],[0.198 0.198],'HeadStyle','plain',...
    'HeadWidth',10,'HeadLength',9,...
    'LineWidth',4,'Color',[0.97,0.69,0.69])
annotation('arrow',[0.304 0.0805],[0.198 0.198],'HeadStyle','plain',...
    'HeadWidth',10,'HeadLength',9,...
    'LineWidth',4,'Color',[0.48,0.77,0.87])
text(0.7151,-7.53,'Positive bias','Fontname','Arial','FontSize',10,'Color',[0.99,0.03,0.03])
text(-1.1537,-7.5,'Negative bias','Fontname','Arial','FontSize',10,'Color',[0,0,1])
annotation('arrow',[0.2758 0.3038],[0.64 0.64],'HeadStyle','vback1',...
    'HeadWidth',10,'HeadLength',7,...
    'LineWidth',2,'Color',[0,0,0])
text(-0.2681,-0.8853,'Bias','Fontname','Arial','FontSize',12,'Color','r')
annotation('arrow',[0.2758 0.3978],[0.55 0.55],'HeadStyle','vback1',...
    'HeadWidth',10,'HeadLength',7,...
    'LineWidth',2,'Color',[0.89,0.01,0.01])
annotation('arrow',[0.2758 0.2008 ],[0.55 0.55],'HeadStyle','vback1',...
    'HeadWidth',10,'HeadLength',7,...
    'LineWidth',2,'Color',[0.89,0.01,0.01])
text(-0.4175,-3.0029,'Variation','Fontname','Arial','FontSize',12)
tt = text(0,-9.72,'$$\frac{model - Obs}{Obs}$$',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',13)  
set(tt,'Interpreter','latex');
plot([-1.2, -1.2],[-8,4],'k-','LineWidth',1)
plot([1.2, 1.2],[-8,4],'k-','LineWidth',1)


%% Panel(b) for CMIP6
cmip56(2) = subtightplot(1,2,2,[0.03 0.02],[0.2 0.015],[0.08 0.01])
hold on
set(gca,'Fontname','Arial','FontSize',11)
set(gca,'linewidth',1,'box','off')
set(gca,'Ycolor','k')
h6_X = raincloud_plot(cLand_bias6,'box_on', 0, ...
    'color', [0.97,0.63,0.59], 'cloud_edge_col', [1.00,0.47,0.41],...
    'bandwidth', 12,'density_type', 'ks','alpha', 0.5);
set(h6_X{2},'Marker','none') 
h6_GPP = raincloud_plot(GPP_bias6,'box_on', 0,...
    'color', [0.80,1.00,0.47], 'cloud_edge_col', [0.56,0.80,0.16],...
    'bandwidth', 12,'density_type', 'ks','alpha', 0.5);
set(h6_GPP{2},'Marker','none')
h6_tauE = raincloud_plot(tauE_bias6,'box_on', 0,...
    'color', [0.36,0.92,0.82], 'cloud_edge_col', [0.36,0.92,0.82],...
    'bandwidth', 12,'density_type', 'ks','alpha', 0.2);
set(h6_tauE{2},'Marker','none')
plot([0 0], [-8 4],'k--','LineWidth',1.5)
plot( [-1.5 1.5],[0 0],'-','LineWidth',1.8,'color',[0.65 0.65 0.65])

avgcLand6 = nanmean(cmip6_3vars.cLand); 
avgGPP6 = nanmean(cmip6_3vars.GPP); 
avgtauE6 = nanmean(cmip6_3vars.tauE); 
for i = 1:11
    
    if cmip6_3vars.Nlimitation(i) == 0 
        leg5_X(i)= plot(cmip6_3vars.cLand(i),-2,'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg5_GPP(i)= plot(cmip6_3vars.GPP(i),-4,'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg5_tauE(i)= plot(cmip6_3vars.tauE(i),-6,'Marker','o', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
    
    else 
        leg5_X(i)= plot(cmip6_3vars.cLand(i),-2,'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg5_GPP(i)= plot(cmip6_3vars.GPP(i),-4,'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
        leg5_tauE(i)= plot(cmip6_3vars.tauE(i),-6,'Marker','>', 'MarkerEdgeColor', [0.6 0.6 0.6],...
        'MarkerSize',9,'LineStyle','none','LineWidth',1.5)
    end
    
end
set(gca, 'YLim',[-8 4],'XLim',[-1.2 1.2]);

H_pa5_cLand = patch([min(cmip6_3vars.cLand),max(cmip6_3vars.cLand),max(cmip6_3vars.cLand),min(cmip6_3vars.cLand)],[-2.5,-2.5,-1.5,-1.5],cb(4,:));
set(H_pa5_cLand,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
line([avgcLand6 avgcLand6],[-2.5 -1.5],'color','k','linewidth',2);

H_pa5_GPP = patch([min(cmip6_3vars.GPP),max(cmip6_3vars.GPP),max(cmip6_3vars.GPP),min(cmip6_3vars.GPP)],[-4.5,-4.5,-3.5,-3.5],[0.56,0.80,0.16]);
set(H_pa5_GPP,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
line([avgGPP6 avgGPP6],[-4.5 -3.5],'color','k','linewidth',2);

H_pa5_tauE = patch([min(cmip6_3vars.tauE),max(cmip6_3vars.tauE),max(cmip6_3vars.tauE),min(cmip6_3vars.tauE)],[-6.5,-6.5,-5.5,-5.5],[0.00,0.77,0.80]);
set(H_pa5_tauE,'EdgeColor','none','EdgeAlpha',0.5,'FaceAlpha',0.3)
line([avgtauE6 avgtauE6],[-6.5 -5.5],'color','k','linewidth',2);

yticks([-6 -4 -2]);
yticklabels([]);
plot([-1.2, -1.2],[-8,4],'k-','LineWidth',1)
plot([1.2, 1.2],[-8,4],'k-','LineWidth',1)
plot([-1.2, 1.2],[4,4],'k-','LineWidth',1)
text(-1.1,3.3,'(b) CMIP6','Fontname','Arial','FontSize',13)

annotation('arrow',[0.7667 0.9936],[0.198 0.198],'HeadStyle','plain',...
    'HeadWidth',10,'HeadLength',9,...
    'LineWidth',4,'Color',[0.97,0.69,0.69])
annotation('arrow',[0.7667 0.5452],[0.198 0.198],'HeadStyle','plain',...
    'HeadWidth',10,'HeadLength',9,...
    'LineWidth',4,'Color',[0.48,0.77,0.87])
text(0.7151,-7.53,'Positive bias','Fontname','Arial','FontSize',10,'Color',[0.99,0.03,0.03])
text(-1.1537,-7.5,'Negative bias','Fontname','Arial','FontSize',10,'Color',[0,0,1])
tt = text(0,-9.72,'$$\frac{model - Obs}{Obs}$$',...
    'HorizontalAlignment','center',...
        'FontName','Arial','FontSize',13)  
set(tt,'Interpreter','latex');





