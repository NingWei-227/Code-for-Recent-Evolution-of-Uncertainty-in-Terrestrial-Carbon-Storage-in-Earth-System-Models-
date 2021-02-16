%% terrestrial C sink simulated by CMIP5 model ensemble
%  period: 1850-2005 except for GFDL-ESM2G(1861-2005) HadGEM2(1860-2005) and MRI-ESM1(1851-2005)
%  unit: PgC
%  CLand = CVeg + CSoil + CLitter expcet for GFDL and HadGEM2 (cLand=cVeg+cSoil): processed in CDO   
%  Csink estimated as the absolute change in CLand(yr) relative to the CLand(1850)
%  variables used in analysis: npp, gpp, cLand, sftlf
%  climate data: precipitation (pr) and surface air temperature (tas)

%% bcc-csm1-1-m
clear
clc
file2 = 'load file path'

pr_BCC = ncread([file2, 'file name of the .nc dataset'],'pr');     % unit: Kg m-2 m-1 = mm
tas_BCC = ncread([file2, 'file name of the .nc dataset'],'tas');  % unit: K

areacell = ncread([file2,'file name of the .nc dataset'],'areacella');              % unit:m2
sftlf_BCC = ncread([file2,'file name of the .nc dataset'],'sftlf');                 % unit: %                
cellarea_BCC = areacell.* sftlf_BCC*10^(-2);                         %unit: m2

nppBCC_3D = ncread([file2,'file name of the .nc dataset'],'npp');     % unit: KgC m-2 s-1
gppBCC_3D = ncread([file2,'file name of the .nc dataset'],'gpp');     % unit: KgC m-2 s-1
cLand_BCC = ncread([file2,'file name of the .nc dataset'],'cLand');   % unit: KgC m-2

nppBCC_3D(nppBCC_3D<0) = NaN; 
gppBCC_3D(gppBCC_3D<0) = NaN; 
cLand_BCC(cLand_BCC<0) = NaN;

nppBCC_3Dpg = nppBCC_3D.* cellarea_BCC.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppBCC_3Dpg = gppBCC_3D.* cellarea_BCC.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLbcc_3Dpg = cLand_BCC.* cellarea_BCC.* 10^(-12);  % unit: PgC

NPPbcc_tmp = nansum(nppBCC_3Dpg,1); NPPbcc_tmp = nansum(NPPbcc_tmp,2); 
NPPbcc_tmp = squeeze(NPPbcc_tmp);

GPPbcc_tmp = nansum(gppBCC_3Dpg,1);GPPbcc_tmp = nansum(GPPbcc_tmp,2);
GPPbcc_tmp = squeeze(GPPbcc_tmp);

CUEbcc_tmp = NPPbcc_tmp./GPPbcc_tmp;

CLbcc_tmp = nansum(CLbcc_3Dpg,1); CLbcc_tmp = nansum(CLbcc_tmp,2);
CLbcc_tmp = squeeze(CLbcc_tmp);

prBCC_tmp = nanmean(pr_BCC_3D,1); prBCC_tmp = nanmean(prBCC_tmp,2);
prBCC_tmp = squeeze(prBCC_tmp);

tasBCC_tmp = nanmean(tas_BCC_3D,1); tasBCC_tmp = nanmean(tasBCC_tmp,2);
tasBCC_tmp = squeeze(tasBCC_tmp);

% estimate terrestrial C residence time
nyear = 156    % from 1850 to 2005
for yr=1:nyear-1
    NetCbcc_tmp(yr+1,1) = CLbcc_tmp(yr+1) - CLbcc_tmp(yr);
    tuaEbcc_tmp(yr+1,1) = CLbcc_tmp(yr+1)./(NPPbcc_tmp(yr+1) - NetCbcc_tmp(yr+1));
end

XcBcc_tmp = tuaEbcc_tmp.*NPPbcc_tmp;
XpBcc_tmp = XcBcc_tmp - CLbcc_tmp;

% save data as .mat
save('file_path\NPP_cmip5_tmp.mat','NPPbcc_tmp')
save('file_path\GPP_cmip5_tmp.mat','GPPbcc_tmp')
save('file_path\CUE_cmip5_tmp.mat','CUEbcc_tmp')
save('file_path\cLand_cmip5_tmp.mat','CLbcc_tmp')
save('file_path\tuaE_cmip5_tmp.mat','tuaEbcc_tmp')
save('file_path\Cnet_cmip5_tmp.mat','NetCbcc_tmp')
save('file_path\Xc_cmip5_tmp.mat','XcBcc_tmp')
save('file_path\Xp_cmip5_tmp.mat','XpBcc_tmp')
save('file_path\pr_cmip5_tmp.mat','prBCC_tmp')
save('file_path\tas_cmip5_tmp.mat','tasBCC_tmp')


%% BNU-ESM
clear
clc 

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                % unit:m2
sftlf_BNU = ncread([file2,'file name of the .nc dataset'],'sftlf');                   % unit: %                
cellarea_BNU = areacell.* sftlf_BNU*10^(-2);                                       % unit: m2

nppBNU_3D = ncread([file2, 'file name of the .nc dataset'],'npp');        % unit: KgC m-2 s-1
gppBNU_3D = ncread([file2, 'file name of the .nc dataset'],'gpp');        % unit: KgC m-2 s-1
cLand_BNU = ncread([file2, 'file name of the .nc dataset'],'cLand');  % unit: KgC m-2
pr_BNU = ncread([file2, 'file name of the .nc dataset'],'pr');        % unit: Kg m-2 m-1 = mm
tas_BNU = ncread([file2, 'file name of the .nc dataset'],'tas');      % unit: K

nppBNU_3D(nppBNU_3D<0) = NaN; 
gppBNU_3D(gppBNU_3D<0) = NaN; 
cLand_BNU(cLand_BNU<0) = NaN; 

nppBNU_3Dpg = nppBNU_3D.* cellarea_BNU.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppBNU_3Dpg = gppBNU_3D.* cellarea_BNU.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLbnu_3Dpg = cLand_BNU.* cellarea_BNU.* 10^(-12);  % unit: PgC

NPPbnu_tmp = nansum(nppBNU_3Dpg,1); NPPbnu_tmp = nansum(NPPbnu_tmp,2); 
NPPbnu_tmp = squeeze(NPPbnu_tmp);

GPPbnu_tmp = nansum(gppBNU_3Dpg,1);GPPbnu_tmp = nansum(GPPbnu_tmp,2);
GPPbnu_tmp = squeeze(GPPbnu_tmp);

CUEbnu_tmp = NPPbnu_tmp./GPPbnu_tmp;

CLbnu_tmp = nansum(CLbnu_3Dpg,1); CLbnu_tmp = nansum(CLbnu_tmp,2);
CLbnu_tmp = squeeze(CLbnu_tmp);

prBNU_tmp = nanmean(pr_BNU_3D,1); prBNU_tmp = nanmean(prBNU_tmp,2);
prBNU_tmp = squeeze(prBNU_tmp);

tasBNU_tmp = nanmean(tas_BNU_3D,1); tasBNU_tmp = nanmean(tasBNU_tmp,2);
tasBNU_tmp = squeeze(tasBNU_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCbnu_tmp(yr+1,1) = CLbnu_tmp(yr+1) - CLbnu_tmp(yr);
    tuaEbnu_tmp(yr+1,1) = CLbnu_tmp(yr+1)./(NPPbnu_tmp(yr+1) - NetCbnu_tmp(yr+1));
end

XcBnu_tmp = tuaEbnu_tmp.*NPPbnu_tmp;
XpBnu_tmp = XcBnu_tmp - CLbnu_tmp;

% save data as .mat
save('file_path\NPP_cmip5_tmp.mat','NPPbnu_tmp','-append')
save('file_path\GPP_cmip5_tmp.mat','GPPbnu_tmp','-append')
save('file_path\CUE_cmip5_tmp.mat','CUEbnu_tmp','-append')
save('file_path\cLand_cmip5_tmp.mat','CLbnu_tmp','-append')
save('file_path\tuaE_cmip5_tmp.mat','tuaEbnu_tmp','-append')
save('file_path\Cnet_cmip5_tmp.mat','NetCbnu_tmp','-append')
save('file_path\Xc_cmip5_tmp.mat','XcBnu_tmp','-append')
save('file_path\Xp_cmip5_tmp.mat','XpBnu_tmp','-append')
save('file_path\pr_cmip5_tmp.mat','prBNU_tmp','-append')
save('file_path\tas_cmip5_tmp.mat','tasBNU_tmp','-append')
%% CanESM2
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');            % unit:m2
sftlf_CAN = ncread([file2,'file name of the .nc dataset'],'sftlf');               % unit: %                
cellarea_CAN = areacell.* sftlf_CAN*10^(-2);                       % unit: m2

nppCAN_3D = ncread([file2, 'file name of the .nc dataset'],'npp');        % unit: KgC m-2 s-1
gppCAN_3D = ncread([file2, 'file name of the .nc dataset'],'gpp');        % unit: KgC m-2 s-1
cLand_CAN = ncread([file2, 'file name of the .nc dataset'],'cLand');      % unit: KgC m-2
pr_CAN = ncread([file2, 'file name of the .nc dataset'],'pr');            % unit: Kg m-2 m-1 = mm
tas_CAN = ncread([file2, 'file name of the .nc dataset'],'tas');          % unit: K

nppCAN_3D(nppCAN_3D<0) = NaN; 
gppCAN_3D(gppCAN_3D<0) = NaN; 
cLand_CAN(cLand_CAN<0) = NaN; 

nppCAN_3Dpg = nppCAN_3D.* cellarea_CAN.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppCAN_3Dpg = gppCAN_3D.* cellarea_CAN.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLcan_3Dpg = cLand_CAN.* cellarea_CAN.* 10^(-12);                % unit: PgC

NPPcan_tmp = nansum(nppCAN_3Dpg,1); NPPcan_tmp = nansum(NPPcan_tmp,2); 
NPPcan_tmp = squeeze(NPPcan_tmp);

GPPcan_tmp = nansum(gppCAN_3Dpg,1);GPPcan_tmp = nansum(GPPcan_tmp,2);
GPPcan_tmp = squeeze(GPPcan_tmp);

CUEcan_tmp = NPPcan_tmp./GPPcan_tmp;

CLcan_tmp = nansum(CLcan_3Dpg,1); CLcan_tmp = nansum(CLcan_tmp,2);
CLcan_tmp = squeeze(CLcan_tmp);

prCAN_tmp = nanmean(pr_CAN_3D,1); prCAN_tmp = nanmean(prCAN_tmp,2);
prCAN_tmp = squeeze(prCAN_tmp);

tasCAN_tmp = nanmean(tas_CAN_3D,1); tasCAN_tmp = nanmean(tasCAN_tmp,2);
tasCAN_tmp = squeeze(tasCAN_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCcan_tmp(yr+1,1) = CLcan_tmp(yr+1) - CLcan_tmp(yr);
    tuaEcan_tmp(yr+1,1) = CLcan_tmp(yr+1)./(NPPcan_tmp(yr+1) - NetCcan_tmp(yr+1));
end

XcCan_tmp = tuaEcan_tmp.*NPPcan_tmp;
XpCan_tmp = XcCan_tmp - CLcan_tmp;

% save temporal data as .mat
save('file_path\NPP_cmip5_tmp.mat','NPPcan_tmp','-append')
save('file_path\GPP_cmip5_tmp.mat','GPPcan_tmp','-append')
save('file_path\CUE_cmip5_tmp.mat','CUEcan_tmp','-append')
save('file_path\cLand_cmip5_tmp.mat','CLcan_tmp','-append')
save('file_path\tuaE_cmip5_tmp.mat','tuaEcan_tmp','-append')
save('file_path\Cnet_cmip5_tmp.mat','NetCcan_tmp','-append')
save('file_path\Xc_cmip5_tmp.mat','XcCan_tmp','-append')
save('file_path\Xp_cmip5_tmp.mat','XpCan_tmp','-append')
save('file_path\pr_cmip5_tmp.mat','prCAN_tmp','-append')
save('file_path\tas_cmip5_tmp.mat','tasCAN_tmp','-append')


%% CCM4
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');             % unit:m2
sftlf_CCSM = ncread([file2,'file name of the .nc dataset'],'sftlf');               % unit: %                
cellarea_CCSM = areacell.* sftlf_CCSM*10^(-2);                                     % unit: m2
nppCCSM_3D = ncread([file2, 'file name of the .nc dataset'],'npp');       % unit: KgC m-2 s-1
gppCCSM_3D = ncread([file2, 'file name of the .nc dataset'],'gpp');       % unit: KgC m-2 s-1
cLand_CCSM = ncread([file2, 'file name of the .nc dataset'],'cLand');     % unit: KgC m-2
pr_CCSM = ncread([file2, 'file name of the .nc dataset'],'pr');           % unit: Kg m-2 m-1 = mm
tas_CCSM = ncread([file2, 'file name of the .nc dataset'],'tas');         % unit: KgC m-2

nppCCSM_3D(nppCCSM_3D<0) = NaN; %nppCCSM_3D(nppCCSM_3D>10000) = NaN;
gppCCSM_3D(gppCCSM_3D<0) = NaN; %gppCCSM_3D(gppCCSM_3D>10000) = NaN;
cLand_CCSM(cLand_CCSM<0) = NaN; %cLand_CCSM(cLand_CCSM>100000) = NaN; 

nppCCSM_3Dpg = nppCCSM_3D.* cellarea_CCSM.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppCCSM_3Dpg = gppCCSM_3D.* cellarea_CCSM.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLccsm_3Dpg = cLand_CCSM.* cellarea_CCSM.* 10^(-12);  % unit: PgC

nppCCSM_3Dpg(nppCCSM_3Dpg>10000) = NaN;
gppCCSM_3Dpg(gppCCSM_3Dpg>10000) = NaN;
CLccsm_3Dpg(CLccsm_3Dpg>=9999) = NaN;

NPPccsm_tmp = nansum(nppCCSM_3Dpg,1); NPPccsm_tmp = nansum(NPPccsm_tmp,2); 
NPPccsm_tmp = squeeze(NPPccsm_tmp);

GPPccsm_tmp = nansum(gppCCSM_3Dpg,1);GPPccsm_tmp = nansum(GPPccsm_tmp,2);
GPPccsm_tmp = squeeze(GPPccsm_tmp);

CUEccsm_tmp = NPPccsm_tmp./GPPccsm_tmp;

CLccsm_tmp = nansum(CLccsm_3Dpg,1); CLccsm_tmp = nansum(CLccsm_tmp,2);
CLccsm_tmp = squeeze(CLccsm_tmp);

prCCSM_tmp = nanmean(pr_CCSM_3D,1); prCCSM_tmp = nanmean(prCCSM_tmp,2);
prCCSM_tmp = squeeze(prCCSM_tmp);

tasCCSM_tmp = nanmean(tas_CCSM_3D,1); tasCCSM_tmp = nanmean(tasCCSM_tmp,2);
tasCCSM_tmp = squeeze(tasCCSM_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCccsm_tmp(yr+1,1) = CLccsm_tmp(yr+1) - CLccsm_tmp(yr);
    tuaEccsm_tmp(yr+1,1) = CLccsm_tmp(yr+1)./(NPPccsm_tmp(yr+1) - NetCccsm_tmp(yr+1));
end

XcCcsm_tmp = tuaEccsm_tmp.*NPPccsm_tmp;
XpCcsm_tmp = XcCcsm_tmp - CLccsm_tmp;

% save temporal data
save('file_path\NPP_cmip5_tmp.mat','NPPccsm_tmp','-append')
save('file_path\GPP_cmip5_tmp.mat','GPPccsm_tmp','-append')
save('file_path\CUE_cmip5_tmp.mat','CUEccsm_tmp','-append')
save('file_path\cLand_cmip5_tmp.mat','CLccsm_tmp','-append')
save('file_path\tuaE_cmip5_tmp.mat','tuaEccsm_tmp','-append')
save('file_path\Cnet_cmip5_tmp.mat','NetCccsm_tmp','-append')
save('file_path\Xc_cmip5_tmp.mat','XcCcsm_tmp','-append')
save('file_path\Xp_cmip5_tmp.mat','XpCcsm_tmp','-append')
save('file_path\pr_cmip5_tmp.mat','prCCSM_tmp','-append')
save('file_path\tas_cmip5_tmp.mat','tasCCSM_tmp','-append')


%% GFDL-ESM2G
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                % unit:m2
sftlf_GF = ncread([file2,'file name of the .nc dataset'],'sftlf');                   % unit: %                
cellarea_GF = areacell.* sftlf_GF*10^(-2);                                       % unit: m2

nppGF_3D = ncread([file2, 'file name of the .nc dataset'],'npp');        % unit: KgC m-2 s-1
gppGF_3D = ncread([file2, 'file name of the .nc dataset'],'gpp');        % unit: KgC m-2 s-1
cLand_GF = ncread([file2, 'file name of the .nc dataset'],'cLand');  % unit: KgC m-2

pr_GF = ncread([file2, 'file name of the .nc dataset'],'pr');     % unit: Kg m-2 m-1 = mm
tas_GF = ncread([file2, 'file name of the .nc dataset'],'tas');  % unit: K

nppGF_3D(nppGF_3D<0) = NaN; 
gppGF_3D(gppGF_3D<0) = NaN; 
cLand_GF(cLand_GF<0) = NaN; 

nppGF_3Dpg = nppGF_3D.* cellarea_GF.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppGF_3Dpg = gppGF_3D.* cellarea_GF.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLgf_3Dpg = cLand_GF.* cellarea_GF.* 10^(-12);  % unit: PgC

NPPgf_tmp = nansum(nppGF_3Dpg,1); NPPgf_tmp = nansum(NPPgf_tmp,2); 
NPPgf_tmp = squeeze(NPPgf_tmp);

GPPgf_tmp = nansum(gppGF_3Dpg,1);GPPgf_tmp = nansum(GPPgf_tmp,2);
GPPgf_tmp = squeeze(GPPgf_tmp);

CUEgf_tmp = NPPgf_tmp./GPPgf_tmp;

CLgf_tmp = nansum(CLgf_3Dpg,1); CLgf_tmp = nansum(CLgf_tmp,2);
CLgf_tmp = squeeze(CLgf_tmp);

prGF_tmp = nanmean(pr_GF_3D,1); prGF_tmp = nanmean(prGF_tmp,2);
prGF_tmp = squeeze(prGF_tmp);

tasGF_tmp = nanmean(tas_GF_3D,1); tasGF_tmp = nanmean(tasGF_tmp,2);
tasGF_tmp = squeeze(tasGF_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCgf_tmp(yr+1,1) = CLgf_tmp(yr+1) - CLgf_tmp(yr);
    tuaEgf_tmp(yr+1,1) = CLgf_tmp(yr+1)./(NPPgf_tmp(yr+1) - NetCgf_tmp(yr+1));
end

XcGf_tmp = tuaEgf_tmp.*NPPgf_tmp;
XpGf_tmp = XcGf_tmp - CLgf_tmp;

% save temporal data
save('file_path\NPP_cmip5_tmp.mat','NPPgf_tmp','-append')
save('file_path\GPP_cmip5_tmp.mat','GPPgf_tmp','-append')
save('file_path\CUE_cmip5_tmp.mat','CUEgf_tmp','-append')
save('file_path\cLand_cmip5_tmp.mat','CLgf_tmp','-append')
save('file_path\tuaE_cmip5_tmp.mat','tuaEgf_tmp','-append')
save('file_path\Cnet_cmip5_tmp.mat','NetCgf_tmp','-append')
save('file_path\Xc_cmip5_tmp.mat','XcGf_tmp','-append')
save('file_path\Xp_cmip5_tmp.mat','XpGf_tmp','-append')
save('file_path\pr_cmip5_tmp.mat','prGF_tmp','-append')
save('file_path\tas_cmip5_tmp.mat','tasGF_tmp','-append')


%% HadGEM2-ES
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                % unit:m2
sftlf_HAD = ncread([file2,'file name of the .nc dataset'],'sftlf');                       % unit: %                
cellarea_HAD = areacell.* sftlf_HAD*10^(-2);                                               % unit: m2

nppHAD_3D = ncread([file2, 'file name of the .nc dataset'],'npp');        % unit: KgC m-2 s-1
gppHAD_3D = ncread([file2, 'file name of the .nc dataset'],'gpp');        % unit: KgC m-2 s-1
cLand_HAD = ncread([file2, 'file name of the .nc dataset'],'cLand');   % unit: KgC m-2
pr_HAD = ncread([file2, 'file name of the .nc dataset'],'pr');     % unit: Kg m-2 m-1 = mm
tas_HAD = ncread([file2, 'file name of the .nc dataset'],'tas');  % unit: K

nppHAD_3D(nppHAD_3D<0) = NaN; 
gppHAD_3D(gppHAD_3D<0) = NaN; 
cLand_HAD(cLand_HAD<0) = NaN; 

nppHAD_3Dpg = nppHAD_3D.* cellarea_HAD.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppHAD_3Dpg = gppHAD_3D.* cellarea_HAD.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLhad_3Dpg = cLand_HAD.* cellarea_HAD.* 10^(-12);  % unit: PgC

NPPhad_tmp = nansum(nppHAD_3Dpg,1); NPPhad_tmp = nansum(NPPhad_tmp,2); 
NPPhad_tmp = squeeze(NPPhad_tmp);

GPPhad_tmp = nansum(gppHAD_3Dpg,1);GPPhad_tmp = nansum(GPPhad_tmp,2);
GPPhad_tmp = squeeze(GPPhad_tmp);

CUEhad_tmp = NPPhad_tmp./GPPhad_tmp;

CLhad_tmp = nansum(CLhad_3Dpg,1); CLhad_tmp = nansum(CLhad_tmp,2);
CLhad_tmp = squeeze(CLhad_tmp);

prHAD_tmp = nanmean(pr_HAD_3D,1); prHAD_tmp = nanmean(prHAD_tmp,2);
prHAD_tmp = squeeze(prHAD_tmp);

tasHAD_tmp = nanmean(tas_HAD_3D,1); tasHAD_tmp = nanmean(tasHAD_tmp,2);
tasHAD_tmp = squeeze(tasHAD_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetChad_tmp(yr+1,1) = CLhad_tmp(yr+1) - CLhad_tmp(yr);
    tuaEhad_tmp(yr+1,1) = CLhad_tmp(yr+1)./(NPPhad_tmp(yr+1) - NetChad_tmp(yr+1));
end

XcHad_tmp = tuaEhad_tmp.*NPPhad_tmp;
XpHad_tmp = XcHad_tmp - CLhad_tmp;

% save temporal data
save('file_path\NPP_cmip5_tmp.mat','NPPhad_tmp','-append')
save('file_path\GPP_cmip5_tmp.mat','GPPhad_tmp','-append')
save('file_path\CUE_cmip5_tmp.mat','CUEhad_tmp','-append')
save('file_path\cLand_cmip5_tmp.mat','CLhad_tmp','-append')
save('file_path\tuaE_cmip5_tmp.mat','tuaEhad_tmp','-append')
save('file_path\Cnet_cmip5_tmp.mat','NetChad_tmp','-append')
save('file_path\Xc_cmip5_tmp.mat','XcHad_tmp','-append')
save('file_path\Xp_cmip5_tmp.mat','XpHad_tmp','-append')
save('file_path\pr_cmip5_tmp.mat','prHAD_tmp','-append')
save('file_path\tas_cmip5_tmp.mat','tasHAD_tmp','-append')


%% IPSL-CM5A-MR
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');  % unit:m2
sftlf_IPSL = ncread([file2,'file name of the .nc dataset'],'sftlf');    % unit: %                
cellarea_IPSL = areacell.* sftlf_IPSL*10^(-2);                          % unit: m2

nppIPSL_3D = ncread([file2, 'file name of the .nc dataset'],'npp');     % unit: KgC m-2 s-1
gppIPSL_3D = ncread([file2, 'file name of the .nc dataset'],'gpp');     % unit: KgC m-2 s-1
cLand_IPSL = ncread([file2, 'file name of the .nc dataset'],'cLand');   % unit: KgC m-2

pr_IPSL = ncread([file2, 'file name of the .nc dataset'],'pr');         % unit: Kg m-2 m-1 = mm
tas_IPSL = ncread([file2, 'file name of the .nc dataset'],'tas');       % unit: Kg

nppIPSL_3D(nppIPSL_3D<0) = NaN; 
gppIPSL_3D(gppIPSL_3D<0) = NaN; 
cLand_IPSL(cLand_IPSL<0) = NaN; 

nppIPSL_3Dpg = nppIPSL_3D.* cellarea_IPSL.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppIPSL_3Dpg = gppIPSL_3D.* cellarea_IPSL.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLipsl_3Dpg = cLand_IPSL.* cellarea_IPSL.* 10^(-12);  % unit: PgC

NPPipsl_tmp = nansum(nppIPSL_3Dpg,1); NPPipsl_tmp = nansum(NPPipsl_tmp,2); 
NPPipsl_tmp = squeeze(NPPipsl_tmp);

GPPipsl_tmp = nansum(gppIPSL_3Dpg,1);GPPipsl_tmp = nansum(GPPipsl_tmp,2);
GPPipsl_tmp = squeeze(GPPipsl_tmp);

CUEipsl_tmp = NPPipsl_tmp./GPPipsl_tmp;

CLipsl_tmp = nansum(CLipsl_3Dpg,1); CLipsl_tmp = nansum(CLipsl_tmp,2);
CLipsl_tmp = squeeze(CLipsl_tmp);

prIPSL_tmp = nanmean(pr_IPSL_3D,1); prIPSL_tmp = nanmean(prIPSL_tmp,2);
prIPSL_tmp = squeeze(prIPSL_tmp);

tasIPSL_tmp = nanmean(tas_IPSL_3D,1); tasIPSL_tmp = nanmean(tasIPSL_tmp,2);
tasIPSL_tmp = squeeze(tasIPSL_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCipsl_tmp(yr+1,1) = CLipsl_tmp(yr+1) - CLipsl_tmp(yr);
    tuaEipsl_tmp(yr+1,1) = CLipsl_tmp(yr+1)./(NPPipsl_tmp(yr+1) - NetCipsl_tmp(yr+1));
end

XcIpsl_tmp = tuaEipsl_tmp.*NPPipsl_tmp;
XpIpsl_tmp = XcIpsl_tmp - CLipsl_tmp;

% save temporal data
save('file_path\NPP_cmip5_tmp.mat','NPPipsl_tmp','-append')
save('file_path\GPP_cmip5_tmp.mat','GPPipsl_tmp','-append')
save('file_path\CUE_cmip5_tmp.mat','CUEipsl_tmp','-append')
save('file_path\cLand_cmip5_tmp.mat','CLipsl_tmp','-append')
save('file_path\tuaE_cmip5_tmp.mat','tuaEipsl_tmp','-append')
save('file_path\Cnet_cmip5_tmp.mat','NetCipsl_tmp','-append')
save('file_path\Xc_cmip5_tmp.mat','XcIpsl_tmp','-append')
save('file_path\Xp_cmip5_tmp.mat','XpIpsl_tmp','-append')
save('file_path\pr_cmip5_tmp.mat','prIPSL_tmp','-append')
save('file_path\tas_cmip5_tmp.mat','tasIPSL_tmp','-append')


%% MIROC-ESM
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');      % unit:m2
sftlf_MIROC = ncread([file2,'file name of the .nc dataset'],'sftlf');              % unit: %                
cellarea_MIROC = areacell.* sftlf_MIROC*10^(-2);                                       % unit: m2

nppMIROC_3D = ncread([file2, 'file name of the .nc dataset'],'npp');        % unit: KgC m-2 s-1
gppMIROC_3D = ncread([file2, 'file name of the .nc dataset'],'gpp');        % unit: KgC m-2 s-1
cLand_MIROC = ncread([file2, 'file name of the .nc dataset'],'cLand');   % unit: KgC m-2
pr_MIROC = ncread([file2, 'file name of the .nc dataset'],'pr');     % unit: Kg m-2 m-1 = mm
tas_MIROC = ncread([file2, 'file name of the .nc dataset'],'tas');  % unit: KgC m-2

nppMIROC_3D(nppMIROC_3D<0) = NaN; 
gppMIROC_3D(gppMIROC_3D<0) = NaN; 
cLand_MIROC(cLand_MIROC<0) = NaN; 

nppMIROC_3Dpg = nppMIROC_3D.* cellarea_MIROC.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppMIROC_3Dpg = gppMIROC_3D.* cellarea_MIROC.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLmiroc_3Dpg = cLand_MIROC.* cellarea_MIROC.* 10^(-12);  % unit: PgC

NPPmiroc_tmp = nansum(nppMIROC_3Dpg,1); NPPmiroc_tmp = nansum(NPPmiroc_tmp,2); 
NPPmiroc_tmp = squeeze(NPPmiroc_tmp);

GPPmiroc_tmp = nansum(gppMIROC_3Dpg,1);GPPmiroc_tmp = nansum(GPPmiroc_tmp,2);
GPPmiroc_tmp = squeeze(GPPmiroc_tmp);

CUEmiroc_tmp = NPPmiroc_tmp./GPPmiroc_tmp;

CLmiroc_tmp = nansum(CLmiroc_3Dpg,1); CLmiroc_tmp = nansum(CLmiroc_tmp,2);
CLmiroc_tmp = squeeze(CLmiroc_tmp);

prMIROC_tmp = nanmean(pr_MIROC_3D,1); prMIROC_tmp = nanmean(prMIROC_tmp,2);
prMIROC_tmp = squeeze(prMIROC_tmp);

tasMIROC_tmp = nanmean(tas_MIROC_3D,1); tasMIROC_tmp = nanmean(tasMIROC_tmp,2);
tasMIROC_tmp = squeeze(tasMIROC_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCmiroc_tmp(yr+1,1) = CLmiroc_tmp(yr+1) - CLmiroc_tmp(yr);
    tuaEmiroc_tmp(yr+1,1) = CLmiroc_tmp(yr+1)./(NPPmiroc_tmp(yr+1) - NetCmiroc_tmp(yr+1));
end

XcMiroc_tmp = tuaEmiroc_tmp.*NPPmiroc_tmp;
XpMiroc_tmp = XcMiroc_tmp - CLmiroc_tmp;

% save temporal data
save('file_path\NPP_cmip5_tmp.mat','NPPmiroc_tmp','-append')
save('file_path\GPP_cmip5_tmp.mat','GPPmiroc_tmp','-append')
save('file_path\CUE_cmip5_tmp.mat','CUEmiroc_tmp','-append')
save('file_path\cLand_cmip5_tmp.mat','CLmiroc_tmp','-append')
save('file_path\tuaE_cmip5_tmp.mat','tuaEmiroc_tmp','-append')
save('file_path\Cnet_cmip5_tmp.mat','NetCmiroc_tmp','-append')
save('file_path\Xc_cmip5_tmp.mat','XcMiroc_tmp','-append')
save('file_path\Xp_cmip5_tmp.mat','XpMiroc_tmp','-append')
save('file_path\pr_cmip5_tmp.mat','prMIROC_tmp','-append')
save('file_path\tas_cmip5_tmp.mat','tasMIROC_tmp','-append')


%% MPI-ESM-MR
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');     % unit:m2
sftlf_MPI = ncread([file2,'file name of the .nc dataset'],'sftlf');            % unit: %                
cellarea_MPI = areacell.* sftlf_MPI*10^(-2);                                               % unit: m2

nppMPI_3D = ncread([file2, 'file name of the .nc dataset'],'npp');        % unit: KgC m-2 s-1
gppMPI_3D = ncread([file2, 'file name of the .nc dataset'],'gpp');        % unit: KgC m-2 s-1
cLand_MPI = ncread([file2, 'file name of the .nc dataset'],'cLand');   % unit: KgC m-2
pr_MPI = ncread([file2, 'file name of the .nc dataset'],'pr');     % unit: Kg m-2 m-1 = mm
tas_MPI = ncread([file2, 'file name of the .nc dataset'],'tas');  % unit: KgC m-2

nppMPI_3D(nppMPI_3D<0) = NaN; 
gppMPI_3D(gppMPI_3D<0) = NaN; 
cLand_MPI(cLand_MPI<0) = NaN; 

nppMPI_3Dpg = nppMPI_3D.* cellarea_MPI.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppMPI_3Dpg = gppMPI_3D.* cellarea_MPI.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLmpi_3Dpg = cLand_MPI.* cellarea_MPI.* 10^(-12);  % unit: PgC

NPPmpi_tmp = nansum(nppMPI_3Dpg,1); NPPmpi_tmp = nansum(NPPmpi_tmp,2); 
NPPmpi_tmp = squeeze(NPPmpi_tmp);

GPPmpi_tmp = nansum(gppMPI_3Dpg,1);GPPmpi_tmp = nansum(GPPmpi_tmp,2);
GPPmpi_tmp = squeeze(GPPmpi_tmp);

CUEmpi_tmp = NPPmpi_tmp./GPPmpi_tmp;

CLmpi_tmp = nansum(CLmpi_3Dpg,1); CLmpi_tmp = nansum(CLmpi_tmp,2);
CLmpi_tmp = squeeze(CLmpi_tmp);

prMPI_tmp = nanmean(pr_MPI_3D,1); prMPI_tmp = nanmean(prMPI_tmp,2);
prMPI_tmp = squeeze(prMPI_tmp);

tasMPI_tmp = nanmean(tas_MPI_3D,1); tasMPI_tmp = nanmean(tasMPI_tmp,2);
tasMPI_tmp = squeeze(tasMPI_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCmpi_tmp(yr+1,1) = CLmpi_tmp(yr+1) - CLmpi_tmp(yr);
    tuaEmpi_tmp(yr+1,1) = CLmpi_tmp(yr+1)./(NPPmpi_tmp(yr+1) - NetCmpi_tmp(yr+1));
end

XcMpi_tmp = tuaEmpi_tmp.*NPPmpi_tmp;
XpMpi_tmp = XcMpi_tmp - CLmpi_tmp;


% save temporal data
save('file_path\NPP_cmip5_tmp.mat','NPPmpi_tmp','-append')
save('file_path\GPP_cmip5_tmp.mat','GPPmpi_tmp','-append')
save('file_path\CUE_cmip5_tmp.mat','CUEmpi_tmp','-append')
save('file_path\cLand_cmip5_tmp.mat','CLmpi_tmp','-append')
save('file_path\tuaE_cmip5_tmp.mat','tuaEmpi_tmp','-append')
save('file_path\Cnet_cmip5_tmp.mat','NetCmpi_tmp','-append')
save('file_path\Xc_cmip5_tmp.mat','XcMpi_tmp','-append')
save('file_path\Xp_cmip5_tmp.mat','XpMpi_tmp','-append')
save('file_path\pr_cmip5_tmp.mat','prMPI_tmp','-append')
save('file_path\tas_cmip5_tmp.mat','tasMPI_tmp','-append')


%% MRI-ESM1
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');     % unit:m2
sftlf_MRI = ncread([file2,'file name of the .nc dataset'],'sftlf');            % unit: %                
cellarea_MRI = areacell.* sftlf_MRI*10^(-2);                                               % unit: m2

nppMRI_3D = ncread([file2, 'file name of the .nc dataset'],'npp');        % unit: KgC m-2 s-1
gppMRI_3D = ncread([file2, 'file name of the .nc dataset'],'gpp');        % unit: KgC m-2 s-1
cLand_MRI = ncread([file2, 'file name of the .nc dataset'],'cLand');   % unit: KgC m-2
pr_MRI = ncread([file2, 'file name of the .nc dataset'],'pr');     % unit: Kg m-2 m-1 = mm
tas_MRI = ncread([file2, 'file name of the .nc dataset'],'tas');  % unit: KgC m-2

nppMRI_3D(nppMRI_3D<0) = NaN; 
gppMRI_3D(gppMRI_3D<0) = NaN; 
cLand_MRI(cLand_MRI<0) = NaN; 

nppMRI_3Dpg = nppMRI_3D.* cellarea_MRI.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppMRI_3Dpg = gppMRI_3D.* cellarea_MRI.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLmri_3Dpg = cLand_MRI.* cellarea_MRI.* 10^(-12);  % unit: PgC

NPPmri_tmp = nansum(nppMRI_3Dpg,1); NPPmri_tmp = nansum(NPPmri_tmp,2); 
NPPmri_tmp = squeeze(NPPmri_tmp);

GPPmri_tmp = nansum(gppMRI_3Dpg,1);GPPmri_tmp = nansum(GPPmri_tmp,2);
GPPmri_tmp = squeeze(GPPmri_tmp);

CUEmri_tmp = NPPmri_tmp./GPPmri_tmp;

CLmri_tmp = nansum(CLmri_3Dpg,1); CLmri_tmp = nansum(CLmri_tmp,2);
CLmri_tmp = squeeze(CLmri_tmp);

prMRI_tmp = nanmean(pr_MRI_3D,1); prMRI_tmp = nanmean(prMRI_tmp,2);
prMRI_tmp = squeeze(prMRI_tmp);

tasMRI_tmp = nanmean(tas_MRI_3D,1); tasMRI_tmp = nanmean(tasMRI_tmp,2);
tasMRI_tmp = squeeze(tasMRI_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCmri_tmp(yr+1,1) = CLmri_tmp(yr+1) - CLmri_tmp(yr);
    tuaEmri_tmp(yr+1,1) = CLmri_tmp(yr+1)./(NPPmri_tmp(yr+1) - NetCmri_tmp(yr+1));
end

XcMri_tmp = tuaEmri_tmp.*NPPmri_tmp;
XpMri_tmp = XcMri_tmp - CLmri_tmp;


% save temporal data
save('file_path\NPP_cmip5_tmp.mat','NPPmri_tmp','-append')
save('file_path\GPP_cmip5_tmp.mat','GPPmri_tmp','-append')
save('file_path\CUE_cmip5_tmp.mat','CUEmri_tmp','-append')
save('file_path\cLand_cmip5_tmp.mat','CLmri_tmp','-append')
save('file_path\tuaE_cmip5_tmp.mat','tuaEmri_tmp','-append')
save('file_path\Cnet_cmip5_tmp.mat','NetCmri_tmp','-append')
save('file_path\Xc_cmip5_tmp.mat','XcMri_tmp','-append')
save('file_path\Xp_cmip5_tmp.mat','XpMri_tmp','-append')
save('file_path\pr_cmip5_tmp.mat','prMRI_tmp','-append')
save('file_path\tas_cmip5_tmp.mat','tasMRI_tmp','-append')


%% NorESM1-M
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');      % unit:m2
sftlf_NOR = ncread([file2,'file name of the .nc dataset'],'sftlf');             % unit: %                
cellarea_NOR = areacell.* sftlf_NOR*10^(-2);                                               % unit: m2

nppNOR_3D = ncread([file2, 'file name of the .nc dataset'],'npp');        % unit: KgC m-2 s-1
gppNOR_3D = ncread([file2, 'file name of the .nc dataset'],'gpp');        % unit: KgC m-2 s-1
cLand_NOR = ncread([file2, 'file name of the .nc dataset'],'cLand3');  % unit: KgC m-2
pr_NOR = ncread([file2, 'file name of the .nc dataset'],'pr');        % unit: Kg m-2 m-1 = mm
tas_NOR = ncread([file2, 'file name of the .nc dataset'],'tas');     % unit: KgC m-2

nppNOR_3D(nppNOR_3D<0) = NaN; 
gppNOR_3D(gppNOR_3D<0) = NaN; 
cLand_NOR(cLand_NOR<0) = NaN; 

nppNOR_3Dpg = nppNOR_3D.* cellarea_NOR.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppNOR_3Dpg = gppNOR_3D.* cellarea_NOR.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLnor_3Dpg = cLand_NOR.* cellarea_NOR.* 10^(-12);                % unit: PgC

NPPnor_tmp = nansum(nppNOR_3Dpg,1); NPPnor_tmp = nansum(NPPnor_tmp,2); 
NPPnor_tmp = squeeze(NPPnor_tmp);

GPPnor_tmp = nansum(gppNOR_3Dpg,1);GPPnor_tmp = nansum(GPPnor_tmp,2);
GPPnor_tmp = squeeze(GPPnor_tmp);

CUEnor_tmp = NPPnor_tmp./GPPnor_tmp;

CLnor_tmp = nansum(CLnor_3Dpg,1); CLnor_tmp = nansum(CLnor_tmp,2);
CLnor_tmp = squeeze(CLnor_tmp);

prNOR_tmp = nanmean(pr_NOR_3D,1); prNOR_tmp = nanmean(prNOR_tmp,2);
prNOR_tmp = squeeze(prNOR_tmp);

tasNOR_tmp = nanmean(tas_NOR_3D,1); tasNOR_tmp = nanmean(tasNOR_tmp,2);
tasNOR_tmp = squeeze(tasNOR_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCnor_tmp(yr+1,1) = CLnor_tmp(yr+1) - CLnor_tmp(yr);
    tuaEnor_tmp(yr+1,1) = CLnor_tmp(yr+1)./(NPPnor_tmp(yr+1) - NetCnor_tmp(yr+1));
end

XcNor_tmp = tuaEnor_tmp.*NPPnor_tmp;
XpNor_tmp = XcNor_tmp - CLnor_tmp;

% save temporal data
save('file_path\NPP_cmip5_tmp.mat','NPPnor_tmp','-append')
save('file_path\GPP_cmip5_tmp.mat','GPPnor_tmp','-append')
save('file_path\CUE_cmip5_tmp.mat','CUEnor_tmp','-append')
save('file_path\cLand_cmip5_tmp.mat','CLnor_tmp','-append')
save('file_path\tuaE_cmip5_tmp.mat','tuaEnor_tmp','-append')
save('file_path\Cnet_cmip5_tmp.mat','NetCnor_tmp','-append')
save('file_path\Xc_cmip5_tmp.mat','XcNor_tmp','-append')
save('file_path\Xp_cmip5_tmp.mat','XpNor_tmp','-append')
save('file_path\pr_cmip5_tmp.mat','prNOR_tmp','-append')
save('file_path\tas_cmip5_tmp.mat','tasNOR_tmp','-append')


