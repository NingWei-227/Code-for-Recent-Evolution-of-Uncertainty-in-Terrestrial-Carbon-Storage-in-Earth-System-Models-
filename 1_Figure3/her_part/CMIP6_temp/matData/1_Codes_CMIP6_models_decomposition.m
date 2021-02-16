%% terrestrial C sink simulated by CMIP6 model ensemble
%  period: 1850-2005
%  unit: PgC
%  CLand = CVeg + CSoil + CLitter
%  Csink estimated as the absolute change in CLand(yr) relative to the CLand(1850)
%  variables used in analysis: npp, gpp, cLand, sftlf
%  climate data: precipitation (pr) and surface air temperature (tas)

%% ACCESS-ESM1-5
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                   % unit:m2
sftlf_ASS = ncread([file2,'file name of the .nc dataset'],'sftlf');                          % unit: %                
cellarea_ASS = areacell.* sftlf_ASS*10^(-2);                                                                    %unit: m2

nppASS_3D = ncread([file2,'file name of the .nc dataset'],'npp');      % unit: KgC m-2 s-1
gppASS_3D = ncread([file2,'file name of the .nc dataset'],'gpp');      % unit: KgC m-2 s-1
cLand_ASS = ncread([file2,'file name of the .nc dataset'],'cLand'); % unit: KgC m-2
pr_ASS = ncread([file, 'file name of the .nc dataset'],'pr');     % unit: Kg m-2 m-1 = mm
tas_ASS = ncread([file, 'file name of the .nc dataset'],'tas');  % unit: K

nppASS_3D(nppASS_3D<0) = NaN; 
gppASS_3D(gppASS_3D<0) = NaN; 
cLand_ASS(cLand_ASS<0) = NaN;

nppASS_3Dpg = nppASS_3D.* cellarea_ASS.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppASS_3Dpg = gppASS_3D.* cellarea_ASS.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLass_3Dpg = cLand_ASS.* cellarea_ASS.* 10^(-12);                % unit: PgC

NPPass_tmp = nansum(nppASS_3Dpg,1); NPPass_tmp = nansum(NPPass_tmp,2); 
NPPass_tmp = squeeze(NPPass_tmp);

GPPass_tmp = nansum(gppASS_3Dpg,1);GPPass_tmp = nansum(GPPass_tmp,2);
GPPass_tmp = squeeze(GPPass_tmp);

CUEass_tmp = NPPass_tmp./GPPass_tmp;

CLass_tmp = nansum(CLass_3Dpg,1); CLass_tmp = nansum(CLass_tmp,2);
CLass_tmp = squeeze(CLass_tmp);

prASS_tmp = nanmean(pr_ASS_3D,1); prASS_tmp = nanmean(prASS_tmp,2);
prASS_tmp = squeeze(prASS_tmp);

tasASS_tmp = nanmean(tas_ASS_3D,1); tasASS_tmp = nanmean(tasASS_tmp,2);
tasASS_tmp = squeeze(tasASS_tmp);


% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCass_tmp(yr+1,1) = CLass_tmp(yr+1) - CLass_tmp(yr);
    tuaEass_tmp(yr+1,1) = CLass_tmp(yr+1)./(NPPass_tmp(yr+1) - NetCass_tmp(yr+1));
end

XcAss_tmp = tuaEass_tmp.*NPPass_tmp;
XpAss_tmp = XcAss_tmp - CLass_tmp;

% save temporal data
save('file_path\2NPP_cmip6_tmp.mat','NPPass_tmp')
save('file_path\2GPP_cmip6_tmp.mat','GPPass_tmp')
save('file_path\2CUE_cmip6_tmp.mat','CUEass_tmp')
save('file_path\2cLand_cmip6_tmp.mat','CLass_tmp')
save('file_path\2tuaE_cmip6_tmp.mat','tuaEass_tmp')
save('file_path\2Cnet_cmip6_tmp.mat','NetCass_tmp')
save('file_path\2Xc_cmip6_tmp.mat','XcAss_tmp')
save('file_path\2Xp_cmip6_tmp.mat','XpAss_tmp')
save('file_path\2pr_cmip6_tmp.mat','prASS_tmp')
save('file_path\2tas_cmip6_tmp.mat','tasASS_tmp')


%% BCC-CSM2-MR
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                   % unit:m2
sftlf_BCC = ncread([file2,'file name of the .nc dataset'],'sftlf');                      % unit: %                
cellarea_BCC = areacell.* sftlf_BCC*10^(-2);                                             %unit: m2

nppBCC_3D = ncread([file2,'file name of the .nc dataset'],'npp');      % unit: KgC m-2 s-1
gppBCC_3D = ncread([file2,'file name of the .nc dataset'],'gpp');      % unit: KgC m-2 s-1
cLand_BCC = ncread([file2,'file name of the .nc dataset'],'cLand');    % unit: KgC m-2
pr_BCC = ncread([file, 'file name of the .nc dataset'],'pr');          % unit: Kg m-2 m-1 = mm
tas_BCC = ncread([file, 'file name of the .nc dataset'],'tas');        % unit: KgC m-2

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
for yr=1:nyear-1
    NetCbcc_tmp(yr+1,1) = CLbcc_tmp(yr+1) - CLbcc_tmp(yr);
    tuaEbcc_tmp(yr+1,1) = CLbcc_tmp(yr+1)./(NPPbcc_tmp(yr+1) - NetCbcc_tmp(yr+1));
end

XcBcc_tmp = tuaEbcc_tmp.*NPPbcc_tmp;
XpBcc_tmp = XcBcc_tmp - CLbcc_tmp;

% save temporal data
save('file_path\2NPP_cmip6_tmp.mat','NPPbcc_tmp','-append')
save('file_path\2GPP_cmip6_tmp.mat','GPPbcc_tmp','-append')
save('file_path\2CUE_cmip6_tmp.mat','CUEbcc_tmp','-append')
save('file_path\2cLand_cmip6_tmp.mat','CLbcc_tmp','-append')
save('file_path\2tuaE_cmip6_tmp.mat','tuaEbcc_tmp','-append')
save('file_path\2Cnet_cmip6_tmp.mat','NetCbcc_tmp','-append')
save('file_path\2Xc_cmip6_tmp.mat','XcBcc_tmp','-append')
save('file_path\2Xp_cmip6_tmp.mat','XpBcc_tmp','-append')
save('file_path\2pr_cmip6_tmp.mat','prBCC_tmp','-append')
save('file_path\2tas_cmip6_tmp.mat','tasBCC_tmp','-append')


%% CanESM5
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                   % unit:m2
sftlf_CAN = ncread([file2,'file name of the .nc dataset'],'sftlf');                      % unit: %                
cellarea_CAN = areacell.* sftlf_CAN*10^(-2);                                             % unit: m2

nppCAN_3D = ncread([file2,'file name of the .nc dataset'],'npp');      % unit: KgC m-2 s-1
gppCAN_3D = ncread([file2,'file name of the .nc dataset'],'gpp');      % unit: KgC m-2 s-1
cLand_CAN = ncread([file2,'file name of the .nc dataset'],'cLand');    % unit: KgC m-2
pr_CAN = ncread([file2, 'file name of the .nc dataset'],'pr');          % unit: Kg m-2 m-1 = mm
tas_CAN = ncread([file2, 'file name of the .nc dataset'],'tas');        % unit: K

nppCAN_3D(nppCAN_3D<0) = NaN; 
gppCAN_3D(gppCAN_3D<0) = NaN; 
cLand_CAN(cLand_CAN<0) = NaN;

nppCAN_3Dpg = nppCAN_3D.* cellarea_CAN.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppCAN_3Dpg = gppCAN_3D.* cellarea_CAN.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLcan_3Dpg = cLand_CAN.* cellarea_CAN.* 10^(-12);  % unit: PgC

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

% save temporal data
save('file_path\2NPP_cmip6_tmp.mat','NPPcan_tmp','-append')
save('file_path\2GPP_cmip6_tmp.mat','GPPcan_tmp','-append')
save('file_path\2CUE_cmip6_tmp.mat','CUEcan_tmp','-append')
save('file_path\2cLand_cmip6_tmp.mat','CLcan_tmp','-append')
save('file_path\2tuaE_cmip6_tmp.mat','tuaEcan_tmp','-append')
save('file_path\2Cnet_cmip6_tmp.mat','NetCcan_tmp','-append')
save('file_path\2Xc_cmip6_tmp.mat','XcCan_tmp','-append')
save('file_path\2Xp_cmip6_tmp.mat','XpCan_tmp','-append')
save('file_path\2pr_cmip6_tmp.mat','prCAN_tmp','-append')
save('file_path\2tas_cmip6_tmp.mat','tasCAN_tmp','-append')

%% CESM2
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                    % unit:m2
sftlf_CESM = ncread([file2,'file name of the .nc dataset'],'sftlf');                          % unit: %                
cellarea_CESM = areacell.* sftlf_CESM*10^(-2);                                                            %unit: m2

nppCESM_3D = ncread([file2,'file name of the .nc dataset'],'npp');      % unit: KgC m-2 s-1
gppCESM_3D = ncread([file2,'file name of the .nc dataset'],'gpp');      % unit: KgC m-2 s-1
cLand_CESM = ncread([file2,'file name of the .nc dataset'],'cLand'); % unit: KgC m-2
pr_CESM = ncread([file, 'file name of the .nc dataset'],'pr');     % unit: Kg m-2 m-1 = mm
tas_CESM = ncread([file, 'file name of the .nc dataset'],'tas');  % unit: K

nppCESM_3D(nppCESM_3D<0) = NaN; 
gppCESM_3D(gppCESM_3D<0) = NaN; 
cLand_CESM(cLand_CESM<0) = NaN;

nppCESM_3Dpg = nppCESM_3D.* cellarea_CESM.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppCESM_3Dpg = gppCESM_3D.* cellarea_CESM.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLcesm_3Dpg = cLand_CESM.* cellarea_CESM.* 10^(-12);  % unit: PgC

NPPcesm_tmp = nansum(nppCESM_3Dpg,1); NPPcesm_tmp = nansum(NPPcesm_tmp,2); 
NPPcesm_tmp = squeeze(NPPcesm_tmp);

GPPcesm_tmp = nansum(gppCESM_3Dpg,1);GPPcesm_tmp = nansum(GPPcesm_tmp,2);
GPPcesm_tmp = squeeze(GPPcesm_tmp);

CUEcesm_tmp = NPPcesm_tmp./GPPcesm_tmp;

CLcesm_tmp = nansum(CLcesm_3Dpg,1); CLcesm_tmp = nansum(CLcesm_tmp,2);
CLcesm_tmp = squeeze(CLcesm_tmp);

prCESM_tmp = nanmean(pr_CESM_3D,1); prCESM_tmp = nanmean(prCESM_tmp,2);
prCESM_tmp = squeeze(prCESM_tmp);

tasCESM_tmp = nanmean(tas_CESM_3D,1); tasCESM_tmp = nanmean(tasCESM_tmp,2);
tasCESM_tmp = squeeze(tasCESM_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCcesm_tmp(yr+1,1) = CLcesm_tmp(yr+1) - CLcesm_tmp(yr);
    tuaEcesm_tmp(yr+1,1) = CLcesm_tmp(yr+1)./(NPPcesm_tmp(yr+1) - NetCcesm_tmp(yr+1));
end

XcCesm_tmp = tuaEcesm_tmp.*NPPcesm_tmp;
XpCesm_tmp = XcCesm_tmp - CLcesm_tmp;


% save temporal data
save('file_path\2NPP_cmip6_tmp.mat','NPPcesm_tmp','-append')
save('file_path\2GPP_cmip6_tmp.mat','GPPcesm_tmp','-append')
save('file_path\2CUE_cmip6_tmp.mat','CUEcesm_tmp','-append')
save('file_path\2cLand_cmip6_tmp.mat','CLcesm_tmp','-append')
save('file_path\2tuaE_cmip6_tmp.mat','tuaEcesm_tmp','-append')
save('file_path\2Cnet_cmip6_tmp.mat','NetCcesm_tmp','-append')
save('file_path\2Xc_cmip6_tmp.mat','XcCesm_tmp','-append')
save('file_path\2Xp_cmip6_tmp.mat','XpCesm_tmp','-append')
save('file_path\2pr_cmip6_tmp.mat','prCESM_tmp','-append')
save('file_path\2tas_cmip6_tmp.mat','tasCESM_tmp','-append')

%% CNRM-ESM2-1
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                    % unit:m2
sftlf_CNRM = ncread([file2,'file name of the .nc dataset'],'sftlf');                          % unit: %                
cellarea_CNRM = areacell.* sftlf_CNRM*10^(-2);                                                            %unit: m2

nppCNRM_3D = ncread([file2,'file name of the .nc dataset'],'npp');      % unit: KgC m-2 s-1
gppCNRM_3D = ncread([file2,'file name of the .nc dataset'],'gpp');      % unit: KgC m-2 s-1
cLand_CNRM = ncread([file2,'file name of the .nc dataset'],'cLand');    % unit: KgC m-2
pr_CNRM = ncread([file2, 'file name of the .nc dataset'],'pr');         % unit: Kg m-2 m-1 = mm
tas_CNRM = ncread([file2, 'file name of the .nc dataset'],'tas');       % unit: K

nppCNRM_3D(nppCNRM_3D<0) = NaN; 
gppCNRM_3D(gppCNRM_3D<0) = NaN; 
cLand_CNRM(cLand_CNRM<0) = NaN;

nppCNRM_3Dpg = nppCNRM_3D.* cellarea_CNRM.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppCNRM_3Dpg = gppCNRM_3D.* cellarea_CNRM.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLcnrm_3Dpg = cLand_CNRM.* cellarea_CNRM.* 10^(-12);  % unit: PgC

NPPcnrm_tmp = nansum(nppCNRM_3Dpg,1); NPPcnrm_tmp = nansum(NPPcnrm_tmp,2); 
NPPcnrm_tmp = squeeze(NPPcnrm_tmp);

GPPcnrm_tmp = nansum(gppCNRM_3Dpg,1);GPPcnrm_tmp = nansum(GPPcnrm_tmp,2);
GPPcnrm_tmp = squeeze(GPPcnrm_tmp);

CUEcnrm_tmp = NPPcnrm_tmp./GPPcnrm_tmp;

CLcnrm_tmp = nansum(CLcnrm_3Dpg,1); CLcnrm_tmp = nansum(CLcnrm_tmp,2);
CLcnrm_tmp = squeeze(CLcnrm_tmp);

prCNRM_tmp = nanmean(pr_CNRM_3D,1); prCNRM_tmp = nanmean(prCNRM_tmp,2);
prCNRM_tmp = squeeze(prCNRM_tmp);

tasCNRM_tmp = nanmean(tas_CNRM_3D,1); tasCNRM_tmp = nanmean(tasCNRM_tmp,2);
tasCNRM_tmp = squeeze(tasCNRM_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCcnrm_tmp(yr+1,1) = CLcnrm_tmp(yr+1) - CLcnrm_tmp(yr);
    tuaEcnrm_tmp(yr+1,1) = CLcnrm_tmp(yr+1)./(NPPcnrm_tmp(yr+1) - NetCcnrm_tmp(yr+1));
end

XcCnrm_tmp = tuaEcnrm_tmp.*NPPcnrm_tmp;
XpCnrm_tmp = XcCnrm_tmp - CLcnrm_tmp;

% save temporal data
save('file_path\2NPP_cmip6_tmp.mat','NPPcnrm_tmp','-append')
save('file_path\2GPP_cmip6_tmp.mat','GPPcnrm_tmp','-append')
save('file_path\2CUE_cmip6_tmp.mat','CUEcnrm_tmp','-append')
save('file_path\2cLand_cmip6_tmp.mat','CLcnrm_tmp','-append')
save('file_path\2tuaE_cmip6_tmp.mat','tuaEcnrm_tmp','-append')
save('file_path\2Cnet_cmip6_tmp.mat','NetCcnrm_tmp','-append')
save('file_path\2Xc_cmip6_tmp.mat','XcCnrm_tmp','-append')
save('file_path\2Xp_cmip6_tmp.mat','XpCnrm_tmp','-append')
save('file_path\2pr_cmip6_tmp.mat','prCNRM_tmp','-append')
save('file_path\2tas_cmip6_tmp.mat','tasCNRM_tmp','-append')

%% EC-Earth3-Veg
clear
clc

% estimate cellarea
file2 = 'load file path'        
lon_bnds = ncread([file2,'file name of the .nc dataset'],'lon_bnds');
lat_bnds = ncread([file2,'file name of the .nc dataset'],'lat_bnds');


E = wgs84Ellipsoid;
[~,m] = size(lat_bnds);
[~,n] = size(lon_bnds);
lat_bnds = double(lat_bnds);
lon_bnds = double(lon_bnds);

for i = 1:m
    for j = 1:n
        cellarea(i,j) = areaquad(lat_bnds(1,i),lon_bnds(1,j),lat_bnds(2,i),lon_bnds(2,j),E); %m2
    end
end

areacell = cellarea';   % unit:m2
clearvars -except file2 areacell

% estimate global values   
sftlf_EC = ncread([file2,'file name of the .nc dataset'],'sftlf');                          % unit: %                
cellarea_EC = areacell.* sftlf_EC*10^(-2);                                                                   %unit: m2

nppEC_3D = ncread([file2,'file name of the .nc dataset'],'npp');      % unit: KgC m-2 s-1
gppEC_3D = ncread([file2,'file name of the .nc dataset'],'gpp');      % unit: KgC m-2 s-1
cLand_EC = ncread([file2,'file name of the .nc dataset'],'cLand'); % unit: KgC m-2
pr_EC = ncread([file2, 'file name of the .nc dataset'],'pr');      % unit: Kg m-2 m-1 = mm
tas_EC = ncread([file2, 'file name of the .nc dataset'],'tas');   % unit: K

nppEC_3D(nppEC_3D<0) = NaN; 
gppEC_3D(gppEC_3D<0) = NaN; 
cLand_EC(cLand_EC<0) = NaN;

nppEC_3Dpg = nppEC_3D.* cellarea_EC.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppEC_3Dpg = gppEC_3D.* cellarea_EC.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLec_3Dpg = cLand_EC.* cellarea_EC.* 10^(-12);  % unit: PgC

NPPec_tmp = nansum(nppEC_3Dpg,1); NPPec_tmp = nansum(NPPec_tmp,2); 
NPPec_tmp = squeeze(NPPec_tmp);

GPPec_tmp = nansum(gppEC_3Dpg,1);GPPec_tmp = nansum(GPPec_tmp,2);
GPPec_tmp = squeeze(GPPec_tmp);

CUEec_tmp = NPPec_tmp./GPPec_tmp;

CLec_tmp = nansum(CLec_3Dpg,1); CLec_tmp = nansum(CLec_tmp,2);
CLec_tmp = squeeze(CLec_tmp);

prEC_tmp = nanmean(pr_EC_3D,1); prEC_tmp = nanmean(prEC_tmp,2);
prEC_tmp = squeeze(prEC_tmp);

tasEC_tmp = nanmean(tas_EC_3D,1); tasEC_tmp = nanmean(tasEC_tmp,2);
tasEC_tmp = squeeze(tasEC_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCec_tmp(yr+1,1) = CLec_tmp(yr+1) - CLec_tmp(yr);
    tuaEec_tmp(yr+1,1) = CLec_tmp(yr+1)./(NPPec_tmp(yr+1) - NetCec_tmp(yr+1));
end

XcEc_tmp = tuaEec_tmp.*NPPec_tmp;
XpEc_tmp = XcEc_tmp - CLec_tmp;

% save temporal data
save('file_path\2NPP_cmip6_tmp.mat','NPPec_tmp','-append')
save('file_path\2GPP_cmip6_tmp.mat','GPPec_tmp','-append')
save('file_path\2CUE_cmip6_tmp.mat','CUEec_tmp','-append')
save('file_path\2cLand_cmip6_tmp.mat','CLec_tmp','-append')
save('file_path\2tuaE_cmip6_tmp.mat','tuaEec_tmp','-append')
save('file_path\2Cnet_cmip6_tmp.mat','NetCec_tmp','-append')
save('file_path\2Xc_cmip6_tmp.mat','XcEc_tmp','-append')
save('file_path\2Xp_cmip6_tmp.mat','XpEc_tmp','-append')
save('file_path\2pr_cmip6_tmp.mat','prEC_tmp','-append')
save('file_path\2tas_cmip6_tmp.mat','tasEC_tmp','-append')

%% IPSL-CM6A-LR
clear
clc

file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                    % unit:m2
sftlf_IPSL = ncread([file2,'file name of the .nc dataset'],'sftlf');                          % unit: %                
cellarea_IPSL = areacell.* sftlf_IPSL*10^(-2);                                                                   %unit: m2

nppIPSL_3D = ncread([file2,'file name of the .nc dataset'],'npp');       % unit: KgC m-2 s-1
gppIPSL_3D = ncread([file2,'file name of the .nc dataset'],'gpp');       % unit: KgC m-2 s-1
cLand_IPSL = ncread([file2,'file name of the .nc dataset'],'cLand3'); % unit: KgC m-2
pr_IPSL = ncread([file, 'file name of the .nc dataset'],'pr');       % unit: Kg m-2 m-1 = mm
tas_IPSL = ncread([file, 'file name of the .nc dataset'],'tas');    % unit: K

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
save('file_path\2NPP_cmip6_tmp.mat','NPPipsl_tmp','-append')
save('file_path\2GPP_cmip6_tmp.mat','GPPipsl_tmp','-append')
save('file_path\2CUE_cmip6_tmp.mat','CUEipsl_tmp','-append')
save('file_path\2cLand_cmip6_tmp.mat','CLipsl_tmp','-append')
save('file_path\2tuaE_cmip6_tmp.mat','tuaEipsl_tmp','-append')
save('file_path\2Cnet_cmip6_tmp.mat','NetCipsl_tmp','-append')
save('file_path\2Xc_cmip6_tmp.mat','XcIpsl_tmp','-append')
save('file_path\2Xp_cmip6_tmp.mat','XpIpsl_tmp','-append')
save('file_path\2pr_cmip6_tmp.mat','prIPSL_tmp','-append')
save('file_path\2tas_cmip6_tmp.mat','tasIPSL_tmp','-append')

%% MIROC-ES2L
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                    % unit:m2
sftlf_MIC = ncread([file2,'file name of the .nc dataset'],'sftlf');                          % unit: %                
cellarea_MIC = areacell.* sftlf_MIC*10^(-2);                                                                   %unit: m2

nppMIC_3D = ncread([file2,'file name of the .nc dataset'],'npp');      % unit: KgC m-2 s-1
gppMIC_3D = ncread([file2,'file name of the .nc dataset'],'gpp');      % unit: KgC m-2 s-1
cLand_MIC = ncread([file2,'file name of the .nc dataset'],'cLand'); % unit: KgC m-2
pr_MIC = ncread([file, 'file name of the .nc dataset'],'pr');      % unit: Kg m-2 m-1 = mm
tas_MIC = ncread([file, 'file name of the .nc datasetc'],'tas');   % unit: K

nppMIC_3D(nppMIC_3D<0) = NaN; 
gppMIC_3D(gppMIC_3D<0) = NaN; 
cLand_MIC(cLand_MIC<0) = NaN;

nppMIC_3Dpg = nppMIC_3D.* cellarea_MIC.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppMIC_3Dpg = gppMIC_3D.* cellarea_MIC.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLmic_3Dpg = cLand_MIC.* cellarea_MIC.* 10^(-12);  % unit: PgC

NPPmic_tmp = nansum(nppMIC_3Dpg,1); NPPmic_tmp = nansum(NPPmic_tmp,2); 
NPPmic_tmp = squeeze(NPPmic_tmp);

GPPmic_tmp = nansum(gppMIC_3Dpg,1);GPPmic_tmp = nansum(GPPmic_tmp,2);
GPPmic_tmp = squeeze(GPPmic_tmp);

CUEmic_tmp = NPPmic_tmp./GPPmic_tmp;

CLmic_tmp = nansum(CLmic_3Dpg,1); CLmic_tmp = nansum(CLmic_tmp,2);
CLmic_tmp = squeeze(CLmic_tmp);

prMIC_tmp = nanmean(pr_MIC_3D,1); prMIC_tmp = nanmean(prMIC_tmp,2);
prMIC_tmp = squeeze(prMIC_tmp);

tasMIC_tmp = nanmean(tas_MIC_3D,1); tasMIC_tmp = nanmean(tasMIC_tmp,2);
tasMIC_tmp = squeeze(tasMIC_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCmic_tmp(yr+1,1) = CLmic_tmp(yr+1) - CLmic_tmp(yr);
    tuaEmic_tmp(yr+1,1) = CLmic_tmp(yr+1)./(NPPmic_tmp(yr+1) - NetCmic_tmp(yr+1));
end

XcMic_tmp = tuaEmic_tmp.*NPPmic_tmp;
XpMic_tmp = XcMic_tmp - CLmic_tmp;

% save temporal data
save('file_path\2NPP_cmip6_tmp.mat','NPPmic_tmp','-append')
save('file_path\2GPP_cmip6_tmp.mat','GPPmic_tmp','-append')
save('file_path\2CUE_cmip6_tmp.mat','CUEmic_tmp','-append')
save('file_path\2cLand_cmip6_tmp.mat','CLmic_tmp','-append')
save('file_path\2tuaE_cmip6_tmp.mat','tuaEmic_tmp','-append')
save('file_path\2Cnet_cmip6_tmp.mat','NetCmic_tmp','-append')
save('file_path\2Xc_cmip6_tmp.mat','XcMic_tmp','-append')
save('file_path\2Xp_cmip6_tmp.mat','XpMic_tmp','-append')
save('file_path\2pr_cmip6_tmp.mat','prMIC_tmp','-append')
save('file_path\2tas_cmip6_tmp.mat','tasMIC_tmp','-append')

%% MPI-ESM1-2-LR
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                    % unit:m2
sftlf_MPI = ncread([file2,'file name of the .nc dataset'],'sftlf');                          % unit: %                
cellarea_MPI = areacell.* sftlf_MPI*10^(-2);                                                                   %unit: m2

nppMPI_3D = ncread([file2,'file name of the .nc dataset'],'npp');      % unit: KgC m-2 s-1
gppMPI_3D = ncread([file2,'file name of the .nc dataset'],'gpp');      % unit: KgC m-2 s-1
cLand_MPI = ncread([file2,'file name of the .nc dataset'],'cLand'); % unit: KgC m-2
pr_MPI = ncread([file2, 'file name of the .nc dataset'],'pr');     % unit: Kg m-2 m-1 = mm
tas_MPI = ncread([file2, 'file name of the .nc dataset'],'tas');  % unit: K

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
save('file_path\2NPP_cmip6_tmp.mat','NPPmpi_tmp','-append')
save('file_path\2GPP_cmip6_tmp.mat','GPPmpi_tmp','-append')
save('file_path\2CUE_cmip6_tmp.mat','CUEmpi_tmp','-append')
save('file_path\2cLand_cmip6_tmp.mat','CLmpi_tmp','-append')
save('file_path\2tuaE_cmip6_tmp.mat','tuaEmpi_tmp','-append')
save('file_path\2Cnet_cmip6_tmp.mat','NetCmpi_tmp','-append')
save('file_path\2Xc_cmip6_tmp.mat','XcMpi_tmp','-append')
save('file_path\2Xp_cmip6_tmp.mat','XpMpi_tmp','-append')
save('file_path\2pr_cmip6_tmp.mat','prMPI_tmp','-append')
save('file_path\2tas_cmip6_tmp.mat','tasMPI_tmp','-append')

%% NorESM2-LM
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                    % unit:m2
sftlf_NOR = ncread([file2,'file name of the .nc dataset'],'sftlf');                          % unit: %                
cellarea_NOR = areacell.* sftlf_NOR*10^(-2);                                                                   %unit: m2

nppNOR_3D = ncread([file2,'file name of the .nc dataset'],'npp');          % unit: KgC m-2 s-1
gppNOR_3D = ncread([file2,'file name of the .nc dataset'],'gpp');          % unit: KgC m-2 s-1
cLand_NOR = ncread([file2,'file name of the .nc dataset'],'cLand');     % unit: KgC m-2
pr_NOR = ncread([file, 'file name of the .nc dataset'],'pr');     % unit: Kg m-2 m-1 = mm
tas_NOR = ncread([file, 'file name of the .nc dataset'],'tas');       % unit: K

nppNOR_3D(nppNOR_3D<0) = NaN; 
gppNOR_3D(gppNOR_3D<0) = NaN; 
cLand_NOR(cLand_NOR<0) = NaN;

nppNOR_3Dpg = nppNOR_3D.* cellarea_NOR.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppNOR_3Dpg = gppNOR_3D.* cellarea_NOR.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLnor_3Dpg = cLand_NOR.* cellarea_NOR.* 10^(-12);  % unit: PgC


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
save('file_path\2NPP_cmip6_tmp.mat','NPPnor_tmp','-append')
save('file_path\2GPP_cmip6_tmp.mat','GPPnor_tmp','-append')
save('file_path\2CUE_cmip6_tmp.mat','CUEnor_tmp','-append')
save('file_path\2cLand_cmip6_tmp.mat','CLnor_tmp','-append')
save('file_path\2tuaE_cmip6_tmp.mat','tuaEnor_tmp','-append')
save('file_path\2Cnet_cmip6_tmp.mat','NetCnor_tmp','-append')
save('file_path\2Xc_cmip6_tmp.mat','XcNor_tmp','-append')
save('file_path\2Xp_cmip6_tmp.mat','XpNor_tmp','-append')
save('file_path\2pr_cmip6_tmp.mat','prNOR_tmp','-append')
save('file_path\2tas_cmip6_tmp.mat','tasNOR_tmp','-append')

%% UKESM1-0-LL
clear
clc

% estimate global values
file2 = 'load file path'
areacell = ncread([file2,'file name of the .nc dataset'],'areacella');                    % unit:m2
sftlf_UK = ncread([file2,'file name of the .nc dataset'],'sftlf');                          % unit: %                
cellarea_UK = areacell.* sftlf_UK*10^(-2);                                                                   %unit: m2

nppUK_3D = ncread([file2,'file name of the .nc dataset'],'npp');      % unit: KgC m-2 s-1
gppUK_3D = ncread([file2,'file name of the .nc dataset'],'gpp');      % unit: KgC m-2 s-1
cLand_UK = ncread([file2,'file name of the .nc dataset'],'cLand'); % unit: KgC m-2
pr_UK = ncread([file, 'file name of the .nc dataset'],'pr');     % unit: Kg m-2 m-1 = mm
tas_UK = ncread([file, 'file name of the .nc dataset'],'tas');  % unit: K

nppUK_3D(nppUK_3D<0) = NaN; 
gppUK_3D(gppUK_3D<0) = NaN; 
cLand_UK(cLand_UK<0) = NaN;

nppUK_3Dpg = nppUK_3D.* cellarea_UK.* 10^(-12).*60*60*24*365; % convert from [KgC m-2 s-1] to [PgC yr-1]; 
gppUK_3Dpg = gppUK_3D.* cellarea_UK.* 10^(-12).*60*60*24*365; % unit: PgC yr-1 
CLuk_3Dpg = cLand_UK.* cellarea_UK.* 10^(-12);  % unit: PgC

NPPuk_tmp = nansum(nppUK_3Dpg,1); NPPuk_tmp = nansum(NPPuk_tmp,2); 
NPPuk_tmp = squeeze(NPPuk_tmp);

GPPuk_tmp = nansum(gppUK_3Dpg,1);GPPuk_tmp = nansum(GPPuk_tmp,2);
GPPuk_tmp = squeeze(GPPuk_tmp);

CUEuk_tmp = NPPuk_tmp./GPPuk_tmp;

CLuk_tmp = nansum(CLuk_3Dpg,1); CLuk_tmp = nansum(CLuk_tmp,2);
CLuk_tmp = squeeze(CLuk_tmp);

prUK_tmp = nanmean(pr_UK_3D,1); prUK_tmp = nanmean(prUK_tmp,2);
prUK_tmp = squeeze(prUK_tmp);

tasUK_tmp = nanmean(tas_UK_3D,1); tasUK_tmp = nanmean(tasUK_tmp,2);
tasUK_tmp = squeeze(tasUK_tmp);

% estimate terrestrial C residence time 
for yr=1:nyear-1
    NetCuk_tmp(yr+1,1) = CLuk_tmp(yr+1) - CLuk_tmp(yr);
    tuaEuk_tmp(yr+1,1) = CLuk_tmp(yr+1)./(NPPuk_tmp(yr+1) - NetCuk_tmp(yr+1));
end

XcUk_tmp = tuaEuk_tmp.*NPPuk_tmp;
XpUk_tmp = XcUk_tmp - CLuk_tmp;

% save temporal data
save('file_path\2NPP_cmip6_tmp.mat','NPPuk_tmp','-append')
save('file_path\2GPP_cmip6_tmp.mat','GPPuk_tmp','-append')
save('file_path\2CUE_cmip6_tmp.mat','CUEuk_tmp','-append')
save('file_path\2cLand_cmip6_tmp.mat','CLuk_tmp','-append')
save('file_path\2tuaE_cmip6_tmp.mat','tuaEuk_tmp','-append')
save('file_path\2Cnet_cmip6_tmp.mat','NetCuk_tmp','-append')
save('file_path\2Xc_cmip6_tmp.mat','XcUk_tmp','-append')
save('file_path\2Xp_cmip6_tmp.mat','XpUk_tmp','-append')
save('file_path\2pr_cmip6_tmp.mat','prUK_tmp','-append')
save('file_path\2tas_cmip6_tmp.mat','tasUK_tmp','-append')



