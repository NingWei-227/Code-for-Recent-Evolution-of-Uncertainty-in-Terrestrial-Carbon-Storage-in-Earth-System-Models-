###### CMIP6 ###########
setwd('File_path/CMIP6_xls')

X_cmip6 = read_excel('CMIP6_temp11.xlsx', sheet = 1, col_names = TRUE, na='Na')
Xc_cmip6 = read_excel('CMIP6_temp11.xlsx', sheet = 2, col_names = TRUE, na='Na')
Xp_cmip6 = read_excel('CMIP6_temp11.xlsx', sheet = 3, col_names = TRUE, na='Na')

NPP_cmip6 = read_excel('CMIP6_temp11.xlsx', sheet = 4, col_names = TRUE, na='Na')
tuaE_cmip6 = read_excel('CMIP6_temp11.xlsx', sheet = 5, col_names = TRUE, na='Na')

GPP_cmip6 = read_excel('CMIP6_temp11.xlsx', sheet = 6, col_names = TRUE, na='Na')
CUE_cmip6 = read_excel('CMIP6_temp11.xlsx', sheet = 7, col_names = TRUE, na='Na')

tas_cmip6 = read_excel('CMIP6_temp11.xlsx', sheet = 8, col_names = TRUE, na='Na')
pr_cmip6 = read_excel('CMIP6_temp11.xlsx', sheet = 9, col_names = TRUE, na='Na')
opTuaE_cmip6 = read_excel('CMIP6_temp11.xlsx', sheet = 10, col_names = TRUE, na='Na')
baseTuaE_cmip6 = read_excel('CMIP6_temp11.xlsx', sheet = 11, col_names = TRUE, na='Na')

X6_ed = X_cmip6[151:155,2:12] 
Xp6_ed = Xp_cmip6[151:155,2:12] 
Xc6_ed = Xc_cmip6[151:155,2:12] 

tuaE6_ed = tuaE_cmip6[151:155,2:12] 
NPP6_ed = NPP_cmip6[151:155,2:12] 

GPP6_ed = GPP_cmip6[151:155,2:12] 
CUE6_ed = CUE_cmip6[151:155,2:12] 

opTuaE6_ed = opTuaE_cmip6[151:155,2:12] 
pr6_ed = pr_cmip6[151:155,2:12] 
tas6_ed = tas_cmip6[151:155,2:12] 
baseTuaE6_ed = baseTuaE_cmip6[151:155,2:12] 

X6_ed = colMeans(X6_ed)
Xp6_ed = colMeans(Xp6_ed)
Xc6_ed = colMeans(Xc6_ed)

tuaE6_ed = colMeans(tuaE6_ed)
NPP6_ed = colMeans(NPP6_ed)

GPP6_ed = colMeans(GPP6_ed)
CUE6_ed = colMeans(CUE6_ed)

opTuaE6_ed = colMeans(opTuaE6_ed)
pr6_ed = colMeans(pr6_ed)
tas6_ed = colMeans(tas6_ed)
baseTuaE6_ed = colMeans(baseTuaE6_ed)

#step1: X~Xp+Xc
X6ed_ft = data.frame(Xc6_ed,Xp6_ed)
step1 = hier.part(X6_ed,X6ed_ft,family = 'gaussian',gof = "Rsqu", barplot = FALSE)
Xcp6_X = step1$I.perc$ind.exp.var
Xc_X = Xcp6_X[1]

#step2: Xc~NPP+tuaE
Xc6ed_ft = data.frame(log(tuaE6_ed),log(NPP6_ed))
step2 = hier.part(log(Xc6_ed),Xc6ed_ft,family = 'gaussian',gof = "Rsqu", barplot = FALSE)
tuaE_npp_2Xc = step2$I.perc$ind.exp.var
tuaE_npp_2X = Xc_X*tuaE_npp_2Xc*0.01
tuaE_2X = tuaE_npp_2X[1]
npp_2X = tuaE_npp_2X[2]

# step3: NPP~gpp+CUE
NPP6ed_ft = data.frame(log(GPP6_ed),log(CUE6_ed))
step3 = hier.part(log(NPP6_ed),NPP6ed_ft,family = 'gaussian',gof = "Rsqu", barplot = FALSE)
gpp_cue_2npp = step3$I.perc$ind.exp.var
gpp_cue_2X = gpp_cue_2npp*npp_2X*0.01

# step4: tuaE~tas+pr+tuaE'
tuaE6ed_ft = data.frame(log(baseTuaE6_ed),log(tas6_ed),log(pr6_ed))
step4 = hier.part(log(opTuaE6_ed),tuaE6ed_ft,family = 'gaussian',gof = "Rsqu", barplot = FALSE)
base_tas_pr_2tuaE = step4$I.perc$ind.exp.var


'Xp_X'
print(Xcp6_X[2])

'NPP_2X'
print(tuaE_npp_2X[2])

'tuaE_2X'
print(tuaE_npp_2X[1])


'GPP_2NPP'
print(gpp_cue_2npp[1])
'CUE_2NPP'
print(gpp_cue_2npp[2])

'basetuaE_2tuaE'
print(base_tas_pr_2tuaE[1])
'tas_2tuaE'
print(base_tas_pr_2tuaE[2])
'pr_2tuaE'
print(base_tas_pr_2tuaE[3])



# values are used for Figure4_TraceGlobal_Pie_circle








