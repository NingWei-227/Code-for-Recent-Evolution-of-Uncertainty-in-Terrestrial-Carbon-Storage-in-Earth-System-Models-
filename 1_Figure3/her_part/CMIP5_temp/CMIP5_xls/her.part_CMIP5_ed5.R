setwd('file path/CMIP5_xls')
library(readxl)
library(vegan)
library(hier.part)
library(gtools)


# CMIP5
X_cmip5 = read_excel('CMIP5_temp11.xlsx', sheet = 1, col_names = TRUE, na='Na')
Xc_cmip5 = read_excel('CMIP5_temp11.xlsx', sheet = 2, col_names = TRUE, na='Na')
Xp_cmip5 = read_excel('CMIP5_temp11.xlsx', sheet = 3, col_names = TRUE, na='Na')

NPP_cmip5 = read_excel('CMIP5_temp11.xlsx', sheet = 4, col_names = TRUE, na='Na')
tuaE_cmip5 = read_excel('CMIP5_temp11.xlsx', sheet = 5, col_names = TRUE, na='Na')

GPP_cmip5 = read_excel('CMIP5_temp11.xlsx', sheet = 6, col_names = TRUE, na='Na')
CUE_cmip5 = read_excel('CMIP5_temp11.xlsx', sheet = 7, col_names = TRUE, na='Na')

tas_cmip5 = read_excel('CMIP5_temp11.xlsx', sheet = 8, col_names = TRUE, na='Na')
pr_cmip5 = read_excel('CMIP5_temp11.xlsx', sheet = 9, col_names = TRUE, na='Na')
opTuaE_cmip5 = read_excel('CMIP5_temp11.xlsx', sheet = 10, col_names = TRUE, na='Na')
baseTuaE_cmip5 = read_excel('CMIP5_temp11.xlsx', sheet = 11, col_names = TRUE, na='Na')

X5_ed = X_cmip5[151:155,2:12] 
Xp5_ed = Xp_cmip5[151:155,2:12]
Xc5_ed = Xc_cmip5[151:155,2:12]

tuaE5_ed = tuaE_cmip5[151:155,2:12]
NPP5_ed = NPP_cmip5[151:155,2:12]

GPP5_ed = GPP_cmip5[151:155,2:12]
CUE5_ed = CUE_cmip5[151:155,2:12]

opTuaE5_ed = opTuaE_cmip5[151:155,2:12]
pr5_ed = pr_cmip5[151:155,2:12]
tas5_ed = tas_cmip5[151:155,2:12]
baseTuaE5_ed = baseTuaE_cmip5[151:155,2:12]

X5_ed = colMeans(X5_ed)
Xp5_ed = colMeans(Xp5_ed)
Xc5_ed = colMeans(Xc5_ed)

tuaE5_ed = colMeans(tuaE5_ed)
NPP5_ed = colMeans(NPP5_ed)

GPP5_ed = colMeans(GPP5_ed)
CUE5_ed = colMeans(CUE5_ed)

opTuaE5_ed = colMeans(opTuaE5_ed)
pr5_ed = colMeans(pr5_ed)
tas5_ed = colMeans(tas5_ed)
baseTuaE5_ed = colMeans(baseTuaE5_ed)

#step1: X~Xp+Xc
X5ed_ft = data.frame(Xc5_ed,Xp5_ed)
step1 = hier.part(X5_ed,X5ed_ft,family = 'gaussian',gof = "Rsqu", barplot = FALSE)
Xcp5_X = step1$I.perc$ind.exp.var
Xc_X = Xcp5_X[1]

#step2: Xc~NPP+tuaE
Xc5ed_ft = data.frame(log(tuaE5_ed),log(NPP5_ed))
step2 = hier.part(log(Xc5_ed),Xc5ed_ft,family = 'gaussian',gof = "Rsqu", barplot = FALSE)
tuaE_npp_2Xc = step2$I.perc$ind.exp.var
tuaE_npp_2X = Xc_X*tuaE_npp_2Xc*0.01
tuaE_2X = tuaE_npp_2X[1]
npp_2X = tuaE_npp_2X[2]

# step3: NPP~gpp+CUE
NPP5ed_ft = data.frame(log(GPP5_ed),log(CUE5_ed))
step3 = hier.part(log(NPP5_ed),NPP5ed_ft,family = 'gaussian',gof = "Rsqu", barplot = FALSE)
gpp_cue_2npp = step3$I.perc$ind.exp.var
gpp_cue_2X = gpp_cue_2npp*npp_2X*0.01

# step4: tuaE~tas+pr+tuaE'
tuaE5ed_ft = data.frame(log(baseTuaE5_ed),log(tas5_ed),log(pr5_ed))
step4 = hier.part(log(opTuaE5_ed),tuaE5ed_ft,family = 'gaussian',gof = "Rsqu", barplot = FALSE)
base_tas_pr_2tuaE = step4$I.perc$ind.exp.var


'Xp_X'
print(Xcp5_X[2])

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

















