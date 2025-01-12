library(dplyr)
library(OmicsPLS)
library(readr)
library(readxl)
library(writexl)
library(openxlsx)
library(readxl)
library(pracma)
library(tidyr)

#integration with less no. of components to reduce overfitting

setwd("~/BINF/o2pls/tr-sec")

rna <- read_excel("~/BINF/BINF_old/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx", 
                  col_types = c("text", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "skip", "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "skip", "skip", 
                                "skip", "skip", "skip"))

secretome <- read_excel("~/BINF/BINF_old/multi-omics/ADAM9 KD HCT116 proteomics/secretomics/perseus output/secretome_results.xlsx")
secretome2 = secretome[,c(1,3,5:14)]
rnan = rna[,-1]
secretomen = secretome2[,c(-1, -2)]

#samples as rows and genes/protein as columns
rnaT = t(rnan)
secretomet = t(secretomen)

#need to center data before cross validation
rnaTS = scale(rnaT, center = TRUE, scale = TRUE)
secrTS = as.matrix(scale(secretomet, center = TRUE, scale = TRUE))

input_checker(secrTS)
#input matrices satisfy the required conditions

#lowering p_thresh and q_thresh coz secretome has < 3k columns
crossval_o2m_adjR2(rnaTS, secrTS, 1:3, 0:3, 0:3, nr_folds = 5,  stripped = TRUE,
                   p_thresh = 2000, q_thresh = 2000, tol = 1e-10,
                   max_iterations = 100)

'
Minimum is at n = 1 
Elapsed time: 11.03 sec
       MSE n nx ny
1 1.869010 1  1  1
2 1.908907 2  2  0
3 1.987989 3  1  0

'

crossval_o2m(rnaTS,secrTS,1:2,0:2,0:2,5,  p_thresh = 2000, q_thresh = 2000)

'
*******************
Elapsed time: 8.7 sec
*******
Minimal 5-CV error is at ax=2 ay=2 a=2 
*******
Minimum MSE is 1.845372 
*******************

'

fit = o2m(rnaTS,secrTS,2,2,2, p_thresh = 2000, q_thresh = 2000)
summary(fit)
setwd("~/BINF/o2pls/tr-sec/output")
saveRDS(fit,"tr-secrfit.rds")

loadings = loadings(fit)
loadings_xjoint = as.data.frame(loadings(fit,"Xjoint"))
write_xlsx(loadings_xjoint, "loading_values_TransJoint.xlsx")
loadings_yjoint = as.data.frame(loadings(fit,"Yjoint"))
write_xlsx(loadings_yjoint,"loading_values_secJoint.xlsx")
loadings_xorth = as.data.frame(loadings(fit,"Xorth"))
write_xlsx(loadings_xorth,"loading_values_Transorth.xlsx")
loadings_yorth =as.data.frame( loadings(fit,"Yorth"))
write_xlsx(loadings_yorth,"loading_values_secOrth.xlsx")

scores_xjoint = as.data.frame(scores(fit,"Xjoint"))
write_xlsx(scores_xjoint, "scores_TransJoint.xlsx")
scores_yjoint =as.data.frame( scores(fit,"Yjoint"))
write_xlsx(scores_yjoint,"scores_secJoint.xlsx")
scores_xorth = as.data.frame(scores(fit,"Xorth"))
write_xlsx(scores_xorth,"scores_Transorth.xlsx")
scores_yorth = as.data.frame(scores(fit,"Yorth"))
write_xlsx(scores_yorth,"scores_secOrth.xlsx")

write_xlsx(as.data.frame(fit$E), "trans_residuals.xlsx")
write_xlsx(as.data.frame(fit$Ff), "sec_residuals.xlsx")


#manually indexed all lvs and scores with gene names and ids

######post o2PLS analysis
#find genes which occur in both datasets
setwd("~/BINF/o2pls/tr-sec/output")

lv_TransJnt <- read_excel("loading_values_TransJoint.xlsx")
lv_secrjnt <- read_excel("loading_values_secJoint.xlsx")
comms = inner_join(lv_TransJnt, lv_secrjnt, by = c("GeneName" ))

#2472 common genes 
common_names = comms[,c(1,4)]
write_xlsx(common_names,"analysis after integration/commonProtsSecTrans.xlsx")

#load transcriptome joint scores
T <- read_xlsx("scores_TransJoint.xlsx")
T = as.matrix(T[-1])

T1 = T[,1]
T2 = T[,2]

#load protein joint scores
U <- read_excel("scores_secJoint.xlsx")
U <- as.matrix(U[,-1])
U1 = U[,1]
U2 = U[,2]

colnames(T)[colnames(T) == "V1"] <- "V1gene"
colnames(T)[colnames(T) == "V2"] <- "V2gene"

colnames(U)[colnames(U) == "V1"] <- "V1sec"
colnames(U)[colnames(U) == "V2"] <- "V2sec"

#load joint loading values
#transcriptome
LX <- read_excel("loading_values_TransJoint.xlsx")
#secretome
LY <- read_excel("loading_values_secJoint.xlsx")

B = diag(c(0.374, 0.417), nrow = 2, ncol = 2)

cov1 = inner_join(LX, LY, by =c("GeneName" ))
#in cov1, X is transcriptome and y is proteome
cov1 = cov1 %>% rename("V1gene" = "V1.x")
cov1 = cov1 %>% rename("V2gene" = "V2.x")

cov1 = cov1 %>% rename("V1sec" = "V1.y")
cov1 = cov1 %>% rename("V2sec" = "V2.y")

C = apply(cov1, 1, function(row){
  
  lvx <- matrix(as.numeric(c(row["V1gene"], row["V2gene"])) , nrow = 2, ncol = 1)
  
  lvy <- matrix( as.numeric(c(row["V1sec"], row["V2sec"])) , nrow = 2, ncol = 1)
  jointx<-T1*lvx[1]+T2*lvx[2]
  jointy<-U1*lvy[1]+U2*lvy[2] 
  dot(jointx, jointy)
})

#create the covariance matrix
cov = cbind(gene = cov1$GeneName,
            ID = cov1$ID,
            covariance = as.numeric(trimws(C)))

Cr = apply(cov1, 1, function(row){
  
  lvx <- matrix(as.numeric(c(row["V1gene"], row["V2gene"])) , nrow = 2, ncol = 1)
  
  lvy <- matrix( as.numeric(c(row["V1sec"], row["V2sec"])) , nrow =2, ncol = 1)
  jointx<-T1*lvx[1]+T2*lvx[2]
  jointy<-U1*lvy[1]+U2*lvy[2] 
  covxy = dot(jointx, jointy)
  varx = dot(jointx, jointx)
  vary = dot(jointy, jointy)
  as.numeric(covxy/sqrt(varx*vary))
})

covcor = as.data.frame(cbind(cov, Cr))
write_xlsx(covcor,"analysis after integration/gene_prot_corrJnt.xlsx")

#find protein level joint correlation with ADAM9

loadings = cov1
lvy1= loadings[,5]
lvy2 = loadings[,6]

#finding joint component of ADAM9
lvy1A9 = as.numeric(cov1[cov1$GeneName== "ADAM9",5])
lvy2A9 = as.numeric(cov1[cov1$GeneName == "ADAM9",6])

jntA9 = U1*lvy1A9 + U2*lvy2A9
varA9 = dot(jntA9,jntA9)

jointsy = matrix(nrow = nrow(lvy1), ncol = 10)
cov_y_a9 = matrix(nrow = nrow(lvy1), ncol = 1)
cor_y_a9 = matrix(nrow = nrow(lvy1), ncol = 1)

for(i in 1:nrow(lvy1)){
  #lvy = c(lvy1[i], lvy2[i])
  jnti =U1*as.numeric(lvy1[i,]) + U2*as.numeric(lvy2[i,]) 
  #print(jnti)
  jointsy[i,] <- jnti
  cova9i = dot(jnti, jntA9)
  cov_y_a9[i,] <- cova9i
  
  vari = dot(jnti,jnti)
  cor_ya9i = cova9i/sqrt(vari*varA9)
  cor_y_a9[i,] <- cor_ya9i
}

cov1 = cbind(cov1,cov_y_a9 )
covcorA9 = cbind(cov1,cor_y_a9 )
write_xlsx(covcorA9,"analysis after integration/covcorA9.xlsx")

##load all sec
secretome <- read_excel("~/BINF/BINF_old/multi-omics/ADAM9 KD HCT116 proteomics/secretomics/perseus output/secretome_results.xlsx")
covcorall <- read_excel("analysis after integration/gene_prot_corrJnt.xlsx")

#find TTs, the proteins with a significant positive correlation between x and y
#one sided fisher, right-handed

n=10

bh1.xy = c()
pp1.xy = c()

cc2 = covcorall

bh1.xy = c()
pp1.xy = c()
for (i in 1:nrow(cc2)){
  r = as.numeric(trimws(cc2$Cr[i]))
  f =  0.5*log((r+1)/(1-r))
  z = sqrt(n-3)*f
  
  #right-handed p
  p = pnorm(z, lower.tail = FALSE)
  pp1.xy = append(pp1.xy,p)
  
}

bh1.xy = p.adjust(pp1.xy, method = "BH")

cc3 = cbind(cc2, bh1.xy)
write_xlsx(cc3,"analysis after integration/alltargets_bh1xy.xlsx")
cc4 = cc3[cc3$bh1.xy < 0.05,]
write_xlsx(cc4, "analysis after integration/tt.xlsx")
cc4 = cc3[cc3$bh1.xy < 0.05,]
write_xlsx(cc4, "analysis after integration/tt.xlsx")

#ptts: no significant correlation between x and y, or negative correlation between x and y

ccneg = cc3[cc3$correlationxy.jnt <0,]
ccposns = cc3[cc3$bh1.xy > 0.1,]
ptt = ccposns
write_xlsx(ptt,"analysis after integration/ptt.xlsx")


#significance of tts and ptts
sigtt = inner_join(cc4, secretome, by ="ID")
sigtt = sigtt[sigtt$FDR < 0.05 ,]
write_xlsx(sigtt,"analysis after integration/ttDEP.xlsx")

sigptt = ptt %>% inner_join(secretome,by = "ID") 
sigptt = sigptt[sigptt$FDR < 0.05 ,]
write_xlsx(sigptt,"analysis after integration/pttDEP.xlsx")

##now find correlation of TTs and PTTs with ADAM9
#load proteome data: the 2 joint lvs, joint scores, orth lv, orth score
setwd("~/BINF/o2pls/tr-sec/output")
Yj_lv <- as.data.frame(read_excel("loading_values_SecJoint.xlsx"))
Yj_sc <- as.data.frame(read_excel("scores_SecJoint.xlsx"))

Yo_lv <- read_excel("loading_values_SecOrth.xlsx")
Yo_lv = as.data.frame(Yo_lv)
Yo_sc <- as.data.frame(read_excel("scores_secOrth.xlsx"))

Yj <- as.matrix(Yj_sc[,2:3]) %*% t(as.matrix(Yj_lv[,3:4]))
Yo <- as.matrix(Yo_sc[,2:3]) %*% t(as.matrix(Yo_lv[,3:4]))  

Yjt = t(Yj)

#bind with names
Yjcomps = cbind(Yj_lv[,c(1,2)], Yjt)
Yocomp = cbind(Yj_sc[,1], Yo)

#Yj and Yo are the joint and orthogonal components
#add them to get predicted values

Ypred = Yj + Yo

#dropping Y from the names since we are going to be working with proteins alone now

Ynames = Yjcomps[,c(1,2)]
Yp2 = rbind(t(Ynames), Ypred)
Yp3 = as.data.frame(t(Yp2))

corr.y.A9 = c()
corr.pvalue = c()
which(Yp3$`GeneName` == "ADAM9") #1
A9 = Yp3[1,]
A9data = as.numeric(A9[, c(-1, -2)])

for(i in 1:nrow(Yp3)){
  idata = as.numeric(Yp3[i, c(-1, -2)])
  cori = cor.test(A9data, idata, method = "pearson")
  corr.y.A9 = append(corr.y.A9, cori$estimate)
  corr.pvalue = append(corr.pvalue, cori$p.value)
  
}

Yp4 = cbind(Yp3, corr.y.A9,corr.pvalue )
#Yp4 gives the correlation between the predicted values for all proteins

#filter for ptts
ptt <- read_excel("analysis after integration/pttDEP.xlsx") 
#ptts that are not DEP filtered
anyDuplicated(ptt)

ptt_pred = inner_join(ptt,Yp4, by = "ID" )
ptt_pred2 = ptt_pred[,-9:-19]
ptt_pred2 = ptt_pred2 %>% rename(corr.pvalue.ya9 = corr.pvalue)

write_xlsx(ptt_pred2,"analysis after integration/ptt_all_predictedcors.xlsx")
pttneg = ptt_pred2[ptt_pred2$corr.y.A9 < 0,]
write_xlsx(pttneg,"analysis after integration/pttneg.xlsx")
pttnegsig = pttneg[pttneg$corr.pvalue.ya9 < 0.05,]
pttneg$corrBH = p.adjust(pttneg$corr.pvalue.ya9, method = "BH")
pttnegsigBH = pttneg[pttneg$corrBH< 0.05,]

write_xlsx(pttnegsigBH,"analysis after integration/pttnegsigBH.xlsx")

pttpos = ptt_pred2[ptt_pred2$corr.y.A9 > 0,]
write_xlsx(pttpos, "analysis after integration/pttpos.xlsx")
pttpossig = pttpos[pttpos$corr.pvalue.ya9 < 0.05,]
pttpos$corrBH = p.adjust(pttpos$corr.pvalue.ya9, method = "BH")
pttpossigBH = pttpos[pttpos$corrBH< 0.05,]

write_xlsx(pttpossigBH,"analysis after integration/pttpossigBH.xlsx")

##tt
tt <- as.data.frame(read_excel("analysis after integration/ttDEP.xlsx"))
tt_pred = inner_join(tt, Yp4, by = c("ID"))
tt_pred2 = tt_pred[,c(-9:-19)]
ttneg = tt_pred2[tt_pred2$corr.y.A9 < 0,]
write_xlsx(ttneg,"analysis after integration/ttneg.xlsx")
ttnegsig = ttneg[ttneg$corr.pvalue < 0.05,]
ttneg$corrBH = p.adjust(ttneg$corr.pvalue, method = "BH")
ttnegsigBH = ttneg[ttneg$corrBH< 0.05,]

write_xlsx(ttnegsigBH, "analysis after integration/ttnegsigBH.xlsx")

ttpos = tt_pred2[tt_pred2$corr.y.A9 > 0,]
write_xlsx(ttpos,"analysis after integration/ttpos.xlsx")
ttpossig = ttpos[ttpos$corr.pvalue < 0.05,]
ttpos$corr.BH = p.adjust(ttpos$corr.pvalue, method = "BH")
ttpossigBH = ttpos[ttpos$corr.BH < 0.05,]
write_xlsx(ttpossigBH,"analysis after integration/ttpossigBH.xlsx")

##find which tt are degs
mRNA_All_Results <- read_excel("~/BINF/BINF_old/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx")
mrna = as.data.frame(mRNA_All_Results[,c(1,14,17)])
ttdegs = inner_join(tt_pred2, mrna, by = c("gene"="GeneName"))

write_xlsx(ttdegs,"tt_deg.xlsx")

###surface annotate pttpossig
pttpos_idmap <- read_excel("analysis after integration/pttpos_idmap.xlsx")
pttpos_idmap = as.data.frame(pttpos_idmap)

pttpos_idmap$`Subcellular location [CC]` <- gsub(
  "\\{.*?\\}"," ",pttpos_idmap$`Subcellular location [CC]`
)

pttpos_idmap$`Subcellular location [CC]` <- 
  gsub("Note.*?\\.", "", pttpos_idmap$`Subcellular location [CC]`)

pttpos_idmap$`Subcellular location [CC]` <- 
  gsub("SUBCELLULAR LOCATION:", "", pttpos_idmap$`Subcellular location [CC]`)

#make subset of 'membrane' and 'secreted' annotated
uniprot_suf_secr <- pttpos_idmap %>% 
  filter(grepl("membrane", pttpos_idmap$`Subcellular location [CC]`,
               ignore.case = TRUE) |
           grepl("secreted", pttpos_idmap$`Subcellular location [CC]`,
                 ignore.case = TRUE, ))

write_xlsx(uniprot_suf_secr,"analysis after integration/surface annotation/uniprot_suf_secr.xlsx")
'''
#FC cut off
setwd("~/BINF/o2pls/tr-sec/output/analysis after integration/fc1.5")
ttpos1.5 = ttpossigBH[abs(ttpossigBH$log2FC) > log2(1.5), ]
write_xlsx(ttpos1.5, "ttpos1.5.xlsx")

ttneg1.5 = ttnegsigBH[abs(ttnegsigBH$log2FC) > log2(1.5), ]
write_xlsx(ttneg1.5, "ttneg1.5.xlsx")

##comparing correlations between rwa data and joint components
#finding correlations based on raw data



pttpos1.5 = pttpossigBH[abs(pttpossigBH$log2FC) > log2(1.5), ]
write_xlsx(pttpos1.5, "pttpos1.5.xlsx")

pttneg1.5 = pttnegsigBH[abs(pttnegsigBH$log2FC) > log2(1.5), ]
write_xlsx(pttneg1.5, "pttneg1.5.xlsx")
''' 



