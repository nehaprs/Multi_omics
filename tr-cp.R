library(dplyr)
library(OmicsPLS)
library(readr)
library(readxl)
library(writexl)
library(openxlsx)
library(readxl)
library(pracma)

setwd("~/BINF/multi-omics/transcriptome-cp redo")

rna <- read_excel("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx", 
                  col_types = c("text", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "skip", "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "skip", "skip", 
                                "skip", "skip", "skip"))

#which(rna$control4 == 0)
'''
> sum(rna$control1 == 0)
[1] 0
> sum(rna$control2 == 0)
[1] 1 at C10orf105
> sum(rna$control3 == 0)
[1] 0
> sum(rna$control4 == 0)
[1] 1 at PIWIL2
> sum(rna$control5 == 0)
[1] 0
> sum(rna$control6 == 0)
[1] 0
> sum(rna$ADAM9KD1 == 0)
[1] 0
> sum(rna$ADAM9KD2 == 0)
[1] 0
> sum(rna$ADAM9KD3 == 0)
[1] 0
> sum(rna$ADAM9KD4 == 0)
[1] 0
> sum(rna$ADAM9KD5 == 0)
[1] 0
> sum(rna$ADAM9KD6 == 0)
[1] 0

'''
sum(rna$control6 == 0)

cp <- read_excel("~/BINF/cellular proteomics/results_volcano_cleaned.xlsx", 
                 col_types = c("text", "skip", "text", 
                               "skip", "skip", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric"))
rnan = rna[,-1]
cpn = cp[,c(-1,-2)]
#samples as rows and genes/protein as columns
rnaT = t(rnan)
cpT = t(cpn)
#which(cp$`Gene Names`=="C10orf105")
#need to center data before cross validation
rnaTS = scale(rnaT, center = TRUE, scale = TRUE)
cpTS = scale(cpT, center = TRUE, scale = TRUE)

crossval_o2m_adjR2(rnaTS,cpTS,1:5,0:5,0:5,nr_folds = 5)
"
  MSE n nx ny
1 1.896183 1  0  5
2 1.787606 2  0  2
3 1.771414 3  1  0
4 1.897528 4  3  0
5 1.838378 5  0  2
"


crossval_o2m(rnaTS,cpTS,2:4,0:2,0:1,5)
"
*******
Minimal 5-CV error is at ax=1 ay=1 a=2 
*******
Minimum MSE is 1.71855 
"
#change rownames so that they match

row.names(rnaTS) <- as.list(row.names(cpTS))
#fit o2pls
fit = o2m(rnaTS, cpTS, 2,1,1)
summary(fit)
write_file(fit, "o2plsfit")
head(fit)

#extract loadings
loadings = loadings(fit)
loadings_xjoint = loadings(fit,"Xjoint")
write.table(loadings_xjoint, "loading_values_TransJoint")
loadings_yjoint = loadings(fit,"Yjoint")
write.table(loadings_yjoint,"loading_values_cpJoint")
loadings_xorth = loadings(fit,"Xorth")
write.table(loadings_xorth,"loading_values_Transorth")
loadings_yorth = loadings(fit,"Yorth")
write.table(loadings_yorth,"loading_values_cpOrth")

scores = scores(fit)

scores_xjoint = scores(fit,"Xjoint")
write.table(scores_xjoint, "scores_TransJoint.xlsx")
scores_yjoint = scores(fit,"Yjoint")
write.table(scores_yjoint,"scores_cpJoint.xlsx")
scores_xorth = scores(fit,"Xorth")
write.table(scores_xorth,"scores_Transorth.xlsx")
scores_yorth = scores(fit,"Yorth")
write.table(scores_yorth,"scores_cpOrth.xlsx")

print(fit)
write_xlsx(as.data.frame(fit$E), "trans_residuals.xlsx")
write_xlsx(as.data.frame(fit$Ff), "cp_residuals.xlsx")
saveRDS(fit,"tr-cpfit.rds")

#find genes which occur in both datasets

lv_TransJnt <- read_excel("lv_TransJnt.xlsx")
lv_cpjnt <- read_excel("lv_cpjnt.xlsx")
comms = inner_join(lv_TransJnt, lv_cpjnt, by = c("GeneName" = "Gene Names"))
common_names = comms[,c(2,6)]

######post o2PLS analysis
#load transcriptome joint scores
T <- read_excel("scores_Transjnt.xlsx")
T = as.matrix(T)

T1= as.numeric(trimws(T[,2]))
T2= as.numeric(trimws(T[,3]))
TT = cbind(T1,T2)

covT=matrix(c(dot(T1,T1),dot(T1,T2),dot(T2,T1),dot(T2,T2)),2,2)

#load protein joint scores
U <- read_excel("scores_cpjnt.xlsx")
U <- as.matrix(U)
rownames(U) <- U[,1]  
U <- U[,-1]  

U1 =  as.numeric(trimws(U[,1]))
U2 =  as.numeric(trimws(U[,2]))
UU = cbind(U1, U2)
covU=matrix(c(dot(U1,U1),dot(U1,U2),dot(U2,U1),dot(U2,U2)),2,2)

colnames(T)[colnames(T) == "V1"] <- "V1gene"
colnames(T)[colnames(T) == "V2"] <- "V2gene"


colnames(U)[colnames(U) == "V1"] <- "V1cp"
colnames(U)[colnames(U) == "V2"] <- "V2cp"
#load joint loading values
#transcriptome
LX <- read_excel("lv_TransJnt.xlsx")
#cellular proteome
LY <- read_excel("lv_cpjnt.xlsx")

B = matrix(c(0.788, 0, 0, 0.727), nrow = 2, ncol = 2)

cov1 = inner_join(LX, LY, by =c("GeneName" = "Gene Names"))
#in cov1, X is transcriptome and y is proteome
cov1 = cov1 %>% rename("V1gene" = "V1.x")
cov1 = cov1 %>% rename("V2gene" = "V2.x")

cov1 = cov1 %>% rename("V1cp" = "V1.y")
cov1 = cov1 %>% rename("V2cp" = "V2.y")
diagT = diag(diag(covT))

Dq=diag(c(0.788*dot(T1,T1),0.727*dot(T2,T2)))
Dn =  diagT %*% B

#calculate the covariance for each row
'''
C = apply(cov1, 1, function(row){
  
  lvx <- matrix(as.numeric(c(row["V1gene"], row["V2gene"])) , nrow = 2, ncol = 1)
  
  lvy <- matrix( as.numeric(c(row["V1cp"], row["V2cp"])) , nrow = 2, ncol = 1)
  as.numeric(t(lvx) %*% diagT %*% B %*% lvy)
})

#create the covariance matrix
cov = cbind(gene = cov1$GeneName,
            ID = cov1$ID,
            covariance = as.numeric(trimws(C)))

Cr = apply(cov1, 1, function(row){
  
  lvx <- matrix(as.numeric(c(row["V1gene"], row["V2gene"])) , nrow = 2, ncol = 1)
  
  lvy <- matrix( as.numeric(c(row["V1cp"], row["V2cp"])) , nrow = 2, ncol = 1)
  cv = as.numeric(t(lvx) %*% diagT %*% B %*% lvy)
  varx = t(lvx) %*% covT %*% lvx
  vary = t(lvy) %*% covU %*% lvy
  
  as.numeric(cv/sqrt(varx*vary))
})

covcor = as.data.frame(cbind(cov, Cr))
write_xlsx(covcor,"analysis after integration/gene_prot_corrJnt.xlsx")
'''
C = apply(cov1, 1, function(row){
  
  lvx <- matrix(as.numeric(c(row["V1gene"], row["V2gene"])) , nrow = 2, ncol = 1)
  
  lvy <- matrix( as.numeric(c(row["V1cp"], row["V2cp"])) , nrow = 2, ncol = 1)
  jointx<-T1*lvx[1]+T2*lvx[2]
  jointy<-U1*lvy[1]+U2*lvy[2]
  dot(jointx, jointy)
})

#create the covariance matrix
cov = cbind(gene = cov1$gene,
            ID = cov1$cp,
            covariance = as.numeric(trimws(C)))

Cr = apply(cov1, 1, function(row){
  
  lvx <- matrix(as.numeric(c(row["V1gene"], row["V2gene"])) , nrow = 2, ncol = 1)
  
  lvy <- matrix( as.numeric(c(row["V1cp"], row["V2cp"])) , nrow = 2, ncol = 1)
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
lvy1 = loadings[,7]
lvy2 = loadings[,8]
length(lvy1)

#finding joint component of ADAM9
lvy1A9 = as.numeric(cov1[cov1$GeneName== "ADAM9",7])
lvy2A9 = as.numeric(cov1[cov1$GeneName == "ADAM9",8])

jntA9 = U1*lvy1A9 + U2*lvy2A9

varA9 = dot(jntA9,jntA9)

jointsy = matrix(nrow = length(lvy1), ncol = 10)
cov_y_a9 = matrix(nrow = length(lvy1), ncol = 1)
cor_y_a9 = matrix(nrow = length(lvy1), ncol = 1)

for(i in 1:length(lvy1)){
  #lvy = c(lvy1[i], lvy2[i])
  jnti = U1*lvy1[i] + U2*lvy2[i]
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


##load all cp
cp1 <- read_excel("~/BINF/cellular proteomics/volcano_cleaned.xlsx")
cp2 = cp1

'''
cp2$pvalue = 10^(- cp2$`-Log(P-value)`)
cp2$FDR = p.adjust(cp2$pvalue, method = "fdr")
write_xlsx(cp2,"~/BINF/cellular proteomics/cp_All_results.xlsx")
'''
cp = cp2

covcorall <- read_excel("analysis after integration/covcorxy_cleaned.xlsx")

# use fisher 2 sided test to find significant corr_y_a9


#find TTs, the proteins with a significant positive correlation between x and y
#one sided fisher, right-handed

bh1.xy = c()
pp1.xy = c()
for (i in 1:nrow(covcorall)){
  r = as.numeric(trimws(covcorall$Corr_xyJnt[i]))
  f =  0.5*log((r+1)/(1-r))
  z = sqrt(n-3)*f
  
  #right-handed p
  p = pnorm(z, lower.tail = FALSE)
  pp1.xy = append(pp1.xy,p)
  
}

bh1.xy = p.adjust(pp1.xy, method = "BH")
cc3 = cbind(covcorall, bh1.xy)
write_xlsx(cc3,"analysis after integration/alltargets_bh1xy.xlsx")

#find TTs, the proteins with a significant positive correlation between x and y
#one sided fisher, right-handed

cc2 = covcorall

bh1.xy = c()
pp1.xy = c()
for (i in 1:nrow(cc2)){
  r = as.numeric(trimws(cc2$Corr_xyJnt[i]))
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

#ptts: no significant correlation between x and y, or negative correlation between x and y
ccneg = cc3[cc3$Corr_xyJnt <0,]
ccposns = cc3[cc3$bh1.xy > 0.1,]
ptt = rbind(ccneg, ccposns)
write_xlsx(ptt,"analysis after integration/ptt.xlsx")

#significance of tts and ptts
sigtt = inner_join(cc4, cp, by ="ID")
sigtt = sigtt[sigtt$FDR < 0.05 ,]
write_xlsx(sigtt,"analysis after integration/ttDEP.xlsx")

sigptt = ptt %>% inner_join(cp,by = "ID") %>% filter(sigptt$FDR < 0.05)
write_xlsx(sigptt,"analysis after integration/pttDEP.xlsx")
#ptt edited in excel to remove double count. upload again
ptt <- as.data.frame(read_excel("analysis after integration/ptt.xlsx"))
sigptt = ptt %>% inner_join(cp,by = "ID") %>% filter(sigptt$FDR < 0.05)
write_xlsx(sigptt,"analysis after integration/pttDEP.xlsx")
##compare with old tts and ptts
pttold <-as.data.frame( read_excel("~/BINF/multi-omics/o2pls_all_data/rnaseq_cp/redo math/diagT/covbyhand/significance of correlation/dep filtered/redo/final results/PTT/ptt.xlsx"))

compareptt = inner_join(pttold, sigptt, by = "ID")
#1791 proteins in common, of 1973 pptold and 1876 new ptt.

##compare tt
ttold <- as.data.frame(read_excel("~/BINF/multi-omics/o2pls_all_data/rnaseq_cp/redo math/diagT/covbyhand/significance of correlation/dep filtered/redo/final results/TT/tt.xlsx"))
comparett = inner_join(ttold, sigtt, by = "ID")
#253 common proteins, of of 323 new and 280 old tts. 


##now find correlation of TTs and PTTs with ADAM9
#load proteome data: the 2 joint lvs, joint scores, orth lv, orth score

Yj_lv <- as.data.frame(read_excel("lv_cpjnt.xlsx"))
Yj_sc <- as.data.frame(read_excel("scores_cpjnt.xlsx"))

Yo_lv <- read_excel("lv_cpOrth.xlsx")
Yo_lv = as.data.frame(Yo_lv)
Yo_sc <- as.data.frame(read_excel("scores_cporth.xlsx"))

Yj <- as.matrix(Yj_sc[,c(2,3)]) %*% t(as.matrix(Yj_lv[,c(4,5)]))
Yo <- as.matrix(Yo_sc[,2]) %*% t(as.matrix(Yo_lv[,c(4)]))                                      

Yjt = t(Yj)

Yjcomps = cbind(Yj_lv[,c(2,3)], Yjt)


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
which(Yp3$`Gene Names` == "ADAM9") #4836
A9 = Yp3[4836,]
A9data = as.numeric(A9[, c(-1, -2)])

#6.10.2024 
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
ptt2 = distinct(ptt)

ptt_pred = inner_join(ptt2,Yp4, by = "ID" )
ptt_pred2 = ptt_pred[,c(-26,-27,-28,-29,-30,-31,-32,-33,-34,-35,-36)]
ptt_pred2 = ptt_pred2 %>% rename(corr.pvalue.ya9 = corr.pvalue)
write_xlsx(ptt_pred2,"analysis after integration/ptt_all_predictedcors.xlsx")
pttneg = ptt_pred2[ptt_pred2$corr.y.A9 < 0,]
write_xlsx(pttneg,"pttneg.xlsx")
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
tt_pred2 = tt_pred[,c(-26,-27,-28,-29,-30,-31,-32,-33,-34,-35,-36)]

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
ttpossigBH = ttpos[ttpos$BH < 0.05,]
write_xlsx(ttpossigBH,"analysis after integration/ttpossigBH.xlsx")

##find which tt are degs
mRNA_All_Results <- read_excel("~/BINF/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx")
mrna = as.data.frame(mRNA_All_Results[,c(1,14,17)])
ttdegs = inner_join(tt_pred2, mrna, by = c("GeneName"))
#all tts are also DEGs


###surface annotate pttnegsig
pttneg_idmap <- as.data.frame(read_excel("analysis after integration/surface annotating/pttneg_idmap.xlsx"))

pttneg_idmap$`Subcellular location [CC]` <- gsub(
  "\\{.*?\\}"," ",pttneg_idmap$`Subcellular location [CC]`
)

pttneg_idmap$`Subcellular location [CC]` <- 
  gsub("Note.*?\\.", "", pttneg_idmap$`Subcellular location [CC]`)

pttneg_idmap$`Subcellular location [CC]` <- 
  gsub("SUBCELLULAR LOCATION:", "", pttneg_idmap$`Subcellular location [CC]`)

#make subset of 'membrane' and 'secreted' annotated
uniprot_suf_secr <- pttneg_idmap %>% 
  filter(grepl("membrane", pttneg_idmap$`Subcellular location [CC]`,
               ignore.case = TRUE) |
           grepl("secreted", pttneg_idmap$`Subcellular location [CC]`,
                 ignore.case = TRUE, ))

write_xlsx(uniprot_suf_secr,"analysis after integration/surface annotating/uniprot_suf_secr.xlsx")

###compare with old results
pttnegold <- read_excel("~/BINF/multi-omics/o2pls_all_data/rnaseq_cp/redo math/diagT/covbyhand/significance of correlation/dep filtered/redo/final results/PTT/pttnegsig.xlsx")

pttcompare = inner_join(pttnegsig, pttnegold, by = "ID")
