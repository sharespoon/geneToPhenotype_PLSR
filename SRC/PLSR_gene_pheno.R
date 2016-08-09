library(pls)
library(chillR)
library(pheatmap)

setwd("/Volumes/Sharon Wang/Data/2015-05 RNAseq/2016-02HSPH_RNAseq_analysis/geneToPhenotype_PLSR/INPUT/")

# input data for model building. Mitochondrial data were from 201508 experiment. CTG from 201506 experiment. RNAseq from 201505 exp. 
alldrugctg=read.table("alldrugctg_forplsr.txt", header = FALSE, as.is = TRUE)
alldrugGenes<-read.csv("ExpconditionvsGeneName_normalizedbyAdamWay.csv", as.is = TRUE, sep = ',')
rtPCRgeneName <- read.csv("/Volumes/Sharon Wang/Data/qRT_PCR/analysis20160228_for2015exp/20151223PCR_analysis/OUTPUT/rtPCR_geneNames.csv",header = FALSE, stringsAsFactors = FALSE)
MMP=read.csv("201508_MitoMP_lowDensity.csv",as.is = TRUE, sep = ',')
mito_Variance=read.csv("201508_mitoVariance_lowDensity.csv",as.is=TRUE, sep = ',')
mito_contrast=read.csv("201508_MitoContrast_lowDensity.csv",as.is = TRUE, sep = ',')
mito_homo=read.csv("201508_MitoHomogeneity_lowDensity.csv",as.is = TRUE, sep = ',')
mito_correlation=read.csv("201508_MitoCorrelation_lowDensity.csv",as.is=TRUE, sep = ',')
nuclei_area=read.csv("201508_nucleiArea_lowDensity.csv",as.is=TRUE, sep = ',')
nuclei_roundness=read.csv("201508_nucleiroundness_lowDensity.csv",as.is=TRUE, sep = ',')

# clean up and label the data
rownames(alldrugGenes)=alldrugGenes[,1];
alldrugGenes2=alldrugGenes[,-1]

rownames(alldrugctg)=rownames(alldrugGenes2);
colnames(alldrugctg)='CTG'

genedata=as.matrix(alldrugGenes2);
alldrugctg=as.matrix(alldrugctg);


# just look at the RT-PCR genes from the RNAseq data
AllgeneNames <- colnames(alldrugGenes2)
rtPCRgeneName2 <-rtPCRgeneName[,2]
rtPCRgeneName_clean <- rtPCRgeneName2[rtPCRgeneName2 %in% AllgeneNames] # remove the genes that was not present in the filtered normalized RNAseq_counts
rtPCRgeneName_notpresent <- rtPCRgeneName2[-which(rtPCRgeneName2 %in% AllgeneNames)]

rtPCRgene_forPLSR <-alldrugGenes2[,rtPCRgeneName_clean] 
rtPCRgene_forPLSR <- as.matrix(rtPCRgene_forPLSR)

# remove the 6h data as this condition was not present in the phenotype data. Here, I used hard-coding. 
Allgene_no6h=alldrugGenes2[-c(2,8,14,20),]
rtPCRgene_no6h=rtPCRgene_forPLSR[-c(2,8,14,20),]

## set up basic crossvalidation PLSR models to regress gene expression against phenotypical measurement

# Wrap pls cross-validation analysis into a function
geneToPhenotypePLSR_Val=function(gene,pheno) {
  df=data.frame(Y=I(pheno),X=I(gene))
  plsr_model=plsr(Y~X,ncomp=10,data=df,validation="LOO")
  summary(plsr_model)
  plot(RMSEP(plsr_model),legendpos="topright",main = NULL)
  title(main = "Root Mean Square Error of PLSR model")
  quartz()
  plot(plsr_model,ncomp=10,asp=1,line=TRUE)
  quartz()
  plot(R2(plsr_model,estimate="all"),legendpos="topright",main = NULL)
  title(main = "R2 and Q2 of PLSR model")
}

# model with cellular ATP values measured in 2015-06 experiments
X1=rtPCRgene_forPLSR    # rtPCR gene expression in RNAseq as regressors
Y1=alldrugctg  # Y is the phenoptye to regress

geneToPhenotypePLSR_Val(X1,Y1)
geneToPhenotypePLSR_Val(genedata,Y1)  # genedata (X) is all the gene expression data. 

# model with mitochondrial membrane potential and texture data from 2015-08 experiments
mito_MP=as.matrix(MMP$fold.change.in.Mitochondrial.membrane.potential)
mito_Var=as.matrix(mito_Variance$fold.change.in.mitochondrial.sum.of.variance)
mito_Con=as.matrix(mito_contrast$mito_contrast)
mito_Homo=as.matrix(mito_homo$mito_homogeneity)
mito_Corr=as.matrix(mito_correlation$mito_correlation)
Y2=as.matrix(data.frame(mito_MP,mito_Var,mito_Con,mito_Homo,mito_Corr))
X2=as.matrix(rtPCRgene_no6h)
X2_2=as.matrix(Allgene_no6h)

geneToPhenotypePLSR_Val(X2,Y2)
geneToPhenotypePLSR_Val(X2_2,Y2)   # X2_2 is all the gene expression data with no 6 hour. 

# model with nuclei data from 201508 live imaging experiment
nuclei_area=as.matrix(nuclei_area$Nuclei_area)
nuclei_round=as.matrix(nuclei_roundness$Nuclei_roundness)
Y3=as.matrix(data.frame(nuclei_area,nuclei_round))  # should I include nuclei number? it didn't change much

geneToPhenotypePLSR_Val(X2_2,Y3)

## after cross-validation, build well-fitted models based on PLSR models with no cross-validation and the right number of component.
## Also, fit all the gene data this time. 

geneToPhenotypePLSR=function(gene,pheno,component) {
  df=data.frame(Y=I(pheno),X=I(gene))
  plsr_model=plsr(Y~X,ncomp=component,data=df,method="oscorespls") # only "oscorespls" method is compatible with VIP function
  summary(plsr_model)
  par(mfrow=c(2,2))
  plot(RMSEP(plsr_model),legendpos="topright",main = NULL)
  title(main = "Root Mean Square Error of PLSR model")
  plot(plsr_model,ncomp=component,asp=1,line=TRUE)
  plot(plsr_model,"loadings",comps=1:2,legendpos="topleft")
  plot(R2(plsr_model,estimate="all"),legendpos="topright",main = NULL)
  title(main = "R2 and Q2 of PLSR model")
  VIP=VIP(plsr_model)
  return(VIP)
}

# From the PLSR cross validation analysis, we found the components to be:
component_MMP=10
component_mitoCon=2
component_mitoHomo=7
component_nuround=5
component_nuArea=2

# only fit the mito membrane potential data and plot in different ways and calculate VIP scores
mitoMP_PLSR=geneToPhenotypePLSR(X2_2,mito_MP,component_MMP)
VIPmitoMP_above4=mitoMP_PLSR[,colSums(mitoMP_PLSR>4)>1]
pheatmap(VIPmitoMP_above4)

# fit mito_contrast data
mitoCon_PLSR=geneToPhenotypePLSR(X2_2,mito_Con,component_mitoCon)
pheatmap(mitoCon_PLSR)

# fit mito Homogeneity data
mitoHomo_PLSR=geneToPhenotypePLSR(X2_2,mito_Homo,component_mitoHomo)
pheatmap(mitoHomo_PLSR)

# fit nuclei roundness data
nuRound_PLSR=geneToPhenotypePLSR(X2_2,nuclei_round,component_nuround)
VIPnuRound_above4=nuRound_PLSR[,colSums(nuRound_PLSR>4)>1]
pheatmap(VIPnuRound_above4)

# fit nuclei area data
nuArea_PLSR=geneToPhenotypePLSR(X2_2,nuclei_area,component_nuArea)
pheatmap(nuArea_PLSR)

# look at the intersect of genes, consider draw a venn diagram
mitoMP_nuRound_intersect=intersect(colnames(VIPmitoMP_above4),colnames(VIPnuRound_above4))






# you can export the gene names with a high VIP by
# write.csv(colnames(VIP)[VIP[1,]>2],"VIP of first component bigger than 2.csv")



### SCRATCH SPACE
#try a random number of ctg to test the model
#ctg_random=alldrugctg
#ctg_random[,]=runif(10,0,1)
#rtPCRgene_ctgrandom_PLSR=data.frame(CTG=I(ctg_random),geneexpression=I(rtPCRgene_forPLSR))

# if you have bigger sample size, you could also consider split the intial dataset into a training set and a test set. This was done in the example of pls package. 

# try reformat X to consider different time points. Only the data at the specific time were non-zero, but others were zero. This didn't produce a better model. 
# the Q2 was bad after component 3. 

# geneToPhenotype=data.frame(Y=I(Y),X=I(X))
#geneToPhenotype_PLSR=plsr(Y ~ X, ncomp=10, data=geneToPhenotype,validation="LOO")
#summary(geneToPhenotype_PLSR)

#plot(RMSEP(geneToPhenotype_PLSR),legendpos="topright")
#quartz()
#geneToPhenotype_PLSR2=plsr(Y ~ X, ncomp=6, data=geneToPhenotype)
#plot(geneToPhenotype_PLSR2,ncomp=6,asp=1,line=TRUE)
# this is R^2
#plot(R2(geneToPhenotype_PLSR, estimate = "train"))
# this is Q2
#plot(R2(geneToPhenotype_PLSR,estimate="all"),legendpos = "topright")

