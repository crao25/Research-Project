library(TCGA2STAT)
library(curatedTCGAData)
## library(MultiAssayExperiment)
## library(RaggedExperiment)
## library(ExperimentHub)
library(devtools)
library(DT)

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("CNTools", version = "3.8")
#cnasnp <- getTCGA(disease="OV", data.type="CNA_SNP")
#mut.ov <- getTCGA(disease="OV", data.type="Mutation", type="somatic")

#retrieving OMICS profiles (RNASeq, Met, Meth)
counts.brca <- getTCGA(disease = "BRCA", data.type = "RNASeq", type = "count", clinical = TRUE)
meth.brca <- getTCGA(disease = "BRCA", data.type = "Methylation", type = "450K", clinical = TRUE)
mut.brca <- getTCGA(disease = "BRCA", data.type = "Mutation", type = "somatic", clinical = TRUE)

#getting matched tumor-normal pofiles
mut.brca.tum.norm <- TumorNormalMatch(mut.brca$dat)
counts.brca.tum.norm <- TumorNormalMatch(counts.brca$dat)



#retrieving CNASNP data
test <- curatedTCGAData(diseaseCode = "COAD",
                        assays = "*GISTIC*",
                        dry.run = FALSE)


data <- datatable(head(assays(test)[[2]]))
dat <- datatable(head(as.data.frame(rowData((test)[[2]]))))

# Get RNA-SeqV2 data for LUSC patients
lusc.rnaseq2 <- getTCGA(disease="LUSC", data.type="RNASeq2")
# tumor-normal matched profiles
lusc.rnaseq2.tum.norm <- TumorNormalMatch(lusc.rnaseq2$dat)

lusc.methyl <- getTCGA(disease="LUSC", data.type="Methylation")
lusc.rnaseq2 <- getTCGA(disease="LUSC", data.type="RNASeq2", clinical=TRUE)

met.var <- apply(lusc.methyl$dat, 1, var)
met.data <- subset(lusc.methyl$dat, met.var >= quantile(met.var, 0.99, na.rm=T) & !is.na(met.var))

rnaseq2.var <- apply(log10(1+lusc.rnaseq2$dat), 1, var)
rnaseq.data <- subset(log10(1+lusc.rnaseq2$dat), rnaseq2.var >= quantile(rnaseq2.var, 0.99, na.rm=T) & !is.na(rnaseq2.var))

met.rnaseq2 <- OMICSBind(dat1 = rnaseq.data, dat2= met.data)

library(CCA)
lusc.cc <- rcc(t(met.rnaseq2$X), t(met.rnaseq2$Y), 0.75025, 0.5005)
lusc.cc2 <- comput(t(met.rnaseq2$X), t(met.rnaseq2$Y), lusc.cc)
plotNice2(lusc.cc2, d1=1, d2=2, XY="X", cex=0.7)

meth.lusc <- getTCGA(disease = "LUSC", data.type = "Methylation", type = "450K")
meth.lusc.var <- apply(meth.OV$dat, 1, var)
meth.lusc.data <- subset(meth.lusc$dat, meth.OV.var >= quantile(meth.OV.var, 0.99, na.rm=T) & !is.na(meth.OV.var))

cna.OV <- curatedTCGAData(diseaseCode = "READ",
                          assays = "*GISTIC*",
                          dry.run = FALSE)
cna.OV.var <- apply((assays(cna.OV)[[2]]), 1, var)
cna.OV.data <- subset((assays(cna.OV)[[2]]), cna.OV.var >= quantile(cna.OV.var, 0.99, na.rm=T) & !is.na(cna.OV.var))
rnaseq.OV <- getTCGA(disease = "READ", data.type = "RNASeq", type = "count", clinical = TRUE)
rnaseq.OV.var <- apply(log10(1+rnaseq.OV$dat), 1, var)
rnaseq.OV.data <- subset(log10(1+rnaseq.OV$dat), rnaseq.OV.var >= quantile(rnaseq.OV.var, 0.99, na.rm=T) & !is.na(rnaseq.OV.var))

mut.OV <- getTCGA(disease = "READ", data.type = "Mutation", type = "somatic")
mut.OV.var <- apply(mut.OV$dat, 1, var)
mut.OV.data <- subset(mut.OV$dat, mut.OV.var >= quantile(mut.OV.var, 0.99, na.rm=T) & !is.na(mut.OV.var))