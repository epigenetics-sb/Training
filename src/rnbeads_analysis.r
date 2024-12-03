

##################################################################################################

#preparing sample sheet and data import
#all samples 
#path to data.files
idat.dir<-file.path("/projects/mnt/deep-cluster/walter/projects/researchers/jil/projects/deNBI/deNBI_WS24/data/idats")

#outputPath
report.dir<-"/projects/researchers/jil/projects/deNBI/20241129_WS24/analysis/rnb"
#dir.create(report.dir)
#samplesheet file path
sample.annotation<-file.path("/projects/researchers/jil/projects/deNBI/20241129_WS24/data/20241202_SampleSheet_SelectedSamples.txt")
data.type<-"idat.dir"


library(RnBeads)
library(data.table)

help("rnb.options") #get info on all available RnB parameters.
rnb.options(assembly = "hg19")
rnb.options(identifiers.column = "SampleID")
rnb.options2xml(pretty=TRUE)
rnb.options(import.table.separator = "\t")

#import rnbead set
rnbSet = rnb.run.import(data.source = c(idat.dir, sample.annotation), 
                        data.type = data.type,
                        dir.reports = report.dir)
rnb = rnbSet$rnb.set

#########################################################

#let's inspect some features of our dataset

#see samplesheet
head(pheno(rnb))

#raw methylation values
mm <- meth(rnb)

#plot raw methylation vbalues histogram
hist(mm[,5], col="steelblue", breaks=50)

#are regions summarized? how?
summarized.regions(rnb)

#inspect available annotation
anno = rnb.annotation2data.frame(rnb.get.annotation("probes450"))
head(anno)

#see annotation available for promoters
annot.promoters <- annotation(rnb, type="promoters")
head(annot.promoters)

#get methylation for promoters (on a subset of samples and promoters just for illustration)
meth(rnb, type="promoters", row.names=TRUE, i=1:5, j=1:3)

#inspect coverage (number of beads per probe)
nbead <- covg(rnb, row.names=TRUE)
nbead[1:5,1:3]

#detection p-values
pvals <- dpval(rnb, row.names=TRUE)

#check control probes
rnb.plot.control.boxplot(rnb)

#some specific control
rnb.plot.control.boxplot(rnb, "BISULFITE CONVERSION I")

#Negative control boxplots are generated with the following command:
rnb.plot.negative.boxplot(rnb)

# barplot of a selected control probe
control.meta.data <- rnb.get.annotation("controls450")
ctrl.probe<-paste0(unique(control.meta.data[["Target"]])[2], ".1")
rnb.plot.control.barplot(rnb, ctrl.probe)

#access to full qc results
qc_data = qc(rnb)

#use some genotyping probes to check for genomic sample similarity
snp.probes = anno[grep("rs", rownames(anno)), ]

rnb.plot.snp.heatmap(rnb)
#rnb.plot.snp.barplot(dataset = rnb, probeID = rownames(snp.probes) [2], writeToFile = TRUE)

#set some options as desired before running the complete quality check module
rnb.options(qc.snp.boxplot=TRUE)
rnb.options(import.sex.prediction = TRUE)
rnb.options(qc.cnv = TRUE) # calculate copy number variations

#see all available options
help("rnb.options") #get info on all available RnB parameters.

####################################################################

#run complete quality report with specific settings

#QC
rnb.options(qc.boxplots = TRUE)
rnb.options(qc.barplots = TRUE)
rnb.options(qc.negative.boxplot = TRUE)
rnb.options(qc.snp.distances = TRUE)
rnb.options(qc.snp.boxplot = TRUE)
rnb.options(qc.snp.barplot = TRUE)
rnb.options(qc.sample.batch.size = 50)
rnb.options(qc.coverage.plots = FALSE)
rnb.options(qc.coverage.threshold.plot = 1:10)
rnb.options(qc.coverage.histograms = FALSE)
rnb.options(qc.coverage.violins = FALSE)

rnb.run.qc(rnb.set = rnb, dir.reports = report.dir) #this will takes some minutes

#####################################################################

#PREPROCESSING
#preprocessing typically includes filtering of undesired datapoints (mostly sites) and data normalization
rnb.unprocessed = rnb #first we store the unprocessed rnbeads object as a back-up
save.rnb.set(rnb.unprocessed, "/projects/researchers/jil/projects/deNBI/20241129_WS24/analysis/rnb/rnbSet_Unprocessed", archive = TRUE)

nrow(meth(rnb.unprocessed)) # the number of sites in the unfiltered object

#some examples of how sites to be removed from dataset
# Remove probes outside of CpG context
rnb.set.filtered <- rnb.execute.context.removal(rnb.set = rnb.unprocessed, contexts = NULL )$dataset
nrow(meth(rnb.set.filtered)) # the number of CpG sites in the unfiltered object
# SNP filtering allowing no SNPs in the probe sequence (range 3 bp)
rnb.set.filtered <- rnb.execute.snp.removal(rnb.set = rnb.set.filtered, snp = "3")$dataset
# Removal of CpG sites in the unfiltered object
# that contain a SNP in the range of 3bp
nrow(meth(rnb.set.filtered))
# Remove CpGs on sex chromosomes
rnb.set.filtered <- rnb.execute.sex.removal(rnb.set = rnb.set.filtered)$dataset
nrow(meth(rnb.set.filtered))
# Remove probes and samples based on a greedycut approach
greedycut.results <- rnb.execute.greedycut(rnb.set = rnb.set.filtered, pval.threshold = 0.05)#$dataset
to_remove = rownames(meth(object = rnb.set.filtered, row.names = TRUE)) [greedycut.results[["sites"]]]
remove.sites(object = rnb.set.filtered, probelist = to_remove)
nrow(meth(rnb.set.filtered))

# Remove probes containing NA for beta values
rnb.set.filtered <- rnb.execute.na.removal(rnb.set.filtered)$dataset
nrow(meth(rnb.set.filtered))
# Remove probes for which the beta values have low standard deviation
rnb.set.filtered <- rnb.execute.variability.removal(rnb.set.filtered, 0.005)$dataset
nrow(meth(rnb.set.filtered))

#we remove sites with any NA or negative values
mm = meth(rnb.set.filtered, row.names = TRUE)
mneg = apply(mm, 1, function(x) any(x <= 0))
mneg = is.na(mneg)
head(mm[mneg, 1:3],10)
rnb.set.filtered = remove.sites(object = rnb.set.filtered, probelist = mneg)
nrow(meth(rnb.set.filtered))
save.rnb.set(rnb.set.filtered, "/projects/researchers/jil/projects/deNBI/20241129_WS24/analysis/rnb/rnbSet_Filtered", archive = TRUE)

####normalization...

#preset some options and then run filtering/normalization with one command
rnb.options(filtering.whitelist = NULL)
rnb.options(filtering.blacklist = NULL)
rnb.options(filtering.snp = "3")
rnb.options(filtering.cross.reactive = FALSE)
rnb.options(filtering.greedycut = TRUE)
rnb.options(filtering.greedycut.pvalue.threshold = 0.05)
rnb.options(filtering.greedycut.rc.ties = "row")
rnb.options(filtering.sex.chromosomes.removal = TRUE)
rnb.options(filtering.missing.value.quantile = 0.8)
rnb.options(filtering.coverage.threshold = 3)
rnb.options(filtering.low.coverage.masking = FALSE)
rnb.options(filtering.high.coverage.outliers = FALSE)
rnb.options(filtering.deviation.threshold = 0)

rnb.options(normalization = NULL)
rnb.options(normalization.method = "bmiq")       #normalization method
rnb.options(normalization.background.method = "methylumi.noob")#background removal yes/no
rnb.options(normalization.plot.shifts = TRUE)

#rnb.set.norm <- rnb.execute.normalization(rnb.set.unfiltered, method="wm.dazen", bgcorr.method="methylumi.noob")
preprocessed = rnb.run.preprocessing(rnb.set = rnb.set.filtered, dir.reports = report.dir) 
rnb = preprocessed$rnb.set

#####################################################################

#inference module

#Tissue deconvolution - reference based celltype estimation using houseman algorithm
ct <- rnb.execute.ct.estimation(rnb, cell.type.column="CellType", test.max.markers=10000, top.markers=500)
rnb.plot.ct.heatmap(ct.obj = ct)

ct$contributions

#immune cell content
immune.content <- rnb.execute.lump(rnb)

#calculate epigenetic age
rnb.execute.age.prediction(object = rnb)

#preset some settings
rnb.options(inference.age.prediction = TRUE)
rnb.options(inference.age.column = "Age")
rnb.options(inference.age.prediction.training = FALSE)
rnb.options(inference.age.prediction.cv = FALSE)
rnb.options(inference.immune.cells = TRUE)
rnb.options(inference.genome.methylation = "Genome-wide methylation")
rnb.options(inference.targets.sva = character())
rnb.options(inference.reference.methylome.column = "CellType")
rnb.options(inference.max.cell.type.markers = 10000)
rnb.options(inference.top.cell.type.markers = 500)
rnb.options(inference.sva.num.method = "leek")

rnb_inference = rnb.run.inference(rnb.set = rnb, dir.reports = report.dir)
rnb = rnb_inference$rnb.set

#let's do the sva analysis as an addendum
sva.obj <- rnb.execute.sva(rnb, cmp.cols = "Group", numSVmethod="be")
s
rnb <- set.covariates.sva(rnb, sva.obj)

#####################################################################

#exploratory module
rnb.run.exploratory(rnb.set = rnb, dir.reports = report.dir)

Sys.setenv("R_ZIPCMD"="/usr/bin/zip")
save.rnb.set(rnb, "/projects/researchers/jil/projects/deNBI/20241129_WS24/analysis/rnb/rnbSet_exploratoryAnalysis", archive = TRUE)

# load.rnb.set

#####################################################################



























