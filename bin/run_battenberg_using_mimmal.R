#!/usr/bin/env Rscript

# Running Battenberg using MiMMAl
# 
# Coded by George Cresswell; 2018-05-03 for use with latest versions of ASCAT and Battenberg

# pick up command line args
args = commandArgs(trailingOnly=TRUE)

# Install packages not available in conda

#We require these libraries
if (!require(ASCAT)) stop("Package 'ASCAT' missing\n.")
if (!require(Battenberg)) stop("Package 'Battenberg' missing\n.")
if (!require(CGHcall)) stop("Package 'CGHcall' missing\n.")
if (!require(mixtools)) stop("Package 'mixtools' missing\n.") #Depends on version 1.0.4 or earlier
if (!require(ggplot2)) stop("Package 'ggplot2' missing\n.")
if (!require(cowplot)) stop("Package 'cowplot' missing\n.")
if (!require(ComplexHeatmap)) stop("Package 'ComplexHeatmap' missing\n.")
if (!require(circlize)) stop("Package 'circlize' missing\n.")
if (!require(MiMMAl)) stop("Package 'MiMMAl' missing\n.")

#File and directory informations
SRC_FILE_PATH         = args[1]
OUTPUT_DIR            = args[2]
BAF_PATH              = args[3]
LRR_PATH              = args[4]
CENTROMERES_PATH      = args[5]
SAMPLENAME            = args[6]

#Now to set a very large number of parameters!
ASCAT_GAMMA           = 1                        #We test a range of gammas dictated by this vector
SEGMENTATION_GAMMA    = 10                       #This segmentation gamma calls pure 2+1s of 44 SNPs 50% of the time with sd = 0.033 (4.72)
KMIN                  = 3                        #Segments must be this number of SNPs minimum
CLONALITY_DIST_METRIC = 0                        #The distance metric that fit.copy.number.ASCAT.only uses
ASCAT_DIST_METRIC     = 2                        #The distance metric that ASCAT uses
MIN_PLOIDY            = 1.6                      #We don't want very low ploidy solutions
MAX_PLOIDY            = 5.5                      #We don't want hexaploid solutions
MIN_RHO               = 0.1                      #We don't want very impure solutions
MIN_GOODNESS_OF_FIT   = 0.63                     #We want good fits
BALANCED_THRESHOLD    = 0.51                     #Used by ASCAT for in the distance metric
PRESET                = FALSE                    #We don't preset
PRESET_RHO            = NA                       #We don't preset
PRESET_PSI            = NA                       #We don't preset
SIGLEVEL              = 0.05                     #The significance level (p-value) used for clonality tests
MAXDIST               = 0.01                     #Max distance used in clonality tests

#George's new Battenberg parameters
CGHCALL.NORMALISATION = TRUE                     #Use CGHcall's Log R Ratio normalisation technique
seed                  = 1                        #Seeding for reproducibility
segmentArms           = TRUE                     #Should we segment the individual chromosome arms to aid segmentation?
# minSplitSegDiff       = -0.5                     #Used when forcing centomere splitting to avoid missing segments due to splitting
sd.width              = 1/3                      #Width either side of standard deviation determined by EM to explore in grid search

#Source the functions that need to be read in
source(SRC_FILE_PATH)

#Read in the arm boundary data
armBoundary = read.table(CENTROMERES_PATH,
                         sep="\t",
                         header=TRUE,
                         stringsAsFactors = FALSE)

#Right how long will this take?
start.time = Sys.time()

#Read in the case data from rds files
case.LRR = readRDS(LRR_PATH)
case.BAF = readRDS(BAF_PATH)

#Rename the columns to these names in the LRR/BAF file just for simplicity
print(paste("sample name:", SAMPLENAME))
colnames(case.BAF)[4] = SAMPLENAME
colnames(case.LRR)[4] = SAMPLENAME

#Chromosome names
chroms = unique(case.LRR$Chr)

#Remove the NAs that are produced
case.LRR = na.omit(case.LRR)
case.BAF = case.BAF[which(case.BAF$Name %in% case.LRR$Name),]
  
#Do we do CGHcall normalisation?
if(CGHCALL.NORMALISATION) {
  
  #Subset what we need for CGHcall
  case.LRR = case.LRR[,c(1:3,3,4:ncol(case.LRR))]
  print(head(case.LRR)) 
  #Use the make_cghRaw function which converts to CGHcall format
  case.LRR = make_cghRaw(case.LRR)
  
  #Normalise the data, smoothing outlier and normalising to median
  case.LRR = normalize(case.LRR,
                       method="median",
                       smoothOutliers=TRUE)

  #If female, reconvert the X chromosome from 23
  lrr.chrs = chromosomes(case.LRR)

  #Make 23s, X
  lrr.chrs[which(lrr.chrs==23)] = "X"
  
  #Recreate the original dataframe
  case.LRR = data.frame(featureNames(case.LRR),
                        lrr.chrs,
                        bpstart(case.LRR),
                        copynumber(case.LRR), stringsAsFactors = FALSE)
  
  #Recreate the original dataframe
  colnames(case.LRR)[1:4] = c("Name", "Chr", "Position", SAMPLENAME)
  
  #Recreate the original dataframe
  rownames(case.LRR)      = NULL
  
}

#Randomly flips ~50% of the SNPs to try to make normal BAF ~0.5
set.seed(seed)
flipped.snps = sample(1:nrow(case.BAF), round(nrow(case.BAF) / 2))
case.BAF[flipped.snps,4:ncol(case.BAF)] = 1 - case.BAF[flipped.snps,4:ncol(case.BAF)]
  
#What's happening?
print(paste0("Analysing sample : ",SAMPLENAME))

#Make a directory for the testing
dir.create(OUTPUT_DIR)
setwd(OUTPUT_DIR)

#Record the parameters now that we are aware of them all
makeSink(TRUE)

#Let's do a plot of the LRR and BAF
plot_gw_LRR(LRR = case.LRR[,SAMPLENAME], CHR = case.LRR$Chr, POS = case.LRR$Position, NAM = SAMPLENAME)
plot_gw_BAF(BAF = case.BAF[,SAMPLENAME], CHR = case.LRR$Chr, POS = case.LRR$Position, NAM = SAMPLENAME)

#Let's write out the file for the sample we want to analysis
write.table(case.BAF[,c("Chr","Position",SAMPLENAME)],
            file = paste0(SAMPLENAME, ".mutantBAF.tab"),
            sep="\t",
            quote=FALSE,
            row.names = FALSE,
            col.names = c("Chromosome", "Position", SAMPLENAME))

#Let's write out the file for the sample we want to analysis
write.table(case.LRR[,c("Chr","Position",SAMPLENAME)],
            file = paste0(SAMPLENAME, ".mutantLogR.tab"),
            sep="\t",
            quote=FALSE,
            row.names = FALSE,
            col.names = c("Chromosome", "Position", SAMPLENAME))

#Let's write the a mimiced heterozygousMutBAFs_haplotyped.txt file
write.table(case.BAF[,c("Chr","Position",SAMPLENAME)],
            file = paste0(SAMPLENAME, "_heterozygousBAFs_raw.txt"),
            sep="\t",
            quote=FALSE,
            row.names = FALSE,
            col.names = c("Chromosome", "Position", SAMPLENAME))

#Import raw data again
inputfile=paste0(SAMPLENAME, "_heterozygousBAFs_raw.txt")

#Read in the raw BAF
BAFraw = read.table(inputfile,sep="\t",header=T, stringsAsFactors=F)

#Make a variable for recording BAFoutput
BAFoutput = NULL

#Run through the chromosomes to segment and model them
for (chr in unique(BAFraw[,1])) {
  
  BAFrawchr = BAFraw[BAFraw[,1]==chr,c(2,3)]
  BAFrawchr = BAFrawchr[!is.na(BAFrawchr[,2]),]
  
  BAF = BAFrawchr[,2]
  pos = BAFrawchr[,1]
  names(BAF) = rownames(BAFrawchr)
  names(pos) = rownames(BAFrawchr)
  
  sdev = Battenberg:::getMad(ifelse(BAF<0.5,BAF,1-BAF),k=25)

  # Standard deviation is not defined for a single value
  if (is.na(sdev)) {
    sdev = 0
  }

  # Have bottom cap at 10,000X
  if(sdev<0.002){
    sdev = 0.002
  }
  
  #Here we will implement our mixture modelling to generate segments with non-truncated means
  #Generate mirrored BAF for first segmention
  mBAF = abs(BAF - 0.5) + 0.5
  
  #Take the mBAF for the p arm
  mBAFp     = mBAF[which(BAFrawchr$Position <= armBoundary[armBoundary$Chr==chr,"Boundary"])]
  
  #Take the mBAF for the p arm
  mBAFq     = mBAF[which(BAFrawchr$Position > armBoundary[armBoundary$Chr==chr,"Boundary"])]
  
  #If we are segmenting the arms separately, apply it here, but don't bother if there isn't enough snps in one arm (>=KMIN)
  if(segmentArms & length(mBAFp)>=KMIN & length(mBAFq)>=KMIN) {
    
    #Now we will segment this arm mBAF
    mBAFpsegm = getBafSegMean(mBAFp, KMIN, SEGMENTATION_GAMMA*sdev, minLen = KMIN)
    
    #Now we will segment this arm mBAF
    mBAFqsegm = getBafSegMean(mBAFq, KMIN, SEGMENTATION_GAMMA*sdev, minLen = KMIN)
    
    #Combine them to get the whole chromosome and call it the broken version
    mBAFsegmBroken   = c(mBAFpsegm, mBAFqsegm)
    
    #The result is split via the centromere
    mBAFsegm = mBAFsegmBroken
    
  } else {
    
    #Now we will segment this mBAF
    mBAFsegm = getBafSegMean(mBAF, KMIN, SEGMENTATION_GAMMA*sdev, minLen = KMIN)
    
  }
  
  #Record the segmentation for the chromosome
  BAFsegsChr = cbind(chr, BAFrawchr, mBAFsegm)
  
  #Add it to the BAF output
  BAFoutput  = rbind(BAFoutput, BAFsegsChr)

}

#Give column names
colnames(BAFoutput) = c("Chromosome", "Position", "BAF", "BAFseg")

#Write out
write.table(BAFoutput,
            file = paste0(SAMPLENAME, "_heterozygousBAFs_segmented.txt"),
            sep="\t",
            quote=FALSE,
            row.names = FALSE)

#Run MiMMAl
runMiMMAl(samplename = SAMPLENAME, 
          inputfile = paste0(SAMPLENAME, "_heterozygousBAFs_segmented.txt"), 
          sd.width = sd.width, plot.lit.plot = FALSE, plot.star.plot = FALSE,
          seed = seed)

#Run ASCAT
fit.copy.number(samplename = SAMPLENAME,
                outputfile.prefix = paste(SAMPLENAME, "_", sep=""),
                inputfile.baf.segmented = paste(SAMPLENAME,".BAFphased.txt", sep=""),
                inputfile.baf = paste(SAMPLENAME,".mutantBAF.tab", sep=""),
                inputfile.logr = paste(SAMPLENAME,".mutantLogR.tab", sep=""),
                dist_choice = CLONALITY_DIST_METRIC,
                ascat_dist_choice = ASCAT_DIST_METRIC,
                min.ploidy = MIN_PLOIDY,
                max.ploidy = MAX_PLOIDY,
                min.rho = MIN_RHO,
                min.goodness = MIN_GOODNESS_OF_FIT,
                uninformative_BAF_threshold = BALANCED_THRESHOLD,
                gamma_param = ASCAT_GAMMA,
                use_preset_rho_psi = PRESET,
                preset_rho = PRESET_RHO,
                preset_psi = PRESET_PSI)

# Go over all segments, determine which segments are a mixture of two states and fit a second CN state
callSubclones(sample.name=SAMPLENAME,
              baf.segmented.file=paste(SAMPLENAME,".BAFphased.txt", sep=""),
              logr.file=paste(SAMPLENAME,".mutantLogR.tab", sep=""),
              rho.psi.file=paste(SAMPLENAME, "_rho_and_psi.txt",sep=""),
              output.file=paste(SAMPLENAME,"_subclones.txt", sep=""),
              output.figures.prefix=paste(SAMPLENAME,"_subclones_chr", sep=""),
              output.gw.figures.prefix=paste(SAMPLENAME,"_BattenbergProfile", sep=""),
              masking_output_file=paste(SAMPLENAME,"_segment_masking_details.txt", sep=""),
              chr_names=chroms,
              gamma=ASCAT_GAMMA,
              segmentation.gamma=NA,
              siglevel=SIGLEVEL,
              maxdist=MAXDIST,
              noperms=1000,
              seed = seed)

#OK, now we've done
end.time   = Sys.time()

#How long did that take?
end.time - start.time

