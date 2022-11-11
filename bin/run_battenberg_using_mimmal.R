#!/usr/bin/env Rscript

# Running Battenberg using MiMMAl
#
# Coded by George Cresswell; 2018-05-03 for use with latest versions of ASCAT and Battenberg

##	Running Battenberg using MiMMAl source file
##
##	This source file was developed to accompany runBattenbergUsingMiMMAl.R
##
##  Coded by: George Cresswell; 2018-02-27

#This is a function to plot LRRs
plot_gw_LRR <- function(LRR, CHR, POS, NAM) {

  CHR[which(CHR=="X")] = 23

  CHR = as.numeric(CHR)

  png(paste0(NAM,".genome_wide_LRR.png"), height=720, width=1784)
  par(mar=c(5, 6, 4, 2) + 0.1)
  chr.total <- table(CHR)
  chrs      <- unique(CHR)

  plot(LRR,
       pch=20,
       cex=0.02,
       ylim=c(-3,3),
       col="red",
       xaxt="n",
       xlab="Chromosomes",
       ylab="LRR",
       cex.axis=2,
       cex.lab=2,
       cex.main=2,
       main=NAM)

  prev.chr <- NULL
  abline(v=0, lty="dotted")

  for (chr in min(chrs):max(chrs)) {

    current.chr <- sum(prev.chr, chr.total[chr])
    text.offset <- chr.total[chr]/2
    abline(v=current.chr, lty="dotted")
    text(current.chr-text.offset, -1, paste0(chr), cex=1)
    prev.chr <- current.chr

  }

  abline(h=0, lty="dotted")

  dev.off()

}

#This is a function to plot BAFs
plot_gw_BAF <- function(BAF, CHR, POS, NAM) {

  CHR[which(CHR=="X")] = 23

  CHR = as.numeric(CHR)

  png(paste0(NAM,".genome_wide_BAF.png"), height=720, width=1784)
  par(mar=c(5, 6, 4, 2) + 0.1)
  chr.total <- table(CHR)
  chrs      <- unique(CHR)

  plot(BAF,
       pch=20,
       cex=0.02,
       ylim=c(0,1),
       col="red",
       xaxt="n",
       yaxt="n",
       xlab="Chromosomes",
       ylab="B-allele Frequency",
       cex.axis=2,
       cex.lab=2,
       cex.main=2,
       main=NAM)

  axis(side = 2, at = seq(0, 1, 0.5), cex.axis = 2)

  prev.chr <- NULL
  abline(v=0, lty="dotted")

  for (chr in min(chrs):max(chrs)) {

    current.chr <- sum(prev.chr, chr.total[chr])
    text.offset <- chr.total[chr]/2
    abline(v=current.chr, lty="dotted")
    text(current.chr-text.offset, 0, paste0(chr), cex=1)
    prev.chr <- current.chr

  }

  abline(h=0, lty="dotted")
  abline(h=0.5, lty="dotted")
  abline(h=1, lty="dotted")

  dev.off()

}

##  A small wrapper for the selectFastPcf function that includes min segment length,
##    deals with duds (length==0) and only gives the mean segment value per snp
##
##  Coded by George Cresswell; 2017-01-05
getBafSegMean = function(Baf, Kmin, Gamma, minLen=50) {

  #Now we will segment the Baf
  if(length(Baf)<minLen){

    #If the chromosome has less than minLen snps, just take the mean
    BafSegm = rep(mean(Baf),length(Baf))

    #If it has no snps make it NULL
    if(length(Baf)==0) {BafSegm = NULL}

  } else {

    #If it has a normal amount of snps, do your magic!
    BafSegm = Battenberg:::selectFastPcf(Baf,Kmin,Gamma,T)$yhat

  }

  #Return the mean Bafs
  return(BafSegm)

}

##  A function to save for each sample the settings that were used to make it
##
##
##  Coded by George Cresswell; 2017-02-15
makeSink = function(perform = TRUE) {

  if(perform == TRUE) {

    #Save our settings
    sink(file = paste0(SAMPLENAME,"_parameters.txt"))
    cat(paste0("Time/date: ",Sys.time(),"\n"))
    cat(paste0("Path to centromeres: ",CENTROMERES_PATH,"\n"))
    cat(paste0("ASCAT gamma: ",ASCAT_GAMMA,"\n"))
    cat(paste0("Segmentation gamma: ",SEGMENTATION_GAMMA,"\n"))
    cat(paste0("Segmentation kmin: ",KMIN,"\n"))
    cat(paste0("Clonal ASCAT distance metric: ",CLONALITY_DIST_METRIC,"\n"))
    cat(paste0("Normal ASCAT distance metric: ",ASCAT_DIST_METRIC,"\n"))
    cat(paste0("Minimum ploidy: ",MIN_PLOIDY,"\n"))
    cat(paste0("Maximum ploidy: ",MAX_PLOIDY,"\n"))
    cat(paste0("Minimum purity: ",MIN_RHO,"\n"))
    cat(paste0("Minimum goodness of fit: ",MIN_GOODNESS_OF_FIT,"\n"))
    cat(paste0("Distance metric imbalance threshold: ",BALANCED_THRESHOLD,"\n"))
    cat(paste0("Did we preset the rho and psi? ",PRESET,"\n"))
    cat(paste0("Preset purity: ",PRESET_RHO,"\n"))
    cat(paste0("Preset ploidy: ",PRESET_PSI,"\n"))
    cat(paste0("Significance level for clonality test: ",SIGLEVEL,"\n"))
    cat(paste0("Maximum distance used in clonality test: ",MAXDIST,"\n"))
    cat(paste0("CGHcall LRR normalisation performed? ",CGHCALL.NORMALISATION,"\n"))
    cat(paste0("Set seed: ",seed,"\n"))
    cat(paste0("Have we used the centromeres to aid segmentation? ",segmentArms,"\n"))
    # cat(paste0("What is the max distance in segmentation fit before splitting is not accepted? ",minSplitSegDiff,"\n"))
    cat(paste0("Standard deviation distance either side of EM sd: ",sd.width,"\n"))
    sink()

  }

}

#######################################################################################################
# Due to shitty coding in ASCAT that they refuse to replace we have to overwrite the sunrise plotting #
#######################################################################################################
ascat.plotSunriseReplacement = function (d, psi_opt1, rho_opt1, minim = T)
{
    par(mar = c(5, 5, 0.5, 0.5), cex = 0.75, cex.lab = 2, cex.axis = 2)
    if (minim) {
        hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(10,
            "RdBu"))(256))
    }
    else {
        hmcol = colorRampPalette(RColorBrewer::brewer.pal(10,
            "RdBu"))(256)
    }
    image(log(d), col = hmcol, axes = F, xlab = "Ploidy", ylab = "Aberrant cell fraction")
    ploidy_min <- as.numeric(rownames(d)[1])
    ploidy_max <- as.numeric(rownames(d)[nrow(d)])
    purity_min <- as.numeric(colnames(d)[1])
    purity_max <- as.numeric(colnames(d)[ncol(d)])
    axis(1, at = seq(0, 1, by = 1/(ploidy_max-1)), labels = seq(ploidy_min, ploidy_max, by = 1))
    axis(2, at = c(0,1), labels = c(purity_min, purity_max))
    if (psi_opt1 > 0 && rho_opt1 > 0) {
        points((psi_opt1 - ploidy_min)/(ploidy_max - 1), (rho_opt1 -
            purity_min)/(1/purity_max), col = "green", pch = "X",
            cex = 2)
    }
}

callSubclonesReplacement = function (sample.name, baf.segmented.file, logr.file, rho.psi.file,
    output.file, output.figures.prefix, output.gw.figures.prefix,
    chr_names, masking_output_file, max_allowed_state = 250,
    sv_breakpoints_file = NULL, gamma = 1, segmentation.gamma = NA,
    siglevel = 0.05, maxdist = 0.01, noperms = 1000, seed = as.integer(Sys.time()),
    calc_seg_baf_option = 1)
{
    set.seed(seed)
    res = Battenberg:::load.rho.psi.file(rho.psi.file)
    rho = res$rho
    psit = res$psit
    psi = rho * psit + 2 * (1 - rho)
    goodness = res$goodness
    BAFvals = read.table(baf.segmented.file, sep = "\t", header = T,
        stringsAsFactors = F)
    if (colnames(BAFvals)[1] == "X") {
        BAFvals = BAFvals[, -1]
    }
    BAF = BAFvals[, 3]
    BAFphased = BAFvals[, 4]
    BAFseg = BAFvals[, 5]
    SNPpos = BAFvals[, c(1, 2)]
    LogRvals = as.data.frame(read_table_generic(logr.file))
    if (colnames(LogRvals)[1] == "X") {
        LogRvals = LogRvals[, -1]
    }
    ctrans = c(1:length(chr_names))
    names(ctrans) = chr_names
    ctrans.logR = c(1:length(chr_names))
    names(ctrans.logR) = chr_names
    BAFpos = as.vector(ctrans[as.vector(BAFvals[, 1])] * 1e+09 +
        BAFvals[, 2])
    res = Battenberg:::determine_copynumber(BAFvals, LogRvals, rho, psi, gamma,
        ctrans, ctrans.logR, maxdist, siglevel, noperms)
    subcloneres = res$subcloneres
    write.table(subcloneres, gsub(".txt", "_1.txt", output.file),
        quote = F, col.names = T, row.names = F, sep = "\t")
    # res = Battenberg:::merge_segments(subcloneres, BAFvals, LogRvals, rho,
    #     psi, gamma, calc_seg_baf_option)
    # BAFvals = res$bafsegmented
    # res = Battenberg:::determine_copynumber(BAFvals, LogRvals, rho, psi, gamma,
    #     ctrans, ctrans.logR, maxdist, siglevel, noperms)
    # subcloneres = res$subcloneres
    BAFpvals = res$BAFpvals
    res = Battenberg:::mask_high_cn_segments(subcloneres, BAFvals, max_allowed_state)
    subcloneres = res$subclones
    BAFvals = res$bafsegmented
    write.table(BAFvals, file = baf.segmented.file, sep = "\t",
        row.names = F, col.names = T, quote = F)
    masking_details = data.frame(samplename = sample.name, masked_count = 0,
        masked_size = 0, max_allowed_state = max_allowed_state)
    write.table(masking_details, file = masking_output_file,
        quote = F, col.names = T, row.names = F, sep = "\t")
    write.table(subcloneres, output.file, quote = F, col.names = T,
        row.names = F, sep = "\t")
    segment_breakpoints = Battenberg:::collapse_bafsegmented_to_segments(BAFvals)
    if (!is.null(sv_breakpoints_file) & !ifelse(is.null(sv_breakpoints_file),
        TRUE, sv_breakpoints_file == "NA") & !ifelse(is.null(sv_breakpoints_file),
        TRUE, is.na(sv_breakpoints_file))) {
        svs = read.table(sv_breakpoints_file, header = T, stringsAsFactors = F)
    }
    for (chr in chr_names) {
        pos = SNPpos[SNPpos[, 1] == chr, 2]
        if (length(pos) == 0) {
            next
        }
        if (!is.null(sv_breakpoints_file) & !ifelse(is.null(sv_breakpoints_file),
            TRUE, sv_breakpoints_file == "NA") & !ifelse(is.null(sv_breakpoints_file),
            TRUE, is.na(sv_breakpoints_file))) {
            svs_pos = svs[svs$chromosome == chr, ]$position/1e+06
        }
        else {
            svs_pos = NULL
        }
        breakpoints_pos = segment_breakpoints[segment_breakpoints$chromosome ==
            chr, ]
        breakpoints_pos = sort(unique(c(breakpoints_pos$start,
            breakpoints_pos$end)/1e+06))
        png(filename = paste(output.figures.prefix, chr, ".png",
            sep = ""), width = 2000, height = 2000, res = 200)
        Battenberg:::create.subclonal.cn.plot(chrom = chr, chrom.position = pos/1e+06,
            LogRposke = LogRvals[LogRvals[, 1] == chr, 2], LogRchr = LogRvals[LogRvals[,
                1] == chr, 3], BAFchr = BAF[SNPpos[, 1] == chr],
            BAFsegchr = BAFseg[SNPpos[, 1] == chr], BAFpvalschr = BAFpvals[SNPpos[,
                1] == chr], subcloneres = subcloneres, breakpoints_pos = breakpoints_pos,
            svs_pos = svs_pos, siglevel = siglevel, x.min = min(pos)/1e+06,
            x.max = max(pos)/1e+06, title = paste(sample.name,
                ", chromosome ", chr, sep = ""), xlab = "Position (Mb)",
            ylab.logr = "LogR", ylab.baf = "BAF (phased)")
        dev.off()
    }
    subclones = as.data.frame(subcloneres)
    subclones[, 2:ncol(subclones)] = sapply(2:ncol(subclones),
        function(x) {
            as.numeric(as.character(subclones[, x]))
        })
    seg_length = floor((subclones$endpos - subclones$startpos)/1000)
    is_subclonal_maj = abs(subclones$nMaj1_A - subclones$nMaj2_A) >
        0
    is_subclonal_min = abs(subclones$nMin1_A - subclones$nMin2_A) >
        0
    is_subclonal_maj[is.na(is_subclonal_maj)] = F
    is_subclonal_min[is.na(is_subclonal_min)] = F
    segment_states_min = subclones$nMin1_A * ifelse(is_subclonal_min,
        subclones$frac1_A, 1) + ifelse(is_subclonal_min, subclones$nMin2_A,
        0) * ifelse(is_subclonal_min, subclones$frac2_A, 0)
    segment_states_maj = subclones$nMaj1_A * ifelse(is_subclonal_maj,
        subclones$frac1_A, 1) + ifelse(is_subclonal_maj, subclones$nMaj2_A,
        0) * ifelse(is_subclonal_maj, subclones$frac2_A, 0)
    ploidy = sum((segment_states_min + segment_states_maj) *
        seg_length, na.rm = T)/sum(seg_length, na.rm = T)
    Battenberg:::plot.gw.subclonal.cn(subclones = subclones, BAFvals = BAFvals,
        rho = rho, ploidy = ploidy, goodness = goodness, output.gw.figures.prefix = output.gw.figures.prefix,
        chr.names = chr_names)
    cellularity_ploidy_output = data.frame(cellularity = c(rho),
        ploidy = c(ploidy), psi = c(psit))
    cellularity_file = gsub("_subclones.txt", "_cellularity_ploidy.txt",
        output.file)
    write.table(cellularity_ploidy_output, cellularity_file,
        quote = F, sep = "\t", row.names = F)
}


# Define some functions
getPQ = function(df, arm_pos) {

  chrs = unique(df[,1])

  psnqs = unlist(lapply(chrs, function(c) {

    cdf = df[df[,1]==c,]

    b = arm_pos[which(arm_pos[,1] == c),2]

    a = cdf[,2] < b

    a = ifelse(a, "p", "q")

    return(a)

  }))

  return(psnqs)

}



#File and directory informations

OUTPUT_DIR <- "mimmal"
BAF_PATH <- args[0]
LRR_PATH <- args[1]
CENTROMERES_PATH <- args[2]
CASENAME <- args[3]
SAMPLENAME <- args[4]

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
print(paste0("Analysing case : ",CASENAME,", sample : ",SAMPLENAME))

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
              seed = seed,
              sv_breakpoints_file = "NA")

#OK, now we've done
end.time   = Sys.time()

#How long did that take?
end.time - start.time
