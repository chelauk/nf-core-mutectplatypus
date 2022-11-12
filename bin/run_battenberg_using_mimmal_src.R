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

