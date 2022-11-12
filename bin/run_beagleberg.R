#!/usr/bin/env Rscript

# Running Battenberg using Beagle and QC'ing with MiMMAl
#
# Coded by George Cresswell; 2019-09-24

# collect command line args

args <- commandArgs(trailingOnly = TRUE)

# We require these libraries
if (!require(ASCAT)) stop("Package 'ASCAT' missing\n.")
if (!require(Battenberg)) stop("Package 'Battenberg' missing\n.")
if (!require(CGHcall)) stop("Package 'CGHcall' missing\n.")
if (!require(mixtools)) stop("Package 'mixtools' missing\n.") # Depends on version 1.0.4 or earlier
if (!require(ggplot2)) stop("Package 'ggplot2' missing\n.")
if (!require(cowplot)) stop("Package 'cowplot' missing\n.")
if (!require(ComplexHeatmap)) stop("Package 'ComplexHeatmap' missing\n.")
if (!require(circlize)) stop("Package 'circlize' missing\n.")
if (!require(MiMMAl)) stop("Package 'MiMMAl' missing\n.")
if (!require(vcfR)) stop("Package 'vcfR' is missing\n.")
if (!require(copynumber)) stop("Package 'copynumber' is missing\n.")

# File and directory informations
SRC_FILE_PATH <- args[1]
OUTPUT_DIR <- args[2]
LRR_PATH <- args[3]
CENTROMERES_PATH <- args[4]
SAMPLENAME <- args[5]
HET_SEQZ <- args[6]
PHASED_VCF_PREFIX <- args[7]

# Now to set a very large number of parameters!
ASCAT_GAMMA <- 1 # Gamma is 1 is sequencing data
SEGMENTATION_GAMMA <- 10 # Segmentation parameter of phased data
KMIN <- 3 # Segments must be this number of SNPs minimum in phased data
PHASING_GAMMA <- 3 # For detecting haplotype blocks this gamma is less than phased gamma
PHASING_KMIN <- 1 # Minimum number of SNPs in a haplotype block
CLONALITY_DIST_METRIC <- 0 # The distance metric that fit.copy.number uses
ASCAT_DIST_METRIC <- 2 # The distance metric that ASCAT uses
MIN_PLOIDY <- 1.6 # We don't want very low ploidy solutions
MAX_PLOIDY <- 5.5 # We don't want hexaploid solutions
MIN_RHO <- 0.95 # This is organoid data so it is pure!
MIN_GOODNESS_OF_FIT <- 0.63 # We want good fits
BALANCED_THRESHOLD <- 0.51 # Used by ASCAT for in the distance metric
PRESET <- FALSE # We don't preset
PRESET_RHO <- NA # We don't preset
PRESET_PSI <- NA # We don't preset
SIGLEVEL <- 0.05 # The significance level (p-value) used for clonality tests
MAXDIST <- 0.01 # Max distance used in clonality tests

# George's new Battenberg parameters
CGHCALL.NORMALISATION <- TRUE # Use CGHcall's Log R Ratio normalisation technique
seed <- 1 # Seeding for reproducibility
segmentArms <- TRUE # Should we segment the individual chromosome arms to aid segmentation?
sd.width <- 1 / 3 # Width either side of standard deviation determined by EM to explore in grid search

# Source the functions that need to be read in
source(SRC_FILE_PATH)

# Replace this function because it's broken
assignInNamespace("ascat.plotSunrise", ascat.plotSunriseReplacement, ns = "ASCAT")
assignInNamespace("callSubclones", callSubclonesReplacement, ns = "Battenberg")

# Read in the arm boundary data
armBoundary <- read.table(CENTROMERES_PATH,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
)

# Edit the chromosome names so they have chr
armBoundary$Chr <- paste0("chr", gsub("chr", "", armBoundary$Chr))

# Right how long will this take?
start.time <- Sys.time()

# Where are we?
ORIGINAL_DIR <- getwd()



#################################################################
##### We start with BAF processing, haplotype based phasing #####
#################################################################



# Read in the data
het_snps_sample <- read.table(HET_SEQZ,
    header = T, stringsAsFactors = F
)

# Add header manually
colnames(het_snps_sample) <- c(
    "chromosome", "position", "base.ref", "depth.normal", "depth.tumor",
    "depth.ratio", "Af", "Bf", "zygosity.normal", "GC.percent",
    "good.reads", "AB.normal", "AB.tumor", "tumor.strand"
)

# Read in the case data from rds files
case.LRR <- readRDS(LRR_PATH)
colnames(case.LRR)[4] <- SAMPLENAME

# Chromosome names
chroms <- unique(case.LRR$Chr)

# Make a MiMMAl output file
mimmal_df <- NULL
# Collect all the phased chromosome data
phased_genome <- NULL

for (chr in paste0("chr", chroms)) {
    # Get only chromosome
    het_snps <- het_snps_sample[het_snps_sample$chromosome == chr, ]

    # Read in the phased data
    phased_snps <- read.vcfR(Sys.glob("*", "phased.", chr, ".vcf.gz"))

    # Immediately remove anything with an indel
    phased_snps <- extract.indels(phased_snps, return.indels = F)

    # Subset the ones in the het snps file
    phased_snps <- phased_snps[getPOS(phased_snps) %in% het_snps$position]

    # Now remove het snps that aren't common in the population
    het_snps <- het_snps[het_snps$position %in% getPOS(phased_snps), ]

    # Now only use ones with LRR data
    het_snps <- het_snps[paste0(het_snps$chromosome, "_", het_snps$position) %in% case.LRR$Name, ]
    phased_snps <- phased_snps[paste0(chr, "_", getPOS(phased_snps)) %in% case.LRR$Name]

    # Make a dataframe of the ref and alt(s) and the phased genotype
    phased_data <- data.frame(
        ref = getREF(phased_snps),
        alt = getALT(phased_snps),
        gt = extract.gt(phased_snps, element = "GT"),
        stringsAsFactors = F
    )

    # Rename columns
    colnames(phased_data) <- c("ref", "alt", "gt")

    # Split out genotypes
    split_gt <- strsplit(phased_data$gt, split = "|")

    # Make chromosome A and B
    phased_data$chromosome_a <- as.numeric(unlist(lapply(split_gt, function(i) i[1])))
    phased_data$chromosome_b <- as.numeric(unlist(lapply(split_gt, function(i) i[3])))

    # Convert numbers to bases
    chromosome_a_bases <- lapply(1:nrow(phased_data), function(i) {
        # What is the number?
        n <- phased_data$chromosome_a[i]

        # If it is zero it is the reference
        if (n == 0) {
            base <- phased_data$ref[i]
        } else {
            # OK, so it is the alternative and it could be multiple things, so we split them out
            alt_bases <- unlist(strsplit(phased_data$alt[i], split = ","))

            # But which alt base is it?
            base <- alt_bases[n]
        }

        # Return the base
        return(base)
    })

    # Convert numbers to bases
    chromosome_b_bases <- lapply(1:nrow(phased_data), function(i) {
        # What is the number?
        n <- phased_data$chromosome_b[i]

        # If it is zero it is the reference
        if (n == 0) {
            base <- phased_data$ref[i]
        } else {
            # OK, so it is the alternative and it could be multiple things, so we split them out
            alt_bases <- unlist(strsplit(phased_data$alt[i], split = ","))

            # But which alt base is it?
            base <- alt_bases[n]
        }

        # Return the base
        return(base)
    })

    # Add to data frame
    phased_data$chromosome_a_base <- unlist(chromosome_a_bases)
    phased_data$chromosome_b_base <- unlist(chromosome_b_bases)

    # Get the BAF for the chromosome
    chromosome_a_baf <- lapply(1:nrow(phased_data), function(i) {
        # What is the base on the chromosome?
        base <- phased_data$chromosome_a_base[i]

        # Split out the two bases in the het.seqz file
        split_seqz_bases <- unlist(strsplit(het_snps[i, "AB.normal"], split = ""))

        # Is it the a allele or b allele by sequenza
        a_b_seqz <- which(split_seqz_bases %in% base)

        # Is it empty and therefore an NA?
        if (length(a_b_seqz) != 1) {
            af <- NA
        } else {
            af <- het_snps[i, a_b_seqz + 6]
        }

        return(af)
    })

    # Get the BAF for the chromosome
    chromosome_b_baf <- lapply(1:nrow(phased_data), function(i) {
        # What is the base on the chromosome?
        base <- phased_data$chromosome_b_base[i]

        # Split out the two bases in the het.seqz file
        split_seqz_bases <- unlist(strsplit(het_snps[i, "AB.normal"], split = ""))

        # Is it the a allele or b allele by sequenza
        a_b_seqz <- which(split_seqz_bases %in% base)

        # Is it empty and therefore an NA?
        if (length(a_b_seqz) != 1) {
            af <- NA
        } else {
            af <- het_snps[i, a_b_seqz + 6]
        }

        return(af)
    })

    # Make a dataframe
    phased_chromosome <- data.frame(
        chr = chr,
        pos = het_snps$position,
        chromosome_a = unlist(chromosome_a_baf),
        chromosome_b = unlist(chromosome_b_baf)
    )

    # Remove the NA
    phased_chromosome <- na.omit(phased_chromosome)

    # For segmentation we need to calculate this sdev value that Battenberg uses to adjust gamma
    sdev <- Battenberg:::getMad(ifelse(phased_chromosome$chromosome_a < 0.5, phased_chromosome$chromosome_a, 1 - phased_chromosome$chromosome_a), k = 25)

    # Have bottom cap at 10,000X
    if (sdev < 0.002) {
        sdev <- 0.002
    }

    print(paste0("sdev: ", sdev, " chr: ", chr))

    # Seg'em
    phased_chromosome$hap_segs <- Battenberg:::selectFastPcf(phased_chromosome$chromosome_a, kmin = PHASING_KMIN, gamma = PHASING_GAMMA * sdev, yest = T)$yhat

    # Flip them
    phased_chromosome$phased <- ifelse(phased_chromosome$hap_segs < 0.5, 1 - phased_chromosome$chromosome_a, phased_chromosome$chromosome_a)


    ####### New addition ######
    chr.LRR <- case.LRR[case.LRR$Chr == gsub("chr", "", chr), ]
    chr.LRR <- chr.LRR[chr.LRR$Position %in% phased_chromosome$pos, ]
    seg_input <- cbind(phased_chromosome[, c("chr", "pos", "phased")], chr.LRR[, 4])
    pq <- getPQ(df = seg_input[, 1:2], armBoundary)
    colnames(seg_input)[3:4] <- c("BAF", "LRR")
    seg_input$chr <- gsub("chr", "", seg_input$chr)
    msseg <- multipcf(seg_input, arms = pq, gamma = SEGMENTATION_GAMMA * 1, fast = T, digits = 8, normalize = F)
    segment <- rep(msseg$BAF, times = msseg$n.probes)
    print(table(segment < 0.45) / length(segment))
    pass <- segment >= 0.45
    seg_input <- seg_input[pass, ]
    pq <- pq[pass]
    phased_chromosome <- phased_chromosome[pass, ]
    msseg <- multipcf(seg_input, arms = pq, gamma = SEGMENTATION_GAMMA * 2, fast = T, digits = 8, normalize = F)
    segment <- rep(msseg$BAF, times = msseg$n.probes)
    ####### Fin addition ######

    # Add segmentation of phasing
    phased_chromosome$phase_segs <- segment

    # Make a dataframe for MiMMAl input
    mimmal_df_chr <- data.frame(chr = chr, pos = phased_chromosome$pos, BAF = phased_chromosome$chromosome_a, BAFseg = phased_chromosome$phase_segs)

    # Collate
    mimmal_df <- rbind(mimmal_df, mimmal_df_chr)

    # Collect the phased data per chromosome too
    phased_genome <- rbind(phased_genome, phased_chromosome)
}

# Make a plot
p <- ggplot(phased_genome, aes(pos, phased)) +
    geom_point(cex = 0.1) +
    geom_step(aes(y = phase_segs, col = "#e50000")) +
    ylim(0, 1) +
    ylab("Phased BAF") +
    xlab("Position") +
    ggtitle(SAMPLENAME) +
    geom_hline(yintercept = c(0.5), linetype = "dashed", color = "black") +
    geom_hline(yintercept = c(0, 1), color = "black") +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5))

# Facet them into chromosomes
p + facet_grid(chr ~ .) + theme(legend.position = "none")

dir.create(OUTPUT_DIR)
ggsave(paste0(OUTPUT_DIR, "/", SAMPLENAME, "_phased_chromosomes.png"), width = 10, height = 40)



############################################
##### Now we will process the LRR data #####
############################################



# Remove the NAs that are produced
case.LRR <- na.omit(case.LRR)

# Match them up
case.LRR <- case.LRR[case.LRR$Name %in% paste0(phased_genome$chr, "_", phased_genome$pos), ]

# Do something about the rownames
rownames(case.LRR) <- NULL

# Do we do CGHcall normalisation?
if (CGHCALL.NORMALISATION) {
    # Subset what we need for CGHcall
    case.LRR <- case.LRR[, c(1:3, 3, 4:ncol(case.LRR))]

    # Use the make_cghRaw function which converts to CGHcall format
    case.LRR <- make_cghRaw(case.LRR)

    # Normalise the data, smoothing outlier and normalising to median
    case.LRR <- normalize(case.LRR,
        method = "median",
        smoothOutliers = TRUE
    )

    # If female, reconvert the X chromosome from 23
    lrr.chrs <- chromosomes(case.LRR)

    # Make 23s, X
    lrr.chrs[which(lrr.chrs == 23)] <- "X"

    # Recreate the original dataframe
    case.LRR <- data.frame(featureNames(case.LRR),
        lrr.chrs,
        bpstart(case.LRR),
        copynumber(case.LRR),
        stringsAsFactors = FALSE
    )

    # Recreate the original dataframe
    colnames(case.LRR)[1:4] <- c("Name", "Chr", "Position", SAMPLENAME)

    # Recreate the original dataframe
    rownames(case.LRR) <- NULL
}

# What's happening?
print(paste0("Analysing case : ", SAMPLENAME))

# Move to our output dir
setwd(OUTPUT_DIR)

# Record the parameters now that we are aware of them all
makeSink(TRUE)

# Let's do a plot of the LRR and BAF
plot_gw_LRR(LRR = case.LRR[, SAMPLENAME], CHR = case.LRR$Chr, POS = case.LRR$Position, NAM = SAMPLENAME)

# Let's write out the file for the sample we want to analysis
write.table(case.LRR[, c("Chr", "Position", SAMPLENAME)],
    file = paste0(SAMPLENAME, ".mutantLogR.tab"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = c("Chromosome", "Position", SAMPLENAME)
)



###################################
##### Time to QC the phasing  #####
###################################



# Write it out
write.table(mimmal_df,
    file = paste0(SAMPLENAME, "_mimmal_input.txt"),
    sep = "\t", quote = F, row.names = F
)

# Run MiMMAl
runMiMMAl(
    samplename = paste0(SAMPLENAME, "_mimmal"), inputfile = paste0(SAMPLENAME, "_mimmal_input.txt"),
    plot.sd.den = F, plot.lit.plot = F, plot.star.plot = F, plot.transformed = F
)

# Get the results
mimmal <- read.table(paste0(SAMPLENAME, "_mimmal.BAFphased.txt"), sep = "\t", header = T)

# Collapse segments
seg_collapse <- rle(phased_genome$phase_segs)

# Segment indices
segment_compare <- data.frame(
    index1 = c(1, c(cumsum(seg_collapse$lengths) + 1)[-length(seg_collapse$lengths)]),
    index2 = cumsum(seg_collapse$lengths),
    phasing_median = seg_collapse$value
)

# Add on the MiMMAl result
segment_compare$mimmal_median <- mimmal$BAFseg[segment_compare$index1]

# How long are segments?
segment_compare$seg_length <- segment_compare$index2 - segment_compare$index1 + 1

# Make a plot of that
ggplot(segment_compare, aes(x = mimmal_median, y = phasing_median, size = seg_length)) +
    geom_point() +
    ylab("Phasing") +
    xlab("Mixture Model") +
    ggtitle(paste0(SAMPLENAME, " - MiMMAl QC")) +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1) +
    geom_hline(yintercept = 0.5, lty = "dotted") +
    geom_vline(xintercept = 0.5, lty = "dotted") +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(SAMPLENAME, "_MiMMAl_QC.pdf"), width = 8, height = 7)



#############################################
##### Now we will run ASCAT-Battenberg  #####
#############################################



# Make Battenberg input
battenberg_input <- data.frame(
    Chromosome = gsub("chr", "", phased_genome$chr),
    Position = phased_genome$pos,
    BAF = phased_genome$chromosome_a,
    BAFphased = phased_genome$phased,
    BAFseg = phased_genome$phase_segs
)

# Write it out
write.table(battenberg_input,
    file = paste0(SAMPLENAME, ".BAFphased.txt"),
    sep = "\t", quote = F, row.names = F
)

# Make a BAF file
case.BAF <- data.frame(
    Chr = gsub("chr", "", phased_genome$chr),
    Position = phased_genome$pos,
    SAMPLENAME = phased_genome$chromosome_a
)

# Let's write out the file for the sample we want to analysis
write.table(case.BAF,
    file = paste0(SAMPLENAME, ".mutantBAF.tab"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = c("Chromosome", "Position", SAMPLENAME)
)

# Run ASCAT
fit.copy.number(
    samplename = SAMPLENAME,
    outputfile.prefix = paste(SAMPLENAME, "_", sep = ""),
    inputfile.baf.segmented = paste(SAMPLENAME, ".BAFphased.txt", sep = ""),
    inputfile.baf = paste(SAMPLENAME, ".mutantBAF.tab", sep = ""),
    inputfile.logr = paste(SAMPLENAME, ".mutantLogR.tab", sep = ""),
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
    preset_psi = PRESET_PSI
)

# Go over all segments, determine which segments are a mixture of two states and fit a second CN state
callSubclonesReplacement(
    sample.name = SAMPLENAME,
    baf.segmented.file = paste(SAMPLENAME, ".BAFphased.txt", sep = ""),
    logr.file = paste(SAMPLENAME, ".mutantLogR.tab", sep = ""),
    rho.psi.file = paste(SAMPLENAME, "_rho_and_psi.txt", sep = ""),
    output.file = paste(SAMPLENAME, "_subclones.txt", sep = ""),
    output.figures.prefix = paste(SAMPLENAME, "_subclones_chr", sep = ""),
    output.gw.figures.prefix = paste(SAMPLENAME, "_BattenbergProfile", sep = ""),
    masking_output_file = paste(SAMPLENAME, "_segment_masking_details.txt", sep = ""), # May need commenting out
    chr_names = chroms,
    gamma = ASCAT_GAMMA,
    segmentation.gamma = NA,
    siglevel = SIGLEVEL,
    maxdist = MAXDIST,
    noperms = 1000,
    seed = seed, # May need commenting out
    sv_breakpoints_file = "NA", # May need commenting out
    calc_seg_baf_option = 1
)
