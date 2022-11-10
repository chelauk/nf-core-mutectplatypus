#!/usr/bin/env Rscript

# Coded by George Cresswell; 2018-01-18

# We want sequenza for GC normalisation
if (!require(sequenza)) stop("Package 'sequenza' missing\n.")

# Get sample information from command line
args = commandArgs(trailingOnly=TRUE)
patient <- args[0]
sample <- args[1]
is.female <- args[2]

# Get min reads parameter
min.norm.reads <- args[3]

# Read in the file
cn.data <- read.table(Sys.glob("*het.seqz"), header = TRUE, stringsAsFactors = FALSE)

# Make chromosome numbers not contain text, i.e. just 1, not chr1
cn.data$chromosome <- as.character(unlist(lapply(strsplit(as.character(cn.data$chromosome), "chr"), function(x) x[length(x)])))

# Calculate mean depth ratio per GC content bin
gc.stats <- gc.norm(
    x = cn.data$depth.ratio,
    gc = cn.data$GC.percent
)

# make a vector of average depth in bins
gc.vect <- setNames(gc.stats$raw.mean, gc.stats$gc.values)

# Normalise by average depth ratio in bin
cn.data$depth.ratio <- cn.data$depth.ratio / gc.vect[as.character(cn.data$GC.percent)]

# How large is this?
n.snps.total <- nrow(cn.data)

# Remove anything where depth in the normal is less than 25 reads
cn.data <- cn.data[which(cn.data$depth.normal >= min.norm.reads), ]

# How many snps are left after filtering for coverage?
n.snps.filtered <- nrow(cn.data)

# Tell me how many SNPs are left after filtering?
print(paste0(round((n.snps.filtered / n.snps.total) * 100), "% of the SNPs remain after coverage filtering (", n.snps.filtered, " in total)"))

# What about coverage stats?
print(paste0("Median coverage in the normal after filtering is ", median(cn.data$depth.normal), "X"))
print(paste0("Median coverage in the tumour after filtering is ", median(cn.data$depth.tumor), "X"))

# Calculate lrr
cn.data$lrr <- log(cn.data$depth.ratio, base = 2) - median(log(cn.data$depth.ratio, base = 2))

# Which chromosomes?
chrs <- 1:22

# Add X if female
if (is.female) {
    chrs <- c(chrs, "X")
}

# Take only the autosomes in males
cn.data <- cn.data[which(cn.data$chromosome %in% chrs), ]

# Sequenza rounds the Bf, to avoid Battenberg merging segments we are going to recalculate BAF and provide more digits
cn.data$Bf <- round(round(cn.data$good.reads * cn.data$Bf) / cn.data$good.reads, digits = 8)

##### Try this to clean up copy number data in future #####

# # Output for mappability assessment
# output = data.frame(chr = paste0("chr",cn.data$chromosome),
# 					start = as.integer(cn.data$position-1),
# 					end = as.integer(cn.data$position),
# 					name = paste0("chr",cn.data$chromosome,"_",cn.data$position))

# # Output this as input
# write.table(output, file = map_input, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# # Calculate mappability
# system(paste0("bigWigAverageOverBed ",map_ref," ",map_input," ",map_output))

# # Read in results
# mapping = read.table(map_output, sep = "\t", stringsAsFactors=FALSE)

# # Show a pie chart
# pdf(pie_chart)
# pie(table(mapping$V6==1), main = "SNPs passing mappability")
# dev.off()

# # Subset only the well mapped ones
# cn.data = cn.data[mapping$V6==1,]

# # Test that they match
# match = all(paste0("chr",cn.data$chromosome,"_",cn.data$position) == mapping[mapping$V6==1,1])

# # Do they match ids?
# if(!match) {stop("ERROR: Mapping filter doesn't subset correctly")}

# Output the BAF
baf <- data.frame(
    Name = paste0("chr", cn.data$chromosome, "_", cn.data$position),
    Chr = cn.data$chromosome,
    Position = cn.data$position,
    BAF = cn.data$Bf
)

# Output the LRR
lrr <- data.frame(
    Name = paste0("chr", cn.data$chromosome, "_", cn.data$position),
    Chr = cn.data$chromosome,
    Position = cn.data$position,
    LRR = cn.data$lrr
)

# Write out the BAF and LRR
saveRDS(baf, file = paste0(patient, "_", id, ".BAF.rds"))
saveRDS(lrr, file = paste0(patient, "_", id, ".LRR.rds"))
