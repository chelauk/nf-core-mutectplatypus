library(vcfR)

# get commandline args
args <- commandArgs(trailingOnly = TRUE)

pan_sides <- as.integer(args[4])
map_ref <- args[3]
map_input <- "intermediate_file"
map_output <- paste0(args[1], "_mappability.out")

# Read vcf file
vcf <- read.vcfR(args[2])

# What is the chromosome and position
chrs <- getCHROM(vcf)
pos <- getPOS(vcf)

# Make file ready for reading by bigWigAverageOverBed
output <- data.frame(chr = chrs, start = pos - 1 - pan_sides, end = pos + pan_sides, name = paste0(chrs, "_", pos))

# Write that file out
write.table(output,
    file = map_input, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
)

# Calculate mappability
system(paste0("bigWigAverageOverBed ", map_ref, " ", map_input, " ", map_output))
