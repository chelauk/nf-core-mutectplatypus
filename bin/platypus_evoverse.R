#!/usr/bin/env Rscript
# Load libs:

library(evoverse)

args <- commandArgs(TRUE)
my_segments <- paste0(args[1], "/", args[2])
my_sample <- args[2]
my_vcf <- args[3]
my_drivers <- args[4]

x <- vcfR::read.vcfR(my_vcf)
normal <- str_replace(x@meta[grep("normal_sample", x@meta)],
"##normal_sample=", "")

# load TCGA drivers
load(my_drivers)

platypus_calls <- evoparse_platypus_mutations(my_vcf)
fit_cnas <- evoparse_Sequenza_CNAs(my_segments)

snvs <- platypus_calls[[my_sample]]$mutations %>%
  dplyr::filter(FILTER == "PASS")

my_colnames <- colnames(snvs)[4:28]
my_colnames <- c(paste0(my_colnames, ".x"))

snvs <- left_join(snvs, mutect_calls[[normal]]$mutations,
    by = c("chr", "from", "to")) %>%
    select("chr", "from", "to", all_of(my_colnames))

my_colnames <-  my_colnames %>% str_replace(".x", "")
my_colnames <- c("chr", "from", "to", my_colnames)
colnames(snvs) <- my_colnames

snv_drivers <- snvs %>%
    dplyr::mutate(effect = sapply(str_split(INFO, "\\|"), function(x) x[4])) %>%
    dplyr::mutate(driver_label = sapply(str_split(INFO, "\\|"), function(x) x[5])) %>% # nolint
    semi_join(PCAWG_TabS3_drivers, by = c("driver_label" = "gene")) %>%
    dplyr::mutate(is_driver = TRUE) %>%
    select(chr, from, to, ref, alt, driver_label, is_driver, effect) %>%
    dplyr::filter(effect %in% c("HIGH", "MODERATE"))
snv_marked <- full_join(snvs, snv_drivers, by = c("chr", "from", "to", "ref", "alt")) # nolint

cna_calls <- fit_cnas$segments
purity <- fit_cnas$purity

my_calls <- list(cna_calls, snv_marked)

my_filename <- paste0(my_sample, ".pdf")
my_rds <- paste0(my_sample, ".rds")
saveRDS(my_calls, file = my_rds)

fit <- pipeline_qc_copynumbercalls(
    mutations = snv_marked,
    cna = cna_calls,
    purity = purity,
    smooth = TRUE,
    reference = "GRCh38",
    description = my_sample,
    ccf_method = "ROUGH",
    peak_method = "closest",
    purity_error = 0.03
)
plot(fit) %>%
    ggsave(filename = my_filename, height = 17, width = 21)
