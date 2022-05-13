devtools::load_all("/scratch/datasets/ukbiobank/ukbkings")
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
prj_dir <- args[1]
psam_file <- args[2]
chrom <- args[3]
greedy <- "/scratch/groups/ukbiobank/Edinburgh_Data/Software/tools/GreedyRelated/GreedyRelated"


psam <- read_delim(
    file = psam_file,
    delim = "\t",
    col_names = c("fid", "iid", "sex"),
    col_types = "iii",
    skip = 1
)

df <- bio_gen_related_remove(
    project_dir = prj_dir,
    greedy_related = greedy,
    keep = psam$fid
)

df %>%
    select(fid = eid, iid = eid) %>%
    write_delim(
        str_interp("data/samples_relatedness_cutoff_chr${chrom}.remove"),
        delim = " ",
        col_names = FALSE
    )