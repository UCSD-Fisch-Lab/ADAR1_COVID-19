# RNA Editing filtering
library(data.table)
library(plyr)

args <- commandArgs(trailingOnly=TRUE)
alu_maf <- args[1]
nonalu_maf <- args[2]
out_maf <- args[3]

cat(out_maf)
cat("\n")

strand_table_path <- "/expanse/lustre/projects/csd691/kfisch/RNA_editing_pipeline/db/gene_strand_hg38.txt"
strand_table <- unique(fread(strand_table_path, stringsAsFactors = FALSE, sep = " ", header=FALSE))
names(strand_table) <- c("strand", "Hugo_Symbol")

readMaf <- function(maf_file, type){
  maf <- fread(maf_file, sep="\t", stringsAsFactors=FALSE, skip=1)
  if(nrow(maf)==0){
    return(NA)
  }
  maf$Alu <- type
  maf
}

alu <- readMaf(alu_maf, "Alu")
nonalu <- readMaf(nonalu_maf, "Nonalu")

maf <- rbind(alu, nonalu)

maf <- subset(maf, is.na(gnomAD_AF))
maf <- subset(maf, dbSNP_RS %in% c("novel", ""))
maf$VAF <- maf$t_alt_count/maf$t_depth
maf <- merge(maf, strand_table, by = "Hugo_Symbol")
maf$mutsig <- paste0(maf$Reference_Allele, ">", maf$Tumor_Seq_Allele2)
maf <- subset(maf, (mutsig == "A>G" & strand=="+") | (mutsig == "T>C" & strand=="-"))
maf$NCBI_Build <- "hg38"
maf$Strand <- NULL
write.table(maf, file=out_maf, row.names = FALSE, quote = FALSE, sep = "\t")
