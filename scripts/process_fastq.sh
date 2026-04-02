#!/bin/bash

# =======================================
# Bash: Download FASTQ Files
# =======================================

# Create directory for storing FASTQ files
mkdir -p /home/user/fastq_files
cd /home/user/fastq_files

# Create a text file with SRR IDs
cat > srr_ids.txt << 'EOF'
SRR7212223
SRR7212224
SRR7212225
SRR7212226
SRR7212227
SRR7212228
SRR7212229
SRR7212230
SRR7212231
SRR7212232
SRR7212233
SRR7212234
SRR7212235
SRR7212236
SRR7212237
SRR7212238
SRR7212239
SRR7212240
SRR7212241
SRR7212242
SRR7212243
SRR7212244
SRR7212245
SRR7212303
SRR7212304
SRR7212305
SRR7212306
SRR7212307
SRR7212308
SRR7212309
SRR7212310
SRR7212311
SRR7212312
SRR7212313
SRR7212314
SRR7212315
SRR7212316
SRR7212317
SRR7212318
SRR7212319
SRR7212320
SRR7212321
SRR7212322
SRR7212323
SRR7212324
SRR7212325
SRR7212326
SRR7212327
SRR7212328
SRR7212329
SRR7212330
SRR7212331
SRR7212332
SRR7212333
SRR7212334
SRR7212335
SRR7212336
EOF

# Download FASTQ files for each SRR ID
while read acc
do
  echo "Downloading $acc"
  wget "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=$acc" -O "${acc}.fastq"
done < srr_ids.txt

# =======================================
# Bash: Check download count
# =======================================
ls *.fastq | wc -l

# =======================================
# R: Start Processing FASTQ Files
# =======================================

Rscript <<'EOF_R'

# Load necessary libraries
library(dada2)

# Define paths
path <- "/home/user/fastq_files"
filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE, recursive = TRUE)

# List raw FASTQ files
fnFs <- sort(list.files(path, pattern = "\\.fastq$", full.names = TRUE))
sample.names <- sub("\\.fastq$", "", basename(fnFs))

# Define filtered output files
filtFs <- file.path(filt_path, paste0(sample.names, "_filt.fastq.gz"))

# Filter and trim the sequences
out <- filterAndTrim(
  fnFs, filtFs,
  truncLen = 240,
  maxN = 0,
  maxEE = 2,
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

# Save filtering summary
write.csv(out, file.path(path, "filtering_summary.csv"))

# Keep only samples that passed filtering
keep <- file.exists(filtFs) & out[, "reads.out"] > 0
fnFs <- fnFs[keep]
filtFs <- filtFs[keep]
sample.names <- sample.names[keep]
out <- out[keep, , drop = FALSE]

# Dereplicate sequences
derepFs <- lapply(filtFs, derepFastq)
names(derepFs) <- sample.names
saveRDS(derepFs, file.path(path, "derepFs_all_samples.rds"))

# Learn error models
errF <- learnErrors(derepFs, multithread = TRUE, verbose = TRUE)
saveRDS(errF, file.path(path, "errF_all_samples.rds"))

# Infer ASVs
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
saveRDS(dadaFs, file.path(path, "dadaFs_all_samples.rds"))

# Build sequence table
seqtab <- makeSequenceTable(dadaFs)
saveRDS(seqtab, file.path(path, "seqtab_all_samples.rds"))
write.csv(seqtab, file.path(path, "ASV_table_prechimera_all_samples.csv"))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

# Save final ASV outputs
saveRDS(seqtab.nochim, file.path(path, "ASV_table_all_samples.rds"))
write.csv(seqtab.nochim, file.path(path, "ASV_table_all_samples.csv"))

# =======================================
# End of R Script
# =======================================
EOF_R

echo "Analysis complete. Results saved in the 'fastq_files' directory."
