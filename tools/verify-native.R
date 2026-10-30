# Run after installing bambamR, with minimap2, samtools and featureCounts on PATH.
# Uses synthetic reads only; no study data or network access is needed.
library(bambamR)
for (tool in c("minimap2", "samtools", "featureCounts")) {
  stopifnot(nzchar(Sys.which(tool)))
}
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else tempfile("bambamR native ; path ")
if (dir.exists(root) && length(list.files(root, all.files = TRUE, no.. = TRUE))) {
  stop("Choose a new or empty native-test output directory.")
}
dir.create(root, recursive = TRUE, showWarnings = FALSE)
set.seed(803)
reference <- paste(sample(c("A", "C", "G", "T"), 6000, TRUE), collapse = "")
fasta <- file.path(root, "reference genome.fa")
writeLines(c(">chr1", reference), fasta)
annotation <- file.path(root, "genes.gtf")
writeLines(c(
  'chr1\tfixture\texon\t100\t2000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
  'chr1\tfixture\texon\t3500\t5500\t.\t+\t.\tgene_id "g2"; transcript_id "t2";'
), annotation)
reads_dir <- file.path(root, "reads")
dir.create(reads_dir)
starts <- list(s1 = c(200, 600, 1000, 3600), s2 = c(400, 3700, 4100, 4500, 4900))
for (sample in names(starts)) {
  records <- lapply(seq_along(starts[[sample]]), function(i) {
    start <- starts[[sample]][i]
    c(paste0("@", sample, "_", i), substr(reference, start, start + 249L),
      "+", strrep("I", 250L))
  })
  writeLines(unlist(records), file.path(reads_dir, paste0(sample, ".fastq")))
}
bams <- vapply(names(starts), function(sample) {
  bb_align(file.path(reads_dir, paste0(sample, ".fastq")), fasta,
           file.path(root, sample), aligner = "minimap2", threads = 1L,
           extra_args = c("-x", "sr"))$bam
}, character(1))
expected <- matrix(c(3, 1, 1, 4), 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
external <- bb_count_reads(bams, annotation, method = "featureCounts", threads = 1L)
stopifnot(identical(dimnames(external), dimnames(expected)), all(external == expected))
internal <- bb_count_reads(bams, annotation, method = "internal")
stopifnot(identical(dimnames(internal), dimnames(expected)), all(internal == expected))
for (aligner in c("STAR", "HISAT2")) {
  required <- if (aligner == "STAR") "STAR" else c("hisat2", "hisat2-build")
  if (!all(nzchar(Sys.which(required)))) {
    if (identical(Sys.getenv("BAMBAMR_REQUIRE_ALL_ALIGNERS"), "true")) {
      stop("Missing required native integration tool: ", aligner)
    }
    cat("NOT QUALIFIED on this host:", aligner, "\n")
    next
  }
  index <- file.path(root, paste0(aligner, " index"))
  if (aligner == "STAR") {
    dir.create(index)
    command <- "STAR"
    arguments <- c("--runMode", "genomeGenerate", "--genomeDir", index,
                   "--genomeFastaFiles", fasta, "--genomeSAindexNbases", "3",
                   "--genomeChrBinNbits", "9", "--runThreadN", "1",
                   "--outFileNamePrefix", file.path(root, "STAR-index-"))
  } else {
    command <- "hisat2-build"
    arguments <- c(fasta, index)
  }
  log <- file.path(root, paste0(aligner, "-index.log"))
  status <- system2(command, shQuote(arguments), stdout = log, stderr = log)
  stopifnot(status == 0L)
  aligned <- vapply(names(starts), function(sample) {
    bb_align(file.path(reads_dir, paste0(sample, ".fastq")), index,
             file.path(root, paste0(aligner, "-", sample)), aligner = aligner,
             threads = 1L)$bam
  }, character(1))
  counted <- bb_count_reads(aligned, annotation, method = "featureCounts", threads = 1L)
  stopifnot(identical(dimnames(counted), dimnames(expected)), all(counted == expected))
  cat("PASS native alignment:", aligner, "\n")
}
metadata <- data.frame(condition = c("a", "b"), row.names = names(starts))
result <- bb_pipeline(fastq_dir = reads_dir, genome_index = fasta, annotation = annotation,
                      sample_info = metadata, aligner = "minimap2", threads = 1L,
                      output_dir = file.path(root, "pipeline"), skip = c("de", "viz"))
stopifnot(all(result$counts == expected), identical(result$status, "complete"))

# One paired fragment must count once, not once per mate.
sam <- file.path(root, "paired.sam")
writeLines(c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:6000",
  paste("pair1", 99, "chr1", 200, 60, "75M", "=", 300, 175,
        substr(reference, 200, 274), strrep("I", 75), sep = "\t"),
  paste("pair1", 147, "chr1", 300, 60, "75M", "=", 200, -175,
        substr(reference, 300, 374), strrep("I", 75), sep = "\t")
), sam)
paired_bam <- file.path(root, "paired.bam")
status <- system2("samtools", shQuote(c("sort", "-o", paired_bam, sam)))
stopifnot(status == 0L)
paired <- bb_count_reads(paired_bam, annotation, method = "featureCounts", paired = TRUE)
stopifnot(paired["g1", 1] == 1L, paired["g2", 1] == 0L)
cat("PASS: native alignment, exact sample counts, internal/external counting, paired fragments\n")
cat("Test products:", root, "\n")
if (!length(args)) cat("Use an output-directory argument to retain products after R exits.\n")
