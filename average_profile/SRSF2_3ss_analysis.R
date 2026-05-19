# ================================
# SRSF2 CLIP-seq Binding at 3′ Splice Sites
# ================================

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicFeatures)

# ---- Load CLIP-seq peaks ----
clip_file <- "peaks_id74671_fdr0.05_group_hSRSF2_genome_hg19_ensembl59.bedgraph"
clip_data <- read.delim(clip_file, header = FALSE, sep = "\t", skip = 1)
colnames(clip_data) <- c("Chrom", "Start", "End", "Strand")
clip_data$Strand <- ifelse(clip_data$Strand > 0, "+", "-")

clip_gr <- GRanges(
  seqnames = clip_data$Chrom,
  ranges = IRanges(start = clip_data$End, width = 1),
  strand = clip_data$Strand,
  seqlengths = seqlengths(Hsapiens)
)

clip_fwd <- coverage(clip_gr[strand(clip_gr) == "+"])
clip_rev <- coverage(clip_gr[strand(clip_gr) == "-"])

# ---- Load gene annotations ----
gtf <- "gencode.v19.annotation.gtf.gz"
annot <- import.gff2(gtf)

pc_ids   <- import.bed("pc_ref.bed")$name
linc_ids <- import.bed("linc_ref.bed")$name

pc_trs   <- annot[!is.na(annot$transcript_id) & annot$transcript_id %in% pc_ids]
linc_trs <- annot[!is.na(annot$transcript_id) & annot$transcript_id %in% linc_ids]

pc_in   <- unlist(intronsByTranscript(makeTxDbFromGRanges(pc_trs)))
linc_in <- unlist(intronsByTranscript(makeTxDbFromGRanges(linc_trs)))

# ---- Define 3′ splice sites (±200 nt around intron end) ----
FLANK <- 200

make_3ss <- function(introns) {
  GRanges(
    seqnames = seqnames(introns),
    ranges = ifelse(strand(introns) == "+",
                    IRanges(end(introns) - FLANK, end(introns) + FLANK),
                    IRanges(start(introns) - FLANK, start(introns) + FLANK)),
    strand = strand(introns)
  )
}

pc_3ss <- make_3ss(pc_in)
linc_3ss <- make_3ss(linc_in)

# ---- Compute coverage ----
senseRle <- function(gr) {
  fwd <- clip_fwd[gr[strand(gr) == "+"]]
  rev <- revElements(clip_rev[gr[strand(gr) == "-"]])
  Reduce("+", fwd) + Reduce("+", rev)
}

antisRle <- function(gr) {
  fwd <- clip_fwd[gr[strand(gr) == "-"]]
  rev <- revElements(clip_rev[gr[strand(gr) == "+"]])
  Reduce("+", fwd) + Reduce("+", rev)
}

sumCov_3ss <- data.frame(
  Pos = seq_len(2 * FLANK + 1),
  SRSF2_sPC = as.vector(senseRle(pc_3ss)),
  SRSF2_aPC = as.vector(antisRle(pc_3ss)),
  SRSF2_sNC = as.vector(senseRle(linc_3ss)),
  SRSF2_aNC = as.vector(antisRle(linc_3ss))
)

# ---- Plotting ----
pdf("SRSF2_3ss_profile.pdf", width = 10, height = 6)
par(mfrow = c(1, 2), bg = "white")

matplot(sumCov_3ss$Pos, sumCov_3ss[, c("SRSF2_sPC", "SRSF2_aPC")],
        type = "l", lty = 1, lwd = 2, col = c("black", "gray70"),
        main = "SRSF2 Binding - Protein-Coding 3'ss",
        xlab = "Position around 3'-ss", ylab = "CLIP-seq signal", xaxt = "n")
axis(1, at = c(1, FLANK + 1, 2 * FLANK + 1), labels = c("-200", "3'ss", "+200"))
abline(v = FLANK + 1, lty = 2)

matplot(sumCov_3ss$Pos, sumCov_3ss[, c("SRSF2_sNC", "SRSF2_aNC")],
        type = "l", lty = 1, lwd = 2, col = c("black", "gray70"),
        main = "SRSF2 Binding - lincRNA 3'ss",
        xlab = "Position around 3'-ss", ylab = "CLIP-seq signal", xaxt = "n")
axis(1, at = c(1, FLANK + 1, 2 * FLANK + 1), labels = c("-200", "3'ss", "+200"))
abline(v = FLANK + 1, lty = 2)

legend("topright", legend = c("Sense", "Antisense"), col = c("black", "gray70"),
       lty = 1, lwd = 2, bty = "n")
dev.off()
