#!/usr/bin/env Rscript
# Splicing Index Analysis Script
library(ggplot2)
library(dplyr)

# Read coverage data
included_pc <- read.table("splicing_coverage/included_pc.csv", header=FALSE, col.names=c("GeneSymbol", "Incscore"))
included_lnc <- read.table("splicing_coverage/included_lnc.csv", header=FALSE, col.names=c("GeneSymbol", "Incscore"))
excluded_pc <- read.table("splicing_coverage/excluded_pc.csv", header=FALSE, col.names=c("GeneSymbol", "Exscore"))
excluded_lnc <- read.table("splicing_coverage/excluded_lnc.csv", header=FALSE, col.names=c("GeneSymbol", "Exscore"))

# Merge data
pc_data <- merge(included_pc, excluded_pc, by="GeneSymbol") %>%
  filter(Incscore > 5 & Exscore > 2) %>%
  mutate(si = Incscore / (Incscore + Exscore),
         type = "Protein-coding")

lnc_data <- merge(included_lnc, excluded_lnc, by="GeneSymbol") %>%
  filter(Incscore > 5 & Exscore > 2) %>%
  mutate(si = Incscore / (Incscore + Exscore),
         type = "lncRNA")

# Combine datasets
combined <- rbind(pc_data, lnc_data)

# 1. Splicing Index Comparison
p1 <- ggplot(combined, aes(x=type, y=si, fill=type)) +
  geom_boxplot(alpha=0.8, outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.3) +
  labs(title="Splicing Index Comparison",
       x="Gene Type", y="Splicing Index") +
  theme_minimal()

ggsave("splicing_index_comparison.pdf", p1, width=8, height=6)

# 2. Density Plot
p2 <- ggplot(combined, aes(x=si, fill=type)) +
  geom_density(alpha=0.6) +
  labs(title="Splicing Index Distribution",
       x="Splicing Index", y="Density") +
  theme_minimal()

ggsave("splicing_index_density.pdf", p2, width=8, height=6)

# 3. Statistical Test
test_result <- wilcox.test(si ~ type, data=combined)
sink("statistical_results.txt")
cat("Wilcoxon Rank Sum Test Results:\n")
print(test_result)
cat("\nSummary Statistics:\n")
cat("Protein-coding genes:\n")
summary(pc_data$si)
cat("\nlncRNAs:\n")
summary(lnc_data$si)
sink()

# Save results
write.csv(combined, "splicing_index_results.csv", row.names=FALSE)