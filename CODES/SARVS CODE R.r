
setwd("C:/Users/Stuart/Downloads")

dir.create("figures", showWarnings = FALSE)
dir.create("analysis", showWarnings = FALSE)

library(ggplot2)

cat("\n╔════════════════════════════════════════════════════════╗\n")
cat("║   SARS-CoV-2 ONT ANALYSIS                             ║\n")
cat("║   Dataset: PRJNA675364                                ║\n")
cat("╚════════════════════════════════════════════════════════╝\n\n")

# ============================================
# ANALYSIS 1: GENOME COVERAGE
# ============================================

cat("[1/3] GENOME COVERAGE\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

consensus_files <- list.files("consensus", pattern = "*.fa$", full.names = TRUE)

genome_stats <- data.frame(
  Sample = gsub(".consensus.fa", "", basename(consensus_files)),
  Length_bp = sapply(consensus_files, function(f) {
    lines <- readLines(f)
    seq_lines <- lines[!grepl("^>", lines)]
    sum(nchar(seq_lines))
  })
)

genome_stats$Coverage_pct <- (genome_stats$Length_bp / 29903) * 100

cat(sprintf("✓ Total samples: %d\n", nrow(genome_stats)))
cat(sprintf("✓ Mean coverage: %.1f%%\n", mean(genome_stats$Coverage_pct)))
cat(sprintf("✓ Range: %.1f%% - %.1f%%\n\n", min(genome_stats$Coverage_pct), max(genome_stats$Coverage_pct)))

write.csv(genome_stats, "analysis/01_genome_coverage.csv", row.names = FALSE)

# Plot 1
p1 <- ggplot(genome_stats, aes(x = reorder(Sample, -Coverage_pct), y = Coverage_pct, fill = Coverage_pct)) +
  geom_col(alpha = 0.8) +
  scale_fill_gradient(low = "orange", high = "green") +
  labs(title = "SARS-CoV-2 Genome Coverage by Sample",
       x = "Sample", y = "Coverage (%)", fill = "Coverage %") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave("figures/01_genome_coverage.png", p1, width = 12, height = 6, dpi = 300)
cat("✓ Figure: figures/01_genome_coverage.png\n")
cat("✓ Data: analysis/01_genome_coverage.csv\n\n")

# ============================================
# ANALYSIS 2: VARIANT ANALYSIS
# ============================================

cat("[2/3] VARIANT ANALYSIS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

variant_files <- list.files("variants", pattern = "*.tsv$", full.names = TRUE)

variant_counts <- data.frame(
  Sample = gsub(".variants.tsv", "", basename(variant_files)),
  Variant_Count = sapply(variant_files, function(f) {
    lines <- readLines(f)
    length(lines) - 1
  })
)

cat(sprintf("✓ Total variants: %d\n", sum(variant_counts$Variant_Count)))
cat(sprintf("✓ Mean per sample: %.1f\n", mean(variant_counts$Variant_Count)))
cat(sprintf("✓ Range: %d - %d\n\n", min(variant_counts$Variant_Count), max(variant_counts$Variant_Count)))

write.csv(variant_counts, "analysis/02_variant_counts.csv", row.names = FALSE)

# Plot 2
p2 <- ggplot(variant_counts, aes(x = reorder(Sample, -Variant_Count), y = Variant_Count, fill = Variant_Count)) +
  geom_col(alpha = 0.8) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Variants per Sample",
       x = "Sample", y = "Variant Count", fill = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave("figures/02_variant_counts.png", p2, width = 12, height = 6, dpi = 300)
cat("✓ Figure: figures/02_variant_counts.png\n")
cat("✓ Data: analysis/02_variant_counts.csv\n\n")

# ============================================
# ANALYSIS 3: QUALITY SUMMARY
# ============================================

cat("[3/3] QUALITY SUMMARY\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

quality <- merge(genome_stats, variant_counts, by = "Sample")
quality <- quality[order(-quality$Coverage_pct), ]

cat(sprintf("✓ Quality summary created\n\n"))

write.csv(quality, "analysis/03_quality_summary.csv", row.names = FALSE)

# Plot 3
p3 <- ggplot(quality, aes(x = Coverage_pct, y = Variant_Count, size = Length_bp, color = Coverage_pct)) +
  geom_point(alpha = 0.7) +
  geom_text(aes(label = Sample), nudge_y = 30, size = 2) +
  scale_color_gradient(low = "red", high = "green") +
  scale_size_continuous(range = c(3, 8)) +
  labs(title = "Sample Quality Assessment",
       x = "Coverage (%)", y = "Variant Count", color = "Coverage %", size = "Length") +
  theme_minimal()

ggsave("figures/03_quality_assessment.png", p3, width = 12, height = 8, dpi = 300)
cat("✓ Figure: figures/03_quality_assessment.png\n")
cat("✓ Data: analysis/03_quality_summary.csv\n\n")

# ============================================
# SUMMARY
# ============================================

cat("╔════════════════════════════════════════════════════════╗\n")
cat("║              ANALYSIS COMPLETE                         ║\n")
cat("╚════════════════════════════════════════════════════════╝\n\n")
cat("✓ 3 figures generated in: figures/\n")
cat("✓ 3 data files in: analysis/\n")
cat("✓ Top samples:\n")
print(head(quality[, c("Sample", "Coverage_pct", "Variant_Count")], 5))
cat("\n")