# ============================================================
# Hi-C LOOP DOWNSTREAM ANALYSIS
# ============================================================
#
# INPUT:
#   chr1_loop_anchor_pair_CTCF_classified.tsv
#
# INPUT STRUCTURE:
#
#   1  Loop_ID
#   2  Anchor1_chr
#   3  Anchor1_start
#   4  Anchor1_end
#   5  Anchor2_chr
#   6  Anchor2_start
#   7  Anchor2_end
#   8  Anchor1_gene_count
#   9  Anchor2_gene_count
#   10 Anchor1_promoter_count
#   11 Anchor2_promoter_count
#   12 Anchor1_cCRE_promoter_count
#   13 Anchor2_cCRE_promoter_count
#   14 Anchor1_enhancer_count
#   15 Anchor2_enhancer_count
#   16 Anchor1_CTCF_count
#   17 Anchor2_CTCF_count
#   18 Anchor1_type
#   19 Anchor2_type
#   20 Anchor_pair_class
#   21 CTCF_status
#
# OUTPUT:
#   Tables
#   Statistical tests
#   Publication-quality figures
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

required_packages <- c(
  "tidyverse",
  "rstatix"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  
}

library(tidyverse)
library(rstatix)


# ============================================================
# 2. DIRECTORY SETUP
# ============================================================

INPUT_FILE <- "D:/hic/hic_project/downstream/loops/chr1_loop_anchor_pair_CTCF_classified.tsv"

OUTPUT_DIR <- "D:/hic/hic_project/downstream/R_downstream_analysis"

FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
TABLE_DIR  <- file.path(OUTPUT_DIR, "tables")
STAT_DIR   <- file.path(OUTPUT_DIR, "statistics")

dir.create(OUTPUT_DIR, showWarnings = FALSE)
dir.create(FIGURE_DIR, showWarnings = FALSE)
dir.create(TABLE_DIR, showWarnings = FALSE)
dir.create(STAT_DIR, showWarnings = FALSE)


# ============================================================
# 3. READ DATA
# ============================================================

loops <- read.delim(
  INPUT_FILE,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 4. BASIC QC
# ============================================================

cat("\n")
cat("============================================================\n")
cat("BASIC DATA QC\n")
cat("============================================================\n")

cat("Number of rows    :", nrow(loops), "\n")
cat("Number of columns :", ncol(loops), "\n")
cat("Unique Loop IDs   :", n_distinct(loops$Loop_ID), "\n")


if (nrow(loops) != n_distinct(loops$Loop_ID)) {
  
  warning("Duplicate Loop IDs detected!")
  
} else {
  
  cat("Loop ID check     : PASS\n")
  
}


if (ncol(loops) == 21) {
  
  cat("Column check      : PASS\n")
  
} else {
  
  warning("Expected 21 columns.")
  
}


# ============================================================
# 5. CALCULATE LOOP-LEVEL REGULATORY BURDEN
# ============================================================

loops <- loops %>%
  
  mutate(
    
    total_genes =
      Anchor1_gene_count +
      Anchor2_gene_count,
    
    total_promoters =
      Anchor1_promoter_count +
      Anchor2_promoter_count,
    
    total_cCRE_promoters =
      Anchor1_cCRE_promoter_count +
      Anchor2_cCRE_promoter_count,
    
    total_enhancers =
      Anchor1_enhancer_count +
      Anchor2_enhancer_count,
    
    total_CTCF =
      Anchor1_CTCF_count +
      Anchor2_CTCF_count,
    
    genomic_distance =
      abs(Anchor2_start - Anchor1_start),
    
    genomic_distance_kb =
      genomic_distance / 1000,
    
    genomic_distance_mb =
      genomic_distance / 1000000
  )


# ============================================================
# 6. BASIC SUMMARY TABLE
# ============================================================

summary_table <- loops %>%
  
  summarise(
    
    loops = n(),
    
    total_genes = sum(total_genes),
    
    mean_genes_per_loop =
      mean(total_genes),
    
    total_promoters =
      sum(total_promoters),
    
    mean_promoters_per_loop =
      mean(total_promoters),
    
    total_cCRE_promoters =
      sum(total_cCRE_promoters),
    
    mean_cCRE_promoters_per_loop =
      mean(total_cCRE_promoters),
    
    total_enhancers =
      sum(total_enhancers),
    
    mean_enhancers_per_loop =
      mean(total_enhancers),
    
    total_CTCF =
      sum(total_CTCF),
    
    mean_CTCF_per_loop =
      mean(total_CTCF),
    
    median_loop_distance_kb =
      median(genomic_distance_kb),
    
    mean_loop_distance_kb =
      mean(genomic_distance_kb)
  )


write.table(
  summary_table,
  file = file.path(
    TABLE_DIR,
    "overall_loop_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 7. ANCHOR FUNCTIONAL STATUS
# ============================================================

anchor_types <- bind_rows(
  
  loops %>%
    count(Anchor1_type) %>%
    rename(type = Anchor1_type),
  
  loops %>%
    count(Anchor2_type) %>%
    rename(type = Anchor2_type)
  
) %>%
  
  group_by(type) %>%
  
  summarise(
    Anchor_Count = sum(n),
    .groups = "drop"
  ) %>%
  
  mutate(
    Percentage =
      100 * Anchor_Count /
      sum(Anchor_Count)
  )


write.table(
  anchor_types,
  file = file.path(
    TABLE_DIR,
    "anchor_functional_status.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 8. FIGURE 1
# ANCHOR FUNCTIONAL STATUS
# ============================================================

p1 <- ggplot(
  anchor_types,
  aes(
    x = type,
    y = Anchor_Count
  )
) +
  
  geom_col() +
  
  labs(
    title = "Functional classification of loop anchors",
    x = "Anchor functional status",
    y = "Number of anchors"
  ) +
  
  theme_classic(base_size = 14)


ggsave(
  file.path(
    FIGURE_DIR,
    "Figure1_anchor_functional_status.png"
  ),
  p1,
  width = 7,
  height = 5,
  dpi = 300
)


ggsave(
  file.path(
    FIGURE_DIR,
    "Figure1_anchor_functional_status.pdf"
  ),
  p1,
  width = 7,
  height = 5
)


# ============================================================
# 9. LOOP-PAIR DISTRIBUTION
# ============================================================

loop_pairs <- loops %>%
  
  count(
    Anchor_pair_class,
    name = "Loops"
  ) %>%
  
  mutate(
    Percentage =
      100 * Loops / sum(Loops)
  ) %>%
  
  arrange(desc(Loops))


write.table(
  loop_pairs,
  file = file.path(
    TABLE_DIR,
    "loop_pair_distribution.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 10. FIGURE 2
# LOOP-PAIR DISTRIBUTION
# ============================================================

p2 <- ggplot(
  loop_pairs,
  aes(
    x = reorder(Anchor_pair_class, -Loops),
    y = Loops
  )
) +
  
  geom_col() +
  
  labs(
    title = "Loop functional architecture",
    x = "Anchor-pair class",
    y = "Number of loops"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


ggsave(
  file.path(
    FIGURE_DIR,
    "Figure2_loop_pair_distribution.png"
  ),
  p2,
  width = 9,
  height = 6,
  dpi = 300
)


ggsave(
  file.path(
    FIGURE_DIR,
    "Figure2_loop_pair_distribution.pdf"
  ),
  p2,
  width = 9,
  height = 6
)


# ============================================================
# 11. CTCF DISTRIBUTION
# ============================================================

ctcf_distribution <- loops %>%
  
  count(
    CTCF_status,
    name = "Loops"
  ) %>%
  
  mutate(
    Percentage =
      100 * Loops / sum(Loops)
  )


write.table(
  ctcf_distribution,
  file = file.path(
    TABLE_DIR,
    "CTCF_distribution.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 12. FIGURE 3
# CTCF DISTRIBUTION
# ============================================================

p3 <- ggplot(
  ctcf_distribution,
  aes(
    x = CTCF_status,
    y = Loops
  )
) +
  
  geom_col() +
  
  labs(
    title = "CTCF-associated loop distribution",
    x = "CTCF status",
    y = "Number of loops"
  ) +
  
  theme_classic(base_size = 14)


ggsave(
  file.path(
    FIGURE_DIR,
    "Figure3_CTCF_distribution.png"
  ),
  p3,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# 13. CTCF × LOOP CLASS
# ============================================================

ctcf_loop_table <- table(
  loops$Anchor_pair_class,
  loops$CTCF_status
)

write.table(
  as.data.frame.matrix(ctcf_loop_table),
  file = file.path(
    TABLE_DIR,
    "CTCF_by_loop_class.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 14. CHI-SQUARE TEST
# ============================================================

chi_test <- chisq.test(ctcf_loop_table)

capture.output(
  chi_test,
  file = file.path(
    STAT_DIR,
    "CTCF_loop_class_chisq_test.txt"
  )
)


# ============================================================
# 15. CTCF PERCENTAGE BY LOOP CLASS
# ============================================================

ctcf_percentage <- loops %>%
  
  count(
    Anchor_pair_class,
    CTCF_status
  ) %>%
  
  group_by(
    Anchor_pair_class
  ) %>%
  
  mutate(
    Percentage =
      100 * n / sum(n)
  ) %>%
  
  ungroup()


write.table(
  ctcf_percentage,
  file = file.path(
    TABLE_DIR,
    "CTCF_percentage_by_loop_class.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 16. FIGURE 4
# CTCF STATUS BY LOOP CLASS
# ============================================================

p4 <- ggplot(
  ctcf_percentage,
  aes(
    x = Anchor_pair_class,
    y = Percentage,
    fill = CTCF_status
  )
) +
  
  geom_col() +
  
  labs(
    title = "CTCF association across loop classes",
    x = "Anchor-pair class",
    y = "Percentage of loops"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


ggsave(
  file.path(
    FIGURE_DIR,
    "Figure4_CTCF_by_loop_class.png"
  ),
  p4,
  width = 10,
  height = 6,
  dpi = 300
)


# ============================================================
# 17. REGULATORY BURDEN BY LOOP CLASS
# ============================================================

burden <- loops %>%
  
  group_by(
    Anchor_pair_class
  ) %>%
  
  summarise(
    
    loops = n(),
    
    mean_genes =
      mean(total_genes),
    
    median_genes =
      median(total_genes),
    
    mean_promoters =
      mean(total_promoters),
    
    median_promoters =
      median(total_promoters),
    
    mean_cCRE_promoters =
      mean(total_cCRE_promoters),
    
    mean_enhancers =
      mean(total_enhancers),
    
    median_enhancers =
      median(total_enhancers),
    
    mean_CTCF =
      mean(total_CTCF),
    
    median_CTCF =
      median(total_CTCF),
    
    .groups = "drop"
  )


write.table(
  burden,
  file = file.path(
    TABLE_DIR,
    "regulatory_burden_by_loop_class.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 18. LONG FORMAT FOR REGULATORY BURDEN
# ============================================================

burden_long <- loops %>%
  
  select(
    Anchor_pair_class,
    total_genes,
    total_promoters,
    total_cCRE_promoters,
    total_enhancers,
    total_CTCF
  ) %>%
  
  pivot_longer(
    cols = -Anchor_pair_class,
    names_to = "Feature",
    values_to = "Count"
  )


# ============================================================
# 19. FIGURE 5
# REGULATORY BURDEN
# ============================================================

p5 <- ggplot(
  burden_long,
  aes(
    x = Anchor_pair_class,
    y = Count
  )
) +
  
  geom_boxplot(
    outlier.size = 0.5
  ) +
  
  facet_wrap(
    ~ Feature,
    scales = "free_y"
  ) +
  
  labs(
    title = "Regulatory burden across loop classes",
    x = "Anchor-pair class",
    y = "Number of annotated features"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


ggsave(
  file.path(
    FIGURE_DIR,
    "Figure5_regulatory_burden.png"
  ),
  p5,
  width = 12,
  height = 8,
  dpi = 300
)


# ============================================================
# 20. LOOP DISTANCE
# ============================================================

distance_summary <- loops %>%
  
  summarise(
    
    Minimum_kb =
      min(genomic_distance_kb),
    
    Q1_kb =
      quantile(
        genomic_distance_kb,
        0.25
      ),
    
    Median_kb =
      median(genomic_distance_kb),
    
    Mean_kb =
      mean(genomic_distance_kb),
    
    Q3_kb =
      quantile(
        genomic_distance_kb,
        0.75
      ),
    
    Maximum_kb =
      max(genomic_distance_kb)
  )


write.table(
  distance_summary,
  file = file.path(
    TABLE_DIR,
    "loop_distance_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 21. FIGURE 6
# LOOP DISTANCE DISTRIBUTION
# ============================================================

p6 <- ggplot(
  loops,
  aes(
    x = genomic_distance_kb
  )
) +
  
  geom_histogram(
    bins = 50
  ) +
  
  scale_x_log10() +
  
  labs(
    title = "Genomic distance distribution of loops",
    x = "Loop genomic distance (kb, log10 scale)",
    y = "Number of loops"
  ) +
  
  theme_classic(base_size = 14)


ggsave(
  file.path(
    FIGURE_DIR,
    "Figure6_loop_distance_distribution.png"
  ),
  p6,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# 22. LOOP DISTANCE BY LOOP CLASS
# ============================================================

p7 <- ggplot(
  loops,
  aes(
    x = Anchor_pair_class,
    y = genomic_distance_kb
  )
) +
  
  geom_boxplot(
    outlier.size = 0.4
  ) +
  
  scale_y_log10() +
  
  labs(
    title = "Genomic distance across loop classes",
    x = "Anchor-pair class",
    y = "Loop distance (kb, log10 scale)"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      )
  )


ggsave(
  file.path(
    FIGURE_DIR,
    "Figure7_loop_distance_by_class.png"
  ),
  p7,
  width = 10,
  height = 6,
  dpi = 300
)


# ============================================================
# 23. DISTANCE STATISTICS
# ============================================================

distance_by_class <- loops %>%
  
  group_by(
    Anchor_pair_class
  ) %>%
  
  summarise(
    n = n(),
    median_kb =
      median(genomic_distance_kb),
    mean_kb =
      mean(genomic_distance_kb),
    .groups = "drop"
  )


write.table(
  distance_by_class,
  file = file.path(
    TABLE_DIR,
    "loop_distance_by_class.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 24. KRUSKAL-WALLIS TEST
# ============================================================

distance_kw <- kruskal_test(
  loops,
  genomic_distance_kb ~ Anchor_pair_class
)


write.table(
  distance_kw,
  file = file.path(
    STAT_DIR,
    "loop_distance_kruskal_wallis.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 25. TOP 20 ENHANCER-RICH LOOPS
# ============================================================

top20_enhancer <- loops %>%
  
  arrange(
    desc(total_enhancers)
  ) %>%
  
  slice_head(n = 20)


write.table(
  top20_enhancer,
  file = file.path(
    TABLE_DIR,
    "R_top20_enhancer_loops.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 26. TOP 20 GENE-RICH LOOPS
# ============================================================

top20_gene <- loops %>%
  
  arrange(
    desc(total_genes)
  ) %>%
  
  slice_head(n = 20)


write.table(
  top20_gene,
  file = file.path(
    TABLE_DIR,
    "R_top20_gene_loops.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 27. PROMOTER-ENHANCER LOOPS
# ============================================================

PE_loops <- loops %>%
  
  filter(
    Anchor_pair_class %in% c(
      "PE",
      "P+E-E",
      "P+E-P",
      "P+E-P+E"
    )
  )


write.table(
  PE_loops,
  file = file.path(
    TABLE_DIR,
    "R_promoter_enhancer_loops.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 28. FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("DOWNSTREAM ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat("Loops analysed :", nrow(loops), "\n")
cat("Columns        :", ncol(loops), "\n")

cat("\nOutput directory:\n")
cat(OUTPUT_DIR, "\n")

cat("\nFigures:\n")
print(list.files(FIGURE_DIR))

cat("\nTables:\n")
print(list.files(TABLE_DIR))

cat("\nStatistics:\n")
print(list.files(STAT_DIR))

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")