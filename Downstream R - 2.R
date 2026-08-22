# ============================================================
# CHR1 3D GENOME DOWNSTREAM ANALYSIS
# ============================================================
#
# PURPOSE:
#   1. TAD / domain statistics
#   2. TAD size distribution
#   3. A/B compartment analysis
#   4. Compartment domain analysis
#   5. Compartment strength analysis
#   6. Pearson correlation analysis
#   7. Publication-quality figures
#
# INPUT FILES:
#
#   tad_domains_25kb.bed
#   AB_compartments_100kb.bed
#   AB_domains_100kb.bed
#   compartment_strength_100kb.txt
#   chr1_pearsons_100kb.txt
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

required_packages <- c(
  "tidyverse"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  
}

library(tidyverse)


# ============================================================
# 2. INPUT / OUTPUT DIRECTORIES
# ============================================================

HIC_DIR <- "D:/hic/hic_project/aligned"
HIC_AB <- "D:/hic/hic_project/downstream/fanc"
              
TAD_FILE <- file.path(
  HIC_DIR,
  "tad_domains_25kb.bed"
)

COMPARTMENT_FILE <- file.path(
  HIC_AB,
  "AB_compartments_100kb.bed"
)

COMPARTMENT_DOMAIN_FILE <- file.path(
  HIC_AB,
  "AB_domains_100kb.bed"
)

STRENGTH_FILE <- file.path(
  HIC_AB,
  "compartment_strength_100kb.txt"
)

PEARSON_FILE <- file.path(
  HIC_DIR,
  "chr1_pearsons_100kb.txt"
)


OUTPUT_DIR <- file.path(
  HIC_DIR,
  "R_3D_genome_analysis"
)

FIGURE_DIR <- file.path(
  OUTPUT_DIR,
  "figures"
)

TABLE_DIR <- file.path(
  OUTPUT_DIR,
  "tables"
)

STAT_DIR <- file.path(
  OUTPUT_DIR,
  "statistics"
)


dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGURE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  STAT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. CHECK INPUT FILES
# ============================================================

input_files <- c(
  TAD_FILE,
  COMPARTMENT_FILE,
  COMPARTMENT_DOMAIN_FILE,
  STRENGTH_FILE,
  PEARSON_FILE
)

cat("\n")
cat("============================================================\n")
cat("INPUT FILE CHECK\n")
cat("============================================================\n")

input_files <- c(
  TAD_FILE,
  COMPARTMENT_FILE,
  COMPARTMENT_DOMAIN_FILE,
  STRENGTH_FILE,
  PEARSON_FILE
)

for (f in input_files) {
  
  if (file.exists(f)) {
    
    cat("FOUND   :", basename(f), "\n")
    
  } else {
    
    cat("MISSING :", basename(f), "\n")
    
  }
  
}

# Stop if important files are missing

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  
  cat("\n")
  cat("WARNING: Some input files are missing.\n")
  cat("The corresponding analysis sections will be skipped.\n")
  
}

# ============================================================
# 4. TAD / DOMAIN ANALYSIS
# ============================================================

if (file.exists(TAD_FILE)) {
  
  cat("\n")
  cat("============================================================\n")
  cat("TAD ANALYSIS\n")
  cat("============================================================\n")
  
  
  # ----------------------------------------------------------
  # Read TAD BED file
  # ----------------------------------------------------------
  
  tads <- read.delim(
    TAD_FILE,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = ""
  )
  
  
  # ----------------------------------------------------------
  # Assign column names
  # ----------------------------------------------------------
  colnames(tads)[1:6] <- c(
    "chr",
    "start",
    "end",
    "TAD_ID",
    "score",
    "strand"
  )
  
  
  # ----------------------------------------------------------
  # Remove comment/header rows
  # ----------------------------------------------------------
  
  tads <- tads %>%
    filter(
      chr != "#",
      chr != "",
      !is.na(chr)
    )
  
  
  # ----------------------------------------------------------
  # Convert coordinates
  # ----------------------------------------------------------
  
  tads <- tads %>%
    mutate(
      start = as.numeric(as.character(start)),
      end   = as.numeric(as.character(end))
    )
  
  
  # ----------------------------------------------------------
  # Count invalid rows BEFORE removing them
  # ----------------------------------------------------------
  
  invalid_tad_rows <- sum(
    is.na(tads$start) |
      is.na(tads$end)
  )
  
  
  # ----------------------------------------------------------
  # Remove invalid coordinates
  # ----------------------------------------------------------
  
  tads <- tads %>%
    filter(
      !is.na(start),
      !is.na(end)
    )
  
  
  # ----------------------------------------------------------
  # Calculate TAD size
  # ----------------------------------------------------------
  tads <- tads %>%
    mutate(
      TAD_size_bp = end - start,
      TAD_size_kb = TAD_size_bp / 1000,
      TAD_size_Mb = TAD_size_bp / 1000000
    )
  
  
  # ----------------------------------------------------------
  # TAD QC
  # ----------------------------------------------------------
  
  cat(
    "Number of valid TADs:",
    nrow(tads),
    "\n"
  )
  
  cat(
    "Invalid TAD rows removed:",
    invalid_tad_rows,
    "\n"
  )
  
  cat(
    "Chromosomes present:\n"
  )
  
  print(
    table(tads$chr)
  )
  
  # ----------------------------------------------------------
  # TAD SUMMARY
  # ----------------------------------------------------------
  
  tad_summary <- tads %>%
    summarise(
      
      Number_of_TADs =
        n(),
      
      Minimum_TAD_size_kb =
        min(TAD_size_kb),
      
      Q1_TAD_size_kb =
        quantile(
          TAD_size_kb,
          0.25
        ),
      
      Median_TAD_size_kb =
        median(TAD_size_kb),
      
      Mean_TAD_size_kb =
        mean(TAD_size_kb),
      
      Q3_TAD_size_kb =
        quantile(
          TAD_size_kb,
          0.75
        ),
      
      Maximum_TAD_size_kb =
        max(TAD_size_kb)
    )
  
  
  print(tad_summary)
  
  
  write.table(
    tad_summary,
    file = file.path(
      TABLE_DIR,
      "TAD_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # ==========================================================
  # FIGURE 1 — TAD SIZE DISTRIBUTION
  # ==========================================================
  # ==========================================================
  # FIGURE 1 — TAD SIZE DISTRIBUTION
  # ==========================================================
  
  p1 <- ggplot(
    tads,
    aes(
      x = TAD_size_kb
    )
  ) +
    
    geom_histogram(
      bins = 50
    ) +
    
    labs(
      title = "TAD size distribution",
      x = "TAD size (kb)",
      y = "Number of TADs"
    ) +
    
    theme_classic(
      base_size = 14
    )
  
  
  ggsave(
    file.path(
      FIGURE_DIR,
      "Figure1_TAD_size_distribution.png"
    ),
    p1,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  ggsave(
    file.path(
      FIGURE_DIR,
      "Figure1_TAD_size_distribution.pdf"
    ),
    p1,
    width = 8,
    height = 6
  )
  
  
  # ==========================================================
  # FIGURE 2 — TAD SIZE BOXPLOT
  # ==========================================================
  
  p2 <- ggplot(
    tads,
    aes(
      y = TAD_size_kb
    )
  ) +
    
    
    geom_boxplot(
      width = 0.3
    ) +
    
    labs(
      title = "Distribution of TAD sizes",
      y = "TAD size (kb)",
      x = NULL
    ) +
    
    theme_classic(
      base_size = 14
    )
  
  ggsave(
    file.path(
      FIGURE_DIR,
      "Figure2_TAD_size_boxplot.png"
    ),
    p2,
    width = 6,
    height = 6,
    dpi = 300
  )
  # ==========================================================
  # TAD SIZE CATEGORIES
  # ==========================================================
  
  tads <- tads %>%
    mutate(
      
      TAD_size_category = case_when(
        
        TAD_size_kb < 100
        ~ "<100 kb",
        
        TAD_size_kb >= 100 &
          TAD_size_kb < 500
        ~ "100–500 kb",
        
        TAD_size_kb >= 500 &
          TAD_size_kb < 1000
        ~ "500 kb–1 Mb",
        
        TAD_size_kb >= 1000
        ~ ">1 Mb"
      )
    )
  
  
  tad_size_categories <- tads %>%
    count(
      TAD_size_category
    ) %>%
    mutate(
      Percentage =
        100 * n / sum(n)
    )
  
  
  write.table(
    tad_size_categories,
    file = file.path(
      TABLE_DIR,
      "TAD_size_categories.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
}

# ============================================================
# 5. A/B COMPARTMENT ANALYSIS
# ============================================================

if (file.exists(COMPARTMENT_FILE)) {
  
  cat("\n")
  cat("============================================================\n")
  cat("A/B COMPARTMENT ANALYSIS\n")
  cat("============================================================\n")
  
  
  # ----------------------------------------------------------
  # Read compartment BED file
  # ----------------------------------------------------------
  
  comp <- read.delim(
    COMPARTMENT_FILE,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = ""
  )
  
  
  # ----------------------------------------------------------
  # Assign columns
  #
  # Actual file structure:
  #
  # chr start end A/B eigenvector_value .
  #
  # Therefore:
  # column 4 = A/B compartment label
  # column 5 = numeric eigenvector value
  # column 6 = strand
  # ----------------------------------------------------------
  
  colnames(comp)[1:6] <- c(
    "chr",
    "start",
    "end",
    "compartment",
    "eigenvector_value",
    "strand"
  )
  
  
  # ----------------------------------------------------------
  # Remove invalid rows
  # ----------------------------------------------------------
  
  comp <- comp %>%
    filter(
      chr != "#",
      chr != "",
      !is.na(chr)
    )
  
  
  # ----------------------------------------------------------
  # Convert numeric columns
  # ----------------------------------------------------------
  
  comp <- comp %>%
    mutate(
      start = as.numeric(as.character(start)),
      end = as.numeric(as.character(end)),
      eigenvector_value =
        as.numeric(as.character(eigenvector_value))
    )
  
  
  # ----------------------------------------------------------
  # Remove invalid numeric rows
  # ----------------------------------------------------------
  
  comp <- comp %>%
    filter(
      !is.na(start),
      !is.na(end),
      !is.na(eigenvector_value)
    )
  
  
  # ----------------------------------------------------------
  # QC
  # ----------------------------------------------------------
  
  cat(
    "Number of compartment bins:",
    nrow(comp),
    "\n"
  )
  
  
  cat(
    "Chromosomes present:\n"
  )
  
  print(
    table(comp$chr)
  )
  
  
  cat(
    "Compartment classes:\n"
  )
  
  print(
    table(comp$compartment)
  )
  
  
  # ==========================================================
  # EIGENVECTOR SUMMARY
  # ==========================================================
  
  eigenvector_summary <- comp %>%
    summarise(
      
      Number_of_bins = n(),
      
      Minimum = min(
        eigenvector_value,
        na.rm = TRUE
      ),
      
      Q1 = quantile(
        eigenvector_value,
        0.25,
        na.rm = TRUE
      ),
      
      Median = median(
        eigenvector_value,
        na.rm = TRUE
      ),
      
      Mean = mean(
        eigenvector_value,
        na.rm = TRUE
      ),
      
      Q3 = quantile(
        eigenvector_value,
        0.75,
        na.rm = TRUE
      ),
      
      Maximum = max(
        eigenvector_value,
        na.rm = TRUE
      )
    )
  
  
  print(eigenvector_summary)
  
  
  write.table(
    eigenvector_summary,
    file = file.path(
      TABLE_DIR,
      "eigenvector_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # ==========================================================
  # A/B COMPARTMENT DISTRIBUTION
  # ==========================================================
  
  compartment_summary <- comp %>%
    count(
      compartment,
      name = "Bins"
    ) %>%
    mutate(
      Percentage = 100 * Bins / sum(Bins)
    )
  
  
  print(compartment_summary)
  
  
  write.table(
    compartment_summary,
    file = file.path(
      TABLE_DIR,
      "AB_compartment_distribution.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # ==========================================================
  # FIGURE 3 — A/B COMPARTMENT DISTRIBUTION
  # ==========================================================
  
  p3 <- ggplot(
    compartment_summary,
    aes(
      x = compartment,
      y = Bins
    )
  ) +
    
    geom_col() +
    
    labs(
      title = "A/B compartment distribution",
      x = "Compartment",
      y = "Number of bins"
    ) +
    
    theme_classic(
      base_size = 14
    )
  
  
  print(p3)
  
  
  ggsave(
    file.path(
      FIGURE_DIR,
      "Figure3_AB_compartment_distribution.png"
    ),
    p3,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  
  ggsave(
    file.path(
      FIGURE_DIR,
      "Figure3_AB_compartment_distribution.pdf"
    ),
    p3,
    width = 7,
    height = 5
  )
  
  
  # ==========================================================
  # FIGURE 4 — EIGENVECTOR DISTRIBUTION
  # ==========================================================
  
  p4 <- ggplot(
    comp,
    aes(
      x = eigenvector_value
    )
  ) +
    
    geom_histogram(
      bins = 50
    ) +
    
    labs(
      title = "Eigenvector distribution",
      x = "Eigenvector value",
      y = "Number of bins"
    ) +
    
    theme_classic(
      base_size = 14
    )
  
  
  print(p4)
  
  
  ggsave(
    file.path(
      FIGURE_DIR,
      "Figure4_eigenvector_distribution.png"
    ),
    p4,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  ggsave(
    file.path(
      FIGURE_DIR,
      "Figure4_eigenvector_distribution.pdf"
    ),
    p4,
    width = 8,
    height = 6
  )
  
  
  # ==========================================================
  # FIGURE 5 — EIGENVECTOR BY COMPARTMENT
  # ==========================================================
  
  p5 <- ggplot(
    comp,
    aes(
      x = compartment,
      y = eigenvector_value
    )
  ) +
    
    geom_boxplot(
      outlier.size = 0.4
    ) +
    
    labs(
      title = "Eigenvector values across A/B compartments",
      x = "Compartment",
      y = "Eigenvector value"
    ) +
    
    theme_classic(
      base_size = 14
    )
  
  
  print(p5)
  
  
  ggsave(
    file.path(
      FIGURE_DIR,
      "Figure5_eigenvector_by_compartment.png"
    ),
    p5,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  
  ggsave(
    file.path(
      FIGURE_DIR,
      "Figure5_eigenvector_by_compartment.pdf"
    ),
    p5,
    width = 7,
    height = 6
  )
  
  
  # ==========================================================
  # COMPLETION MESSAGE
  # ==========================================================
  
  cat("\n")
  cat("A/B compartment analysis completed successfully.\n")
  cat("Figures 3, 4 and 5 generated.\n")
  cat("Tables generated:\n")
  cat("  - eigenvector_summary.tsv\n")
  cat("  - AB_compartment_distribution.tsv\n")
  
  
} else {
  
  cat("\n")
  cat(
    "WARNING: Compartment file not found:\n",
    COMPARTMENT_FILE,
    "\n"
  )
  
}

# ============================================================
# 6. A/B COMPARTMENT DOMAIN ANALYSIS
# ============================================================

if (file.exists(COMPARTMENT_DOMAIN_FILE)) {
  
  cat("\n")
  cat("============================================================\n")
  cat("A/B COMPARTMENT DOMAIN ANALYSIS\n")
  cat("============================================================\n")
  
  
  # ----------------------------------------------------------
  # Read domain file
  # ----------------------------------------------------------
  
  ab_domains <- read.delim(
    COMPARTMENT_DOMAIN_FILE,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    comment.char = ""
  )
  
  
  # ----------------------------------------------------------
  # Inspect number of columns
  # ----------------------------------------------------------
  
  cat(
    "Number of columns in AB domain file:",
    ncol(ab_domains),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Assign standard BED names
  # ----------------------------------------------------------
  
  if (ncol(ab_domains) >= 4) {
    
    colnames(ab_domains)[1:4] <- c(
      "chr",
      "start",
      "end",
      "compartment"
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Remove comment rows
  # ----------------------------------------------------------
  
  ab_domains <- ab_domains %>%
    filter(
      chr != "#",
      chr != "",
      !is.na(chr)
    )
  
  
  # ----------------------------------------------------------
  # Convert coordinates
  # ----------------------------------------------------------
  
  ab_domains <- ab_domains %>%
    mutate(
      start = as.numeric(as.character(start)),
      end = as.numeric(as.character(end))
    ) %>%
    filter(
      !is.na(start),
      !is.na(end)
    ) %>%
    mutate(
      domain_size_bp =
        end - start,
      
      domain_size_kb =
        domain_size_bp / 1000,
      
      domain_size_Mb =
        domain_size_bp / 1000000
    )
  
  
  cat(
    "Number of A/B domains:",
    nrow(ab_domains),
    "\n"
  )
  
  
  # ==========================================================
  # DOMAIN SUMMARY
  # ==========================================================
  
  ab_domain_summary <- ab_domains %>%
    group_by(
      compartment
    ) %>%
    summarise(
      
      Number_of_domains =
        n(),
      
      Minimum_size_kb =
        min(
          domain_size_kb,
          na.rm = TRUE
        ),
      
      Median_size_kb =
        median(
          domain_size_kb,
          na.rm = TRUE
        ),
      
      Mean_size_kb =
        mean(
          domain_size_kb,
          na.rm = TRUE
        ),
      
      Maximum_size_kb =
        max(
          domain_size_kb,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
  
  
  print(ab_domain_summary)
  
  
  write.table(
    ab_domain_summary,
    file = file.path(
      TABLE_DIR,
      "AB_domain_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # ==========================================================
  # FIGURE 6 — A/B DOMAIN SIZE
  # ==========================================================
  
  p6 <- ggplot(
    ab_domains,
    aes(
      x = compartment,
      y = domain_size_kb
    )
  ) +
    
    geom_boxplot(
      outlier.size = 0.4
    ) +
    
    labs(
      title = "A/B compartment domain sizes",
      x = "Compartment",
      y = "Domain size (kb)"
    ) +
    
    theme_classic(
      base_size = 14
    )
  
  
  ggsave(
    file.path(
      FIGURE_DIR,
      "Figure6_AB_domain_size.png"
    ),
    p6,
    width = 7,
    height = 6,
    dpi = 300
  )
  
}


# ============================================================
# 7. COMPARTMENT STRENGTH
# ============================================================

if (file.exists(STRENGTH_FILE)) {
  
  cat("\n")
  cat("============================================================\n")
  cat("COMPARTMENT STRENGTH\n")
  cat("============================================================\n")
  
  
  # ----------------------------------------------------------
  # Read file
  # ----------------------------------------------------------
  
  strength <- read.delim(
    STRENGTH_FILE,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
  
  
  cat(
    "Columns in strength file:\n"
  )
  
  print(
    colnames(strength)
  )
  
  
  # ----------------------------------------------------------
  # Convert numeric columns where possible
  # ----------------------------------------------------------
  
  strength <- strength %>%
    mutate(
      across(
        where(is.character),
        ~ suppressWarnings(
          as.numeric(.x)
        )
      )
    )
  
  
  # ----------------------------------------------------------
  # Summary of numeric columns
  # ----------------------------------------------------------
  
  numeric_strength <- strength %>%
    select(
      where(is.numeric)
    )
  
  
  if (ncol(numeric_strength) > 0) {
    
    strength_summary <- numeric_strength %>%
      summarise(
        across(
          everything(),
          list(
            mean = ~ mean(
              .x,
              na.rm = TRUE
            ),
            median = ~ median(
              .x,
              na.rm = TRUE
            )
          )
        )
      )
    
    
    print(strength_summary)
    
    
    write.table(
      strength_summary,
      file = file.path(
        TABLE_DIR,
        "compartment_strength_summary.tsv"
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
  }
  
}


# ============================================================
# PEARSON CORRELATION ANALYSIS
# ============================================================
    
    cat("\n")
    cat("============================================================\n")
    cat("PEARSON CORRELATION ANALYSIS\n")
    cat("============================================================\n")
    
    
    # ------------------------------------------------------------
    # CHECK PEARSON FILE
    # ------------------------------------------------------------
    
    if (!file.exists(PEARSON_FILE)) {
      
      stop(
        paste(
          "Pearson file not found:",
          PEARSON_FILE
        )
      )
      
    }
    
    cat(
      "Pearson file:",
      PEARSON_FILE,
      "\n"
    )
    
    cat(
      "File size:",
      file.info(PEARSON_FILE)$size,
      "bytes\n"
    )
    
    
    # ------------------------------------------------------------
    # READ PEARSON MATRIX
    # ------------------------------------------------------------
    #
    # IMPORTANT:
    # The FAN-C Pearson file is whitespace-separated,
    # NOT tab-separated.
    #
    # Therefore:
    # sep = ""
    #
    # This automatically handles spaces/tabs.
    
    cat("\nReading Pearson matrix...\n")
    
    pearson_matrix <- as.matrix(
      read.table(
        PEARSON_FILE,
        header = FALSE,
        sep = "",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        comment.char = "",
        fill = TRUE
      )
    )
    
    
    # ------------------------------------------------------------
    # CHECK MATRIX DIMENSIONS
    # ------------------------------------------------------------
    
    cat(
      "Pearson matrix dimensions:",
      nrow(pearson_matrix),
      "x",
      ncol(pearson_matrix),
      "\n"
    )
    
    
    if (
      nrow(pearson_matrix) !=
      ncol(pearson_matrix)
    ) {
      
      stop(
        paste(
          "Pearson file is not a square matrix.",
          "Rows:",
          nrow(pearson_matrix),
          "Columns:",
          ncol(pearson_matrix)
        )
      )
      
    }
    
    
    # ------------------------------------------------------------
    # CONVERT TO NUMERIC
    # ------------------------------------------------------------
    
    pearson_matrix <- matrix(
      suppressWarnings(
        as.numeric(pearson_matrix)
      ),
      nrow = nrow(pearson_matrix),
      ncol = ncol(pearson_matrix)
    )
    
    
    # ------------------------------------------------------------
    # CHECK MATRIX
    # ------------------------------------------------------------
    
    cat(
      "Total matrix elements:",
      length(pearson_matrix),
      "\n"
    )
    
    cat(
      "Finite values before filtering:",
      sum(is.finite(pearson_matrix)),
      "\n"
    )
    
    cat(
      "NaN / Inf values:",
      sum(!is.finite(pearson_matrix)),
      "\n"
    )
    
    
    # ------------------------------------------------------------
    # REMOVE DIAGONAL
    # ------------------------------------------------------------
    #
    # Diagonal values are self-correlations:
    #
    #   bin 1 vs bin 1 = 1
    #   bin 2 vs bin 2 = 1
    #   etc.
    #
    # These should not be included in the distribution.
    
    diag(
      pearson_matrix
    ) <- NA
    
    
    # ------------------------------------------------------------
    # EXTRACT PEARSON VALUES
    # ------------------------------------------------------------
    
    pearson_values <- as.numeric(
      pearson_matrix
    )
    
    
    # Keep only finite correlations
    
    pearson_values <- pearson_values[
      is.finite(pearson_values)
    ]
    
    
    cat(
      "Valid Pearson values after filtering:",
      length(pearson_values),
      "\n"
    )
    
    
    # ------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------
    
    if (
      length(pearson_values) > 0
    ) {
      
      cat("\nPearson summary:\n")
      
      print(
        summary(
          pearson_values
        )
      )
      
      cat(
        "\nMinimum:",
        min(pearson_values),
        "\n"
      )
      
      cat(
        "Maximum:",
        max(pearson_values),
        "\n"
      )
      
      cat(
        "Mean:",
        mean(pearson_values),
        "\n"
      )
      
      cat(
        "Median:",
        median(pearson_values),
        "\n"
      )
      
    }
    
    
    # ------------------------------------------------------------
    # SAVE PEARSON SUMMARY TABLE
    # ------------------------------------------------------------
    
    if (
      length(pearson_values) > 0
    ) {
      
      pearson_summary <- data.frame(
        
        Number_of_correlations =
          length(pearson_values),
        
        Minimum =
          min(pearson_values),
        
        Q1 =
          quantile(
            pearson_values,
            0.25
          ),
        
        Median =
          median(
            pearson_values
          ),
        
        Mean =
          mean(
            pearson_values
          ),
        
        Q3 =
          quantile(
            pearson_values,
            0.75
          ),
        
        Maximum =
          max(
            pearson_values
          )
        
      )
      
      
      write.table(
        pearson_summary,
        file = file.path(
          TABLE_DIR,
          "Pearson_summary.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
      
    }
    
    
    # ============================================================
    # FIGURE 7 — PEARSON DISTRIBUTION
    # ============================================================
    
    if (
      length(pearson_values) > 0
    ) {
      
      pearson_df <- data.frame(
        Pearson = pearson_values
      )
      
      
      p7 <- ggplot(
        pearson_df,
        aes(
          x = Pearson
        )
      ) +
        
        geom_histogram(
          bins = 50,
          na.rm = TRUE
        ) +
        
        labs(
          title = "Pearson correlation distribution",
          x = "Pearson correlation",
          y = "Number of bin pairs"
        ) +
        
        theme_classic(
          base_size = 14
        )
      
      
      # Display plot
      print(p7)
      
      
      # Save PNG
      
      ggsave(
        filename = file.path(
          FIGURE_DIR,
          "Figure7_Pearson_distribution.png"
        ),
        plot = p7,
        width = 8,
        height = 6,
        dpi = 300
      )
      
      
      # Save PDF
      
      ggsave(
        filename = file.path(
          FIGURE_DIR,
          "Figure7_Pearson_distribution.pdf"
        ),
        plot = p7,
        width = 8,
        height = 6
      )
      
      
      cat(
        "\nPearson distribution plots generated successfully.\n"
      )
      
      
    } else {
      
      warning(
        "No valid Pearson values were found."
      )
      
    }
    
    
    # ============================================================
    # FINAL PEARSON CHECK
    # ============================================================
    
    cat("\n")
    cat("============================================================\n")
    cat("PEARSON ANALYSIS COMPLETE\n")
    cat("============================================================\n")
    
    cat(
      "Matrix size:",
      nrow(pearson_matrix),
      "x",
      ncol(pearson_matrix),
      "\n"
    )
    
    cat(
      "Valid Pearson correlations:",
      length(pearson_values),
      "\n"
    )
    
    cat(
      "Pearson summary file:",
      file.path(
        TABLE_DIR,
        "Pearson_summary.tsv"
      ),
      "\n"
    )
    
    cat(
      "Pearson figure:",
      file.path(
        FIGURE_DIR,
        "Figure7_Pearson_distribution.png"
      ),
      "\n"
    )
    
    cat("============================================================\n")
    cat("DONE\n")
    cat("============================================================\n")
    
    # ------------------------------------------------------------
    # 8. FINAL CHECK
    # ------------------------------------------------------------
    
    cat("\n")
    cat("============================================================\n")
    cat("PEARSON ANALYSIS COMPLETE\n")
    cat("============================================================\n")
    
    cat(
      "Valid Pearson values:",
      length(
        pearson_values
      ),
      "\n"
    )
    
    cat(
      "Pearson summary table:",
      file.path(
        TABLE_DIR,
        "Pearson_summary.tsv"
      ),
      "\n"
    )
    
    cat(
      "Pearson PNG:",
      file.path(
        FIGURE_DIR,
        "Figure7_Pearson_distribution.png"
      ),
      "\n"
    )
    
    cat(
      "Pearson PDF:",
      file.path(
        FIGURE_DIR,
        "Figure7_Pearson_distribution.pdf"
      ),
      "\n"
    )
# ============================================================
# 9. FINAL OUTPUT SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("3D GENOME DOWNSTREAM ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat(
  "Output directory:\n",
  OUTPUT_DIR,
  "\n"
)

cat("\n")

cat("Figures generated:\n")

print(
  list.files(
    FIGURE_DIR
  )
)

cat("\n")

cat("Tables generated:\n")

print(
  list.files(
    TABLE_DIR
  )
)

cat("\n")

cat("Statistics generated:\n")

print(
  list.files(
    STAT_DIR
  )
)

cat("\n")
cat("============================================================\n")
cat("DONE\n")
cat("============================================================\n")