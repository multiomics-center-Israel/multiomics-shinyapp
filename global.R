library(shiny)
library(shinythemes)
library(ggplot2)
library(plotly)
library(shinycssloaders)
library(bslib)
library(heatmaply)
library(RColorBrewer)
library(readxl)
library(writexl)
library(openxlsx)
library(DT)
library(dplyr)
library(pheatmap)
library(grid)
library(bcrypt)


# ---------- Appstream authentication ----------

# Shared settings: the same for every researcher.
APPSTREAM_SITE_DIR  <- "/srv/shiny-server"
APPSTREAM_USERS_DIR <- file.path(APPSTREAM_SITE_DIR, "Appstream", "Users")
APPSTREAM_LOGIN_URL <- "http://genomics22.bgu.ac.il:8080/Appstream/"

# Authentication is required only on Shiny Server (local RStudio runs are not checked).
# Shiny Server sets SHINY_PORT; the app folder check is a second, independent signal.
APPSTREAM_REQUIRED <- nzchar(Sys.getenv("SHINY_PORT")) ||
                      startsWith(getwd(), paste0(APPSTREAM_SITE_DIR, "/"))

# This app's identity = its folder under the site directory, e.g.
# /srv/shiny-server/Researcher_Name/Project_App -> "Researcher_Name/Project_App"
APPSTREAM_APP_PATH <- if (startsWith(getwd(), paste0(APPSTREAM_SITE_DIR, "/"))) {
  substring(getwd(), nchar(APPSTREAM_SITE_DIR) + 2)
} else NA_character_

message("Appstream authentication: ",
        if (APPSTREAM_REQUIRED) paste0("ON for ", APPSTREAM_APP_PATH, " (getwd: ", getwd(), ")")
        else "OFF (local run)")

# הגבלת כמות הגנים בקלאסטרים של partitionכי הפלוטלי לא מצליח להציג כלכך הרבה גנים  
MAX_HEATMAP_GENES <- 10000



source("functions_shiny.R")

# ---------- tab panel 1 - Home ----------

home_info <- read.delim("data/home_info.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
RESEARCHER_NAME <- home_info[home_info$V1 == "Researcher", "V2"]
PROJECT_NAME    <- home_info[home_info$V1 == "Project", "V2"]
raw_description <- home_info[home_info$V1 == "Description", "V2"]
EXPERIMENT_DESCRIPTION_HTML <- gsub("\\\\n", "<br>", raw_description)

# Parse Bioinformatics_Team field (comma-separated lowercase keys)
raw_team <- home_info[home_info$V1 == "Bioinformatics_Team", "V2"]
BIOINFO_TEAM <- if (length(raw_team) == 0 || is.na(raw_team) || !nzchar(trimws(raw_team))) {
  character(0)
} else {
  trimws(strsplit(raw_team, ",")[[1]])
}

# Registry of all team member UI cards, keyed by lowercase name.
# Each entry is a function() returning a flex-item div (no Bootstrap grid).
.team_card_style <- "flex:1 1 0; min-width:140px; text-align:center; padding:0 10px;"

TEAM_REGISTRY <- list(
  vered = function() shiny::tags$div(style = .team_card_style,
    shiny::tags$img(src = "vered.jpg",
      style = "width:120px; height:120px; object-fit:cover; border-radius:50%; margin-bottom:10px; border: 2px solid #ccc;"),
    shiny::tags$br(),
    shiny::tags$strong("Vered Chalifa-Caspi, Ph.D."), shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-graduation-cap", style = "margin-right: 5px;"), " Ben-Gurion University", shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-envelope",       style = "margin-right: 5px;"), " veredcc@bgu.ac.il",    shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-phone",          style = "margin-right: 5px;"), " 08-6479034"
  ),
  michal = function() shiny::tags$div(style = .team_card_style,
    shiny::tags$img(src = "michal.jpg",
      style = "width:120px; height:120px; object-fit:cover; border-radius:50%; margin-bottom:10px; border: 2px solid #ccc;"),
    shiny::tags$br(),
    shiny::tags$strong("Michal Sharabi-Schwager, Ph.D."), shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-graduation-cap", style = "margin-right: 5px;"), " Ben-Gurion University", shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-envelope",       style = "margin-right: 5px;"), " sharabmi@bgu.ac.il",   shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-phone",          style = "margin-right: 5px;"), " 054-2284280"
  ),
  oz = function() shiny::tags$div(style = .team_card_style,
    shiny::tags$img(src = "oz.jpg",
      style = "width:120px; height:120px; object-fit:cover; border-radius:50%; margin-bottom:10px; border: 2px solid #ccc;"),
    shiny::tags$br(),
    shiny::tags$strong("Oz Solomon, Ph.D."), shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-graduation-cap", style = "margin-right: 5px;"), " Technion University", shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-envelope",       style = "margin-right: 5px;"), " bioinfo.multiomics@technion.ac.il", shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-phone",          style = "margin-right: 5px;"), " 073-378-1979"
  ),
  shelly = function() shiny::tags$div(style = .team_card_style,
    shiny::tags$img(src = "shelly.jpg",
      style = "width:120px; height:120px; object-fit:cover; border-radius:50%; margin-bottom:10px; border: 2px solid #ccc;"),
    shiny::tags$br(),
    shiny::tags$strong("Shelly Haimlich, M.Sc."), shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-graduation-cap", style = "margin-right: 5px;"), " Ben-Gurion University", shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-envelope",       style = "margin-right: 5px;"), " haimlich@bgu.ac.il",   shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-phone",          style = "margin-right: 5px;"), " 052-6373652"
  ),
  or = function() shiny::tags$div(style = .team_card_style,
    shiny::tags$img(src = "or.jpg",
      style = "width:120px; height:120px; object-fit:cover; border-radius:50%; margin-bottom:10px; border: 2px solid #ccc;"),
    shiny::tags$br(),
    shiny::tags$strong("Or Lazarescu, Ph.D."), shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-graduation-cap", style = "margin-right: 5px;"), " Ben-Gurion University", shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-envelope",       style = "margin-right: 5px;"), " lazaresc@post.bgu.ac.il"
  ),
  dor = function() shiny::tags$div(style = .team_card_style,
    shiny::tags$img(src = "dor.png",
      style = "width:120px; height:120px; object-fit:cover; border-radius:50%; margin-bottom:10px; border: 2px solid #ccc;"),
    shiny::tags$br(),
    shiny::tags$strong("Dor Shalev, Ph.D."), shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-graduation-cap", style = "margin-right: 5px;"), " Technion University",       shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-envelope",       style = "margin-right: 5px;"), " dorshalev@technion.ac.il",  shiny::tags$br(),
    shiny::tags$i(class = "fa-solid fa-phone",          style = "margin-right: 5px;"), " 052-7568598"
  )
)

# Build the team UI: select members in the requested order, evenly spaced.
# Uses flexbox so any number of cards fills the row without Bootstrap grid gaps.
# Falls back to all members if Bioinformatics_Team is empty/missing.
build_team_ui <- function() {
  default_order <- c("vered", "michal", "oz", "shelly", "or", "dor")
  keys <- if (length(BIOINFO_TEAM) > 0) {
    BIOINFO_TEAM[BIOINFO_TEAM %in% names(TEAM_REGISTRY)]
  } else {
    default_order
  }
  if (length(keys) == 0) keys <- default_order
  cards <- lapply(keys, function(k) TEAM_REGISTRY[[k]]())
  do.call(
    shiny::tags$div,
    c(list(style = "display:flex; justify-content:center; flex-wrap:wrap; gap:10px;"), cards)
  )
}

# Researcher-specific Home content (researcher, project, description, team).
# Inserted into #home_info_ui by insertUI() inside app_server() (server.R), i.e. only after a
# successful Appstream check, so it is never part of the static HTML (ui.R).
# Defined here because Shiny loads ui.R and server.R in separate environments that
# share only global.R; ui.R must keep `ui` as its last expression.
home_info_ui <- function() {
  tagList(

    # Researcher name banner
    tags$div(
      style = "background-color: #e8e8e8; padding: 25px; border-radius: 10px; margin-bottom: 20px; text-align: center;",
      h1(RESEARCHER_NAME, style = "margin: 0;")
    ),

    # Project name banner
    tags$div(
      style = "background-color: #f5f5f5; padding: 15px; border-radius: 10px; margin-bottom: 20px; text-align: center;",
      h4(PROJECT_NAME, style = "margin: 0;")
    ),

    # Experiment description block
    tags$div(
      style = paste(
        "padding:15px; border-radius:10px;",
        "margin:20px auto; margin-bottom:20px;",
        "text-align:left; max-width:800px;"
      ),
      p(HTML(EXPERIMENT_DESCRIPTION_HTML), style = "margin:0;")
    ),

    # Bioinformatics Team section (members selected from home_info.txt)
    tags$div(
      style = "background-color: #f5f5f5; padding: 15px; border-radius: 10px; margin-bottom: 20px; text-align: center;",
      tags$h3("Bioinformatics Team", style = "text-align: center; margin-bottom: 30px;"),
      build_team_ui()
    )
  )
}

# ---------- tab panel 3 - Data & Files ----------

DOC_DIR <- "data/Data_&_Files"
DOC_MANIFEST  <- "data/Data_&_Files/Data_&_Files_manifest.xlsx"


# ---------- Omics payloads: standardized filenames and tab visibility ----------
# The researcher-specific omics data are these files in data/ (written by the
# upstream pipelines). No other (legacy) filenames are supported.
PAYLOAD_FILES <- c(
  rnaseq       = "data/shiny_payload_rnaseq.rds",
  proteomics   = "data/shiny_payload_proteomics.rds",
  metabolomics = "data/shiny_payload_metabolomics.rds",
  lipidomics   = "data/shiny_payload_lipidomics.rds",   # tab is still a placeholder; not loaded
  multiomics   = "data/shiny_payload_multiomics.rds"    # part of the contract; not loaded yet
)

# Top-level omics tabs (used by ui.R): a tab is shown when its payload file exists.
# File contents are NOT inspected here; missing analyses inside a payload are
# handled by the existing per-tab "not available" logic.
INDIVIDUAL_OMICS <- c("rnaseq", "proteomics", "metabolomics", "lipidomics")
SHOW_TAB <- as.list(setNames(file.exists(PAYLOAD_FILES[INDIVIDUAL_OMICS]), INDIVIDUAL_OMICS))
# MultiOmics: shown when at least two individual omics payloads exist
# (shiny_payload_multiomics.rds is neither required nor counted).
SHOW_TAB$multiomics <- sum(unlist(SHOW_TAB)) >= 2

.shown <- names(SHOW_TAB)[unlist(SHOW_TAB)]
message("Omics tabs shown: ", if (length(.shown)) paste(.shown, collapse = ", ") else "none")
rm(.shown)


# ---------- tab panel 4 - RNAseq ----------

rds_path_rnaseq <- PAYLOAD_FILES[["rnaseq"]]
rnaseq_env <- new.env(parent = emptyenv())

if (file.exists(rds_path_rnaseq)) {
  list2env(readRDS(rds_path_rnaseq), envir = rnaseq_env)
  # Detect sample_meta columns (sample_id / color / shape / group) so the app
  # works regardless of the names the upstream pipeline used.
  apply_sample_meta_aes(rnaseq_env)
}


# ---------- tab panel 5 - Proteomics ----------

rds_path_proteomics <- PAYLOAD_FILES[["proteomics"]]
proteomics_env <- new.env(parent = emptyenv())


if (file.exists(rds_path_proteomics)) {
  list2env(readRDS(rds_path_proteomics), envir = proteomics_env)
  apply_sample_meta_aes(proteomics_env)
}


# ---------- tab panel 6 - Metabolomics ----------

rds_path_metabolomics <- PAYLOAD_FILES[["metabolomics"]]
metabolomics_env <- new.env(parent = emptyenv())

if (file.exists(rds_path_metabolomics)) {
  list2env(readRDS(rds_path_metabolomics), envir = metabolomics_env)
  # Resolve sample-meta columns BEFORE the QC filter below: the filter checks
  # for an expr_long$SampleID column, and the resolver mirrors whichever
  # sample-id column the payload actually uses to "SampleID".
  apply_sample_meta_aes(metabolomics_env)

  # Remove QC samples from metabolomics objects used for display.
  # QC samples are identified by "QC" appearing in their sample name.
  # This affects: Expression Profile (expr_raw, expr_norm), Clustering
  # (de_expr_norm), and Boxplot/Density (expr_long).
  #
  # NOTE: sample_meta is intentionally NOT filtered — PCA uses pca_object$x
  # (precomputed on all 22 samples) joined with sample_meta, so filtering
  # sample_meta would cause a row-count mismatch. The Expression Profile code
  # (get_explorer_df_from_mat) internally subsets sample_meta to the columns
  # present in the expression matrix, so keeping QC rows in sample_meta is
  # safe and does not cause QC to appear in the Expression Profile.
  .drop_qc_cols <- function(mat) {
    if (is.null(mat) || is.null(colnames(mat))) return(mat)
    keep <- !grepl("QC", colnames(mat), ignore.case = TRUE)
    mat[, keep, drop = FALSE]
  }

  for (.nm in c("expr_raw", "expr_norm", "de_expr_norm")) {
    if (exists(.nm, envir = metabolomics_env, inherits = FALSE)) {
      assign(.nm, .drop_qc_cols(get(.nm, envir = metabolomics_env)),
             envir = metabolomics_env)
    }
  }

  if (exists("expr_long", envir = metabolomics_env, inherits = FALSE)) {
    .el <- get("expr_long", envir = metabolomics_env)
    if ("SampleID" %in% colnames(.el)) {
      .keep <- !grepl("QC", .el$SampleID, ignore.case = TRUE)
      assign("expr_long", .el[.keep, , drop = FALSE], envir = metabolomics_env)
    }
  }

  rm(.drop_qc_cols)
}


# ---------- tab panel 8 - MultiOmics ----------

multi_matching_path <- "data/Matching_genes_proteins.xlsx"
multi_matching_tbl  <- NULL

if (file.exists(multi_matching_path)) {
  .raw_tbl <- readxl::read_excel(multi_matching_path)

  # Validate structure: must have 'genes' and 'proteins' columns
  if (all(c("genes", "proteins") %in% colnames(.raw_tbl)) && nrow(.raw_tbl) > 0) {

    # Validate content: at least some gene IDs must exist in RNAseq data,
    # and at least some protein IDs must match Proteomics data.
    .rna_ids <- character(0)
    if (exists("expr_raw",  envir = rnaseq_env, inherits = FALSE)) {
      .rna_ids <- rownames(get("expr_raw",  envir = rnaseq_env, inherits = FALSE))
    } else if (exists("expr_norm", envir = rnaseq_env, inherits = FALSE)) {
      .rna_ids <- rownames(get("expr_norm", envir = rnaseq_env, inherits = FALSE))
    }

    .prot_ids <- character(0)
    if (exists("expr_raw",  envir = proteomics_env, inherits = FALSE)) {
      .prot_ids <- rownames(get("expr_raw",  envir = proteomics_env, inherits = FALSE))
    } else if (exists("expr_norm", envir = proteomics_env, inherits = FALSE)) {
      .prot_ids <- rownames(get("expr_norm", envir = proteomics_env, inherits = FALSE))
    }

    .gene_overlap <- any(.raw_tbl$genes %in% .rna_ids)
    .prot_stripped_tbl  <- sub("\\.\\d+$", "", as.character(.raw_tbl$proteins))
    .prot_stripped_data <- sub("\\.\\d+$", "", .prot_ids)
    .prot_overlap <- any(.prot_stripped_tbl %in% .prot_stripped_data)

    if (.gene_overlap && .prot_overlap) {
      multi_matching_tbl <- .raw_tbl
    } else {
      warning("Matching_genes_proteins.xlsx loaded but IDs do not overlap with ",
              "the current RNAseq/Proteomics data; MultiOmics Explorer will not be available.")
    }
  } else {
    warning("Matching_genes_proteins.xlsx is missing required columns ('genes', 'proteins') ",
            "or is empty; MultiOmics Explorer will not be available.")
  }
  rm(.raw_tbl)
} else {
  warning("Matching_genes_proteins.xlsx not found; MultiOmics Explorer will not be available.")
}



