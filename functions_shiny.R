
# ---------- tab panel 3 - Data & Files ----------

# Read the documentation manifest Excel file and return a cleaned data.frame for UI rendering
read_Data_and_Files_manifest <- function(manifest_path) {
  
  # If manifest is missing -> return empty df (soft fail)
  if (!file.exists(manifest_path)) {
    warning("Data_&_Files_manifest manifest not found: ", manifest_path)
    return(data.frame(
      title       = character(),
      file_path   = character(),
      description = character(),
      category    = character(),
      updated_at  = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Safe read (Excel)
  df <- tryCatch(
    readxl::read_excel(manifest_path),
    error = function(e) {
      warning("Failed to read Data & Files manifest: ", e$message)
      return(NULL)
    }
  )
  
  # On read error -> empty df
  if (is.null(df)) {
    return(data.frame(
      title       = character(),
      file_path   = character(),
      description = character(),
      category    = character(),
      updated_at  = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Normalize column names
  names(df) <- tolower(names(df))
  
  # Ensure required columns exist (soft fail: add missing as NA)
  required <- c("title", "file_path", "description", "category", "updated_at")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    warning("Data & Files manifest is missing columns: ", paste(missing, collapse = ", "))
    for (m in missing) df[[m]] <- NA
  }
  
  # Clean text fields
  df$title       <- trimws(as.character(df$title))
  df$file_path   <- trimws(as.character(df$file_path))
  df$description <- trimws(as.character(df$description))
  df$category    <- trimws(as.character(df$category))
  
  # description: NA -> "" (so it won't be printed)
  df$description[is.na(df$description)] <- ""
  
  # category: empty -> NA (so it will show without category header)
  df$category[df$category == "" | is.na(df$category)] <- NA
  
  # updated_at: accept Excel Date OR text DD/MM/YYYY; otherwise NA
  raw_date <- df$updated_at
  d <- suppressWarnings(as.Date(raw_date))  # Excel Date often works here
  if (all(is.na(d)) && any(!is.na(raw_date))) {
    d <- suppressWarnings(as.Date(as.character(raw_date), format = "%d/%m/%Y"))
  }
  df$updated_at <- ifelse(is.na(d), NA, format(d, "%d/%m/%Y"))
  
  # title empty -> drop row (nothing to display)
  df <- df[!is.na(df$title) & df$title != "", , drop = FALSE]
  
  # file_path empty -> NA (will be marked missing)
  df$file_path[df$file_path == "" | is.na(df$file_path)] <- NA
  
  # Build full file path from DOC_DIR (global)
  df$.full_path <- ifelse(is.na(df$file_path), NA, file.path(DOC_DIR, df$file_path))
  
  # File exists flag
  df$.exists <- ifelse(is.na(df$.full_path), FALSE, file.exists(df$.full_path))
  
  # Internal id (stable within a read)
  df$.id <- paste0("doc_", seq_len(nrow(df)))
  
  df
}

# Build a single documentation row UI (title/link + optional description + optional date)
doc_file_row_ui <- function(item, output_id_download) {
  
  # Format date string (optional)
  date_txt <- ifelse(
    is.na(item$updated_at) || item$updated_at == "",
    "",
    item$updated_at
  )
  
  # Title as downloadLink if file exists; otherwise grey text
  file_link <- if (isTRUE(item$.exists)) {
    downloadLink(
      output_id_download,
      label = item$title,
      style = "font-weight:600; font-size:16px; text-decoration: underline;"
    )
  } else {
    tags$span(item$title, style = "font-weight:600; font-size:16px; color:#777;")
  }
  
  tags$div(
    style = "margin-bottom: 22px;",
    
    # Optional description
    if (!is.na(item$description) && nzchar(item$description)) {
      tags$div(item$description, style = "color:#444; margin-bottom: 4px;")
    },
    
    # Title + date + missing flag
    tags$div(
      style = "display:flex; align-items:baseline; gap:10px; flex-wrap:wrap;",
      
      file_link,
      
      if (nzchar(date_txt)) tags$span(date_txt, style = "color:#777; font-size: 13px;") else NULL,
      
      if (!isTRUE(item$.exists)) {
        tags$span(
          "Missing file",
          style = "color:#b00020; font-size: 13px; font-weight:600; margin-left:10px;"
        )
      } else NULL
    )
  )
}


# ---------- tab panel 4-5 - RNAseq & Proteomics----------


## ---- Subtab: Exploratory Analyses ----

## ---- Subtab: Data distributions ----

# Get required object from an environment (throws a clean error if missing)
env_get <- function(env, name, label = "Data", inherits = FALSE) {
  if (!exists(name, envir = env, inherits = inherits)) {
    stop(sprintf("%s not available (%s missing).", label, name), call. = FALSE)
  }
  get(name, envir = env, inherits = inherits)
}

# Get optional object from an environment (returns default if missing)
env_get_optional <- function(env, name, default = NULL, inherits = FALSE) {
  if (!exists(name, envir = env, inherits = inherits)) return(default)
  get(name, envir = env, inherits = inherits)
}

# Normalize a config value (trim, take first element, return NULL for "none"/empty)
normalize_opt_var <- function(x) {
  v <- trimws(as.character(x)[1])
  if (!nzchar(v) || tolower(v) %in% c("none", "null", "na", "n/a")) return(NULL)
  v
}

# ============================================================================
# SHARED METADATA-COLUMN RESOLVER (works across all omics)
# ----------------------------------------------------------------------------
# Different payloads use different sample_meta column names (Condition, Stress,
# Treatment, group, time, ...). resolve_sample_meta_aes() inspects the env's
# sample_meta and decides which column should be used for sample-id, color,
# optional shape, and the patterns-plot grouping. It honors a pre-set
# env$color / env$shape / env$group / env$sample_id only when those values
# point to a real column; otherwise it infers from the data.
#
# apply_sample_meta_aes() writes the resolved values back to the env so the
# rest of the app (which already reads env$color / env$shape / env$group)
# keeps working unchanged across payloads. It also mirrors the resolved
# sample-id column to "SampleID" inside sample_meta and expr_long, since a
# few plot helpers historically assume that name.
# ============================================================================

# Sample-id column candidates, in priority order.
.AES_ID_CANDIDATES <- c(
  "SampleID", "sample_id", "Sample_ID", "Sample.ID",
  "Sample", "sample", "sampleID",
  "sample_name", "SampleName",
  "origin_id"
)

# Common biological grouping columns, in priority order.
.AES_GROUP_CANDIDATES <- c(
  "Condition", "condition",
  "Stress", "stress",
  "Treatment", "treatment",
  "Group", "group",
  "Time", "time", "timepoint", "Timepoint", "day", "Day", "hour", "Hour",
  "Genotype", "genotype",
  "Phenotype", "phenotype",
  "celltype", "cell_type", "CellType",
  "tissue", "Tissue", "tissue_type",
  "Batch", "batch"
)

# Column names that look like technical identifiers and should be skipped
# when picking a grouping column.
.AES_ID_LIKE_REGEX <-
  "(^|_)id$|^id$|^id_|^file$|filename|barcode|run_id|library_id|^lane$"

.is_id_like_colname <- function(col) {
  if (col %in% .AES_ID_CANDIDATES) return(TRUE)
  grepl(.AES_ID_LIKE_REGEX, col, ignore.case = TRUE)
}

.is_categorical_vec <- function(x) {
  if (is.factor(x) || is.character(x) || is.logical(x)) return(TRUE)
  if (is.numeric(x)) {
    nz <- x[!is.na(x)]
    if (!length(nz)) return(FALSE)
    n_unique <- length(unique(nz))
    # Treat numeric as categorical only if it has very few distinct values
    return(n_unique <= max(2L, ceiling(length(nz) / 2L)))
  }
  FALSE
}

.has_repeats <- function(x) {
  nz <- x[!is.na(x)]
  if (!length(nz)) return(FALSE)
  length(unique(nz)) < length(nz)
}

# Pick a grouping/color/shape column from sample_meta.
.resolve_aes_one <- function(sm, role = c("color", "shape", "group"),
                             preferred = NULL, exclude = character(0)) {
  role <- match.arg(role)
  cols <- colnames(sm)
  if (!length(cols)) return(NA_character_)

  # 1) Honor an explicit pre-set if it refers to a real column.
  if (!is.null(preferred) && length(preferred) == 1L &&
      is.character(preferred) && !is.na(preferred) && nzchar(preferred) &&
      preferred %in% cols) {
    return(preferred)
  }

  # 2) Prefer well-known biological columns (in declared order).
  for (col in intersect(.AES_GROUP_CANDIDATES, cols)) {
    if (col %in% exclude) next
    if (.is_id_like_colname(col)) next
    if (.is_categorical_vec(sm[[col]]) && .has_repeats(sm[[col]])) return(col)
  }

  # 3) Any non-ID categorical column with repeated values.
  for (col in cols) {
    if (col %in% exclude) next
    if (.is_id_like_colname(col)) next
    if (.is_categorical_vec(sm[[col]]) && .has_repeats(sm[[col]])) return(col)
  }

  # Shape is allowed to remain unresolved (it is genuinely optional).
  if (role == "shape") return(NA_character_)

  # 4) Last resort: first non-ID column, even if numeric / unique.
  for (col in cols) {
    if (col %in% exclude) next
    if (.is_id_like_colname(col)) next
    return(col)
  }
  NA_character_
}

# Pick a sample-id column from sample_meta.
.resolve_aes_id <- function(sm, preferred = NULL) {
  cols <- colnames(sm)
  if (!length(cols)) return(NA_character_)

  if (!is.null(preferred) && length(preferred) == 1L &&
      is.character(preferred) && !is.na(preferred) && nzchar(preferred) &&
      preferred %in% cols) {
    return(preferred)
  }

  hit <- intersect(.AES_ID_CANDIDATES, cols)
  if (length(hit)) return(hit[1L])

  # Fallback: any column whose values are all unique (one per sample).
  for (col in cols) {
    v <- sm[[col]]
    if (length(v) == nrow(sm) && length(unique(v)) == nrow(sm)) return(col)
  }
  NA_character_
}

# Public: returns list(sample_id, color, shape, group). Caches in env.
resolve_sample_meta_aes <- function(env, refresh = FALSE) {
  if (!refresh && exists(".aes_resolved", envir = env, inherits = FALSE)) {
    return(get(".aes_resolved", envir = env, inherits = FALSE))
  }

  sm <- if (exists("sample_meta", envir = env, inherits = FALSE)) {
    get("sample_meta", envir = env, inherits = FALSE)
  } else NULL

  default <- list(sample_id = NA_character_,
                  color     = NA_character_,
                  shape     = NULL,
                  group     = NA_character_)

  if (is.null(sm) || !is.data.frame(sm) || nrow(sm) == 0L || ncol(sm) == 0L) {
    assign(".aes_resolved", default, envir = env)
    return(default)
  }

  pref_id    <- if (exists("sample_id", envir = env, inherits = FALSE))
                  get("sample_id", envir = env, inherits = FALSE) else NULL
  pref_color <- if (exists("color",     envir = env, inherits = FALSE))
                  get("color",     envir = env, inherits = FALSE) else NULL
  pref_shape <- if (exists("shape",     envir = env, inherits = FALSE))
                  get("shape",     envir = env, inherits = FALSE) else NULL
  pref_group <- if (exists("group",     envir = env, inherits = FALSE))
                  get("group",     envir = env, inherits = FALSE) else NULL

  sample_id_col <- .resolve_aes_id(sm, preferred = pref_id)
  excl <- if (!is.na(sample_id_col)) sample_id_col else character(0)

  color_col <- .resolve_aes_one(sm, role = "color",
                                preferred = pref_color, exclude = excl)

  shape_col <- .resolve_aes_one(sm, role = "shape",
                                preferred = pref_shape,
                                exclude = c(excl, color_col))
  if (is.na(shape_col) || !nzchar(shape_col)) shape_col <- NULL

  group_col <- .resolve_aes_one(sm, role = "group",
                                preferred = pref_group, exclude = excl)
  # If group can't be resolved, fall back to color so the patterns plot has
  # a valid X-axis grouping.
  if (is.na(group_col) || !nzchar(group_col)) group_col <- color_col

  out <- list(sample_id = sample_id_col,
              color     = color_col,
              shape     = shape_col,
              group     = group_col)
  assign(".aes_resolved", out, envir = env)
  out
}

# Apply the resolver to env: write env$color / env$shape / env$group /
# env$sample_id, and mirror the resolved sample-id column to "SampleID" in
# both sample_meta and expr_long when missing. Returns the resolved aes list.
apply_sample_meta_aes <- function(env) {
  aes <- resolve_sample_meta_aes(env, refresh = TRUE)

  # Write resolved values back. Use NULL for shape if not resolvable so the
  # existing `if (!is.null(shape))` guards continue to work.
  if (!is.na(aes$color))     assign("color", aes$color, envir = env)
  if (!is.null(aes$shape))   assign("shape", aes$shape, envir = env) else
    if (exists("shape", envir = env, inherits = FALSE))
      assign("shape", NULL, envir = env)
  if (!is.na(aes$group))     assign("group", aes$group, envir = env)
  if (!is.na(aes$sample_id)) assign("sample_id", aes$sample_id, envir = env)

  # Mirror sample-id to "SampleID" in sample_meta when missing — several
  # downstream helpers read sample_meta$SampleID directly.
  if (exists("sample_meta", envir = env, inherits = FALSE)) {
    sm <- get("sample_meta", envir = env, inherits = FALSE)
    if (is.data.frame(sm) && !"SampleID" %in% colnames(sm) &&
        !is.na(aes$sample_id) && aes$sample_id %in% colnames(sm)) {
      sm$SampleID <- sm[[aes$sample_id]]
      assign("sample_meta", sm, envir = env)
    }
  }

  # Mirror sample-id to "SampleID" in expr_long when missing — the boxplot /
  # density helpers reference expr_long$SampleID by that name.
  if (exists("expr_long", envir = env, inherits = FALSE)) {
    el <- get("expr_long", envir = env, inherits = FALSE)
    if (is.data.frame(el) && !"SampleID" %in% colnames(el)) {
      candidate <- NULL
      if (!is.na(aes$sample_id) && aes$sample_id %in% colnames(el)) {
        candidate <- aes$sample_id
      } else if ("sample_id" %in% colnames(el)) {
        candidate <- "sample_id"
      } else {
        for (cand in .AES_ID_CANDIDATES) {
          if (cand %in% colnames(el)) { candidate <- cand; break }
        }
      }
      if (!is.null(candidate)) {
        el$SampleID <- el[[candidate]]
        assign("expr_long", el, envir = env)
      }
    }
  }

  invisible(aes)
}

# Build UI controls for the boxplot "Color by" selector
build_norm_boxplot_controls <- function(expr_long, color_var, shape_var = NULL,
                                         color_by_input_id = "prot_norm_boxplot_color_by") {
  stopifnot(color_var %in% colnames(expr_long))

  choices <- color_var
  if (!is.null(shape_var) && shape_var %in% colnames(expr_long) && shape_var != color_var) {
    choices <- c(color_var, shape_var)
  }

  shiny::div(
    style = "padding: 0 0 10px 0;",
    shiny::selectInput(
      inputId  = color_by_input_id,
      label    = "Color by:",
      choices  = choices,
      selected = color_var,
      width    = "280px"
    )
  )
}

# Build a static ggplot boxplot from long-format expression data
build_norm_boxplot_static <- function(expr_long, map_var) {
  stopifnot(all(c("SampleID", "value") %in% colnames(expr_long)))
  stopifnot(map_var %in% colnames(expr_long))

  # Keep sample order stable
  expr_long$SampleID <- factor(expr_long$SampleID, levels = unique(expr_long$SampleID))

  ggplot2::ggplot(
    expr_long,
    ggplot2::aes(
      x      = SampleID,
      y      = value,
      colour = .data[[map_var]]
    )
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::labs(
      x      = "Sample",
      y      = "log2 (normalized intensity)",
      colour = map_var
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1,
                                           size = 12, colour = "#2c3e50"),
      axis.text.y  = ggplot2::element_text(size = 12, colour = "#2c3e50"),
      axis.title.x = ggplot2::element_text(size = 14, face = "bold",
                                           margin = ggplot2::margin(t = 12)),
      axis.title.y = ggplot2::element_text(size = 14, face = "bold",
                                           margin = ggplot2::margin(r = 10)),
      legend.title = ggplot2::element_text(size = 13, face = "bold"),
      legend.text  = ggplot2::element_text(size = 12),
      plot.margin  = ggplot2::margin(t = 10, r = 14, b = 14, l = 10)
    )
}

#  Compute per-sample density curves into one long data.frame for Plotly
compute_density_df <- function(expr_long, sample_col = "SampleID", value_col = "value") {
  if (!all(c(sample_col, value_col) %in% colnames(expr_long))) return(NULL)
  
  expr_long <- expr_long[is.finite(expr_long[[value_col]]), , drop = FALSE]
  if (nrow(expr_long) == 0) return(NULL)
  
  sample_levels <- unique(expr_long[[sample_col]])
  expr_long[[sample_col]] <- factor(expr_long[[sample_col]], levels = sample_levels)
  
  dens_list <- lapply(levels(expr_long[[sample_col]]), function(sid) {
    v <- expr_long[[value_col]][expr_long[[sample_col]] == sid]
    v <- v[is.finite(v)]
    if (length(v) < 2 || length(unique(v)) < 2) return(NULL)
    
    d <- stats::density(v, na.rm = TRUE, n = 512)
    data.frame(SampleID = sid, x = d$x, y = d$y, stringsAsFactors = FALSE)
  })
  
  dens_df <- do.call(rbind, dens_list)
  if (is.null(dens_df) || nrow(dens_df) == 0) return(NULL)
  dens_df
}

# Build a Plotly multi-line density plot from the density data.frame
build_density_plotly <- function(dens_df) {
  plotly::plot_ly(
    data = dens_df,
    x = ~x,
    y = ~y,
    split = ~SampleID,
    type = "scatter",
    mode = "lines",
    text = ~SampleID,
    hovertemplate = "%{text}<extra></extra>"
  ) |>
    plotly::layout(
      hovermode = "closest",
      xaxis = list(title = "Intensity"),
      yaxis = list(title = "Density"),
      legend = list(orientation = "v")
    )
}


## ---- Subtab: Data structure ----

## PCA + Sample correlation 

# Compute PCA (prcomp) from normalized matrix (features x samples)
compute_pca <- function(expr_norm) {
  stats::prcomp(t(expr_norm))
}

# PCA scatter (PCx vs PCy) from prcomp + sample_meta
# shape can be NULL; color is required
pca_plot_from_prcomp <- function(pr, sample_meta, pc_x = 1, pc_y = 2, shape = NULL, color) {
  stopifnot(is.character(color), length(color) == 1, color %in% names(sample_meta))
  if (!is.null(shape)) {
    stopifnot(is.character(shape), length(shape) == 1, shape %in% names(sample_meta))
  }
  
  axis1 <- paste0("PC", pc_x)
  axis2 <- paste0("PC", pc_y)
  
  percentVar <- pr$sdev^2 / sum(pr$sdev^2)
  names(percentVar) <- colnames(pr$x)
  
  df <- data.frame(
    pr$x[, c(axis1, axis2), drop = FALSE],
    Color  = sample_meta[[color]],
    Sample = if ("SampleID" %in% names(sample_meta)) sample_meta$SampleID else rownames(sample_meta)
  )
  
  if (!is.null(shape)) {
    df$Shape <- sample_meta[[shape]]
    aes_map <- ggplot2::aes(x = .data[[axis1]], y = .data[[axis2]], color = Color, shape = Shape)
  } else {
    aes_map <- ggplot2::aes(x = .data[[axis1]], y = .data[[axis2]], color = Color)
  }
  
  ggplot2::ggplot(df, aes_map) +
    ggplot2::geom_point(size = 4) +
    ggplot2::xlab(sprintf("%s: %1.2f%% of variance", axis1, percentVar[axis1] * 100)) +
    ggplot2::ylab(sprintf("%s: %1.2f%% of variance", axis2, percentVar[axis2] * 100)) +
    ggplot2::theme_minimal()
}

# Sample correlation heatmap (native plotly heatmap)
make_sample_correlation_heatmap_plotly <- function(expr_norm) {
  sampleDists <- stats::dist(t(expr_norm))
  sampleDistMatrix <- as.matrix(sampleDists)

  # Hierarchical clustering for row/column ordering
  hc <- stats::hclust(sampleDists)
  ord <- hc$order
  sampleDistMatrix <- sampleDistMatrix[ord, ord]

  labels <- rownames(sampleDistMatrix)
  pal <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  colorscale <- lapply(seq_along(pal), function(i) list((i - 1) / (length(pal) - 1), pal[i]))

  plotly::plot_ly(
    x = labels, y = labels, z = sampleDistMatrix,
    type = "heatmap",
    colorscale = colorscale,
    hovertemplate = "Row: %{y}<br>Col: %{x}<br>Distance: %{z:.2f}<extra></extra>"
  ) |>
    plotly::layout(
      xaxis = list(tickangle = 45, automargin = TRUE),
      yaxis = list(autorange = "reversed", automargin = TRUE)
    )
}

# Shiny wrappers (env-based) — NEW names expected:
# pca_object, sample_meta, expr_norm, plus optional color/shape
get_pca_aesthetics <- function(env) {
  shape_col <- if (exists("shape", envir = env, inherits = FALSE)) get("shape", envir = env, inherits = FALSE) else NULL
  color_col <- if (exists("color", envir = env, inherits = FALSE)) get("color", envir = env, inherits = FALSE) else NULL
  list(shape = shape_col, color = color_col)
}

make_pca_plotly <- function(env, pc_x, pc_y) {
  shiny::validate(
    shiny::need(exists("pca_object",  envir = env, inherits = FALSE), "pca_object is not available."),
    shiny::need(exists("sample_meta", envir = env, inherits = FALSE), "sample_meta is not available.")
  )

  pr          <- get("pca_object",  envir = env, inherits = FALSE)
  sample_meta <- get("sample_meta", envir = env, inherits = FALSE)
  aes         <- get_pca_aesthetics(env)

  shiny::validate(
    shiny::need(!is.null(aes$color) && !is.na(aes$color) &&
                  nzchar(aes$color) && aes$color %in% colnames(sample_meta),
                "No grouping/color column could be resolved from sample_meta.")
  )

  axis1 <- paste0("PC", pc_x)
  axis2 <- paste0("PC", pc_y)

  percentVar <- pr$sdev^2 / sum(pr$sdev^2)
  names(percentVar) <- colnames(pr$x)

  df <- data.frame(
    x_val  = pr$x[, axis1],
    y_val  = pr$x[, axis2],
    Color  = sample_meta[[aes$color]],
    Sample = if ("SampleID" %in% names(sample_meta)) sample_meta$SampleID else rownames(sample_meta),
    stringsAsFactors = FALSE
  )
  if (!is.null(aes$shape)) df$Shape <- sample_meta[[aes$shape]]

  # Dynamic high-contrast color map: keyed to the actual group levels in the

  # data so each condition gets a stable, distinct color regardless of how many
  # traces plotly creates (e.g. when both color and symbol are mapped).
  # The palette is recycled if there are more groups than base colors, and
  # works for any dataset without hard-coding group names.
  palette_base <- c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e",
                    "#9467bd", "#8c564b", "#e377c2", "#17becf",
                    "#aec7e8", "#ff9896", "#98df8a", "#ffbb78")
  group_levels <- unique(as.character(df$Color))
  n <- length(group_levels)
  pal <- if (n <= length(palette_base)) {
    palette_base[seq_len(n)]
  } else {
    grDevices::colorRampPalette(palette_base)(n)
  }
  color_map <- setNames(pal, group_levels)

  p <- if (!is.null(aes$shape)) {
    plotly::plot_ly(
      df, x = ~x_val, y = ~y_val,
      color  = ~Color, symbol = ~Shape,
      colors = color_map,
      text   = ~paste0("Sample: ", Sample, "<br>", aes$color, ": ", Color,
                        "<br>", aes$shape, ": ", Shape),
      type = "scatter", mode = "markers",
      marker = list(size = 10),
      hovertemplate = "%{text}<extra></extra>"
    )
  } else {
    plotly::plot_ly(
      df, x = ~x_val, y = ~y_val,
      color  = ~Color,
      colors = color_map,
      text   = ~paste0("Sample: ", Sample, "<br>", aes$color, ": ", Color),
      type = "scatter", mode = "markers",
      marker = list(size = 10),
      hovertemplate = "%{text}<extra></extra>"
    )
  }

  plotly::layout(
    p,
    xaxis = list(
      title    = sprintf("%s: %1.2f%% of variance", axis1, percentVar[axis1] * 100),
      zeroline = FALSE
    ),
    yaxis = list(
      title    = sprintf("%s: %1.2f%% of variance", axis2, percentVar[axis2] * 100),
      zeroline = FALSE
    ),
    legend   = list(title = list(text = aes$color)),
    autosize = TRUE
  )
}

make_sample_corr_plotly <- function(env) {
  shiny::validate(
    shiny::need(exists("expr_norm", envir = env, inherits = FALSE), "expr_norm is not available.")
  )
  make_sample_correlation_heatmap_plotly(get("expr_norm", envir = env))
}

## ---- Subtab: DE Genes/Protein ----

# ---- Get contrasts table (expects: contrasts)
get_contrasts_tbl <- function(env) {
  if (!exists("contrasts", envir = env, inherits = FALSE)) return(NULL)
  x <- get("contrasts", envir = env, inherits = FALSE)
  if (!is.data.frame(x) || nrow(x) == 0) return(NULL)
  x
}

# ---- Update contrast dropdown from env
update_de_contrast_choices_from_env <- function(session, env, input_id) {
  ct <- get_contrasts_tbl(env)
  if (is.null(ct)) {
    shiny::updateSelectInput(session, input_id, choices = NULL)
    return(invisible(NULL))
  }
  ch <- setNames(ct$Contrast_name, ct$Contrast_name)
  shiny::updateSelectInput(session, input_id, choices = ch, selected = ch[1])
}

# ---- Build DESeq2 contrast vector c(Factor, Numerator, Denominator)
get_contrast_vec <- function(contrasts_tbl, chosen_name) {
  row <- dplyr::filter(contrasts_tbl, Contrast_name == chosen_name) |>
    dplyr::slice(1)
  
  c(
    as.character(row$Factor),
    as.character(row$Numerator),
    as.character(row$Denominator)
  )
}

# ---- Get contrast vector from env
get_contrast_vec_from_env <- function(env, chosen_name) {
  ct <- get_contrasts_tbl(env)
  shiny::validate(shiny::need(!is.null(ct), "Contrasts table not available."))
  get_contrast_vec(ct, chosen_name)
}

# ---- Get DE results by backend (DESeq2 or limma)
get_de_results_from_env <- function(env, chosen_contrast, backend = c("deseq2", "limma")) {
  backend <- match.arg(backend)
  
  if (backend == "deseq2") {
    shiny::validate(
      shiny::need(exists("de_model", envir = env, inherits = FALSE), "de_model not available."),
      shiny::need(exists("padj_cutoff", envir = env, inherits = FALSE), "padj_cutoff not available.")
    )

    dds <- get("de_model", envir = env, inherits = FALSE)
    shiny::validate(shiny::need(inherits(dds, "DESeqDataSet"), "de_model must be a DESeqDataSet."))

    ct <- get_contrasts_tbl(env)
    shiny::validate(shiny::need(!is.null(ct), "Contrasts table not available."))

    row <- dplyr::filter(ct, Contrast_name == chosen_contrast) |>
      dplyr::slice(1)
    shiny::validate(shiny::need(
      nrow(row) == 1,
      paste0("Chosen contrast '", chosen_contrast, "' not found in contrasts table.")
    ))

    factor_str <- as.character(row$Factor)
    num_str    <- as.character(row$Numerator)
    den_str    <- as.character(row$Denominator)
    shiny::validate(
      shiny::need(length(factor_str) == 1L && !is.na(factor_str) && nzchar(factor_str),
                  "Contrast row is missing 'Factor'."),
      shiny::need(length(num_str)    == 1L && !is.na(num_str)    && nzchar(num_str),
                  "Contrast row is missing 'Numerator'."),
      shiny::need(length(den_str)    == 1L && !is.na(den_str)    && nzchar(den_str),
                  "Contrast row is missing 'Denominator'.")
    )

    alpha <- get("padj_cutoff", envir = env, inherits = FALSE)

    # Use DESeq2's official 3-element character contrast form. It works for
    # both reference-level coding (single coefficient '<Factor>_<Lvl>_vs_<Ref>',
    # as used by current RNA-seq payloads) and cell-means coding
    # ('<Factor><Lvl>' per level, as used by older payloads). DESeq2 also
    # handles level-name mangling (e.g. 'Oxidative stress' -> 'Oxidative.stress'
    # inside resultsNames) internally for this form, so we no longer need to
    # rebuild coefficient names with paste0(Factor, Numerator).
    res <- tryCatch(
      DESeq2::results(dds,
                      contrast = c(factor_str, num_str, den_str),
                      alpha    = alpha),
      error = function(e) e
    )
    shiny::validate(shiny::need(
      !inherits(res, "error"),
      paste0("DESeq2::results failed for contrast '", chosen_contrast, "': ",
             if (inherits(res, "error")) conditionMessage(res) else "")
    ))
    return(res)
  }
  
  # limma backend
  shiny::validate(
    shiny::need(exists("de_model", envir = env, inherits = FALSE), "de_model not available.")
  )
  
  fit <- get("de_model", envir = env, inherits = FALSE)
  shiny::validate(shiny::need(inherits(fit, "MArrayLM"), "de_model must be MArrayLM."))
  
  shiny::validate(
    shiny::need(!is.null(colnames(fit$coefficients)), "limma fit has no coefficient names."),
    shiny::need(chosen_contrast %in% colnames(fit$coefficients),
                paste0("Chosen contrast '", chosen_contrast, "' not found in limma coefficients."))
  )
  
  tt <- limma::topTable(
    fit,
    coef = chosen_contrast,
    number = Inf,
    sort.by = "none"
  )
  
  tt2 <- tt[, c("logFC", "AveExpr", "P.Value", "adj.P.Val"), drop = FALSE]
  
  data.frame(
    log2FoldChange = tt2$logFC,
    baseMean       = 2^(tt2$AveExpr),
    pvalue         = tt2$P.Value,
    padj           = tt2$adj.P.Val,
    row.names      = rownames(tt2),
    check.names    = FALSE
  )
}

# ---- Build volcano dataframe
make_volcano_df <- function(res, log_fc_cutoff, padj_cutoff) {
  tibble::as_tibble(res, rownames = "gene") |>
    dplyr::filter(!is.na(pvalue)) |>
    dplyr::mutate(
      neglog10 = -log10(pvalue),
      highFC   = abs(log2FoldChange) > log_fc_cutoff,
      lowPADJ  = !is.na(padj) & padj <= padj_cutoff,
      signif   = highFC & lowPADJ,
      sig_lab  = ifelse(signif, "Significant", "Not significant")
    )
}

# ---- Build MA dataframe
make_ma_df <- function(res, padj_cutoff) {
  tibble::as_tibble(res, rownames = "gene") |>
    dplyr::filter(!is.na(baseMean), !is.na(log2FoldChange)) |>
    dplyr::mutate(
      sig_lab = ifelse(!is.na(padj) & padj <= padj_cutoff,
                       "Significant", "Not significant")
    )
}

# ---- Volcano df from env
make_volcano_df_from_env <- function(env, chosen_contrast, backend = c("deseq2", "limma")) {
  backend <- match.arg(backend)
  
  shiny::validate(
    shiny::need(exists("log_fc_cutoff", envir = env, inherits = FALSE), "log_fc_cutoff not available."),
    shiny::need(exists("padj_cutoff",   envir = env, inherits = FALSE), "padj_cutoff not available.")
  )
  
  res <- get_de_results_from_env(env, chosen_contrast, backend = backend)
  
  df <- make_volcano_df(
    res,
    log_fc_cutoff = get("log_fc_cutoff", envir = env, inherits = FALSE),
    padj_cutoff   = get("padj_cutoff",   envir = env, inherits = FALSE)
  )
  
  shiny::validate(shiny::need(nrow(df) > 0, "No rows to plot."))
  df |> dplyr::select(gene, log2FoldChange, neglog10, sig_lab)
}

# ---- MA df from env
make_ma_df_from_env <- function(env, chosen_contrast, backend = c("deseq2", "limma")) {
  backend <- match.arg(backend)
  
  shiny::validate(
    shiny::need(exists("padj_cutoff", envir = env, inherits = FALSE), "padj_cutoff not available.")
  )
  
  res <- get_de_results_from_env(env, chosen_contrast, backend = backend)
  
  df <- make_ma_df(
    res,
    padj_cutoff = get("padj_cutoff", envir = env, inherits = FALSE)
  )
  
  shiny::validate(shiny::need(nrow(df) > 0, "No data for MA plot."))
  df |> dplyr::select(gene, baseMean, log2FoldChange, sig_lab)
}

# ---- Volcano plotly from env
plot_volcano_plotly_from_env <- function(env, chosen_contrast, backend = c("deseq2", "limma"), source_id = "de_volcano") {
  backend <- match.arg(backend)
  df <- make_volcano_df_from_env(env, chosen_contrast, backend = backend)
  
  plotly::plot_ly(
    data   = df,
    x      = ~log2FoldChange,
    y      = ~neglog10,
    type   = "scattergl",
    mode   = "markers",
    key    = ~gene,
    text   = ~gene,
    color  = ~sig_lab,
    colors = c("Not significant" = "black", "Significant" = "red"),
    source = source_id,
    marker     = list(size = 6, opacity = 0.85),
    # Keep both selected and unselected points visually identical to the default
    # state — selection should not dim/recolor unselected points.
    selected   = list(marker = list(opacity = 0.85)),
    unselected = list(marker = list(opacity = 0.85)),
    hovertemplate = paste(
      "<b>%{text}</b><br>",
      "log2FC: %{x}<br>",
      "-log10(p): %{y}<extra></extra>"
    )
  ) %>%
    plotly::event_register("plotly_click") %>%
    plotly::event_register("plotly_selected") %>%
    plotly::layout(
      clickmode = "event",
      xaxis = list(title = "log2 fold change"),
      yaxis = list(title = "-log10 p-value")
    )
}

# ---- MA plotly from env
plot_ma_plotly_from_env <- function(env, chosen_contrast, backend = c("deseq2", "limma"), source_id = "de_ma") {
  backend <- match.arg(backend)
  df <- make_ma_df_from_env(env, chosen_contrast, backend = backend)
  
  plotly::plot_ly(
    data   = df,
    x      = ~baseMean,
    y      = ~log2FoldChange,
    type   = "scattergl",
    mode   = "markers",
    key    = ~gene,
    text   = ~gene,
    color  = ~sig_lab,
    colors = c("Not significant" = "gray70", "Significant" = "red"),
    source = source_id,
    marker     = list(size = 6, opacity = 0.9),
    # Keep both selected and unselected points visually identical to the default
    # state — selection should not dim/recolor unselected points.
    selected   = list(marker = list(opacity = 0.9)),
    unselected = list(marker = list(opacity = 0.9)),
    hovertemplate = paste(
      "<b>%{text}</b><br>",
      "baseMean: %{x}<br>",
      "log2FC: %{y}<extra></extra>"
    )
  ) %>%
    plotly::event_register("plotly_click") %>%
    plotly::event_register("plotly_selected") %>%
    plotly::layout(
      clickmode = "event",
      xaxis = list(title = "mean of normalized counts (log10)", type = "log"),
      yaxis = list(title = "log fold change")
    )
}

## ---- Subtab: Stats results ----
wire_stats_results <- function(output,
                               env,
                               data_dir,
                               id_prefix,
                               obj_de_stats = "de_summary",
                               obj_de_mat   = "de_expr_norm",
                               obj_all_mats = list(raw = "expr_norm"),
                               out_de_table_id,
                               out_results_buttons_id,
                               btn_results_all_id,
                               btn_results_de_id,
                               dl_all_mat_btn_ids = list(raw = NULL),
                               dl_de_mat_btn_id,
                               obj_results_all = "de_stats",
                               obj_results_de  = "de_final_results",
                               label_results_all = "Result tables - All genes",
                               label_results_de  = "Result tables - DE genes",
                               file_label_all = "All_genes",
                               file_label_de  = "DE_genes") {
  
  stopifnot(is.character(data_dir), nzchar(data_dir))
  stopifnot(is.environment(env))
  
  # Summary DE statistics table
  output[[out_de_table_id]] <- shiny::renderTable({
    shiny::validate(
      shiny::need(exists(obj_de_stats, envir = env, inherits = FALSE),
                  paste0(obj_de_stats, " not available."))
    )
    get(obj_de_stats, envir = env, inherits = FALSE)
  }, striped = TRUE, hover = TRUE, bordered = TRUE, rownames = TRUE)
  
  # Render Excel download buttons (from objects in env)
  output[[out_results_buttons_id]] <- shiny::renderUI({
    shiny::tagList(
      if (exists(obj_results_all, envir = env, inherits = FALSE)) {
        shiny::div(
          shiny::downloadButton(
            btn_results_all_id,
            label = label_results_all,
            class = "btn btn-success btn-sm",
            style = "margin:4px;"
          ),
          style = "margin-bottom:8px;"
        )
      } else {
        shiny::tags$p(paste0(obj_results_all, " not available."), style = "color:#888;")
      },
      
      if (exists(obj_results_de, envir = env, inherits = FALSE)) {
        shiny::div(
          shiny::downloadButton(
            btn_results_de_id,
            label = label_results_de,
            class = "btn btn-success btn-sm",
            style = "margin:4px;"
          ),
          style = "margin-bottom:8px;"
        )
      } else {
        shiny::tags$p(paste0(obj_results_de, " not available."), style = "color:#888;")
      }
    )
  })
  
  # Excel: All results + DE results.
  # The upstream pipeline stores the preformatted Excel files as raw byte
  # vectors directly in the RDS payload (verified: all six *_final_xlsx
  # objects across the three omics envs are `raw` vectors starting with
  # the ZIP magic `50 4b 03 04` = `PK\x03\x04`). Writing the bytes
  # verbatim with writeBin() reproduces the exact .xlsx file the pipeline
  # produced — preserving any formulas, styles, and sheets — without any
  # openxlsx (re-)serialization that would hang on a raw vector.
  .write_results_handler <- function(obj_name) {
    function(file) {
      if (!exists(obj_name, envir = env, inherits = FALSE)) {
        msg <- paste0(obj_name, " not available.")
        shiny::showNotification(msg, type = "error", duration = 6)
        stop(msg)
      }
      x <- get(obj_name, envir = env, inherits = FALSE)
      if (!is.raw(x)) {
        msg <- sprintf("%s has unexpected class [%s]; expected raw bytes of an .xlsx file.",
                       obj_name, paste(class(x), collapse = "/"))
        message("Stats Results download: ", msg)
        shiny::showNotification(msg, type = "error", duration = 8)
        stop(msg)
      }
      writeBin(x, file)
    }
  }

  output[[btn_results_all_id]] <- shiny::downloadHandler(
    filename = function() paste0(id_prefix, "_RESULTS_ALL_", file_label_all, ".xlsx"),
    content  = .write_results_handler(obj_results_all)
  )

  output[[btn_results_de_id]] <- shiny::downloadHandler(
    filename = function() paste0(id_prefix, "_RESULTS_DE_", file_label_de, ".xlsx"),
    content  = .write_results_handler(obj_results_de)
  )
  
  # Download: All matrices
  for (nm in names(obj_all_mats)) {
    local({
      nm_i <- nm
      obj_name_i <- obj_all_mats[[nm_i]]
      btn_id_i <- dl_all_mat_btn_ids[[nm_i]]
      
      if (is.null(btn_id_i) || !nzchar(btn_id_i)) return()
      
      output[[btn_id_i]] <- shiny::downloadHandler(
        filename = function() {
          paste0(id_prefix, "_", nm_i, "_matrix_", file_label_all, ".txt")
        },
        content = function(file) {
          shiny::validate(
            shiny::need(exists(obj_name_i, envir = env, inherits = FALSE),
                        paste0(obj_name_i, " not available."))
          )
          mat <- get(obj_name_i, envir = env, inherits = FALSE)
          
          df <- data.frame(
            Name = rownames(mat),
            as.data.frame(mat),
            check.names = FALSE
          )
          
          utils::write.table(df, file, sep = "\t",
                             quote = FALSE, row.names = FALSE, col.names = TRUE)
        }
      )
    })
  }
  
  # Download: DE matrix
  output[[dl_de_mat_btn_id]] <- shiny::downloadHandler(
    filename = function() {
      paste0(id_prefix, "_matrix_", file_label_de, ".txt")
    },
    content = function(file) {
      shiny::validate(
        shiny::need(exists(obj_de_mat, envir = env, inherits = FALSE),
                    paste0(obj_de_mat, " not available."))
      )
      mat <- get(obj_de_mat, envir = env, inherits = FALSE)
      
      df <- data.frame(
        Name = rownames(mat),
        as.data.frame(mat),
        check.names = FALSE
      )
      
      utils::write.table(df, file, sep = "\t",
                         quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
  )
}


## ---- Subtab: Clustering ----

# Check pheatmap row dendrogram integrity
.check_dend <- function(p) {
  !is.null(p$tree_row) && !is.null(p$tree_row$order) && !is.null(p$tree_row$labels)
}

# Build hierarchical clustering heatmap inputs (row-scaled + dendrogram order)
prepare_de_heatmap_inputs <- function(de_expr_norm, clust_heatmap_hier,
                                      palette = "RdYlBu", n_colors = 100) {
  if (missing(de_expr_norm)) stop("de_expr_norm not found.")
  if (missing(clust_heatmap_hier)) stop("clust_heatmap_hier not found.")
  if (!.check_dend(clust_heatmap_hier)) stop("Row dendrogram missing or incomplete.")
  
  mat <- as.matrix(de_expr_norm)
  mat <- t(scale(t(mat)))
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  if (nrow(mat) == 0 || ncol(mat) == 0) stop("Matrix is empty after scaling/cleaning.")
  
  ord_rows <- clust_heatmap_hier$tree_row$labels[clust_heatmap_hier$tree_row$order]
  ord_rows <- ord_rows[ord_rows %in% rownames(mat)]
  if (!length(ord_rows)) stop("No overlapping rows between dendrogram and matrix.")
  if (length(miss <- setdiff(rownames(mat), ord_rows))) ord_rows <- c(ord_rows, miss)
  mat <- mat[ord_rows, , drop = FALSE]
  
  rn <- rownames(mat); if (is.null(rn)) rn <- as.character(seq_len(nrow(mat)))
  cn <- colnames(mat); if (is.null(cn)) cn <- as.character(seq_len(ncol(mat)))
  
  hmcol_rev <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, palette))(n_colors))
  colorscale <- Map(function(pp, cc) list(pp, cc),
                    seq(0, 1, length.out = length(hmcol_rev)), hmcol_rev)
  
  vmax <- max(abs(mat), na.rm = TRUE)
  list(mat = mat, rn = rn, cn = cn, colorscale = colorscale, vmax = vmax)
}

# Plot hierarchical clustering heatmap (Plotly)
plot_de_heatmap_plotly <- function(prepped, source_id = NULL) {
  p <- plotly::plot_ly(
    x = prepped$cn, y = prepped$rn, z = prepped$mat, type = "heatmap",
    colorscale = prepped$colorscale,
    zmin = -prepped$vmax, zmax = prepped$vmax, zmid = 0,
    hovertemplate = "Feature: %{y}<br>Sample: %{x}<br>Z: %{z:.3f}<extra></extra>",
    colorbar = list(title = "Z-score"),
    source = source_id
  )
  p <- plotly::event_register(p, "plotly_click")
  plotly::layout(
    p,
    xaxis = list(title = "Samples", tickangle = 45),
    yaxis = list(title = "Features", autorange = "reversed", showticklabels = FALSE, ticks = ""),
    margin = list(l = 90, b = 100)
  )
}

# End-to-end hierarchical clustering heatmap builder
build_de_heatmap_plot <- function(de_expr_norm, clust_heatmap_hier,
                                  palette = "RdYlBu", n_colors = 100,
                                  source_id = NULL) {
  prepped <- prepare_de_heatmap_inputs(de_expr_norm, clust_heatmap_hier, palette, n_colors)
  plot_de_heatmap_plotly(prepped, source_id = source_id)
}

# ---------------------------------------------------------------------------
# Selected-features hierarchical heatmap (used by the Multi-* Explorers).
#
# Built dynamically from env$expr_norm restricted to the user's current
# selection — so non-DE features are included alongside DE features. The
# pipeline is:
#   1. Subset expr_norm to selected_ids that exist as rownames.
#   2. Row-wise Z-score the subset:  t(scale(t(...))).
#   3. Drop rows that became NA/NaN (zero-variance) or Inf.
#   4. Recompute hierarchical clustering (dist + hclust) ON THIS SUBSET
#      and reorder rows by the resulting dendrogram. This is independent of
#      the precomputed clust_heatmap_hier, which only covers DE-significant
#      features.
#   5. Build the plotly inputs (matrix + row/col labels + colorscale +
#      symmetric color limits).
# Returns NULL if < 2 valid features remain after cleaning; callers display
# a friendly "Select at least 2 valid features..." message.
# ---------------------------------------------------------------------------
prepare_selected_features_heatmap_inputs <- function(env, selected_ids,
                                                     palette = "RdYlBu",
                                                     n_colors = 100) {
  selected_ids <- as.character(selected_ids)
  selected_ids <- selected_ids[!is.na(selected_ids) & nzchar(selected_ids)]
  if (length(selected_ids) < 2L) return(NULL)

  # Build from env$expr_norm (NOT de_expr_norm) so non-DE features are
  # included. Z-scores are computed on the selected-feature subset only,
  # and the row order is determined by a fresh hierarchical clustering on
  # that subset — independent of the precomputed clust_heatmap_hier object
  # which only covers DE-significant features.
  if (!exists("expr_norm", envir = env, inherits = FALSE)) return(NULL)
  en <- get("expr_norm", envir = env, inherits = FALSE)
  if (is.null(en)) return(NULL)

  # Resolve the row-membership intersection BEFORE materializing the full
  # expr_norm as a matrix. For RNAseq this matrix is ~10–30k rows; calling
  # as.matrix() on the whole thing every time the user adds/removes a
  # feature is wasted work — we only need the <=100 selected rows. By
  # subsetting first via [present, ] (which preserves rownames on both
  # matrices and data frames), the final as.matrix() runs on a tiny
  # matrix, regardless of how expr_norm is stored upstream.
  en_rn <- rownames(en)
  if (is.null(en_rn)) return(NULL)

  # Keep only selected IDs that exist in expr_norm; missing IDs are dropped
  # silently per spec.
  present <- intersect(selected_ids, en_rn)
  if (length(present) < 2L) return(NULL)

  sub <- tryCatch(as.matrix(en[present, , drop = FALSE]),
                  error = function(e) NULL)
  if (is.null(sub) || ncol(sub) == 0L) return(NULL)

  # Row-wise Z-score on the subset. scale() yields NaN for zero-variance
  # rows; we replace any Inf with NA and drop rows that became all-NA so
  # downstream distance computation stays well-defined.
  z_mat <- t(scale(t(sub)))
  is_inf <- is.infinite(z_mat)
  if (any(is_inf)) z_mat[is_inf] <- NA
  z_mat <- z_mat[stats::complete.cases(z_mat), , drop = FALSE]
  if (nrow(z_mat) < 2L || ncol(z_mat) == 0L) return(NULL)

  # Recompute hierarchical clustering on the selected-feature Z-score
  # matrix. Default Euclidean distance + complete linkage matches typical
  # heatmap conventions. On failure (single-feature, degenerate distance)
  # fall back to selection order.
  row_ord <- tryCatch({
    d <- stats::dist(z_mat)
    h <- stats::hclust(d)
    h$order
  }, error = function(e) seq_len(nrow(z_mat)))
  z_mat <- z_mat[row_ord, , drop = FALSE]

  rn <- rownames(z_mat); if (is.null(rn)) rn <- as.character(seq_len(nrow(z_mat)))
  cn <- colnames(z_mat); if (is.null(cn)) cn <- as.character(seq_len(ncol(z_mat)))

  hmcol_rev <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, palette))(n_colors))
  colorscale <- Map(function(pp, cc) list(pp, cc),
                    seq(0, 1, length.out = length(hmcol_rev)), hmcol_rev)

  vmax <- max(abs(z_mat), na.rm = TRUE)
  list(mat = z_mat, rn = rn, cn = cn, colorscale = colorscale, vmax = vmax)
}

# Render the selected-features heatmap. Same visual style as
# plot_de_heatmap_plotly(). Y-axis tick labels are hidden to keep the plot
# clean when many features are selected; feature IDs remain visible on hover
# via the hovertemplate.
plot_selected_features_heatmap_plotly <- function(prepped, source_id = NULL) {
  p <- plotly::plot_ly(
    x = prepped$cn, y = prepped$rn, z = prepped$mat, type = "heatmap",
    colorscale = prepped$colorscale,
    zmin = -prepped$vmax, zmax = prepped$vmax, zmid = 0,
    hovertemplate = "Feature: %{y}<br>Sample: %{x}<br>Z: %{z:.3f}<extra></extra>",
    colorbar = list(title = "Z-score"),
    source = source_id
  )
  plotly::layout(
    p,
    xaxis = list(title = "Samples", tickangle = 45),
    yaxis = list(title = "Features", autorange = "reversed",
                 showticklabels = FALSE, ticks = "", automargin = TRUE),
    margin = list(l = 120, b = 100)
  )
}

#' §17 S4: prepare a GROUPED features heatmap — the selected-features Z-score
#' matrix (reusing prepare_selected_features_heatmap_inputs) with rows ORDERED by
#' group (group-major; within-group order kept from the global clustering) and a
#' per-row `groups` vector aligned to the final row order. No gene-count cap.
#' `gene_groups` is a named vector gene -> group label. Returns NULL if < 2 genes
#' resolve. Generic + omics-agnostic (reads only env$expr_norm).
#' @noRd
prepare_grouped_heatmap_inputs <- function(env, gene_groups,
                                           palette = "RdYlBu", n_colors = 100) {
  gene_groups <- gene_groups[!is.na(names(gene_groups)) & nzchar(names(gene_groups))]
  if (length(gene_groups) < 2L) return(NULL)
  base <- prepare_selected_features_heatmap_inputs(env, names(gene_groups), palette, n_colors)
  if (is.null(base)) return(NULL)
  g <- unname(gene_groups[base$rn]); g[is.na(g)] <- "Other"
  ord <- order(factor(g, levels = sort(unique(g))))   # stable: group-major, keep within-group order
  base$mat    <- base$mat[ord, , drop = FALSE]
  base$rn     <- base$rn[ord]
  base$groups <- g[ord]
  base
}

#' §17 S4: render a grouped features heatmap (D13) — rows split by group with
#' visual WHITESPACE gaps between groups (blank NA rows, no separator line) plus
#' left-margin group labels, from existing expr_norm. Y tick labels hidden (many
#' rows); group is the visible organizing annotation. Native plotly.
#' @noRd
plot_grouped_heatmap_plotly <- function(prepped, source_id = NULL, gap = 2L) {
  if (is.null(prepped) || is.null(prepped$groups)) return(NULL)
  g <- prepped$groups; mat <- prepped$mat; nc <- ncol(mat)
  tbl <- rle(g); k <- length(tbl$lengths)
  # Stack each group's rows with `gap` blank (NA) rows between them; NA cells draw
  # nothing, so the groups read as separated by clean whitespace (item 6).
  parts <- list(); centers <- numeric(k); rc <- 0L; from <- 1L
  for (i in seq_len(k)) {
    len <- tbl$lengths[i]
    parts[[length(parts) + 1]] <- mat[from:(from + len - 1L), , drop = FALSE]
    centers[i] <- rc + (len + 1L) / 2                     # group label y-centre
    rc <- rc + len; from <- from + len
    if (i < k) { parts[[length(parts) + 1]] <- matrix(NA_real_, nrow = gap, ncol = nc); rc <- rc + gap }
  }
  zg <- do.call(rbind, parts); n <- nrow(zg)
  p <- plotly::plot_ly(
    x = prepped$cn, y = seq_len(n), z = zg, type = "heatmap",
    colorscale = prepped$colorscale, zmin = -prepped$vmax, zmax = prepped$vmax, zmid = 0,
    hovertemplate = "Sample: %{x}<br>Z: %{z:.3f}<extra></extra>",
    colorbar = list(title = "Z-score"), source = source_id)
  anns <- Map(function(m, lab) list(
    x = -0.008, xref = "paper", y = m, yref = "y", text = lab, showarrow = FALSE,
    xanchor = "right", font = list(size = 12)), centers, tbl$values)
  plotly::layout(p, annotations = anns,
    xaxis = list(title = "Samples", tickangle = 45, tickfont = list(size = 12)),
    yaxis = list(title = "", autorange = "reversed", showticklabels = FALSE,
                 ticks = "", automargin = TRUE),
    margin = list(l = 140, b = 100))
}

# Sort cluster ID labels numerically when all are numeric, lexically otherwise.
# Why: cluster IDs are character vectors so "10" sorts before "2" by default;
# users expect 1, 2, ..., 10, 11 in dropdowns, legends, and heatmap row order.
.sort_cluster_ids <- function(x) {
  x <- as.character(x)
  if (!length(x)) return(x)
  nums <- suppressWarnings(as.numeric(x))
  if (!any(is.na(nums))) x[order(nums)] else sort(x)
}

# Subset sample_meta to the columns of an expression matrix, using sample IDs.
# Why: Metabolomics keeps QC rows in sample_meta (PCA depends on the unfiltered
# 22-row table) while de_expr_norm has QCs stripped (16 cols). Any code that
# aligns sample_meta with expression columns by row position must subset by
# sample ID first. For RNAseq/Proteomics the row count already matches and this
# is effectively a no-op reorder.
align_sample_meta_to_cols <- function(sample_meta, cols) {
  cols <- as.character(cols)
  if (is.null(rownames(sample_meta)) || !all(cols %in% rownames(sample_meta))) {
    return(sample_meta)
  }
  sample_meta[cols, , drop = FALSE]
}

# Build Plotly colorscale for heatmaps
.make_colorscale <- function(palette = "RdYlBu", n_colors = 64, reverse = TRUE) {
  cols <- NULL
  if (requireNamespace("RColorBrewer", quietly = TRUE) &&
      palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    maxn <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    cols <- RColorBrewer::brewer.pal(min(maxn, max(3, n_colors)), palette)
    if (n_colors > length(cols)) cols <- grDevices::colorRampPalette(cols)(n_colors)
  } else {
    base <- c("#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8",
              "#FFFFBF","#FEE090","#FDAE61","#F46D43","#D73027","#A50026")
    cols <- grDevices::colorRampPalette(base)(n_colors)
  }
  if (reverse) cols <- rev(cols)
  pos <- seq(0, 1, length.out = length(cols))
  lapply(seq_along(cols), function(i) list(pos[i], cols[i]))
}

# Order rows within clusters using correlation + hclust
.order_rows_within_clusters <- function(mat, clv, cor_method = "pearson", hclust_method = "ward.D2") {
  stopifnot(nrow(mat) == length(clv))
  lvls <- levels(clv)
  new_order_idx <- integer(0)
  
  for (lv in lvls) {
    idx <- which(clv == lv)
    sub <- mat[idx, , drop = FALSE]
    if (nrow(sub) >= 2) {
      C <- stats::cor(t(sub), method = cor_method, use = "pairwise.complete.obs")
      D <- stats::as.dist(1 - C)
      hc <- stats::hclust(D, method = hclust_method)
      new_order_idx <- c(new_order_idx, idx[hc$order])
    } else {
      new_order_idx <- c(new_order_idx, idx)
    }
  }
  
  list(mat = mat[new_order_idx, , drop = FALSE], clv = clv[new_order_idx], order_idx = new_order_idx)
}

# Prepare partition clustering inputs (row-scaled + cluster order)
prepare_partition_inputs <- function(de_expr_norm, clust_partition,
                                     palette = "RdYlBu", n_colors = 64,
                                     order_within = c("hclust", "mean"),
                                     hclust_method = "ward.D2",
                                     cor_method = "pearson") {
  if (missing(de_expr_norm)) stop("de_expr_norm not found.")
  if (missing(clust_partition)) stop("clust_partition not found.")
  order_within <- match.arg(order_within)
  
  mat <- as.matrix(de_expr_norm)
  if (is.null(rownames(mat))) stop("de_expr_norm must have rownames (feature IDs).")
  if (is.null(names(clust_partition))) stop("clust_partition must be a named vector with feature names.")
  
  common <- names(clust_partition)[names(clust_partition) %in% rownames(mat)]
  if (!length(common)) stop("No overlap between rownames(de_expr_norm) and names(clust_partition).")
  
  mat <- mat[common, , drop = FALSE]
  clv_chr <- as.character(clust_partition[common])
  clv <- factor(clv_chr, levels = .sort_cluster_ids(unique(clv_chr)))
  names(clv) <- common  # factor() drops names; restore so name-based subset below works

  mat <- t(scale(t(mat)))
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  clv <- droplevels(clv[rownames(mat)])
  
  if (order_within == "hclust") {
    ord0 <- order(clv)
    mat0 <- mat[ord0, , drop = FALSE]
    clv0 <- droplevels(clv[ord0])
    oo <- .order_rows_within_clusters(mat0, clv0, cor_method = cor_method, hclust_method = hclust_method)
    mat <- oo$mat
    clv <- oo$clv
  } else {
    row_means <- rowMeans(mat, na.rm = TRUE)
    ord <- order(clv, -row_means)
    mat <- mat[ord, , drop = FALSE]
    clv <- clv[ord]
  }
  
  rn <- rownames(mat); if (is.null(rn)) rn <- as.character(seq_len(nrow(mat)))
  cn <- colnames(mat); if (is.null(cn)) cn <- as.character(seq_len(ncol(mat)))
  
  colorscale <- .make_colorscale(palette, n_colors, reverse = TRUE)
  vmax <- max(abs(mat), na.rm = TRUE)
  
  list(mat = mat, rn = rn, cn = cn, colorscale = colorscale, vmax = vmax, clusters = clv)
}

# Plot partition clustering heatmap (Plotly)
plot_partition_heatmap_plotly <- function(prepped, source_id = NULL) {
  use_gl <- (prod(dim(prepped$mat)) > 2e5)
  hm_type <- if (use_gl) "heatmapgl" else "heatmap"
  
  p <- plotly::plot_ly(
    type = hm_type,
    x = prepped$cn, y = prepped$rn, z = prepped$mat,
    colorscale = prepped$colorscale,
    zmin = -prepped$vmax, zmax = prepped$vmax, zmid = 0,
    hovertemplate = "Feature: %{y}<br>Sample: %{x}<br>Z: %{z:.3f}<extra></extra>",
    colorbar = list(title = "Z-score"),
    source = source_id
  )
  p <- plotly::event_register(p, "plotly_click")
  plotly::layout(
    p,
    xaxis = list(title = "Samples", tickangle = 45),
    yaxis = list(title = "Features", autorange = "reversed", showticklabels = FALSE, ticks = ""),
    margin = list(l = 90, b = 100)
  )
}

# Get latest results Excel file by prefix (RESULTS_DE / RESULTS_ALL)
get_latest_results_file <- function(data_dir, prefix = "RESULTS_DE") {
  paths <- list.files(
    data_dir,
    pattern = paste0("^", prefix, "_.*\\.(xlsx|xls)$"),
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(paths)) return(NULL)
  paths[which.max(file.info(paths)$mtime)]
}

# Binary pattern: available patterns that have heatmaps
bp_available_patterns <- function(clust_patterns, clust_heatmaps_by_pattern) {
  
  # patterns from clust_patterns (data.frame with best_pattern)
  if (is.data.frame(clust_patterns)) {
    pat <- unique(na.omit(as.character(clust_patterns$best_pattern)))
    pat <- trimws(pat)
    pat <- pat[nzchar(pat)]
  } else {
    pat <- unique(na.omit(as.character(clust_patterns)))
    pat <- trimws(pat)
    pat <- pat[nzchar(pat)]
  }
  
  hm_names <- names(clust_heatmaps_by_pattern)
  if (is.null(hm_names)) hm_names <- character()
  
  pat0 <- gsub("^pattern_", "", pat)
  hm0  <- gsub("^pattern_", "", hm_names)
  
  hm_names[hm0 %in% pat0]
}

# Binary pattern: TRUE if pattern objects exist and at least one pattern is valid
bp_data_available_env <- function(env) {
  if (!exists("clust_patterns", envir = env, inherits = FALSE)) return(FALSE)
  if (!exists("clust_heatmaps_by_pattern", envir = env, inherits = FALSE)) return(FALSE)

  pat <- get("clust_patterns", envir = env)
  hm  <- get("clust_heatmaps_by_pattern", envir = env)
  if (is.null(pat) || is.null(hm)) return(FALSE)

  length(bp_available_patterns(pat, hm)) > 0
}

# ============================================================================
# CLUSTERING AVAILABILITY PREDICATES
# ----------------------------------------------------------------------------
# Cheap defensive checks the wire-layer uses BEFORE calling the heavy
# prepare_*_inputs() functions. The prepare_* functions are strict (they
# stop() on contract violations); these predicates let the renderers fall
# back to a friendly UI message instead of letting that stop() bubble up
# as a red Shiny error.
# ============================================================================

# TRUE if x is a non-empty vector with at least one non-empty name.
.is_named_nonempty_vec <- function(x) {
  if (is.null(x)) return(FALSE)
  if (!length(x)) return(FALSE)
  nm <- names(x)
  if (is.null(nm)) return(FALSE)
  any(nzchar(nm))
}

# Is the object usable as a renderable static figure (pheatmap gtable or grob)?
.is_renderable_grob <- function(x) {
  if (is.null(x)) return(FALSE)
  inherits(x, c("gtable", "grob", "gTree", "recordedplot", "ggplot"))
}

# TRUE if hierarchical-clustering heatmap inputs are present and well-formed.
hc_data_available_env <- function(env) {
  if (!exists("de_expr_norm",       envir = env, inherits = FALSE)) return(FALSE)
  if (!exists("clust_heatmap_hier", envir = env, inherits = FALSE)) return(FALSE)
  hier <- get("clust_heatmap_hier", envir = env, inherits = FALSE)
  if (is.null(hier)) return(FALSE)
  den <- get("de_expr_norm", envir = env, inherits = FALSE)
  if (is.null(den)) return(FALSE)
  .check_dend(hier)
}

# TRUE if partition-clustering inputs are present and well-formed.
partition_data_available_env <- function(env) {
  if (!exists("de_expr_norm",   envir = env, inherits = FALSE)) return(FALSE)
  if (!exists("clust_partition", envir = env, inherits = FALSE)) return(FALSE)
  cp <- get("clust_partition", envir = env, inherits = FALSE)
  if (!.is_named_nonempty_vec(cp)) return(FALSE)
  den <- get("de_expr_norm", envir = env, inherits = FALSE)
  if (is.null(den)) return(FALSE)
  rn <- tryCatch(rownames(as.matrix(den)), error = function(e) NULL)
  if (is.null(rn) || !length(rn)) return(FALSE)
  any(names(cp) %in% rn)
}

# TRUE if a static clustering figure object exists and is renderable.
clust_static_fig_available_env <- function(env, name) {
  if (!exists(name, envir = env, inherits = FALSE)) return(FALSE)
  .is_renderable_grob(get(name, envir = env, inherits = FALSE))
}

# Standard "not available" message for clustering UI.
.clust_unavailable_msg <- function(text = "Clustering results are not available for this dataset.") {
  shiny::tags$div(text, style = "color:#888; padding:8px 0;")
}

# Binary pattern: draw static pheatmap gtable
bp_draw_static <- function(hm) {
  grid::grid.newpage()
  if (!is.null(hm$gtable)) {
    grid::grid.draw(hm$gtable)
  } else if (!is.null(hm$pheatmap$gtable)) {
    grid::grid.draw(hm$pheatmap$gtable)
  } else {
    stop("No gtable found in hm (expected hm$gtable or hm$pheatmap$gtable).")
  }
}

# Binary pattern: build interactive Plotly heatmap safely
bp_build_plotly_safe <- function(hm) {
  stopifnot(is.list(hm), !is.null(hm$mat))
  mat <- hm$mat
  
  # Match pheatmap ordering (supports character names or numeric indices)
  # Rows
  if (!is.null(hm$row_order)) {
    ro <- hm$row_order
    if (is.character(ro)) {
      ro <- ro[ro %in% rownames(mat)]
      if (length(ro)) mat <- mat[ro, , drop = FALSE]
    } else {
      ro <- ro[ro >= 1 & ro <= nrow(mat)]
      if (length(ro)) mat <- mat[ro, , drop = FALSE]
    }
  } else if (!is.null(hm$tree_row) && !is.null(hm$tree_row$order)) {
    ro <- hm$tree_row$order
    ro <- ro[ro >= 1 & ro <= nrow(mat)]
    if (length(ro)) mat <- mat[ro, , drop = FALSE]
  }
  
  # Cols
  if (!is.null(hm$col_order)) {
    co <- hm$col_order
    if (is.character(co)) {
      co <- co[co %in% colnames(mat)]
      if (length(co)) mat <- mat[, co, drop = FALSE]
    } else {
      co <- co[co >= 1 & co <= ncol(mat)]
      if (length(co)) mat <- mat[, co, drop = FALSE]
    }
  }
  
  # round for display
  mat <- round(mat, 2)
  
  x_labs <- colnames(mat)
  y_labs <- rownames(mat)
  
  x_labs[is.na(x_labs) | x_labs == ""] <- paste0("S", seq_along(x_labs))
  x_labs <- make.unique(x_labs, sep = "_")
  
  y_labs[is.na(y_labs) | y_labs == ""] <- paste0("Feature_", seq_along(y_labs))
  
  txt <- outer(y_labs, x_labs, function(yg, xg) paste0(
    "<b>Feature</b>: ", yg, "<br>",
    "<b>Sample</b>: ", xg
  ))
  
  r <- max(abs(range(mat, na.rm = TRUE)))
  
  # Proper Plotly colorscale (list of [pos,color])
  pal <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "RdYlBu"))(100))
  colorscale <- lapply(seq_along(pal), function(i) {
    list((i - 1) / (length(pal) - 1), pal[i])
  })
  
  # Always use standard heatmap trace (more robust than heatmapgl)
  plt <- plotly::plot_ly(
    x = x_labs,
    y = y_labs,
    z = mat,
    type = "heatmap",
    colorscale = colorscale,
    zmin = -r, zmax = r, zmid = 0,
    zsmooth = FALSE,
    text = txt,
    hovertemplate = "%{text}<br><b>Z-score</b>: %{z:.2f}<extra></extra>"
  )
  
  plt <- plotly::layout(
    plt,
    xaxis = list(tickangle = 45, automargin = TRUE),
    yaxis = list(showticklabels = FALSE, ticks = "", autorange = "reversed", automargin = TRUE)
  )
  
  plotly::config(plt, displayModeBar = TRUE)
}

# Detect feature column in RESULTS_DE (gene/protein).
# When `feats` is supplied and standard column-name detection fails, fall back
# to picking the column whose values overlap most with the expected IDs. This
# rescues schemas like the Proteomics de_final_table, where Excel-style header
# rows are embedded as data so the ID column ends up named "Sample_ID" with
# noise rows ("Condition", "Protein", NAs) above the real protein IDs.
.detect_feature_col <- function(df, feats = NULL) {
  cn <- colnames(df)

  # Priority exact matches
  priority <- c("FeatureID", "Protein")
  hit <- intersect(priority, cn)
  if (length(hit)) return(hit[1])

  # Other common fallbacks
  cand <- c("gene", "Gene", "protein", "Protein", "ID", "id", "feature", "Feature", "feature_id", "gene_id")
  hit2 <- intersect(cand, cn)
  if (length(hit2)) return(hit2[1])

  # Case-insensitive match as last resort
  cn_l <- tolower(trimws(cn))
  if ("featureid" %in% cn_l) return(cn[which(cn_l == "featureid")[1]])
  if ("protein" %in% cn_l)   return(cn[which(cn_l == "protein")[1]])

  # Content-based fallback: column with the most matches against expected IDs.
  if (!is.null(feats) && length(feats)) {
    feats_chr <- unique(as.character(feats))
    counts <- vapply(cn, function(c) {
      v <- df[[c]]
      if (!(is.character(v) || is.factor(v))) return(0L)
      sum(as.character(v) %in% feats_chr, na.rm = TRUE)
    }, integer(1))
    if (max(counts) > 0) return(cn[which.max(counts)])
  }

  NULL
}

# ---------------------------------------------------------------------------
# CLUSTER-PROFILE FACTOR RESOLUTION (Partition Clustering -> Cluster Profiles)
# ---------------------------------------------------------------------------
# Decides which experimental factor goes on the X-axis and which one (if any)
# splits the profile into separate lines. The ONLY source of truth is the
# already-resolved aesthetic pair env$color / env$shape written by
# apply_sample_meta_aes() — the same pair PCA (get_pca_aesthetics()), the
# Expression Profile "Color bars by" selectors and the normalization boxplots
# read. No factor name and no level name is ever hard-coded, so the same code
# serves any experimental design.
#
# Rule (fixed — do NOT invert):
#   two valid, distinct, varying factors -> X = env$shape, series = env$color
#   anything else                        -> one line, X = the single factor
#
# `sm` MUST be the sample_meta already aligned to the clustering columns via
# align_sample_meta_to_cols(). Both the usable/varies checks and the level
# order are read from it, so QC or otherwise unused levels that are absent
# from de_expr_norm can never reach the plot.
#
# The two-factor branch requires BOTH columns to be present in `sm` and to
# carry at least two distinct non-missing values. Every other case (shape
# absent, shape == color, column missing from the aligned metadata, constant
# column) degrades to the pre-existing one-factor plot instead of erroring.
#
# One-factor X column: env$group when usable, else env$color, else env$shape.
# env$group is read ONLY here and ONLY as a scalar — it stays the historical
# single-factor setting and is deliberately NOT extended to carry two values.
.resolve_cluster_profile_factors <- function(env, sm) {
  get_slot <- function(nm) {
    if (!exists(nm, envir = env, inherits = FALSE)) return(NULL)
    v <- get(nm, envir = env, inherits = FALSE)
    if (is.null(v) || !length(v)) return(NULL)
    v <- as.character(v)[1L]
    if (is.na(v) || !nzchar(v)) return(NULL)
    v
  }
  # Usable = the column exists in the aligned metadata.
  usable <- function(col) !is.null(col) && col %in% colnames(sm)
  # Varies = usable AND carries >= 2 distinct non-missing values.
  varies <- function(col) {
    if (!usable(col)) return(FALSE)
    v <- as.character(sm[[col]])
    length(unique(v[!is.na(v)])) >= 2L
  }

  color_col <- get_slot("color")
  shape_col <- get_slot("shape")
  group_col <- get_slot("group")

  if (varies(color_col) && varies(shape_col) &&
      !identical(color_col, shape_col)) {
    return(list(x_col = shape_col, series_col = color_col, n_factors = 2L))
  }

  x_col <- NULL
  for (cand in c(group_col, color_col, shape_col)) {
    if (usable(cand)) { x_col <- cand; break }
  }
  list(x_col      = x_col,
       series_col = NULL,
       n_factors  = if (is.null(x_col)) 0L else 1L)
}

# Level order for a profile axis or series: order of first appearance in the
# aligned metadata (the convention the single-factor plot has always used),
# with missing values dropped so NA can never become a phantom level. Reading
# the values rather than levels() also drops unused factor levels.
.cluster_profile_levels <- function(v) {
  v <- as.character(v)
  unique(v[!is.na(v)])
}

# Wire the full clustering tab (HC + Partition + Binary Pattern) from NEW object names
wire_clustering <- function(session, input, output, env, data_dir,
                            id, source,
                            label = list(feature = "genes", results_de_obj = "de_final_table")) {
  
  # Render hierarchical clustering heatmap
  output[[id$heatmap_hc]] <- plotly::renderPlotly({
    shiny::validate(shiny::need(
      hc_data_available_env(env),
      "Clustering results are not available for this dataset."
    ))
    tryCatch({
      build_de_heatmap_plot(
        de_expr_norm = get("de_expr_norm", envir = env),
        clust_heatmap_hier = get("clust_heatmap_hier", envir = env),
        palette = "RdYlBu",
        n_colors = 100,
        source_id = source$hc
      )
    }, error = function(e) {
      # Last-resort fallback: surface a friendly inline message via validate
      # (no red toast notification — clustering being absent is expected for
      # some payloads and shouldn't look like an error).
      shiny::validate(shiny::need(FALSE, "Clustering heatmap could not be rendered."))
    })
  })
  
  # Partition title — uses the predicate so any missing/NULL/malformed
  # clust_partition yields the friendly message instead of a crash.
  output[[id$pc_title_ui]] <- shiny::renderUI({
    if (!partition_data_available_env(env)) return(.clust_unavailable_msg())
    tryCatch({
      cl <- get("clust_partition", envir = env)
      k  <- length(unique(na.omit(cl)))
      ng <- sum(rownames(get("de_expr_norm", envir = env)) %in% names(cl))
      shiny::h4(
        sprintf("Partition clustering of DE %s (k=%d; n=%d)", label$feature, k, ng),
        class = "plot-title full"
      )
    }, error = function(e) .clust_unavailable_msg())
  })

  # Cluster choices
  shiny::observe({
    if (partition_data_available_env(env)) {
      cls <- unique(as.character(na.omit(get("clust_partition", envir = env))))
      cls <- .sort_cluster_ids(cls)
      if (length(cls)) {
        shiny::updateSelectInput(session, id$pc_cluster_pick, choices = cls, selected = cls[1])
      } else {
        shiny::updateSelectInput(session, id$pc_cluster_pick, choices = NULL)
      }
    } else {
      shiny::updateSelectInput(session, id$pc_cluster_pick, choices = NULL)
    }
  })

  # Precompute partition inputs. Returns NULL (instead of stop()ing) when
  # clust_partition is missing/NULL/malformed, so downstream req()s halt
  # silently and the UI shows the friendly message via the renderers below.
  prepped_all <- shiny::reactive({
    if (!partition_data_available_env(env)) return(NULL)
    tryCatch(
      prepare_partition_inputs(
        de_expr_norm     = get("de_expr_norm", envir = env),
        clust_partition  = get("clust_partition", envir = env),
        palette          = "RdYlBu",
        n_colors         = 64,
        order_within     = "hclust"
      ),
      error = function(e) NULL
    )
  })
  
  # Debounced cluster pick
  cluster_pick_d <- shiny::debounce(shiny::reactive(input[[id$pc_cluster_pick]]), 150)
  
  # Cluster title
  output[[id$pc_cluster_title_ui]] <- shiny::renderUI({
    if (!partition_data_available_env(env)) return(NULL)
    if (is.null(cluster_pick_d()) || !nzchar(cluster_pick_d())) return(NULL)
    all_clusters <- get("clust_partition", envir = env, inherits = FALSE)
    cid <- as.character(cluster_pick_d())
    n <- sum(as.character(all_clusters) == cid, na.rm = TRUE)
    shiny::h4(sprintf("Interactive Cluster %s (n=%d)", cid, n), class = "plot-title full")
  })
  
  # Cluster note (large clusters)
  output[[id$pc_cluster_note_ui]] <- shiny::renderUI({
    shiny::req(prepped_all(), cluster_pick_d())
    P <- prepped_all()
    cid <- as.character(cluster_pick_d())
    n_g <- sum(as.character(P$clusters) == cid, na.rm = TRUE)
    
    max_n <- if (exists("MAX_HEATMAP_GENES", inherits = TRUE)) MAX_HEATMAP_GENES else 10000
    if (n_g > max_n) {
      shiny::tags$div(
        shiny::tags$p(
          shiny::tags$strong("Note: "),
          paste0("Only the top ", format(max_n, big.mark = ","), " most variable features are displayed due to cluster size")
        )
      )
    }
  })
  
  # Per-cluster prepared inputs (P2) — extracted so the gene order can be
  # reused for the heatmap render AND exposed for downstream consumers
  # (e.g. zoom-based Multi-Gene Explorer). Returns NULL if not ready.
  pc_prepped_cluster <- shiny::reactive({
    if (is.null(prepped_all()) || is.null(cluster_pick_d())) return(NULL)
    P <- prepped_all()
    cid <- as.character(cluster_pick_d())
    idx <- which(as.character(P$clusters) == cid)
    if (length(idx) <= 1) return(NULL)

    mat_sub <- P$mat[idx, , drop = FALSE]
    rn_sub  <- P$rn[idx]

    max_n <- if (exists("MAX_HEATMAP_GENES", inherits = TRUE)) MAX_HEATMAP_GENES else 10000
    if (nrow(mat_sub) > max_n) {
      v <- apply(mat_sub, 1, stats::var, na.rm = TRUE)
      keep <- order(v, decreasing = TRUE)[seq_len(max_n)]
      mat_sub <- mat_sub[keep, , drop = FALSE]
      rn_sub  <- rn_sub[keep]
      ord <- order(match(rn_sub, P$rn))
      rn_sub  <- rn_sub[ord]
      mat_sub <- mat_sub[ord, , drop = FALSE]
    }

    list(
      mat        = mat_sub,
      rn         = rn_sub,
      cn         = unname(as.character(P$cn)),
      colorscale = P$colorscale,
      vmax       = P$vmax
    )
  })

  # Public reactive: ordered gene IDs as they appear on the partition heatmap's
  # y-axis for the currently selected cluster (or character(0)).
  pc_gene_order <- shiny::reactive({
    P2 <- pc_prepped_cluster()
    if (is.null(P2) || is.null(P2$rn)) return(character(0))
    as.character(P2$rn)
  })

  # Partition heatmap for selected cluster
  output[[id$pc_heatmap_cluster]] <- plotly::renderPlotly({
    shiny::validate(shiny::need(
      partition_data_available_env(env),
      "Clustering results are not available for this dataset."
    ))
    shiny::req(prepped_all(), cluster_pick_d())
    P2 <- pc_prepped_cluster()
    shiny::validate(shiny::need(!is.null(P2), "Selected cluster has no features or only a single feature."))

    tryCatch({
      p <- plot_partition_heatmap_plotly(P2, source_id = source$pc)
      plotly::config(p, displayModeBar = TRUE)
    }, error = function(e) {
      shiny::validate(shiny::need(FALSE, "Cluster heatmap could not be rendered."))
    })
  })
  
  # Cluster mean profile plot. Returns NULL when partition data is missing
  # so the downstream renderer can halt gracefully via req().
  cluster_means <- shiny::reactive({
    if (!exists("sample_meta", envir = env, inherits = FALSE)) return(NULL)
    P <- prepped_all()
    if (is.null(P) || is.null(P$mat) || is.null(P$clusters)) return(NULL)
    clv <- as.character(P$clusters)
    tryCatch({
      agg <- rowsum(P$mat, group = clv, reorder = FALSE) / as.vector(table(clv))
      list(agg = agg, cn = P$cn)
    }, error = function(e) NULL)
  })

  output[[id$pc_cluster_subgraph]] <- shiny::renderPlot({
    shiny::validate(shiny::need(
      partition_data_available_env(env),
      "Clustering results are not available for this dataset."
    ))
    shiny::validate(
      shiny::need(exists("sample_meta", envir = env, inherits = FALSE), "sample_meta not available.")
    )
    shiny::req(cluster_means(), cluster_pick_d(), prepped_all())

    P <- prepped_all()
    cid <- as.character(cluster_pick_d())
    n_g <- sum(as.character(P$clusters) == cid, na.rm = TRUE)

    cm <- cluster_means()
    shiny::validate(shiny::need(cid %in% rownames(cm$agg), "Cluster not found in precomputed means."))

    y  <- as.numeric(cm$agg[cid, ])
    sm <- get("sample_meta", envir = env)

    # Align sample_meta rows to the columns used in clustering (cm$cn). Required
    # for Metabolomics where sample_meta retains QC samples but de_expr_norm
    # does not — without this the data.frame() below errors out with mismatched
    # row counts. No-op for RNAseq/Proteomics. Factor resolution and level order
    # both read the ALIGNED table, so QC / unused levels never reach the plot.
    sm <- align_sample_meta_to_cols(sm, cm$cn)

    # X-axis + optional series, resolved generically from env$shape / env$color.
    # See .resolve_cluster_profile_factors() for the rule and the fallbacks.
    fac <- .resolve_cluster_profile_factors(env, sm)
    shiny::validate(shiny::need(
      !is.null(fac$x_col),
      "No experimental factor could be resolved from sample_meta."
    ))
    two_factor <- identical(fac$n_factors, 2L)

    df   <- data.frame(X = sm[[fac$x_col]], EXP = y)
    df$X <- factor(as.character(df$X),
                   levels = .cluster_profile_levels(sm[[fac$x_col]]))
    if (two_factor) {
      df$Series <- factor(as.character(sm[[fac$series_col]]),
                          levels = .cluster_profile_levels(sm[[fac$series_col]]))
    }

    # A sample with no value for a plotted factor carries no profile
    # information, so drop it rather than let NA become a phantom X position or
    # a phantom series. droplevels() then removes any level left with no
    # samples; both are no-ops when the metadata is complete, which is the case
    # for every current payload.
    keep <- !is.na(df$X)
    if (two_factor) keep <- keep & !is.na(df$Series)
    df <- droplevels(df[keep, , drop = FALSE])
    shiny::validate(shiny::need(
      nrow(df) > 0,
      "No samples with usable experimental-factor values for this cluster."
    ))

    if (two_factor) {
      # Same high-contrast palette make_pca_plotly() uses, so a given level
      # keeps one colour across the app. Recycled when there are more levels
      # than base colours; no level name is assumed.
      palette_base <- c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e",
                        "#9467bd", "#8c564b", "#e377c2", "#17becf",
                        "#aec7e8", "#ff9896", "#98df8a", "#ffbb78")
      lv  <- levels(df$Series)
      pal <- if (length(lv) <= length(palette_base)) {
        palette_base[seq_along(lv)]
      } else {
        grDevices::colorRampPalette(palette_base)(length(lv))
      }
      names(pal) <- lv

      # group = Series makes each summary layer aggregate within one X × series
      # cell, so means and standard errors are computed per combination. The one
      # dodge object is shared by all three layers so line, point and errorbar
      # stay aligned while the series remain visually distinguishable.
      pd <- ggplot2::position_dodge(width = 0.15)

      p <- ggplot2::ggplot(df, ggplot2::aes(X, EXP, colour = Series, group = Series)) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 2, position = pd) +
        ggplot2::stat_summary(fun = mean, geom = "point", size = 3.5, position = pd) +
        ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar",
                              width = 0.15, linewidth = 1.1, position = pd) +
        ggplot2::scale_colour_manual(values = pal, name = fac$series_col)
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes(X, EXP, group = 1)) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::stat_summary(fun = mean, geom = "line", color = "#F26B5B", linewidth = 2) +
        ggplot2::stat_summary(fun = mean, geom = "point", color = "#F26B5B", size = 3.5) +
        ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar",
                              color = "#F26B5B", width = 0.15, linewidth = 1.1)
    }

    p +
      ggplot2::labs(
        title = paste0("Cluster ", cid, "\n", n_g),
        x     = fac$x_col,
        y     = "Scaled expression (Z)"
      ) +
      ggplot2::theme(
        legend.position  = if (two_factor) "bottom" else "none",
        plot.title       = ggplot2::element_text(face = "bold", size = 18),
        axis.text.x      = ggplot2::element_text(size = 11, angle = 30, hjust = 1, vjust = 1),
        axis.text.y      = ggplot2::element_text(size = 11),
        axis.title       = ggplot2::element_text(size = 13, face = "bold")
      )
  })
  
  # RESULTS_DE latest + read header/data
  results_de_data <- shiny::reactive({
    obj_nm <- label$results_de_obj
    shiny::validate(
      shiny::need(exists(obj_nm, envir = env, inherits = FALSE),
                  paste0(obj_nm, " not available in env."))
    )
    as.data.frame(get(obj_nm, envir = env, inherits = FALSE), check.names = FALSE)
  })
  
  # Partition: feature IDs in the selected cluster (from the heatmap rows)
  partition_cluster_features <- shiny::reactive({
    shiny::req(prepped_all(), cluster_pick_d())
    P <- prepped_all()
    cid <- as.character(cluster_pick_d())
    
    feats <- P$rn[as.character(P$clusters) == cid]
    feats <- feats[!is.na(feats) & nzchar(feats)]
    unique(feats)
  })
  
  # Partition: filtered RESULTS_DE table by feature IDs
  results_de_filtered <- shiny::reactive({
    df <- results_de_data()
    feats <- partition_cluster_features()

    fcol <- .detect_feature_col(df, feats = feats)
    shiny::validate(shiny::need(!is.null(fcol), "Could not detect feature column in RESULTS_DE (gene/protein)."))

    df[df[[fcol]] %in% feats, , drop = FALSE]
  })

  
  output[[id$pc_results_title_ui]] <- shiny::renderUI({
    shiny::req(cluster_pick_d())
    n  <- tryCatch(nrow(results_de_filtered()), error = function(...) 0)
    cl <- cluster_pick_d()
    shiny::tags$h4(
      sprintf("Download result tables Cluster %s (n=%s)", cl, format(n, big.mark = ",")),
      class = "plot-title full",
      style = "margin-top:20px;margin-bottom:10px;"
    )
  })
  
  output[[id$dl_results_de_ui]] <- shiny::renderUI({
    shiny::req(cluster_pick_d())
    cl <- cluster_pick_d()
    n  <- tryCatch(nrow(results_de_filtered()), error = function(...) 0)
    if (n > 0) {
      shiny::downloadButton(
        outputId = id$dl_results_de_btn,
        label    = sprintf("Result tables Cluster %s", cl),
        class    = "btn btn-primary btn-sm"
      )
    } else {
      shiny::tags$p(sprintf("No rows found for Cluster %s in RESULTS_DE file.", cl),
                    style = "color:#888;margin-top:6px;")
    }
  })
  
  output[[id$dl_results_de_btn]] <- shiny::downloadHandler(
    filename = function() paste0("RESULTS_DE_cluster_", cluster_pick_d(), ".xlsx"),
    content = function(file) {
      df_flt <- results_de_filtered()
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Results")
      openxlsx::writeData(wb, "Results", x = df_flt, withFilter = TRUE)
      openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
    }
  )
  
  # Binary pattern titles + selector
  output[[id$bp_static_title_ui]] <- shiny::renderUI({
    if (!bp_data_available_env(env)) return(shiny::tags$div("Binary pattern data is not available.", style = "color:#888;"))
    shiny::req(input[[paste0(id$bp_pattern_ui, "_pattern")]])
  })

  output[[id$bp_interactive_title_ui]] <- shiny::renderUI({
    if (!bp_data_available_env(env)) return(NULL)
    shiny::req(input[[paste0(id$bp_pattern_ui, "_pattern")]])
    shiny::h4("Interactive Binary Pattern", class = "plot-title full")
  })
  
  output[[id$bp_pattern_ui]] <- shiny::renderUI({
    if (!bp_data_available_env(env)) return(NULL)
    patterns_local <- get("clust_patterns", envir = env)
    hm_list_local  <- get("clust_heatmaps_by_pattern", envir = env)
    ch <- bp_available_patterns(patterns_local, hm_list_local)
    if (!length(ch)) return(NULL)
    
    shiny::selectInput(
      inputId = paste0(id$bp_pattern_ui, "_pattern"),
      label   = "Choose pattern:",
      choices = ch,
      selected = ch[1]
    )
  })
  
  # Binary pattern plot outputs
  output[[id$bp_static]] <- shiny::renderPlot({
    if (!bp_data_available_env(env)) return(invisible(NULL))
    hm_list <- get("clust_heatmaps_by_pattern", envir = env)
    pat_id  <- input[[paste0(id$bp_pattern_ui, "_pattern")]]
    shiny::req(pat_id)
    hm <- hm_list[[pat_id]]
    shiny::req(!is.null(hm))
    bp_draw_static(hm)
  })
  
  output[[id$bp_interactive]] <- plotly::renderPlotly({
    if (!bp_data_available_env(env)) return(plotly::plot_ly(type = "scatter", mode = "markers"))
    hm_list <- get("clust_heatmaps_by_pattern", envir = env)
    pat_id  <- input[[paste0(id$bp_pattern_ui, "_pattern")]]
    shiny::req(pat_id)
    
    hm <- hm_list[[pat_id]]
    shiny::req(!is.null(hm), !is.null(hm$mat), is.matrix(hm$mat),
               nrow(hm$mat) > 0, ncol(hm$mat) > 0)
    
    p <- bp_build_plotly_safe(hm)
    p$x$source <- source$bp
    plotly::event_register(p, "plotly_click")
  })
  
  # Binary pattern filtered RESULTS_DE
  bp_pattern_genes <- shiny::reactive({
    pat_id <- input[[paste0(id$bp_pattern_ui, "_pattern")]]
    shiny::req(pat_id)
    hm_list <- get("clust_heatmaps_by_pattern", envir = env)
    hm <- hm_list[[pat_id]]
    shiny::req(!is.null(hm), !is.null(hm$mat), is.matrix(hm$mat))

    feats <- rownames(hm$mat)
    feats <- feats[!is.na(feats) & nzchar(feats)]
    unique(feats)
  })

  # Public reactive: gene IDs in the order they appear on the BP heatmap's
  # y-axis for the currently selected pattern (replicates the row reordering
  # logic in bp_build_plotly_safe). Returns character(0) when not ready.
  bp_gene_order <- shiny::reactive({
    if (!bp_data_available_env(env)) return(character(0))
    pat_id <- input[[paste0(id$bp_pattern_ui, "_pattern")]]
    if (is.null(pat_id) || !nzchar(pat_id)) return(character(0))
    hm_list <- get("clust_heatmaps_by_pattern", envir = env)
    hm <- hm_list[[pat_id]]
    if (is.null(hm) || is.null(hm$mat) || !is.matrix(hm$mat)) return(character(0))

    rn <- rownames(hm$mat)
    if (is.null(rn)) return(character(0))

    # Mirror bp_build_plotly_safe row reordering
    if (!is.null(hm$row_order)) {
      ro <- hm$row_order
      if (is.character(ro)) {
        ro <- ro[ro %in% rn]
        if (length(ro)) rn <- ro
      } else {
        ro <- ro[ro >= 1 & ro <= length(rn)]
        if (length(ro)) rn <- rn[ro]
      }
    } else if (!is.null(hm$tree_row) && !is.null(hm$tree_row$order)) {
      ro <- hm$tree_row$order
      ro <- ro[ro >= 1 & ro <= length(rn)]
      if (length(ro)) rn <- rn[ro]
    }
    as.character(rn)
  })
  
  results_de_bp_filtered <- shiny::reactive({
    pat_id <- input[[paste0(id$bp_pattern_ui, "_pattern")]]
    shiny::req(pat_id)
    df <- results_de_data()
    feats <- bp_pattern_genes()

    fcol <- .detect_feature_col(df, feats = feats)
    shiny::validate(shiny::need(!is.null(fcol), "Could not detect feature column in RESULTS_DE (gene/protein)."))
    df[df[[fcol]] %in% feats, , drop = FALSE]
  })
  
  output[[id$bp_results_title_ui]] <- shiny::renderUI({
    pat_id <- input[[paste0(id$bp_pattern_ui, "_pattern")]]
    shiny::req(pat_id)
    n <- tryCatch(nrow(results_de_bp_filtered()), error = function(...) 0)
    shiny::h4(
      sprintf("Download result tables Binary Pattern %s (n=%s)", pat_id, format(n, big.mark = ",")),
      class = "plot-title full",
      style = "margin-top:20px;margin-bottom:10px;"
    )
  })
  
  output[[id$dl_results_de_bp_ui]] <- shiny::renderUI({
    pat_id <- input[[paste0(id$bp_pattern_ui, "_pattern")]]
    shiny::req(pat_id)
    n <- tryCatch(nrow(results_de_bp_filtered()), error = function(...) 0)
    
    if (n > 0) {
      shiny::downloadButton(
        outputId = id$dl_results_de_bp_btn,
        label    = sprintf("Result tables Pattern %s", pat_id),
        class    = "btn btn-primary btn-sm"
      )
    } else {
      shiny::tags$p(sprintf("No features from Pattern %s were found in RESULTS_DE.", pat_id),
                    style = "color:#888;margin-top:6px;")
    }
  })
  
  output[[id$dl_results_de_bp_btn]] <- shiny::downloadHandler(
    filename = function() {
      pat_id <- input[[paste0(id$bp_pattern_ui, "_pattern")]]
      paste0("RESULTS_DE_binary_pattern_", pat_id, ".xlsx")
    },
    content = function(file) {
      df_flt <- results_de_bp_filtered()
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Results")
      openxlsx::writeData(wb, "Results", x = df_flt, withFilter = TRUE)
      openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
    }
  )

  # Public controller: gene-order reactives for the per-cluster (PC) and
  # per-pattern (BP) heatmaps. Each returns a character vector matching the
  # y-axis order of the currently displayed heatmap, or character(0) when not
  # ready. Used by wire_heatmap_zoom_to_multi_gene_explorer().
  invisible(list(
    pc_gene_order = pc_gene_order,
    bp_gene_order = bp_gene_order
  ))
}


## ---- Subtab: Enrichment ----
## ---- Subtab: GO Enrichmen ----
## ---- Subtab: KEGG Enrichment ----



## ---- Subtab: Gene Explorer ----

# Safe filenames for downloads (avoid '/', ':', etc.)
sanitize_filename <- function(x) {
  x <- as.character(x)
  x <- gsub("[/\\\\:*?\"<>|]+", "_", x)
  x <- gsub("\\s+", "_", x)
  trimws(x)
}

# Format small numeric values in scientific notation (below threshold)
format_scientific_threshold <- function(x, threshold = 0.01) {
  sapply(x, function(v) {
    if (is.na(v)) return(NA)
    if (abs(v) >= threshold) return(format(v, scientific = FALSE, trim = TRUE))
    format(v, scientific = TRUE, digits = 3)
  })
}

# DE table legend (shared)
de_table_legend_html <- function() {
  shiny::HTML("
    <span style='display:inline-block;width:14px;height:14px;background-color:#f5d49a;border:1px solid #aaa;margin-right:6px;'></span>
    Rows highlighted in <b>orange</b> indicate significant <b>up-regulation</b><br>
    <span style='display:inline-block;width:14px;height:14px;background-color:#cedef2;border:1px solid #aaa;margin-right:6px;'></span>
    Rows highlighted in <b>blue</b> indicate significant <b>down-regulation</b>
  ")
}

# §3: one concise, site-consistent instruction shown ABOVE a selectable gene table,
# explaining the Explorer hand-off (1 row -> single Explorer; 2+ -> Multi-Gene
# Explorer). SHARED by Keyword Search and the Enrichment Member-Genes / Core-Genes
# tables so the two read identically. `feature` is the row noun ("gene"/"protein"/…).
.explorer_instruction <- function(single_tab_label = "Gene Explorer",
                                  multi_tab_label  = "Multi-Gene Explorer",
                                  feature = "gene") {
  plural <- paste0(feature, "s")
  # §1: modestly more prominent than a plain muted hint (larger font, stronger
  # weight, the site accent colour, a bit more spacing) — but NOT a heading or an
  # alert. ONE shared helper, so Enrichment + Keyword Search change together.
  shiny::div(style = paste("font-size:14px; font-weight:600; color:#2c3e50;",
                           "margin:6px 0 10px;"),
    sprintf(paste0("Select one row to open the %s in %s, or select multiple rows ",
                   "to open the selected %s in %s."),
            feature, single_tab_label, plural, multi_tab_label))
}

# Shared DE cutoffs UI: a single significance cutoff input. Its "label" is a
# compact inline dropdown that switches between "padj cutoff" and
# "pvalue cutoff"; the numeric input itself is reused. Used by both Gene and
# Protein explorers (Gene / Protein / Metabolite / Multiomics tabs).
build_de_cutoffs_ui <- function(id_sig_mode, id_sig_cutoff, id_fc_cutoff,
                                sig0 = 0.05, fc0 = 1, sig_mode0 = "padj") {
  shiny::tags$div(
    style = "display:flex; gap:14px; align-items:flex-end; flex-wrap:wrap;",

    # Significance cutoff: dropdown-as-label + numeric input (reused across modes)
    shiny::tags$div(
      style = "width:130px;",
      shiny::tags$div(
        class = "sig-mode-dropdown",
        shiny::selectInput(
          inputId   = id_sig_mode,
          label     = NULL,
          choices   = c("padj cutoff" = "padj", "pvalue cutoff" = "pvalue"),
          selected  = sig_mode0,
          width     = "120px",
          selectize = FALSE
        )
      ),
      shiny::numericInput(id_sig_cutoff, label = NULL,
                          value = sig0, min = 0, max = 1, step = 0.01,
                          width = "100%")
    ),

    # linearFC cutoff — unchanged
    shiny::tags$div(
      style = "width:140px;",
      shiny::numericInput(id_fc_cutoff, label = "linearFC cutoff",
                          value = fc0, min = 0, max = 10, step = 0.1,
                          width = "100%")
    )
  )
}

# Render DE stats DT with row coloring.
# Direction (up/down) is always determined by linearFC sign and fc_thr.
# A single significance cutoff (`sig_thr`) is applied to either padj or pvalue
# depending on `sig_mode`.
render_de_stats_dt <- function(df,
                               sig_thr,
                               fc_thr,
                               sig_default = 0.05,
                               fc_default  = 1,
                               sig_mode    = c("padj", "pvalue"),
                               paginate    = FALSE,
                               page_length = 10L,
                               searching   = FALSE,
                               ordering    = FALSE,
                               selection   = "multiple",
                               colour_rows = TRUE) {

  sig_mode <- match.arg(sig_mode)

  if (is.null(sig_thr) || is.na(sig_thr)) sig_thr <- sig_default
  if (is.null(fc_thr)  || is.na(fc_thr))  fc_thr  <- fc_default

  sig_thr <- as.numeric(sig_thr)
  fc_thr  <- as.numeric(fc_thr)

  # Choose significance column based on mode; threshold is the same input
  sig_val <- if (sig_mode == "pvalue") df$pvalue else df$padj

  color_flag <- ifelse(
    !is.na(sig_val) & !is.na(df$linearFC) & (sig_val < sig_thr) & (df$linearFC >= fc_thr),
    "up",
    ifelse(
      !is.na(sig_val) & !is.na(df$linearFC) & (sig_val < sig_thr) & (df$linearFC <= -fc_thr),
      "down",
      "none"
    )
  )
  
  df_display <- df
  df_display$linearFC <- format_scientific_threshold(df$linearFC)
  df_display$pvalue   <- format_scientific_threshold(df$pvalue)
  df_display$padj     <- format_scientific_threshold(df$padj)
  
  df_display$color_flag <- color_flag
  n_cols <- ncol(df_display)
  
  tips <- list(
    linearFC = "Fold change is presented on a linear scale, where the sign indicates the direction of change",
    padj     = "FDR - Adjusted p-value"
  )
  tip_json <- jsonlite::toJSON(tips, auto_unbox = TRUE)
  
  dt <- DT::datatable(
    df_display,
    rownames = FALSE,
    class    = "cell-border",
    selection = selection,
    options  = list(
      # Multi Explorer default (paginate = FALSE): dom "t", no paging / search /
      # sort — UNCHANGED. Enrichment passes paginate = TRUE for a 20-row paged,
      # searchable, sortable table with the SAME formatting + Up/Down colouring.
      dom        = if (paginate) "ftip" else "t",
      paging     = paginate,
      pageLength = page_length,
      searching  = searching,
      ordering   = ordering,
      order      = list(),          # keep the provided row order until the user sorts
      info       = paginate,
      stripeClasses = c("", ""),
      columnDefs = list(
        list(targets = n_cols - 1, visible = FALSE, searchable = FALSE)
      ),
      initComplete = DT::JS(sprintf(
        "function(settings, json) {
           var tips = %s;
           var api = this.api();
           var header = $(api.table().header());

           header.find('th').each(function() {
             var th  = $(this);
             var txt = th.text().trim();

             if (tips[txt]) {
               th.empty();
               th.append(document.createTextNode(txt + ' '));

               var icon = $('<i class=\"fa-solid fa-circle-info dt-tip-info\"></i>');
               icon.attr('title', tips[txt]);
               icon.attr('data-toggle', 'tooltip');
               icon.attr('data-placement', 'top');

               th.append(icon);
             }
           });

           if ($.fn.tooltip) {
             header.find('.dt-tip-info[data-toggle=\"tooltip\"]').tooltip({container: 'body'});
           }
         }",
        tip_json
      ))
    )
  )
  # §1: significance-based Up/Down row colouring is OPT-OUT. The Multi/Gene
  # Explorer keeps it (colour_rows = TRUE default); the Enrichment selected-gene
  # tables pass colour_rows = FALSE (presentation-only — stats unchanged).
  if (isTRUE(colour_rows)) {
    dt <- dt %>% DT::formatStyle(
      "color_flag", target = "row",
      backgroundColor = DT::styleEqual(c("none", "up", "down"),
                                       c(NA, "#f5d49a", "#cedef2")))
  }
  dt
}

# Prepare selectize choices from rownames(matrix)
prepare_id_choices_from_matrix <- function(mat) {
  if (is.null(mat) || is.null(rownames(mat))) return(NULL)
  ids <- trimws(as.character(rownames(mat)))
  sort(unique(ids))
}

# Shared formatter for the log2-scale option label across all omics. Takes
# a method name (e.g. "VST", "pqn", "median"), lowercases it, and pastes
# "<method> normalized counts (log2 scale)". On any degenerate input
# (NULL, NA, empty string, "none"/"na"/"null", zero-length) falls back to
# the generic "normalized counts (log2 scale)" — no method prefix.
.format_log2_norm_label <- function(method) {
  m <- tolower(as.character(method))
  if (length(m) == 0L || is.na(m) || !nzchar(m) ||
      m %in% c("none", "na", "null")) {
    return("normalized counts (log2 scale)")
  }
  paste0(m, " normalized counts (log2 scale)")
}

# Build the RNAseq log2 dropdown / y-axis label from env$norm_method.
# Used by both the RNAseq Normalization Method dropdown (log2 option) and
# the y-axis label of the log2 Expression Profile plot, so they stay in sync.
get_norm_method_label_from_env <- function(env) {
  nm <- if (exists("norm_method", envir = env, inherits = FALSE)) {
    get("norm_method", envir = env, inherits = FALSE)
  } else NULL
  .format_log2_norm_label(nm)
}

# Build the Metabolomics log2 dropdown label from env$chosen_norm.
# Same formatting rules as the RNAseq helper, just sourced from a
# metabolomics-specific slot.
get_chosen_norm_label_from_env <- function(env) {
  nm <- if (exists("chosen_norm", envir = env, inherits = FALSE)) {
    get("chosen_norm", envir = env, inherits = FALSE)
  } else NULL
  .format_log2_norm_label(nm)
}

# Scale selector labels (RNAseq: linear vs log2). Linear label is fixed;
# log2 label is dynamic via get_norm_method_label_from_env().
get_scale_choices_from_env <- function(env) {
  vsd_label <- get_norm_method_label_from_env(env)
  setNames(
    object = c("norm", "log2"),
    nm     = c("normalized counts (linear scale)", vsd_label)
  )
}

# Build long df for a single expression matrix (features x samples) using sample_meta
get_explorer_df_from_mat <- function(feature_id, mat, color_by, sample_meta) {
  if (is.null(mat) || is.null(rownames(mat))) stop("Expression matrix has no rownames.")
  if (is.null(colnames(mat))) stop("Expression matrix has no colnames (samples).")
  
  idx <- which(rownames(mat) == feature_id)
  if (length(idx) == 0) stop(paste("No row found for:", feature_id))
  if (length(idx) > 1) idx <- idx[1]
  
  counts  <- as.numeric(mat[idx, ])
  samples <- colnames(mat)
  
  if (!all(samples %in% rownames(sample_meta))) {
    stop("Some samples from expression matrix not found in sample_meta rownames.")
  }
  sm <- sample_meta[samples, , drop = FALSE]
  
  if (!color_by %in% colnames(sm)) stop(paste("Column", color_by, "not found in sample_meta."))
  
  data.frame(
    Sample     = factor(samples, levels = samples),
    norm_count = counts,
    Effect     = factor(sm[[color_by]])
  )
}

# Wrapper for RNA: choose matrix by scale (linear vs log2)
get_explorer_df_from_mats <- function(feature_id,
                                      scale_sel = c("norm", "log2"),
                                      color_by,
                                      expr_raw,
                                      expr_norm = NULL,
                                      sample_meta) {
  scale_sel <- match.arg(scale_sel)
  
  mat <- expr_raw
  if (scale_sel == "log2") {
    if (is.null(expr_norm)) stop("expr_norm is NULL but scale_sel is 'log2'.")
    mat <- expr_norm
  }
  
  get_explorer_df_from_mat(
    feature_id  = feature_id,
    mat         = mat,
    color_by    = color_by,
    sample_meta = sample_meta
  )
}

# Build expression profile plot (shared)
make_explorer_plot <- function(df,
                               color_by,
                               scale_sel = c("norm", "log2"),
                               ylab_log2,
                               theme_variant = c("gene", "protein")) {
  scale_sel <- match.arg(scale_sel)
  theme_variant <- match.arg(theme_variant)
  
  ylab_txt <- if (scale_sel == "log2") ylab_log2 else "Normalized counts (linear scale)"
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = norm_count, fill = Effect)) +
    ggplot2::geom_col() +
    ggplot2::xlab("Sample") +
    ggplot2::ylab(ylab_txt) +
    ggplot2::labs(fill = color_by)
  
  if (theme_variant == "gene") {
    p + ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(size = 14),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        axis.text.y = ggplot2::element_text(size = 12),
        legend.text  = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 13)
      )
  } else {
    p + ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11),
        axis.text.y = ggplot2::element_text(size = 11),
        legend.position = "right",
        legend.text  = ggplot2::element_text(size = 11),
        legend.title = ggplot2::element_text(size = 12)
      )
  }
}

# ============================================================================
# HEATMAP-ZOOM MULTI-SELECT (revertible add-on)
# ----------------------------------------------------------------------------
# Self-contained helper. Listens to plotly_relayout on a categorical heatmap
# whose y-axis is gene IDs, computes the visible gene range from
# yaxis.range[0]/[1], renders a small status message + "Open visible genes in
# Multi-Gene Explorer" button. The button REPLACES the current Multi-Gene
# Explorer selection (capped at max_n) and switches tabs.
#
# To revert: remove this function definition, the matching server invocation,
# the matching pre-registration of "plotly_relayout-<source>" in
# session$userData$plotlyShinyEventIDs, and the matching UI block (all
# banner-marked).
# ============================================================================
wire_heatmap_zoom_to_multi_gene_explorer <- function(session, input, output,
                                                     source_id,
                                                     gene_order_reactive,
                                                     mge_controller,
                                                     tabs_id,
                                                     multi_tab_label = "Multi-Gene Explorer",
                                                     id_prefix,
                                                     empty_msg = "Zoom into the heatmap to detect visible genes.") {
  id_msg <- paste0(id_prefix, "_msg_ui")
  id_btn <- paste0(id_prefix, "_open_btn")

  visible_genes <- shiny::reactiveVal(character(0))

  # When the displayed gene set changes (e.g. user picks a different cluster
  # or pattern), clear any previous zoom-derived selection so the count UI
  # resets to the empty-state message.
  shiny::observeEvent(gene_order_reactive(), {
    visible_genes(character(0))
  }, ignoreInit = TRUE)

  shiny::observeEvent(
    suppressWarnings(plotly::event_data("plotly_relayout", source = source_id)),
    {
      ed <- plotly::event_data("plotly_relayout", source = source_id)
      shiny::req(!is.null(ed))

      gene_order <- gene_order_reactive()
      n <- length(gene_order)
      if (n == 0) { visible_genes(character(0)); return() }

      # Reset / autoscale on y -> clear
      yauto <- ed[["yaxis.autorange"]]
      if (!is.null(yauto) && (isTRUE(yauto) || identical(yauto, "true"))) {
        visible_genes(character(0)); return()
      }

      r0 <- ed[["yaxis.range[0]"]]
      r1 <- ed[["yaxis.range[1]"]]
      if (is.null(r0) || is.null(r1)) {
        # Not a y-axis change (e.g. x-only zoom/pan). Leave selection untouched.
        return()
      }

      r0 <- suppressWarnings(as.numeric(r0))
      r1 <- suppressWarnings(as.numeric(r1))
      if (any(is.na(c(r0, r1)))) return()

      lo <- max(0,     floor(min(r0, r1)))
      hi <- min(n - 1, ceiling(max(r0, r1)))
      if (hi < lo) { visible_genes(character(0)); return() }

      visible_genes(gene_order[(lo + 1):(hi + 1)])
    },
    ignoreInit = TRUE
  )

  output[[id_msg]] <- shiny::renderUI({
    g <- visible_genes()
    if (length(g) == 0) {
      return(shiny::tags$div(
        style = "color:#888; font-size:13px;",
        empty_msg
      ))
    }
    max_n <- mge_controller$max_n
    suffix <- if (length(g) > max_n)
      sprintf(" (only the first %d will be opened)", max_n) else ""
    shiny::tags$div(
      style = "font-size:13px; color:#2c3e50; font-weight:600;",
      sprintf("Visible genes in zoom: %d%s", length(g), suffix)
    )
  })

  shiny::observeEvent(input[[id_btn]], {
    g <- visible_genes()
    if (length(g) == 0) {
      shiny::showNotification(
        "Zoom into the heatmap first to select visible genes.",
        type = "message", duration = 4
      )
      return()
    }
    max_n <- mge_controller$max_n
    if (length(g) > max_n) {
      shiny::showNotification(
        sprintf("The current zoom contains %d genes. Only the first %d will be opened.",
                length(g), max_n),
        type = "warning", duration = 5
      )
      g <- g[seq_len(max_n)]
    }
    mge_controller$set(g)
    bslib::nav_select(id = tabs_id, selected = multi_tab_label)
  }, ignoreInit = TRUE)
}
# ============================================================================
# END HEATMAP-ZOOM MULTI-SELECT (revertible add-on)
# ============================================================================

# ---- Multi-Gene Explorer (RNAseq only) -------------------------------------

# Build a long-format data frame for multiple features (genes / proteins /
# metabolites). Each row: Gene, Sample, norm_count, Effect.
# (The "Gene" column name is kept internally regardless of the omics type;
# downstream rendering re-labels it for display.)
#
# `mode` chooses the matrix-selection strategy:
#   "rnaseq"     -> dual-matrix logic via get_explorer_df_from_mats(scale_sel)
#                   where scale_sel is one of "norm" / "log2".
#   "proteomics" -> single-matrix logic: scale_sel="raw" picks expr_raw,
#                   anything else (e.g. "imp1") picks expr_norm.

# ---------------------------------------------------------------------------
# Shared Expression-Data table helpers used by both single-feature and
# multi-feature Explorers. Pivots a long expression data frame into a wide
# 1-feature-per-row table and renders it with dynamic decimals (0 if all
# values are whole numbers, 2 otherwise), no search/page/order/info chrome.
# ---------------------------------------------------------------------------

# Pivot a long expression data frame into wide format.
# `long_df` must have columns `Sample` and `norm_count`. If it also has a
# `Gene` column (multi-feature case), every distinct value becomes a row in
# the order it appears (or its factor levels if Gene is a factor).
# Otherwise (single-feature case) `feature_id` is used as the single row id.
build_explorer_wide_df <- function(long_df,
                                   feature_id        = NULL,
                                   feature_col_name  = "Gene") {
  if (is.null(long_df) || nrow(long_df) == 0L) return(NULL)

  if (!"Gene" %in% colnames(long_df)) {
    if (is.null(feature_id) || !nzchar(as.character(feature_id))) return(NULL)
    long_df$Gene <- as.character(feature_id)
  }

  feat_levels   <- if (is.factor(long_df$Gene))   levels(long_df$Gene)
                   else unique(as.character(long_df$Gene))
  sample_levels <- if (is.factor(long_df$Sample)) levels(long_df$Sample)
                   else unique(as.character(long_df$Sample))

  wide <- matrix(NA_real_,
                 nrow = length(feat_levels),
                 ncol = length(sample_levels),
                 dimnames = list(feat_levels, sample_levels))
  g_chr <- as.character(long_df$Gene)
  s_chr <- as.character(long_df$Sample)
  for (i in seq_len(nrow(long_df))) {
    wide[g_chr[i], s_chr[i]] <- long_df$norm_count[i]
  }

  out <- cbind(
    setNames(data.frame(rownames(wide), stringsAsFactors = FALSE),
             feature_col_name),
    as.data.frame(wide, stringsAsFactors = FALSE, check.names = FALSE)
  )
  rownames(out) <- NULL
  out
}

# Render a wide expression-data table (output of build_explorer_wide_df) as a
# compact DT::datatable. Used by every Explorer (single + multi) so styling is
# byte-identical across them.
#
# Optional clickable feature names: when both `click_input_id` and
# `clickable_col` are supplied AND `clickable_col` is a column in `wide_df`,
# each value in that column is wrapped in an inline `<a>` whose onclick pushes
# the raw feature ID into the Shiny input named `click_input_id` (with
# {priority:'event'} so re-clicks always fire). All other cells stay plain
# text. Used by the Multi-Feature Explorers to navigate to the single
# Explorer; single-feature Explorers and MultiOmics omit these args, so their
# tables continue to render verbatim text.
render_explorer_wide_dt <- function(wide_df, feature_col_name = "Gene",
                                    click_input_id = NULL,
                                    clickable_col  = NULL) {
  shiny::validate(shiny::need(!is.null(wide_df) && nrow(wide_df) > 0, ""))
  sample_cols <- setdiff(colnames(wide_df), feature_col_name)

  sample_vals <- unlist(wide_df[, sample_cols, drop = FALSE], use.names = FALSE)
  sample_vals <- sample_vals[!is.na(sample_vals)]
  digits_to_use <- if (length(sample_vals) == 0L) {
    2L
  } else if (all(abs(sample_vals - round(sample_vals)) < 1e-8)) {
    0L
  } else {
    2L
  }

  # Default: every column is HTML-escaped by DT (i.e. rendered as plain text).
  escape_arg <- TRUE

  if (!is.null(click_input_id) && length(click_input_id) == 1L &&
      nzchar(as.character(click_input_id)) &&
      !is.null(clickable_col) && length(clickable_col) == 1L &&
      clickable_col %in% colnames(wide_df)) {
    raw_vals  <- as.character(wide_df[[clickable_col]])
    safe_attr <- htmltools::htmlEscape(raw_vals, attribute = TRUE)
    safe_text <- htmltools::htmlEscape(raw_vals)
    # Anchor pushes the raw feature ID (via data-feature-id) into the
    # Shiny input. Inline style keeps the visual consistent with Bootstrap
    # link styling without requiring a separate CSS rule.
    wide_df[[clickable_col]] <- sprintf(
      '<a href="#" class="mge-feature-link" data-feature-id="%s" onclick="Shiny.setInputValue(\'%s\', this.dataset.featureId, {priority:\'event\'}); return false;" style="color:#2c7be5; text-decoration:underline; cursor:pointer;">%s</a>',
      safe_attr, click_input_id, safe_text
    )
    # Escape every column EXCEPT the clickable one (1-based indices).
    escape_arg <- which(colnames(wide_df) != clickable_col)
  }

  DT::datatable(
    wide_df,
    rownames = FALSE,
    escape   = escape_arg,
    class    = "cell-border stripe compact",
    options  = list(
      dom          = "t",
      paging       = FALSE,
      searching    = FALSE,
      ordering     = FALSE,
      info         = FALSE,
      lengthChange = FALSE,
      scrollX      = TRUE
    )
  ) |>
    DT::formatRound(columns = sample_cols, digits = digits_to_use)
}

build_multi_gene_df_from_env <- function(env, gene_ids,
                                          scale_sel,
                                          color_by,
                                          expression_data_mode = c("rnaseq", "single_matrix", "proteomics")) {
  expression_data_mode <- match.arg(expression_data_mode)
  # Legacy alias: "proteomics" used to mean "single_matrix" before the rename.
  if (expression_data_mode == "proteomics") expression_data_mode <- "single_matrix"
  if (length(gene_ids) == 0) return(NULL)

  if (!exists("expr_raw",    envir = env, inherits = FALSE)) return(NULL)
  if (!exists("sample_meta", envir = env, inherits = FALSE)) return(NULL)

  expr_raw    <- get("expr_raw",    envir = env, inherits = FALSE)
  sample_meta <- get("sample_meta", envir = env, inherits = FALSE)
  expr_norm   <- if (exists("expr_norm", envir = env, inherits = FALSE)) {
    get("expr_norm", envir = env, inherits = FALSE)
  } else NULL

  one_df <- function(g) {
    tryCatch({
      if (expression_data_mode == "rnaseq") {
        # Validate scale_sel for RNAseq path
        sc <- if (is.null(scale_sel)) "norm" else scale_sel
        if (!sc %in% c("norm", "log2")) sc <- "norm"
        get_explorer_df_from_mats(
          feature_id  = g,
          scale_sel   = sc,
          color_by    = color_by,
          expr_raw    = expr_raw,
          expr_norm   = expr_norm,
          sample_meta = sample_meta
        )
      } else {
        # proteomics-style: scale_sel is "raw" or "imp1"
        mat <- if (identical(scale_sel, "raw")) expr_raw else expr_norm
        if (is.null(mat)) return(NULL)
        get_explorer_df_from_mat(
          feature_id  = g,
          mat         = mat,
          color_by    = color_by,
          sample_meta = sample_meta
        )
      }
    }, error = function(e) NULL)
  }

  parts <- lapply(gene_ids, function(g) {
    df <- one_df(g)
    if (is.null(df)) return(NULL)
    df$Gene <- g
    df
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0) return(NULL)

  out <- do.call(rbind, parts)
  out$Gene <- factor(out$Gene, levels = gene_ids[gene_ids %in% out$Gene])
  out
}

# Faceted multi-gene expression plot (one panel per gene, small multiples).
make_multi_gene_facet_plot <- function(df, color_by,
                                       scale_sel = c("norm", "log2"),
                                       ylab_log2) {
  scale_sel <- match.arg(scale_sel)
  ylab_txt <- if (scale_sel == "log2") ylab_log2 else "Normalized counts (linear scale)"

  ngenes <- length(unique(df$Gene))
  ncol_facet <- if (ngenes <= 2) ngenes else if (ngenes <= 6) 3 else 4

  ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = norm_count, fill = Effect)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ Gene, ncol = ncol_facet, scales = "free_y") +
    ggplot2::xlab("Sample") +
    ggplot2::ylab(ylab_txt) +
    ggplot2::labs(fill = color_by) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.text       = ggplot2::element_text(face = "bold", size = 11),
      strip.background = ggplot2::element_rect(fill = "#EDF0F1", colour = NA),
      axis.text.x      = ggplot2::element_text(angle = 90, vjust = 0.5,
                                               hjust = 1, size = 9),
      axis.text.y      = ggplot2::element_text(size = 9),
      legend.position  = "bottom",
      panel.spacing    = ggplot2::unit(0.6, "lines")
    )
}

# Wire a Multi-Feature Explorer (RNAseq / Proteomics / Metabolomics).
# Generalizes the original RNAseq-only Multi-Gene Explorer:
#   - `feature_label` controls singular/plural wording in the UI
#       (e.g. "gene" -> "Selected genes", "protein" -> "Selected proteins")
#   - `feature_col_name` controls the column header in the wide expression
#     table and the DE statistics table ("Gene" / "Protein" / "Metabolite")
#   - `mode` ("rnaseq" | "proteomics") switches the matrix-selection strategy
#     used to build the long expression data frame
#   - `de_stats_fn(stats, contrasts, id)` builds DE stats for one feature
#     (RNA: get_gene_de_stats_for_table; prot/metab: get_protein_de_stats_for_table)
#   - `scale_choices_fn(env)`, `default_scale`, and `ylab_log2_fn(env, scale_sel)`
#     control the Normalization Method selector and y-axis label
wire_multi_gene_explorer <- function(session, input, output, env,
                                     id_prefix = "mge", max_n = 100,
                                     page_size = 20L,
                                     feature_label        = "gene",
                                     feature_label_plural = NULL,
                                     feature_col_name     = "Gene",
                                     # expression_data_mode selects the shape of the
                                     # underlying expression matrices and the scale-choice
                                     # vocabulary. NOT an omics label.
                                     #   "rnaseq"       — two matrices (expr_raw + expr_norm),
                                     #                    user picks raw / norm / log2.
                                     #   "single_matrix"— one matrix per choice (expr_raw=raw,
                                     #                    expr_norm=imp1). Used by Proteomics
                                     #                    AND Metabolomics, which both ship
                                     #                    log2-already data with optional imputation.
                                     # (Legacy alias "proteomics" still accepted for backwards
                                     # compatibility with existing call sites.)
                                     expression_data_mode = c("rnaseq", "single_matrix", "proteomics"),
                                     de_stats_fn          = NULL,
                                     scale_choices_fn     = NULL,
                                     default_scale        = NULL,
                                     ylab_log2_fn         = NULL,
                                     theme_variant        = c("gene", "protein"),
                                     plot_scale_arg_fn    = NULL,
                                     # Force a fixed scale value instead of
                                     # reading input[[id_scale]]. When set
                                     # (e.g. "imp1"), the Normalization
                                     # selectInput is NOT updated and every
                                     # `input[[id_scale]]` read inside this
                                     # helper falls back to fixed_scale.
                                     # Used by Metabolomics to lock the
                                     # data path to the normalized matrix
                                     # without rendering a dropdown.
                                     fixed_scale          = NULL,
                                     # Single-Explorer navigation (clickable
                                     # feature names in the Expression Data
                                     # table). When all three are non-NULL we
                                     # wire the same modal flow used by
                                     # volcano/MA/clustering single-feature
                                     # picks; when any is NULL we skip wiring
                                     # and the table renders plain text.
                                     single_tab_label     = NULL,
                                     tabs_id              = NULL,
                                     single_selectize_id  = NULL) {
  # ============================================================
  # Internal organization
  # ------------------------------------------------------------
  # This function is large because it has to orchestrate many
  # independent UI sections (selection, sort, pagination, paste,
  # plot, data, DE stats, heatmap, single-explorer navigation)
  # that all share the same selection state. To keep navigation
  # tractable, the body is split into:
  #   1) Parameter normalization + ID strings + simple constants.
  #   2) Reactive state + closure controllers + shared reactives
  #      held at OUTER scope (so all helpers below see them).
  #   3) `.wire_mfe_*` inner helper closures, one per section.
  #   4) A short orchestration block that calls each helper in
  #      the original observer-registration order.
  #   5) The public controller list returned to the caller.
  # The helpers are pure side-effect closures: each registers a
  # group of `output[[…]]` + `shiny::observeEvent(…)` for a single
  # logical section. They take NO arguments — they read state via
  # lexical scope from the enclosing function.
  # ============================================================

  expression_data_mode <- match.arg(expression_data_mode)
  # Map the legacy "proteomics" value to the new "single_matrix" name so the
  # downstream branches don't need to know about both spellings.
  if (expression_data_mode == "proteomics") expression_data_mode <- "single_matrix"
  theme_variant <- match.arg(theme_variant)
  if (is.null(feature_label_plural)) feature_label_plural <- paste0(feature_label, "s")
  if (is.null(de_stats_fn))      de_stats_fn      <- get_gene_de_stats_for_table
  if (is.null(scale_choices_fn)) scale_choices_fn <- get_scale_choices_from_env
  if (is.null(default_scale))    default_scale    <- if (expression_data_mode == "rnaseq") "norm" else "raw"
  if (is.null(ylab_log2_fn))     ylab_log2_fn     <- function(env, scale_sel) get_norm_method_label_from_env(env)
  # plot_scale_arg_fn returns "norm" or "log2" — controls whether
  # make_explorer_plot uses the linear or log2 y-axis label.
  if (is.null(plot_scale_arg_fn)) {
    plot_scale_arg_fn <- if (expression_data_mode == "rnaseq") {
      function(scale_sel) scale_sel
    } else {
      function(scale_sel) "log2"  # protein/metab data is already log2
    }
  }

  id_input    <- paste0(id_prefix, "_genes_input")
  id_add      <- paste0(id_prefix, "_add_btn")
  id_clear    <- paste0(id_prefix, "_clear_btn")
  id_summary  <- paste0(id_prefix, "_summary_ui")
  id_count    <- paste0(id_prefix, "_count_ui")
  id_plot     <- paste0(id_prefix, "_plot")
  id_dl_plot  <- paste0("download_", id_prefix, "_plot")
  id_scale    <- paste0(id_prefix, "_scale")
  id_color_by <- paste0(id_prefix, "_color_by")
  id_table        <- paste0(id_prefix, "_data_table")
  id_dl_table     <- paste0("download_", id_prefix, "_table")
  id_remove_gene  <- paste0(id_prefix, "_remove_gene")
  # Clickable feature-name → single Explorer navigation. id_explorer_click is
  # the Shiny input set from the table cell's onclick; id_explorer_view is
  # the "View in <Explorer>" button rendered inside the confirmation modal.
  id_explorer_click <- paste0(id_prefix, "_explorer_click")
  id_explorer_view  <- paste0(id_prefix, "_explorer_view_btn")
  # DE statistics: single DT (id_de_table) showing one contrast at a time.
  # The contrast picker (id_de_contrast_pick) is rendered into id_de_contrast_ui
  # only when there is more than one contrast in env$contrasts; otherwise the
  # subtitle alone is shown. The cutoff inputs and legend are unchanged.
  id_de_table     <- paste0(id_prefix, "_de_table")
  id_de_contrast_ui   <- paste0(id_prefix, "_de_contrast_ui")
  id_de_contrast_pick <- paste0(id_prefix, "_de_contrast_pick")
  id_de_legend_ui <- paste0(id_prefix, "_de_legend_ui")
  id_de_cutoffs_ui<- paste0(id_prefix, "_de_cutoffs_ui")
  id_de_padj_cut  <- paste0(id_prefix, "_de_padj_cutoff")
  id_de_fc_cut    <- paste0(id_prefix, "_de_fc_cutoff")
  id_de_sig_mode  <- paste0(id_prefix, "_de_sig_mode")
  id_has_sel      <- paste0(id_prefix, "_has_selection")
  id_sel_hm       <- paste0(id_prefix, "_selected_heatmap")
  id_sel_hm_title <- paste0(id_prefix, "_selected_heatmap_title")
  id_sel_hm_msg   <- paste0(id_prefix, "_selected_heatmap_msg")
  # Global pagination UI: text + Previous/Next buttons. Applies to Expression
  # Profile, Expression Data, and DE Statistics (the heatmap and the selected-
  # tag list show ALL selected features regardless of page).
  id_pagination_ui <- paste0(id_prefix, "_pagination_ui")
  id_page_prev     <- paste0(id_prefix, "_page_prev")
  id_page_next     <- paste0(id_prefix, "_page_next")
  # Paste-list workflow: small "Paste feature list" button in the manual
  # selection row opens a modal dialog. The modal contains the textarea
  # (id_paste_input) and the Load button (id_paste_btn). On a successful
  # load the modal is dismissed via removeModal(); on a parse/validation
  # error the modal stays open so the user can fix the input.
  id_paste_open  <- paste0(id_prefix, "_paste_open")
  id_paste_input <- paste0(id_prefix, "_paste_input")
  id_paste_btn   <- paste0(id_prefix, "_paste_btn")
  # Feature-ordering controls: a "Sort by" dropdown plus a conditional
  # contrast selector that is only rendered when sort mode is "significance"
  # AND the env has more than one contrast.
  id_sort_by             <- paste0(id_prefix, "_sort_by")
  id_sort_contrast_ui    <- paste0(id_prefix, "_sort_contrast_ui")
  id_sort_contrast_pick  <- paste0(id_prefix, "_sort_contrast_pick")

  page_size <- as.integer(page_size)
  if (is.na(page_size) || page_size < 1L) page_size <- 20L

  selected_genes <- shiny::reactiveVal(character(0))
  current_page   <- shiny::reactiveVal(1L)

  # Auto-clamp the page when the selection shrinks below it (e.g. user removes
  # a feature that was on the last page and that page is now empty).
  shiny::observe({
    n <- length(selected_genes())
    n_pages <- if (n == 0L) 1L else as.integer(ceiling(n / page_size))
    if (current_page() > n_pages) current_page(n_pages)
    if (current_page() < 1L)      current_page(1L)
  })

  # Capture contrast names at wire time — used by both the DE-stats picker and
  # the significance-ordering mode.
  contrasts_at_wire <- if (exists("contrasts", envir = env, inherits = FALSE))
    get("contrasts", envir = env, inherits = FALSE) else NULL
  contrast_names_all <- if (is.data.frame(contrasts_at_wire) &&
                            "Contrast_name" %in% colnames(contrasts_at_wire)) {
    as.character(contrasts_at_wire$Contrast_name)
  } else character(0)

  # Ordering layer: selected_genes() -> ordered_selected_genes() -> paged.
  # Only display order changes; the underlying selection is untouched.
  # Modes:
  #   - "selection"    : insertion order (default)
  #   - "alpha"        : alphabetical A -> Z
  #   - "significance" : pvalue   ascending in the currently picked contrast;
  #                      features with missing pvalue go LAST.
  #   - "linearfc"     : linearFC ascending in the currently picked contrast;
  #                      features with missing linearFC go LAST.
  # The heatmap continues to read selected_genes() directly, so its row
  # order is unaffected by sorting; it is driven by a fresh hierarchical
  # clustering computed on the selected-feature Z-score matrix from
  # env$expr_norm (see prepare_selected_features_heatmap_inputs).
  ordered_selected_genes <- shiny::reactive({
    sel <- selected_genes()
    if (length(sel) == 0L) return(character(0))

    mode <- input[[id_sort_by]]
    if (is.null(mode) || !nzchar(mode)) mode <- "selection"

    if (mode == "alpha") {
      return(sel[order(sel, method = "radix")])
    }

    # "significance" and "linearfc" share the same contrast-resolution and
    # per-feature DE-stats lookup; only the column they sort on differs.
    if (mode %in% c("significance", "linearfc")) {
      # If there are no contrasts, or de_stats is missing, fall back to
      # selection order silently.
      if (length(contrast_names_all) == 0L) return(sel)
      if (!exists("de_stats", envir = env, inherits = FALSE)) return(sel)

      # Pick the contrast: explicit dropdown value if set & valid, else the
      # first contrast in env$contrasts (also used when only 1 contrast).
      cur <- input[[id_sort_contrast_pick]]
      if (is.null(cur) || !nzchar(as.character(cur)) ||
          !cur %in% contrast_names_all) {
        cur <- contrast_names_all[1L]
      }

      de_stats       <- get("de_stats",  envir = env, inherits = FALSE)
      contrasts_data <- get("contrasts", envir = env, inherits = FALSE)
      if (is.null(de_stats) || is.null(contrasts_data)) return(sel)

      sort_field <- if (mode == "linearfc") "linearFC" else "pvalue"

      # Reuse the same per-feature DE helper the table uses; we just look at
      # one column for the chosen contrast. NA when the helper returns NULL
      # or the contrast row is absent for that feature -> sorts last.
      vals <- vapply(sel, function(g) {
        one <- tryCatch(
          de_stats_fn(de_stats, contrasts_data, g),
          error = function(e) NULL
        )
        if (is.null(one) || nrow(one) == 0L) return(NA_real_)
        row <- one[as.character(one$Contrast) == cur, , drop = FALSE]
        if (nrow(row) == 0L) return(NA_real_)
        v <- row[[sort_field]]
        if (is.null(v)) NA_real_ else as.numeric(v[1L])
      }, numeric(1L), USE.NAMES = FALSE)

      ord <- order(vals, na.last = TRUE)
      return(sel[ord])
    }

    # Default mode: preserve insertion order.
    sel
  })

  # The feature subset visible on the current page. Reads the ORDERED list so
  # pagination + Expression Profile + Expression Data + DE Statistics all
  # respect the chosen sort order.
  paged_selected_genes <- shiny::reactive({
    sel <- ordered_selected_genes()
    if (length(sel) == 0L) return(character(0))
    pg <- current_page()
    lo <- (pg - 1L) * page_size + 1L
    hi <- min(length(sel), pg * page_size)
    if (lo > length(sel)) return(character(0))
    sel[lo:hi]
  })

  # Populate the selectize input choices server-side
  gene_choices <- NULL
  if (exists("expr_raw", envir = env, inherits = FALSE)) {
    gene_choices <- prepare_id_choices_from_matrix(get("expr_raw", envir = env, inherits = FALSE))
  } else if (exists("expr_norm", envir = env, inherits = FALSE)) {
    gene_choices <- prepare_id_choices_from_matrix(get("expr_norm", envir = env, inherits = FALSE))
  }

  session$onFlushed(function() {
    shiny::updateSelectizeInput(session, id_input,
                                choices = gene_choices,
                                selected = character(0),
                                server = TRUE)
  }, once = TRUE)

  # Populate scale + color-by inputs
  shiny::observe({
    color_col <- if (exists("color", envir = env, inherits = FALSE))
      get("color", envir = env, inherits = FALSE) else NULL
    shape_col <- if (exists("shape", envir = env, inherits = FALSE))
      get("shape", envir = env, inherits = FALSE) else NULL

    choices_vec <- unique(stats::na.omit(c(color_col, shape_col)))
    if (length(choices_vec) > 0) {
      shiny::updateSelectInput(session, id_color_by,
                               choices = choices_vec,
                               selected = color_col)
    }
    # When the Normalization control is hidden (fixed_scale set), skip the
    # update — there is no widget to populate and the data path uses
    # fixed_scale directly.
    if (is.null(fixed_scale)) {
      shiny::updateSelectInput(session, id_scale,
                               choices  = scale_choices_fn(env),
                               selected = default_scale)
    }
  })

  # Empty-state signal for the UI: "yes" when there is at least one selected
  # feature, "no" otherwise. The Multi-* Explorer UI uses this via a
  # conditionalPanel(condition = "output.<prefix>_has_selection == 'yes'") to
  # hide the Expression Profile / data / DE sections until something is
  # selected.
  output[[id_has_sel]] <- shiny::reactive({
    if (length(selected_genes()) > 0) "yes" else "no"
  })
  shiny::outputOptions(output, id_has_sel, suspendWhenHidden = FALSE)

  # Public helper: add genes to the current selection, dedupe, enforce max_n.
  # Returns invisibly the new selection length.
  add_genes <- function(ids, notify = TRUE) {
    if (length(ids) == 0) return(invisible(length(selected_genes())))
    ids <- trimws(as.character(ids))
    ids <- unique(ids[nzchar(ids) & !is.na(ids)])
    if (!is.null(gene_choices)) ids <- ids[ids %in% gene_choices]
    if (length(ids) == 0) return(invisible(length(selected_genes())))

    cur <- selected_genes()
    combined <- unique(c(cur, ids))
    if (length(combined) > max_n) {
      combined <- combined[seq_len(max_n)]
      if (notify) {
        shiny::showNotification(
          sprintf("Selection limit reached: only the first %d genes are kept.", max_n),
          type = "warning", duration = 5
        )
      }
    }
    selected_genes(combined)
    invisible(length(combined))
  }

  # set_genes() and clear_genes() are bulk replacements (volcano/MA box-select,
  # clustering zoom, "Clear" button). Reset the page so the user always lands
  # on page 1 after such a bulk change.
  set_genes   <- function(ids) {
    selected_genes(character(0))
    current_page(1L)
    add_genes(ids)
  }
  clear_genes <- function() {
    selected_genes(character(0))
    current_page(1L)
  }

  # Public helper: remove a single gene from the current selection.
  # Does not reset the page — the clamp observer above handles cases where
  # the current page becomes invalid.
  remove_gene <- function(id) {
    if (is.null(id) || !nzchar(trimws(as.character(id)))) return(invisible())
    cur <- selected_genes()
    selected_genes(cur[cur != as.character(id)])
    invisible()
  }

  # Pagination buttons
  shiny::observeEvent(input[[id_page_prev]], {
    pg <- current_page()
    if (pg > 1L) current_page(pg - 1L)
  }, ignoreInit = TRUE)

  shiny::observeEvent(input[[id_page_next]], {
    pg <- current_page()
    n  <- length(selected_genes())
    n_pages <- if (n == 0L) 1L else as.integer(ceiling(n / page_size))
    if (pg < n_pages) current_page(pg + 1L)
  }, ignoreInit = TRUE)

  # Any sort-related change reshuffles the feature order, so a "page 4 of 73"
  # view becomes meaningless after the switch. Snap back to page 1 so the user
  # always sees the start of the freshly-ordered list. ignoreInit = TRUE keeps
  # the initial page-1 default unchanged on app load.
  shiny::observeEvent(input[[id_sort_by]], {
    current_page(1L)
  }, ignoreInit = TRUE)
  shiny::observeEvent(input[[id_sort_contrast_pick]], {
    current_page(1L)
  }, ignoreInit = TRUE)

  # Sort-contrast -> DE-contrast synchronization.
  # When the user picks a contrast in the Sort-by area (only meaningful for
  # the "significance" / "linearfc" modes), mirror that choice into the DE
  # Statistics contrast dropdown so the two are aligned by default. The DE
  # picker remains independently editable afterward — this is a one-way
  # push (sort -> DE), never the reverse, so no reactive loop is possible.
  # updateSelectInput is a no-op when the DE picker isn't currently rendered
  # (e.g. single-contrast envs), so the call is safe to make unconditionally
  # once we've already verified the sort mode requires a contrast.
  shiny::observeEvent(input[[id_sort_contrast_pick]], {
    mode <- input[[id_sort_by]]
    if (is.null(mode) || !mode %in% c("significance", "linearfc")) return()
    new_val <- input[[id_sort_contrast_pick]]
    if (is.null(new_val) || !nzchar(as.character(new_val))) return()
    if (length(contrast_names_all) <= 1L) return()
    if (!new_val %in% contrast_names_all) return()
    if (isTRUE(input[[id_de_contrast_pick]] == new_val)) return()
    shiny::updateSelectInput(session, id_de_contrast_pick, selected = new_val)
  }, ignoreInit = TRUE)

  # Click on a gene tag's × triggers Shiny.setInputValue with the gene ID.
  shiny::observeEvent(input[[id_remove_gene]], {
    remove_gene(input[[id_remove_gene]])
  }, ignoreInit = TRUE)

  # Add genes button: union with current, dedupe, enforce max_n
  shiny::observeEvent(input[[id_add]], {
    new_ids <- input[[id_input]]
    if (is.null(new_ids) || length(new_ids) == 0) {
      shiny::showNotification("No genes selected to add.",
                              type = "message", duration = 3)
      return()
    }
    add_genes(new_ids)

    shiny::updateSelectizeInput(session, id_input,
                                choices  = gene_choices,
                                selected = character(0),
                                server   = TRUE)
  }, ignoreInit = TRUE)

  # Clear button: empty the selection
  shiny::observeEvent(input[[id_clear]], {
    clear_genes()
  }, ignoreInit = TRUE)

  # ---------- Paste-list workflow ----------
  # Opens a modal containing the textarea + Load button. The Load handler
  # parses, validates against gene_choices, caps at max_n, calls set_genes()
  # (which replaces the selection and resets pagination to page 1), and emits
  # the spec'd toasts. Defined as a closure for navigability — see the
  # "Internal organization" block at the top of this function.
  .wire_mfe_paste_list <- function() {
    shiny::observeEvent(input[[id_paste_open]], {
      shiny::showModal(shiny::modalDialog(
        title = paste0("Paste ", feature_label_plural),
        size  = "m",
        easyClose = TRUE,
        shiny::textAreaInput(
          inputId     = id_paste_input,
          label       = NULL,
          placeholder = paste0("Paste one ", feature_label, " ID per line"),
          rows        = 8,
          width       = "100%",
          resize      = "vertical"
        ),
        shiny::tags$small(
          style = "color:#888;",
          "One ID per line, or comma / space / tab separated. ",
          sprintf("Selection will be replaced (max %d).", max_n)
        ),
        footer = shiny::tagList(
          shiny::modalButton("Cancel"),
          shiny::actionButton(
            id_paste_btn,
            paste0("Load ", feature_label, " list"),
            icon  = shiny::icon("clipboard-list"),
            class = "btn btn-primary"
          )
        )
      ))
    }, ignoreInit = TRUE)

    # Load handler. Pasted text is split on commas and any whitespace, deduped,
    # validated against the same gene_choices vector used by manual Add, capped
    # at max_n, then loaded via set_genes() so the selection is REPLACED (not
    # appended) and pagination resets to page 1. On success the modal closes.
    # Notifications follow the spec:
    #   - "Loaded N valid <plural>."
    #   - "M <plural> were not found and were skipped."  (when any invalid)
    #   - "Only the first max_n <plural> were loaded."   (when over cap)
    #   - "No valid <plural> were found."                (when none matched)
    shiny::observeEvent(input[[id_paste_btn]], {
      raw <- input[[id_paste_input]]
      if (is.null(raw) || !nzchar(trimws(as.character(raw)))) {
        shiny::showNotification(
          sprintf("Paste at least one %s ID, then click Load.", feature_label),
          type = "message", duration = 3
        )
        return()
      }

      # Split on commas / any whitespace (spaces, tabs, newlines), then trim,
      # drop empties, dedupe — preserving first-seen order.
      tokens <- strsplit(as.character(raw), "[,\\s]+", perl = TRUE)[[1]]
      tokens <- trimws(tokens)
      tokens <- tokens[nzchar(tokens)]
      tokens <- unique(tokens)

      if (length(tokens) == 0L) {
        shiny::showNotification(
          sprintf("Paste at least one %s ID, then click Load.", feature_label),
          type = "message", duration = 3
        )
        return()
      }

      # Validate against the same choices set the manual Add path uses.
      valid   <- if (is.null(gene_choices)) character(0)
                 else tokens[tokens %in% gene_choices]
      invalid <- setdiff(tokens, valid)

      if (length(valid) == 0L) {
        shiny::showNotification(
          sprintf("No valid %s were found.", feature_label_plural),
          type = "warning", duration = 5
        )
        return()  # modal stays open so user can fix
      }

      capped <- length(valid) > max_n
      if (capped) valid <- valid[seq_len(max_n)]

      # Replace the current selection. set_genes() also resets current_page to
      # 1; the underlying add_genes() will not emit its own cap notification
      # because the list is already trimmed to <= max_n.
      set_genes(valid)

      # Successful load -> close the modal.
      shiny::removeModal()

      one_valid <- length(valid) == 1L
      shiny::showNotification(
        sprintf("Loaded %d valid %s.",
                length(valid),
                if (one_valid) feature_label else feature_label_plural),
        type = "default", duration = 4
      )
      if (length(invalid) > 0L) {
        one_inv <- length(invalid) == 1L
        shiny::showNotification(
          sprintf("%d %s %s not found and %s skipped.",
                  length(invalid),
                  if (one_inv) feature_label else feature_label_plural,
                  if (one_inv) "was" else "were",
                  if (one_inv) "was" else "were"),
          type = "warning", duration = 6
        )
      }
      if (capped) {
        shiny::showNotification(
          sprintf("Only the first %d %s were loaded.", max_n, feature_label_plural),
          type = "warning", duration = 6
        )
      }
    }, ignoreInit = TRUE)
  }
  .wire_mfe_paste_list()

  # Selection summary as removable tags. Each tag's × triggers
  # Shiny.setInputValue(id_remove_gene, <gene>, {priority:'event'}) so the
  # same gene can be re-removed after re-adding (priority:'event' makes the
  # input fire even when the value is identical to the previous one).
  # Tags reflect the CURRENT sort order so they line up with the paginated
  # Expression Profile / Data / DE rows below.
  output[[id_summary]] <- shiny::renderUI({
    sel <- ordered_selected_genes()
    if (length(sel) == 0) return(NULL)
    js_input_id <- id_remove_gene
    shiny::tags$div(
      style = "display:flex; flex-wrap:wrap; gap:6px; margin-top:8px;",
      lapply(sel, function(g) {
        # JSON-escape the gene ID for safe embedding in the onclick string
        g_js <- jsonlite::toJSON(as.character(g), auto_unbox = TRUE)
        shiny::tags$span(
          class = "mge-tag",
          shiny::tags$span(g),
          shiny::tags$span(
            class = "mge-tag-x",
            title = paste0("Remove ", g),
            onclick = sprintf(
              "Shiny.setInputValue('%s', %s, {priority:'event'});",
              js_input_id, g_js
            ),
            HTML("&times;")
          )
        )
      })
    )
  })

  # Selection count. The cap is still enforced inside add_genes() (max_n);
  # we just don't show "/ N" in the display anymore.
  output[[id_count]] <- shiny::renderUI({
    sel <- selected_genes()
    shiny::tags$div(
      style = "font-size:13px; color:#555; margin-top:6px;",
      sprintf("Selected %s: %d", feature_label_plural, length(sel))
    )
  })

  # Pagination UI: "Showing <plural> X-Y of N" + Previous/Next buttons.
  # Buttons are hidden when everything fits on one page. Disabled at boundaries.
  output[[id_pagination_ui]] <- shiny::renderUI({
    sel <- selected_genes()
    n   <- length(sel)
    if (n == 0L) return(NULL)

    pg      <- current_page()
    n_pages <- max(1L, as.integer(ceiling(n / page_size)))
    lo      <- (pg - 1L) * page_size + 1L
    hi      <- min(n, pg * page_size)

    text <- sprintf("Showing %s %d–%d of %d",
                    feature_label_plural, lo, hi, n)

    buttons <- if (n <= page_size) {
      NULL
    } else {
      prev_btn <- shiny::actionButton(
        id_page_prev,
        sprintf("Previous %d", page_size),
        class = "btn btn-default btn-sm"
      )
      next_btn <- shiny::actionButton(
        id_page_next,
        sprintf("Next %d", page_size),
        class = "btn btn-default btn-sm"
      )
      if (pg <= 1L)        prev_btn$attribs$disabled <- "disabled"
      if (pg >= n_pages)   next_btn$attribs$disabled <- "disabled"
      shiny::tagList(prev_btn, next_btn)
    }

    shiny::div(
      style = "display:flex; align-items:center; gap:12px; margin:8px 0;",
      shiny::span(text, style = "font-size:14px; color:#555;"),
      buttons
    )
  })

  # Long-format expression data + plot — paginated. Uses the current-page
  # subset; the full selection is shown only in the heatmap and the tag list.
  mge_data <- shiny::reactive({
    sel <- paged_selected_genes()
    shiny::validate(
      shiny::need(length(sel) > 0,
                  sprintf("Select up to %d %s to view their expression profiles.",
                          max_n, feature_label_plural))
    )
    color_by  <- shiny::req(input[[id_color_by]])
    scale_sel <- if (!is.null(fixed_scale)) fixed_scale else shiny::req(input[[id_scale]])

    df <- build_multi_gene_df_from_env(env, sel,
                                       scale_sel            = scale_sel,
                                       color_by             = color_by,
                                       expression_data_mode = expression_data_mode)
    shiny::validate(
      shiny::need(!is.null(df) && nrow(df) > 0,
                  sprintf("Expression data is not available for one or more selected %s.",
                          feature_label_plural))
    )
    list(df = df, color_by = color_by, scale_sel = scale_sel)
  })

  mge_plot <- shiny::reactive({
    d <- mge_data()
    make_multi_gene_facet_plot(
      df        = d$df,
      color_by  = d$color_by,
      scale_sel = plot_scale_arg_fn(d$scale_sel),
      ylab_log2 = ylab_log2_fn(env, d$scale_sel)
    )
  })

  # Dynamic height: ~240px per facet row (up to 4 cols). Sized to the current
  # page so the plot area doesn't grow unboundedly with the full selection.
  output[[id_plot]] <- shiny::renderPlot({
    mge_plot()
  }, height = function() {
    n <- length(paged_selected_genes())
    if (n == 0) return(320)
    cols <- if (n <= 2) n else if (n <= 6) 3 else 4
    rows <- ceiling(n / cols)
    max(320, rows * 240)
  }, res = 100)

  # Download plot — uses the current-page subset (matches what is on screen).
  dl_prefix <- paste0("Multi", tools::toTitleCase(feature_label), "Explorer")
  output[[id_dl_plot]] <- shiny::downloadHandler(
    filename = function() paste0(dl_prefix, "_", Sys.Date(), ".png"),
    content = function(file) {
      n <- length(paged_selected_genes())
      cols <- if (n <= 2) max(1, n) else if (n <= 6) 3 else 4
      rows <- ceiling(max(n, 1) / cols)
      ggplot2::ggsave(file, plot = mge_plot(),
                      width  = max(8, cols * 3.2),
                      height = max(4, rows * 2.6),
                      dpi    = 150)
    }
  )

  # Wide-format expression data view (one row per gene, one column per sample).
  # Reuses the shared build_explorer_wide_df() / render_explorer_wide_dt()
  # helpers, identical to what the single-feature Explorers render.
  mge_table_df <- shiny::reactive({
    d <- mge_data()
    build_explorer_wide_df(d$df, feature_col_name = feature_col_name)
  })

  # If we have navigation wiring (tabs_id + single_selectize_id +
  # single_tab_label all set), make the feature-name column clickable. The
  # render helper falls back to plain-text rendering when click_input_id is
  # NULL, so this is a no-op for callers that don't supply the nav params.
  nav_enabled <- !is.null(single_tab_label) && nzchar(as.character(single_tab_label)) &&
                 !is.null(tabs_id)          && nzchar(as.character(tabs_id)) &&
                 !is.null(single_selectize_id) && nzchar(as.character(single_selectize_id))
  click_input_for_table <- if (nav_enabled) id_explorer_click else NULL

  output[[id_table]] <- DT::renderDT({
    render_explorer_wide_dt(
      mge_table_df(),
      feature_col_name = feature_col_name,
      click_input_id   = click_input_for_table,
      clickable_col    = feature_col_name
    )
  })

  # Download keeps the RAW (un-rounded) numeric values so analysts can re-do
  # any precision-sensitive calculations downstream. The on-screen DT only
  # affects display.
  output[[id_dl_table]] <- shiny::downloadHandler(
    filename = function() paste0(dl_prefix, "_data_", Sys.Date(), ".csv"),
    content  = function(file) {
      df <- mge_table_df()
      utils::write.csv(if (is.null(df)) data.frame() else df,
                       file = file, row.names = FALSE)
    }
  )

  # ---- Clickable feature → single Explorer navigation -------------------
  # Mirrors the modal flow used by volcano / MA / clustering single-point
  # clicks (see wire_click_to_explorer):
  #   1. The Expression Data table's anchor onclick fires
  #      Shiny.setInputValue(id_explorer_click, <feature ID>, {priority:'event'}).
  #   2. We capture it in pending_explorer_nav and pop a confirmation modal
  #      titled "<Feat> selected" with a "View in <Tab>" primary button.
  #   3. On confirm we ONLY switch tabs and keep the value in
  #      pending_explorer_nav. We do NOT update the single Explorer's
  #      selectize here, because that races with the single Explorer's own
  #      tabs_id observer (ge_sync_choices / prot_sync_choices / etc.)
  #      which resets selected = character(0) using an isolate() read of
  #      the selectize input value — which the client hasn't round-tripped
  #      back yet, so our selection gets wiped.
  #   4. A separate observeEvent(input[[tabs_id]]) with priority = -100 reads
  #      the pending value once the tab is active and writes it into the
  #      selectize. The explicit negative priority is what guarantees our
  #      message is queued AFTER the single Explorer's reset; registration
  #      order alone is NOT a reliable guarantee (see the comment on that
  #      observer below).
  # priority:'event' on the JS side guarantees re-clicks of the same
  # feature re-fire the click observer, so the modal can be reopened.
  # Gated on nav_enabled so a Multi-Explorer wired without navigation params
  # silently renders plain text in the table and skips these observers.
  if (nav_enabled) {
    feat_cap <- paste0(toupper(substr(feature_label, 1L, 1L)),
                       substr(feature_label, 2L, nchar(feature_label)))
    pending_explorer_nav <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input[[id_explorer_click]], {
      val <- trimws(as.character(input[[id_explorer_click]]))
      shiny::req(nzchar(val))
      pending_explorer_nav(val)
      shiny::showModal(shiny::modalDialog(
        title = paste0(feat_cap, " selected"),
        shiny::tags$p(shiny::tags$b(val)),
        footer = shiny::tagList(
          shiny::modalButton("Cancel"),
          shiny::actionButton(id_explorer_view,
                              paste0("View in ", single_tab_label),
                              class = "btn-primary")
        ),
        easyClose = TRUE
      ))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input[[id_explorer_view]], {
      shiny::removeModal()
      val <- pending_explorer_nav()
      shiny::req(!is.null(val), nzchar(val))
      if (!is.null(gene_choices) && !(val %in% gene_choices)) {
        shiny::showNotification(
          paste0("Selected ", feature_label, " not found: ", val),
          type = "warning", duration = 3
        )
        pending_explorer_nav(NULL)
        return()
      }
      # Switch tabs. The selectize update lands in the post-nav observer
      # below, which fires AFTER the single Explorer's ge_sync_choices /
      # prot_sync_choices reset.
      bslib::nav_select(id = tabs_id, selected = single_tab_label)
    }, ignoreInit = TRUE)

    # Post-navigation selectize update. Fires whenever input[[tabs_id]]
    # changes; only acts when the user is on the target single Explorer tab
    # AND we have a pending feature to apply.
    #
    # priority = -100 is load-bearing. Relying on registration order (this
    # helper being called after wire_gene_explorer / wire_protein_explorer)
    # is NOT a reliable guarantee of observer execution order: when the
    # Explorer's own ge_sync_choices / prot_sync_choices resync happened to
    # run last it overwrote this write with character(0) and the selection
    # was silently lost. Shiny runs higher-priority observers first, so a
    # negative priority puts this write deterministically after every
    # default-priority (0) observer on the same input. Do not remove.
    shiny::observeEvent(input[[tabs_id]], {
      if (!identical(input[[tabs_id]], single_tab_label)) return()
      val <- pending_explorer_nav()
      if (is.null(val) || !nzchar(as.character(val))) return()
      shiny::updateSelectizeInput(session, single_selectize_id,
                                  choices  = gene_choices,
                                  selected = val,
                                  server   = TRUE)
      pending_explorer_nav(NULL)
    }, ignoreInit = TRUE, priority = -100)
  }

  # ---- Differential Expression Statistics (mirrors single Explorer) ----
  # Calls de_stats_fn() per feature, prepends a feature column under the
  # configured name (Gene / Protein / Metabolite), then renders with the same
  # render_de_stats_dt() helper that drives the single-feature DE table so
  # styling, formatting, and row coloring are identical.
  # Paginated: only the current-page features contribute rows.
  mge_de_stats_df <- shiny::reactive({
    sel <- paged_selected_genes()
    if (length(sel) == 0) return(NULL)
    if (!exists("de_stats",  envir = env, inherits = FALSE)) return(NULL)
    if (!exists("contrasts", envir = env, inherits = FALSE)) return(NULL)

    de_stats  <- get("de_stats",  envir = env, inherits = FALSE)
    contrasts <- get("contrasts", envir = env, inherits = FALSE)
    feat_col  <- feature_col_name

    parts <- lapply(sel, function(g) {
      one <- tryCatch(
        de_stats_fn(de_stats, contrasts, g),
        error = function(e) NULL
      )
      if (is.null(one) || nrow(one) == 0) return(NULL)
      lead <- setNames(data.frame(g, stringsAsFactors = FALSE), feat_col)
      cbind(lead, one)
    })
    parts <- parts[!vapply(parts, is.null, logical(1))]
    if (length(parts) == 0) return(NULL)

    out <- do.call(rbind, parts)
    out[[feat_col]] <- factor(out[[feat_col]],
                              levels = sel[sel %in% out[[feat_col]]])
    out
  })

  # Conditional contrast selector. The contrast picker is needed by ANY
  # sort mode that scores against a contrast — currently "significance"
  # (pvalue ascending) and "linearfc" (linearFC ascending). Renders:
  #   * nothing when sort mode doesn't need a contrast
  #   * nothing when env has <=1 contrast (auto-selected without a dropdown)
  #   * a small Contrast: dropdown alongside the Sort-by selector otherwise
  output[[id_sort_contrast_ui]] <- shiny::renderUI({
    mode <- input[[id_sort_by]]
    if (is.null(mode) || !mode %in% c("significance", "linearfc")) return(NULL)
    if (length(contrast_names_all) <= 1L) return(NULL)
    # Read the picker value via isolate so this renderUI only re-fires when
    # the sort mode (its framing condition) changes — not on every contrast
    # change. The picker's live value flows through input[[id_sort_contrast_pick]]
    # directly; downstream sort + sync observers consume it without needing
    # the dropdown widget to be re-mounted.
    cur <- shiny::isolate(input[[id_sort_contrast_pick]])
    if (is.null(cur) || !cur %in% contrast_names_all) cur <- contrast_names_all[1L]
    shiny::div(
      style = "min-width:200px;",
      shiny::selectInput(
        inputId  = id_sort_contrast_pick,
        label    = "Contrast:",
        choices  = contrast_names_all,
        selected = cur,
        width    = "100%"
      )
    )
  })

  # Cutoffs UI — same shared helper used by the single Gene Explorer.
  output[[id_de_cutoffs_ui]] <- shiny::renderUI({
    padj0 <- if (exists("padj_cutoff", envir = env, inherits = FALSE)) {
      get("padj_cutoff", envir = env, inherits = FALSE)
    } else 0.05
    fc0 <- if (exists("fc_cutoff", envir = env, inherits = FALSE)) {
      get("fc_cutoff", envir = env, inherits = FALSE)
    } else 1
    build_de_cutoffs_ui(
      id_sig_mode   = id_de_sig_mode,
      id_sig_cutoff = id_de_padj_cut,
      id_fc_cutoff  = id_de_fc_cut,
      sig0          = padj0,
      fc0           = fc0
    )
  })

  # No need to update the numericInput's label — build_de_cutoffs_ui()'s
  # selectInput (selectize = FALSE) is itself the visible label, styled via
  # the .sig-mode-dropdown CSS to look like a small label with a ▼.
  # The numeric input below it is rendered with label = NULL.

  # Same legend HTML used by the single Gene Explorer
  output[[id_de_legend_ui]] <- shiny::renderUI({
    shiny::req(length(selected_genes()) > 0)
    de_table_legend_html()
  })

  # ---- DE statistics, one contrast at a time --------------------------
  # mge_de_stats_df() returns a long data frame with a `Contrast` column.
  # This block renders a single DT for the user-picked contrast, with a
  # selectInput dropdown above it (only when there is more than one
  # contrast). The cutoff inputs and legend are unchanged.
  # contrast_names_all was computed earlier (used by the sort-ordering layer).

  # Default cutoffs from env (used as fallbacks before the inputs render)
  padj0_default <- if (exists("padj_cutoff", envir = env, inherits = FALSE))
    get("padj_cutoff", envir = env, inherits = FALSE) else 0.05
  fc0_default   <- if (exists("fc_cutoff", envir = env, inherits = FALSE))
    get("fc_cutoff", envir = env, inherits = FALSE) else 1

  # Resolve the currently selected contrast, defaulting to the first one
  # so the DT can render before the dropdown ever fires.
  selected_contrast <- shiny::reactive({
    if (length(contrast_names_all) == 0L) return(NA_character_)
    cur <- input[[id_de_contrast_pick]]
    if (is.null(cur) || !nzchar(as.character(cur)) ||
        !cur %in% contrast_names_all) {
      return(contrast_names_all[1L])
    }
    as.character(cur)
  })

  # Contrast picker UI rendered into the slot above the DT in ui.R.
  #   * 0 contrasts -> friendly inline message.
  #   * 1 contrast  -> show only the contrast name as a subtitle (no
  #                    "Selected contrast:" prefix, no dropdown).
  #   * >1 contrasts -> show only the dropdown selector. No additional
  #                     subtitle: the dropdown's current value is itself
  #                     the indicator of which contrast is displayed.
  output[[id_de_contrast_ui]] <- shiny::renderUI({
    shiny::req(length(selected_genes()) > 0)
    if (length(contrast_names_all) == 0L) {
      return(shiny::tags$div(
        "No DE contrasts available for this dataset.",
        style = "color:#888; padding:8px 0;"
      ))
    }
    if (length(contrast_names_all) == 1L) {
      return(shiny::tags$div(
        contrast_names_all[1L],
        style = "font-size:15px; font-weight:600; color:#2c3e50; margin-top:6px; margin-bottom:8px;"
      ))
    }
    # selected_contrast() reads input[[id_de_contrast_pick]]; calling it
    # here directly would make this renderUI re-fire on every contrast
    # change and re-mount the dropdown widget unnecessarily. The DT below
    # consumes selected_contrast() reactively and re-renders correctly on
    # its own. isolate() keeps this renderUI tied only to selection
    # presence + contrast count (its framing conditions).
    shiny::div(
      style = "max-width:420px; margin-top:6px; margin-bottom:8px;",
      shiny::selectInput(
        inputId  = id_de_contrast_pick,
        label    = "Contrast:",
        choices  = contrast_names_all,
        selected = shiny::isolate(selected_contrast()),
        width    = "100%"
      )
    )
  })

  # Single DT for the currently selected contrast. Updates whenever the
  # selected features, the contrast pick, or the cutoffs change.
  output[[id_de_table]] <- DT::renderDataTable({
    shiny::req(length(selected_genes()) > 0)
    if (length(contrast_names_all) == 0L) {
      shiny::validate(shiny::need(FALSE, "No DE contrasts available for this dataset."))
    }
    cur <- selected_contrast()
    df_all <- mge_de_stats_df()
    shiny::validate(shiny::need(
      !is.null(df_all) && nrow(df_all) > 0,
      sprintf("Select %s to view their differential expression statistics.",
              feature_label_plural)
    ))
    sub <- df_all[as.character(df_all$Contrast) == cur, , drop = FALSE]
    shiny::validate(shiny::need(
      nrow(sub) > 0,
      sprintf("No statistics for selected %s in contrast '%s'.",
              feature_label_plural, cur)
    ))
    # Drop the redundant Contrast column — the subtitle above already names it.
    sub$Contrast <- NULL
    render_de_stats_dt(
      df          = sub,
      sig_thr     = input[[id_de_padj_cut]],
      fc_thr      = input[[id_de_fc_cut]],
      sig_default = padj0_default,
      fc_default  = fc0_default,
      sig_mode    = if (!is.null(input[[id_de_sig_mode]])) input[[id_de_sig_mode]] else "padj"
    )
  })

  # ---- Selected-features hierarchical heatmap -------------------------
  # Built dynamically from env$expr_norm restricted to the FULL current
  # selection (not paged, not sorted). Z-scores and hierarchical clustering
  # are recomputed on the subset, so non-DE features are included alongside
  # DE features and the row order reflects clustering of THIS selection.
  # See prepare_selected_features_heatmap_inputs() for the build pipeline.
  output[[id_sel_hm_title]] <- shiny::renderUI({
    shiny::req(length(selected_genes()) > 0)
    shiny::h4(
      paste0("Selected ", tools::toTitleCase(feature_label_plural),
             " Hierarchical Heatmap"),
      class = "plot-title full"
    )
  })

  output[[id_sel_hm_msg]] <- shiny::renderUI({
    sel <- selected_genes()
    if (length(sel) == 0L) return(NULL)
    if (length(sel) < 2L) {
      return(shiny::tags$div(
        "Select at least 2 features to display a heatmap.",
        style = "color:#888; padding:8px 0;"
      ))
    }
    # The selected-features heatmap is now driven by env$expr_norm so it
    # can include non-DE features. Show a friendly message only when
    # expr_norm itself is unavailable on the env.
    if (!exists("expr_norm", envir = env, inherits = FALSE) ||
        is.null(get("expr_norm", envir = env, inherits = FALSE))) {
      return(shiny::tags$div(
        "Heatmap data is not available for this dataset.",
        style = "color:#888; padding:8px 0;"
      ))
    }
    NULL
  })

  output[[id_sel_hm]] <- plotly::renderPlotly({
    sel <- selected_genes()
    shiny::validate(shiny::need(
      length(sel) >= 2L,
      "Select at least 2 features to display a heatmap."
    ))
    prepped <- tryCatch(
      prepare_selected_features_heatmap_inputs(env, sel),
      error = function(e) NULL
    )
    # prepped is NULL when fewer than 2 of the selected IDs survive the
    # cleaning pass (missing from expr_norm, all-NA after scaling, zero
    # variance, etc.). Surface that with the spec'd message rather than
    # the generic "could not be built" string.
    shiny::validate(shiny::need(
      !is.null(prepped),
      "Select at least 2 valid features to display a heatmap."
    ))
    plot_selected_features_heatmap_plotly(prepped, source_id = NULL)
  })

  # Annotation section (canonical feature_annot from this env) for the
  # full current selection — not paged/sorted, mirrors the heatmap's use of the
  # complete selection.
  wire_feature_annotation(
    output           = output,
    env              = env,
    id_prefix        = id_prefix,
    ids_reactive     = selected_genes,
    feature_col_name = feature_col_name,
    multiple         = TRUE
  )

  # Public controller — lets other handlers (e.g. plotly box-select) update
  # the multi-gene selection without reaching into private state.
  list(
    add     = add_genes,
    set     = set_genes,
    remove  = remove_gene,
    clear   = clear_genes,
    current = function() selected_genes(),
    max_n   = max_n
  )
}

# Wire plot/data download handlers (shared)
wire_downloads <- function(output,
                           id_plot_btn,
                           id_data_btn,
                           get_plot,
                           get_df,
                           prefix,
                           get_selected_id,
                           get_scale_sel,
                           norm_method_label = NULL,
                           log2_levels = c("log2")) {
  
  output[[id_plot_btn]] <- shiny::downloadHandler(
    filename = function() {
      scale_sel <- shiny::req(get_scale_sel())
      feat <- sanitize_filename(shiny::req(get_selected_id()))
      
      is_log2 <- scale_sel %in% log2_levels
      
      suffix <- if (is_log2) {
        nm <- if (scale_sel == "log2") {
          if (is.null(norm_method_label) || !nzchar(norm_method_label)) "LOG2" else toupper(gsub("\\s+.*$", "", norm_method_label))
        } else {
          toupper(as.character(scale_sel))
        }
        paste0("_", nm, "_log2.png")
      } else {
        "_norm_counts.png"
      }
      
      paste0(prefix, "_", feat, suffix)
    },
    content = function(file) {
      g <- get_plot()
      g_small <- g + ggplot2::theme(
        text = ggplot2::element_text(size = 10),
        axis.text.x = ggplot2::element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 8),
        legend.text = ggplot2::element_text(size = 8),
        legend.title = ggplot2::element_text(size = 9)
      )
      ggplot2::ggsave(filename = file, plot = g_small, width = 8, height = 4, dpi = 300)
    }
  )
  
  output[[id_data_btn]] <- shiny::downloadHandler(
    filename = function() {
      scale_sel <- shiny::req(get_scale_sel())
      feat <- sanitize_filename(shiny::req(get_selected_id()))
      
      is_log2 <- scale_sel %in% log2_levels
      
      suffix <- if (is_log2) {
        nm <- if (scale_sel == "log2") {
          if (is.null(norm_method_label) || !nzchar(norm_method_label)) "LOG2" else toupper(gsub("\\s+.*$", "", norm_method_label))
        } else {
          toupper(as.character(scale_sel))
        }
        paste0("_", nm, "_log2.txt")
      } else {
        "_norm_counts.txt"
      }
      
      paste0(prefix, "_", feat, suffix)
    },
    content = function(file) {
      df <- get_df()
      utils::write.table(df, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )
}

# Wire plotly click events to jump into Explorer tab and preselect the clicked ID (shared)
wire_click_to_explorer <- function(session, input, tabs_id,
                                   selected_tab_label,
                                   selectize_input_id,
                                   choices_vec,
                                   choices_ready_rv,
                                   pending_rv,
                                   clicked_rv,
                                   view_btn_id,
                                   modal_title,
                                   click_sources) {
  
  for (src in click_sources) {
    local({
      src_i <- src
      shiny::observeEvent(plotly::event_data("plotly_click", source = src_i), {
        ed <- plotly::event_data("plotly_click", source = src_i)
        shiny::req(ed)
        
        val <- if (!is.null(ed$key)) ed$key else if (!is.null(ed$y)) ed$y else NA
        val <- trimws(as.character(val))
        shiny::req(nzchar(val))
        
        clicked_rv(val)
        
        shiny::showModal(shiny::modalDialog(
          title = modal_title,
          shiny::tags$p(shiny::tags$b(val)),
          footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(view_btn_id, paste0("View in ", selected_tab_label), class = "btn-primary")
          ),
          easyClose = TRUE
        ))
      }, ignoreInit = TRUE)
    })
  }
  
  shiny::observeEvent(input[[view_btn_id]], {
    shiny::removeModal()
    val <- trimws(as.character(clicked_rv()))
    shiny::req(nzchar(val))
    pending_rv(val)
    bslib::nav_select(id = tabs_id, selected = selected_tab_label)
  }, ignoreInit = TRUE)
  
  # Deferred phase 2: write the pending feature into the Explorer's selectize
  # once the target tab is active and its choices are ready.
  #
  # priority = -100 is load-bearing. This write must land AFTER the single
  # Explorer's own ge_sync_choices / prot_sync_choices resync, which fires on
  # the same tab-change flush and resets selected to character(0). Registration
  # order alone is NOT a reliable guarantee of observer execution order — when
  # the resync ran last it overwrote this write and the selection was silently
  # lost. Shiny runs higher-priority observers first, so a negative priority
  # puts this write deterministically after every default-priority (0)
  # observer. Do not remove.
  shiny::observe({
    shiny::req(identical(input[[tabs_id]], selected_tab_label))
    shiny::req(isTRUE(choices_ready_rv()))

    val <- pending_rv()
    shiny::req(!is.null(val), nzchar(val))

    if (!is.null(choices_vec) && val %in% choices_vec) {
      shiny::updateSelectizeInput(
        session,
        inputId  = selectize_input_id,
        choices  = choices_vec,
        selected = val,
        server   = TRUE
      )
      shiny::invalidateLater(150, session)
      pending_rv(NULL)
    } else {
      shiny::showNotification(paste0("Selected item not found in list: ", val), type = "warning", duration = 3)
      pending_rv(NULL)
    }
  }, priority = -100)
}

# Wire plotly box/lasso selection events to the Multi-Gene Explorer.
# - length(selected) == 1 -> show single-gene modal that, on confirm, sets
#   the single Gene Explorer selectize and switches to that tab.
# - length(selected) >  1 -> show multi-gene modal that, on confirm, adds
#   the IDs to the multi-gene controller and switches to that tab.
wire_multi_select_to_explorer <- function(session, input, output, tabs_id,
                                          single_tab_label   = "Gene Explorer",
                                          multi_tab_label    = "Multi-Gene Explorer",
                                          single_selectize_id,           # e.g. "gene_explorer_gene"
                                          single_choices_vec  = NULL,
                                          mge_controller,                # from wire_multi_gene_explorer()
                                          select_sources,
                                          id_prefix = "mge_select",
                                          feature_label_singular = "gene",
                                          feature_label_plural   = "genes") {
  feat_cap <- paste0(toupper(substr(feature_label_singular, 1, 1)),
                     substr(feature_label_singular, 2, nchar(feature_label_singular)))
  id_pending_single <- paste0(id_prefix, "_pending_single")
  id_pending_multi  <- paste0(id_prefix, "_pending_multi")
  id_open_single    <- paste0(id_prefix, "_open_single_btn")
  id_open_multi     <- paste0(id_prefix, "_open_multi_btn")

  pending_single <- shiny::reactiveVal(NULL)
  pending_multi  <- shiny::reactiveVal(NULL)

  for (src in select_sources) {
    local({
      src_i <- src
      shiny::observeEvent(
        suppressWarnings(plotly::event_data("plotly_selected", source = src_i)),
        {
          ed <- plotly::event_data("plotly_selected", source = src_i)
          shiny::req(!is.null(ed))

          keys <- if (!is.null(ed$key)) ed$key else NULL
          keys <- trimws(as.character(keys))
          keys <- unique(keys[nzchar(keys) & !is.na(keys)])
          if (length(keys) == 0) return()

          if (length(keys) == 1) {
            # Single-feature flow — keep existing UX (modal + view button)
            pending_single(keys[1])
            shiny::showModal(shiny::modalDialog(
              title = paste0(feat_cap, " selected"),
              shiny::tags$p(shiny::tags$b(keys[1])),
              footer = shiny::tagList(
                shiny::modalButton("Cancel"),
                shiny::actionButton(id_open_single,
                                    paste0("View in ", single_tab_label),
                                    class = "btn-primary")
              ),
              easyClose = TRUE
            ))
          } else {
            # Multi-feature flow
            max_n <- mge_controller$max_n
            if (length(keys) > max_n) {
              shiny::showNotification(
                sprintf("You selected %d %s; only the first %d will be added.",
                        length(keys), feature_label_plural, max_n),
                type = "warning", duration = 5
              )
              keys <- keys[seq_len(max_n)]
            }
            pending_multi(keys)
            shiny::showModal(shiny::modalDialog(
              title = paste0("Open ", multi_tab_label, "?"),
              shiny::tags$p(sprintf(
                "You selected multiple %s. Would you like to open them in the %s?",
                feature_label_plural, multi_tab_label
              )),
              shiny::tags$p(shiny::tags$small(
                sprintf("Selected: %d %s", length(keys), feature_label_plural)
              )),
              footer = shiny::tagList(
                shiny::modalButton("Cancel"),
                shiny::actionButton(id_open_multi, paste0("Open ", multi_tab_label),
                                    class = "btn-primary")
              ),
              easyClose = TRUE
            ))
          }
        },
        ignoreInit = TRUE
      )
    })
  }

  # Confirm — single-gene path: validate, then ONLY switch tabs. The selectize
  # write is deferred to the post-navigation observer below, which runs after
  # the single Explorer's own choices resync. Writing it here (before
  # nav_select) let the resync clobber the selection with character(0).
  shiny::observeEvent(input[[id_open_single]], {
    shiny::removeModal()
    val <- pending_single()
    shiny::req(!is.null(val), nzchar(val))

    if (!is.null(single_choices_vec) && !(val %in% single_choices_vec)) {
      shiny::showNotification(
        sprintf("Selected %s not found in list: %s", feature_label_singular, val),
        type = "warning", duration = 3
      )
      pending_single(NULL)
      return()
    }
    # Keep pending_single set; phase 2 consumes and clears it.
    bslib::nav_select(id = tabs_id, selected = single_tab_label)
  }, ignoreInit = TRUE)

  # Deferred phase 2 of the single-gene path.
  #
  # priority = -100 is load-bearing. This write must land AFTER the single
  # Explorer's own ge_sync_choices / prot_sync_choices resync, which fires on
  # the same tab-change flush and resets selected to character(0). Registration
  # order alone is NOT a reliable guarantee of observer execution order. Shiny
  # runs higher-priority observers first, so a negative priority puts this
  # write deterministically after every default-priority (0) observer on the
  # same input. Do not remove.
  shiny::observeEvent(input[[tabs_id]], {
    if (!identical(input[[tabs_id]], single_tab_label)) return()
    val <- pending_single()
    if (is.null(val) || !nzchar(as.character(val))) return()
    shiny::updateSelectizeInput(session, single_selectize_id,
                                choices  = single_choices_vec,
                                selected = val,
                                server   = TRUE)
    pending_single(NULL)
  }, ignoreInit = TRUE, priority = -100)

  # Confirm — multi-gene path: REPLACE the current selection with the newly
  # selected genes (still deduped, still capped at max_n inside set/add).
  # This matches the user's expectation that a fresh box-select starts a new
  # multi-gene exploration rather than appending to a stale list.
  shiny::observeEvent(input[[id_open_multi]], {
    shiny::removeModal()
    ids <- pending_multi()
    shiny::req(!is.null(ids), length(ids) > 0)
    mge_controller$set(ids)
    bslib::nav_select(id = tabs_id, selected = multi_tab_label)
    pending_multi(NULL)
  }, ignoreInit = TRUE)
}

# Build Gene Explorer DE stats table from env-format (gene + linearFC.<contrast>, pvalue.<contrast>, padj.<contrast>)
get_gene_de_stats_for_table <- function(de_stats, contrasts, gene_id) {
  if (is.null(de_stats) || is.null(contrasts)) return(NULL)
  
  # ID column in RNA DE stats
  if (!"FeatureID" %in% colnames(de_stats)) return(NULL)
  
  # Contrast names table
  if (!"Contrast_name" %in% colnames(contrasts)) return(NULL)
  
  # Find the selected gene row
  row <- de_stats[de_stats$FeatureID == gene_id, , drop = FALSE]
  if (nrow(row) == 0) return(NULL)
  
  # Build one row per contrast
  cnames <- contrasts$Contrast_name
  out_list <- list()
  
  for (cname in cnames) {
    col_fc   <- paste0("linearFC.", cname)
    col_p    <- paste0("pvalue.",  cname)
    col_padj <- paste0("padj.",    cname)
    
    # Skip contrasts that don't exist in the stats table
    if (!all(c(col_fc, col_p, col_padj) %in% colnames(row))) next
    
    out_list[[cname]] <- data.frame(
      Contrast = cname,
      linearFC = as.numeric(row[[col_fc]]),
      pvalue   = as.numeric(row[[col_p]]),
      padj     = as.numeric(row[[col_padj]]),
      stringsAsFactors = FALSE
    )
  }
  
  # Return NULL if no contrasts were found
  if (length(out_list) == 0) return(NULL)
  
  do.call(rbind, out_list)
}

# ---------------------------------------------------------------------------
# Canonical feature annotation (feature_annot) helpers
# ---------------------------------------------------------------------------
# Each omics env carries a canonical `feature_annot` data.frame:
#   * rownames(feature_annot) are the feature IDs,
#   * columns are annotation / metadata fields,
#   * it may be missing, NULL, have zero rows, or zero annotation columns.
# These helpers read it defensively and never assume rownames are perfectly
# synchronized with the selected feature IDs. They fully replace the previous
# global `annot` / `trinotate_main` logic and the USE_EXTERNAL_ANNOTATION /
# USE_TRINOTATE_ANNOTATION flags.

# Pull the canonical feature_annot data.frame from an omics env, or NULL.
# preserves rownames (the feature IDs) and column names verbatim.
get_feature_annot_df <- function(env) {
  if (is.null(env) || !is.environment(env)) return(NULL)
  if (!exists("feature_annot", envir = env, inherits = FALSE)) return(NULL)
  fa <- get("feature_annot", envir = env, inherits = FALSE)
  if (is.null(fa)) return(NULL)
  fa <- tryCatch(
    as.data.frame(fa, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL)
  if (is.null(fa) || !is.data.frame(fa)) return(NULL)
  fa
}

# Classify a feature_annot table:
#   "missing"    -> NULL / not a data.frame / zero rows (no features at all)
#   "no_columns" -> a valid table that carries feature IDs but no annotation
#                   columns (a valid-but-empty annotation table)
#   "ok"         -> usable with >= 1 annotation column and >= 1 row
feature_annot_status <- function(fa) {
  if (is.null(fa) || !is.data.frame(fa)) return("missing")
  if (ncol(fa) == 0L) return("no_columns")
  if (nrow(fa) == 0L) return("missing")
  "ok"
}

# Extract annotation rows for one or more feature IDs into a clean data.frame
# with an explicit leading feature-ID column (feature_col_name) followed by all
# annotation columns. Feature IDs absent from rownames(fa) are skipped, so a
# multi-feature selection returns only the rows that exist. Returns a 0-row
# data.frame (with the correct columns) when none of the IDs are present.
build_feature_annot_table <- function(fa, feature_ids,
                                       feature_col_name = "Feature") {
  ids <- unique(as.character(feature_ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  ann_cols <- if (is.null(fa)) character(0) else colnames(fa)

  present <- ids[ids %in% rownames(fa)]

  out <- data.frame(.fid = present, stringsAsFactors = FALSE,
                    check.names = FALSE)
  names(out)[1] <- feature_col_name

  if (length(ann_cols) > 0L) {
    sub <- fa[present, ann_cols, drop = FALSE]
    rownames(sub) <- NULL
    out <- cbind(out, sub, stringsAsFactors = FALSE)
  }
  rownames(out) <- NULL
  out
}

# Detect "empty / invalid" annotation values. A value is blank if, after
# trimming, it is NA / empty, a known null-token (NA, N/A, null, NULL), or
# composed solely of separator / whitespace characters (e.g. ".", "-", "|",
# "||", "| |"). Vectorized; NA-safe. Used to drop columns that carry no useful
# information (e.g. a column where every selected row is just "|").
.is_blank_annot_value <- function(x) {
  x <- trimws(as.character(x))
  out <- is.na(x) | x == ""
  out[is.na(out)] <- TRUE
  nn <- !out
  if (any(nn)) {
    tok <- grepl("^(na|n/?a|null)$", x[nn], ignore.case = TRUE)
    sep <- grepl("^[|;,./_[:space:]-]+$", x[nn])
    out[nn] <- tok | sep
  }
  out
}

# Drop annotation columns whose values are ALL blank (per
# .is_blank_annot_value), always keeping the leading id column(s) named in
# keep_cols. For multi-feature tables this runs AFTER subsetting to the
# selected features, so a column that is populated globally but empty for the
# current selection is removed from the displayed table.
drop_blank_annot_columns <- function(df, keep_cols) {
  ann_cols <- setdiff(colnames(df), keep_cols)
  if (length(ann_cols) == 0L) return(df)
  drop <- vapply(ann_cols,
                 function(cn) all(.is_blank_annot_value(df[[cn]])),
                 logical(1))
  df[, c(keep_cols, ann_cols[!drop]), drop = FALSE]
}

# Interactive DT for an annotation table: explicit Feature-ID column, a global
# search box, and horizontal scrolling for many columns. No DT export/column-
# visibility buttons (Excel export is a separate custom download button below
# the table). Styling matches the app's other DT tables (compact/stripe/hover).
#
# selection: leave NULL (the default) for every read-only annotation table --
#   the argument is then NOT passed to DT::datatable() at all, so the existing
#   Annotation sections keep DT's own default behaviour unchanged. The Keyword
#   Search results table passes "multiple" to enable native DT row selection.
#   (Defaulting this to "none" instead of NULL would silently change all the
#   pre-existing annotation tables, since DT's own default is "multiple".)
# searching: leave NULL (the default) to keep DT's built-in search field, as all
#   the pre-existing annotation tables do. The Keyword Search results table
#   passes FALSE to drop it -- the panel already has its own keyword box, and a
#   second search field there is confusing. Only the `f` is removed from `dom`;
#   sorting, paging, scrollX, pageLength, row selection and styling are
#   unaffected.
# escape: leave NULL (the default) so DT escapes every cell, as all the
#   pre-existing annotation tables require. ONLY the Keyword Search results
#   table passes FALSE, and only after build_keyword_highlight_df() has already
#   HTML-escaped every cell itself and injected the controlled <mark> markup.
#   Never pass FALSE for a table whose values have not been escaped upstream.
render_feature_annot_dt <- function(table_df, selection = NULL, searching = NULL,
                                    escape = NULL, length_menu = NULL) {
  show_search <- if (is.null(searching)) TRUE else isTRUE(searching)
  esc         <- if (is.null(escape)) TRUE else escape
  show_len    <- !is.null(length_menu)
  # dom parts: l=page-length selector, f=search box, t=table, i=info, p=paging.
  # Default (both NULL) is "ftip" -- byte-identical to the pre-existing tables.
  opts <- list(
    dom        = paste0(if (show_len) "l" else "",
                        if (show_search) "f" else "", "tip"),
    scrollX    = TRUE,
    searching  = show_search,
    ordering   = TRUE,
    paging     = TRUE,
    pageLength = 10,
    autoWidth  = FALSE
  )
  if (show_len) opts$lengthMenu <- length_menu

  if (is.null(selection)) {
    DT::datatable(
      table_df,
      rownames = FALSE,
      escape   = esc,
      options  = opts,
      class    = "compact stripe hover"
    )
  } else {
    DT::datatable(
      table_df,
      rownames  = FALSE,
      escape    = esc,
      options   = opts,
      class     = "compact stripe hover",
      selection = selection
    )
  }
}

# Wire the "Annotation" section for any explorer (single or multi).
#   output           : Shiny output object
#   env              : omics env holding the canonical feature_annot
#   id_prefix        : output-id prefix; renders <prefix>_annot_title_ui,
#                      <prefix>_annot_dt, <prefix>_annot_dl_ui (download button)
#                      and the <prefix>_annot_download handler.
#   ids_reactive     : reactive() returning the currently selected feature ID(s)
#   feature_col_name : leading ID-column label ("Gene"/"Protein"/"Metabolite")
#   multiple         : TRUE for multi-feature explorers (plural fallback wording)
# Never crashes when feature_annot is missing / NULL / empty / zero-columns,
# when the selected IDs are absent, or when every annotation column is blank for
# the current selection — a clear user-facing message is shown instead of an
# empty table, and the download button is hidden in that case.
wire_feature_annotation <- function(output, env, id_prefix, ids_reactive,
                                     feature_col_name = "Feature",
                                     multiple = FALSE) {
  id_title    <- paste0(id_prefix, "_annot_title_ui")
  id_dt       <- paste0(id_prefix, "_annot_dt")
  id_dl_ui    <- paste0(id_prefix, "_annot_dl_ui")
  id_download <- paste0(id_prefix, "_annot_download")

  msg_none <- if (multiple)
    "No annotation information is available for the selected features."
  else
    "No annotation information is available for the selected feature."

  msg_nocols <- if (multiple)
    "No annotation information is currently available for these features."
  else
    "No annotation information is currently available for this feature."

  # Single source of truth for this section: returns list(table, msg). `table`
  # is the cleaned, display-ready data.frame (NULL when a message should show);
  # `msg` is the fallback text (NULL when a table is available). Shared by the
  # DT, the download handler, and the download-button visibility.
  annot_state <- shiny::reactive({
    ids <- as.character(ids_reactive())
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (length(ids) == 0L) return(list(table = NULL, msg = NULL))

    fa     <- get_feature_annot_df(env)
    status <- feature_annot_status(fa)
    if (identical(status, "missing"))    return(list(table = NULL, msg = msg_none))
    if (identical(status, "no_columns")) return(list(table = NULL, msg = msg_nocols))

    tbl <- build_feature_annot_table(fa, ids,
                                     feature_col_name = feature_col_name)
    if (nrow(tbl) == 0L) return(list(table = NULL, msg = msg_none))

    # Hide annotation columns that are entirely blank for the current selection.
    tbl <- drop_blank_annot_columns(tbl, keep_cols = feature_col_name)
    if (ncol(tbl) <= 1L) return(list(table = NULL, msg = msg_nocols))

    list(table = tbl, msg = NULL)
  })

  output[[id_title]] <- shiny::renderUI({
    shiny::req(length(ids_reactive()) > 0)
    shiny::h4("Annotation", class = "plot-title full",
              style = "margin-top:10px; margin-bottom:10px;")
  })

  output[[id_dt]] <- DT::renderDT({
    st <- annot_state()
    if (!is.null(st$msg)) shiny::validate(shiny::need(FALSE, st$msg))
    shiny::req(!is.null(st$table))
    render_feature_annot_dt(st$table)
  })

  # Excel-only export of exactly what is displayed (cleaned columns).
  output[[id_download]] <- shiny::downloadHandler(
    filename = function() paste0(id_prefix, "_annotation_", Sys.Date(), ".xlsx"),
    content  = function(file) {
      st <- annot_state()
      df <- if (is.null(st$table)) data.frame() else st$table
      writexl::write_xlsx(df, path = file)
    }
  )

  # Single custom download button, shown below the table on the right, only
  # when there is a table to download. Matches the app's other download buttons.
  output[[id_dl_ui]] <- shiny::renderUI({
    st <- annot_state()
    shiny::req(!is.null(st$table))
    shiny::div(
      style = "display:flex; justify-content:flex-end; margin-top:8px;",
      shiny::downloadButton(id_download, "Download Annotation")
    )
  })
}

# ---------------------------------------------------------------------------
# KEYWORD-MATCH HIGHLIGHTING (display only)
# ---------------------------------------------------------------------------
# Produces a DISPLAY-ONLY copy of a Keyword Search result table in which each
# occurrence of the searched keyword is wrapped in <mark class="kws-hit">.
#
# Hard rules:
#   * The original result data frame is never modified. It stays the source of
#     truth for row->ID mapping, Explorer navigation, Download Results and all
#     counts. Only the copy handed to DT carries markup.
#   * Every cell that reaches DT is HTML-escaped HERE, so the results table can
#     safely be rendered with escape = FALSE. Escaping is applied per SEGMENT
#     (before the markup is inserted), so escaping can never shift the match
#     offsets and a value containing "<script>" stays inert.
#   * Matching mirrors the search itself: both sides lowercased, fixed = TRUE,
#     so it is literal and case-insensitive while the ORIGINAL capitalisation
#     is what gets displayed. Every occurrence in a cell is highlighted.
#   * The leading technical-ID column is escaped but never highlighted -- it is
#     the navigation key, not annotation. An ordinary annotation column whose
#     values happen to look like IDs (e.g. NameID) IS highlighted, because it is
#     part of the searchable annotation data.
#   * Columns with no match are left untouched when they are numeric/logical
#     (safe as-is, and this preserves DT's numeric sorting and costs nothing).
#     Character/factor columns are always escaped even with no match, because
#     escape = FALSE applies to the whole table.

.KWS_HL_OPEN  <- '<mark class="kws-hit">'
.KWS_HL_CLOSE <- '</mark>'

# Highlight one already-non-NA string given its gregexpr() match object.
# `needs_escape` is FALSE when the string provably contains no &, < or >, in
# which case escaping is a no-op and is skipped entirely -- that is the case for
# the overwhelming majority of annotation values and is what keeps the hot loop
# cheap. k == 1 (one match in the cell) is the common case and is handled
# without allocating the segment buffer.
.kws_hl_one <- function(s, mi, needs_escape) {
  st <- as.integer(mi)
  ml <- attr(mi, "match.length")
  k  <- length(st)
  n  <- nchar(s)
  esc <- if (needs_escape) htmltools::htmlEscape else identity

  if (k == 1L) {
    p <- st[1L]; len <- ml[1L]
    return(paste0(
      if (p > 1L) esc(substr(s, 1L, p - 1L)) else "",
      .KWS_HL_OPEN, esc(substr(s, p, p + len - 1L)), .KWS_HL_CLOSE,
      if (p + len <= n) esc(substr(s, p + len, n)) else ""))
  }

  out <- character(4L * k + 1L)
  j <- 0L; pos <- 1L
  for (i in seq_len(k)) {
    p <- st[i]; len <- ml[i]
    out[j + 1L] <- if (p > pos) esc(substr(s, pos, p - 1L)) else ""
    out[j + 2L] <- .KWS_HL_OPEN
    out[j + 3L] <- esc(substr(s, p, p + len - 1L))
    out[j + 4L] <- .KWS_HL_CLOSE
    j   <- j + 4L
    pos <- p + len
  }
  out[j + 1L] <- if (pos <= n) esc(substr(s, pos, n)) else ""
  paste0(out, collapse = "")
}

# Escape + highlight one column. Returns NULL when the column needs no change
# at all, so the caller can leave it untouched (preserving numeric sorting).
#
# Cost control, in order of importance:
#   1. grepl() (cheap) decides which cells match; the far more expensive
#      gregexpr() then runs ONLY on those cells, not on the whole column.
#   2. Cells containing no &, < or > need no escaping at all, so htmlEscape()
#      is kept out of the hot path for essentially all real annotation values.
#   3. A column with no match and no escapable character is returned as-is.
.kws_highlight_col <- function(v, needle_lc, must_escape) {
  x  <- as.character(v)
  na <- is.na(x)
  if (any(na)) x[na] <- ""

  has  <- grepl(needle_lc, tolower(x), fixed = TRUE)
  spec <- grepl("[&<>]", x, useBytes = TRUE)
  any_has <- any(has); any_spec <- any(spec)

  # Nothing to mark and nothing to escape.
  if (!any_has && !any_spec) return(if (must_escape) x else NULL)

  out <- x
  # Cells that need escaping but carry no match.
  esc_only <- spec & !has
  if (any(esc_only)) out[esc_only] <- htmltools::htmlEscape(x[esc_only])

  if (any_has) {
    idx <- which(has)
    mm  <- gregexpr(needle_lc, tolower(x[idx]), fixed = TRUE)
    sp  <- spec[idx]
    xs  <- x[idx]
    out[idx] <- vapply(seq_along(idx),
                       function(i) .kws_hl_one(xs[i], mm[[i]], sp[i]),
                       character(1))
  }
  if (any(na)) out[na] <- ""
  out
}

# Public helper: display-only highlighted copy of `df`.
#   df       : the ORIGINAL result data frame (unmodified by this function)
#   keyword  : the searched keyword
#   id_cols  : leading technical-ID column name(s) -- escaped, never highlighted
# Returns `df` unchanged when there is nothing to do.
build_keyword_highlight_df <- function(df, keyword, id_cols = character(0)) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L || ncol(df) == 0L) return(df)
  kw <- trimws(as.character(keyword)[1L])
  if (is.na(kw) || !nzchar(kw)) return(df)
  needle_lc <- tolower(kw)

  id_cols  <- as.character(id_cols)
  ann_cols <- setdiff(colnames(df), id_cols)
  out <- df

  for (cn in ann_cols) {
    v  <- df[[cn]]
    hl <- .kws_highlight_col(v, needle_lc,
                             must_escape = is.character(v) || is.factor(v))
    if (!is.null(hl)) out[[cn]] <- hl
  }
  # ID column(s): escaped so escape = FALSE is safe, but never highlighted.
  for (cn in intersect(colnames(out), id_cols))
    out[[cn]] <- htmltools::htmlEscape(as.character(out[[cn]]))

  out
}


# ---------------------------------------------------------------------------
# KEYWORD SEARCH over the canonical feature_annot
# ---------------------------------------------------------------------------
# Reverse lookup: annotation VALUE -> feature ID. Every other annotation helper
# above goes the other way (ID -> annotation rows).
#
# Design rules (do not relax):
#   * The searchable haystack is built from ALL columns of feature_annot. No
#     column is ever excluded — not by name, not by value overlap, and NOT
#     because its values currently happen to equal the feature IDs (e.g. the
#     metabolomics `NameID` column). A future payload may put different values
#     in the same column, so any such exclusion would be a latent bug.
#   * rownames(feature_annot) — the technical feature IDs — are NOT searchable.
#     They are absent from the haystack simply because they are not a column.
#     Users who already know an ID use the Explorer's selectize instead.
#   * Matching is case-insensitive LITERAL substring: both sides are lowercased
#     and grepl() is called with fixed = TRUE. Note that R *ignores*
#     ignore.case = TRUE when fixed = TRUE, so pre-lowercasing is the only
#     correct construction. fixed = TRUE is also what makes "[M+H]+1",
#     "C18 H35 N O5", "3-[(3-Hydroxyundecanoyl)oxy]" and KEGG/GO IDs match
#     literally instead of being parsed as regular expressions.
#   * useBytes stays at its default FALSE so UTF-8 works (beta-Estradiol, F1alpha,
#     (3beta,24R,24'R)-... are real values in the metabolomics payload).

# Coerce one feature_annot column to a plain character vector of length n.
# Handles factors (labels, not integer codes), numerics, logicals and list
# columns (flattened so multi-value cells stay searchable). Returns NULL when
# the column cannot be coerced or has the wrong length, in which case the
# caller skips it rather than failing the whole search.
.annot_col_to_chr <- function(x, n) {
  out <- tryCatch({
    if (is.list(x) && !is.data.frame(x)) {
      vapply(x, function(z) paste(as.character(unlist(z)), collapse = " "),
             character(1))
    } else {
      as.character(x)
    }
  }, error = function(e) NULL)
  if (is.null(out) || length(out) != n) return(NULL)
  out
}

# Build a lowercase search index over every column of a feature_annot table.
# Returns list(ids, hay, n_searchable_cols, searchable_cols) where `hay` holds
# one collapsed lowercase string per feature (columns joined with the unit-
# separator control character, which no keyword can contain, so a match can
# never span two columns).
#
# Blank-ish values (NA, "", "NA"/"N/A"/"null", separator-only strings such as
# "|" or "-") are masked to "" via .is_blank_annot_value() BEFORE indexing.
# Without this, searching "-" or "na" would match thousands of rows — the
# active proteomics payload has two columns that are literally " | " for all
# 4103 features. A column left entirely blank contributes nothing and is not
# counted in n_searchable_cols.
build_feature_annot_search_index <- function(fa) {
  ids <- if (is.null(fa)) NULL else rownames(fa)
  ids <- if (is.null(ids)) character(0) else as.character(ids)
  n   <- length(ids)

  no_cols <- function() list(ids = ids,
                             hay = if (n > 0L) rep("", n) else character(0),
                             n_searchable_cols = 0L,
                             searchable_cols   = character(0))

  if (is.null(fa) || !is.data.frame(fa) || n == 0L || ncol(fa) == 0L) {
    return(no_cols())
  }

  cols     <- list()
  ok_names <- character(0)
  for (cn in colnames(fa)) {
    v <- .annot_col_to_chr(fa[[cn]], n)
    if (is.null(v)) next
    v <- trimws(v)
    v[.is_blank_annot_value(v)] <- ""
    if (!any(nzchar(v))) next
    cols[[length(cols) + 1L]] <- v
    ok_names <- c(ok_names, cn)
  }

  if (length(cols) == 0L) return(no_cols())

  hay <- if (length(cols) == 1L) cols[[1L]]
         else do.call(paste, c(cols, list(sep = "\u001f")))  # unit separator

  list(ids               = ids,
       hay               = tolower(hay),
       n_searchable_cols = length(ok_names),
       searchable_cols   = ok_names)
}

# Case-insensitive literal substring search against a prebuilt index.
# Returns the matching feature IDs in rownames(feature_annot) order (which is
# what build_feature_annot_table() then preserves). Never errors.
search_feature_annot_index <- function(index, keyword) {
  if (is.null(index) || is.null(index$hay) || length(index$hay) == 0L) {
    return(character(0))
  }
  kw <- trimws(as.character(keyword)[1L])
  if (is.na(kw) || !nzchar(kw)) return(character(0))
  hit <- tryCatch(grepl(tolower(kw), index$hay, fixed = TRUE),
                  error = function(e) rep(FALSE, length(index$hay)))
  hit[is.na(hit)] <- FALSE
  index$ids[which(hit)]
}

# Availability predicate for the Keyword Search sub-tab, in the style of the
# hc_data_available_env() / partition_data_available_env() family. Bundles the
# annotation table, the search index and the Explorer ID universe so the wiring
# layer can cache all of it in a single lazy reactive(). Never errors.
#
# status:
#   "no_env"               -> no usable omics environment
#   "no_annot"             -> feature_annot missing / NULL / not a df / zero rows
#   "no_columns"           -> valid table carrying IDs but zero annotation columns
#                             (the active RNAseq payload: 4928 x 0)
#   "no_searchable_values" -> columns exist but every value is blank/invalid
#   "no_id_universe"       -> no expr_raw/expr_norm rownames to open features with
#   "unlinked"             -> annotation rows share no IDs with the expression data
#   "ok"                   -> searchable and features can be opened in the Explorers
keyword_search_status <- function(env) {
  out <- list(status = "no_env", fa = NULL, index = NULL, id_universe = NULL,
              n_linked = 0L, n_searchable_cols = 0L)
  if (is.null(env) || !is.environment(env)) return(out)

  fa <- get_feature_annot_df(env)
  st <- feature_annot_status(fa)
  if (identical(st, "missing")) { out$status <- "no_annot"; return(out) }
  out$fa <- fa
  if (identical(st, "no_columns")) { out$status <- "no_columns"; return(out) }

  idx <- tryCatch(build_feature_annot_search_index(fa), error = function(e) NULL)
  if (is.null(idx)) { out$status <- "no_annot"; return(out) }
  out$index             <- idx
  out$n_searchable_cols <- idx$n_searchable_cols
  if (idx$n_searchable_cols == 0L) {
    out$status <- "no_searchable_values"; return(out)
  }

  # Explorer ID universe: same expr_raw -> expr_norm fallback the Explorers use.
  ids <- NULL
  if (exists("expr_raw", envir = env, inherits = FALSE)) {
    ids <- prepare_id_choices_from_matrix(
      get("expr_raw", envir = env, inherits = FALSE))
  }
  if ((is.null(ids) || length(ids) == 0L) &&
      exists("expr_norm", envir = env, inherits = FALSE)) {
    ids <- prepare_id_choices_from_matrix(
      get("expr_norm", envir = env, inherits = FALSE))
  }
  out$id_universe <- ids
  if (is.null(ids) || length(ids) == 0L) {
    out$status <- "no_id_universe"; return(out)
  }

  out$n_linked <- length(intersect(rownames(fa), ids))
  if (out$n_linked == 0L) { out$status <- "unlinked"; return(out) }

  out$status <- "ok"
  out
}

# Output / input IDs owned by one Keyword Search panel. Single source of truth
# shared by keyword_search_panel_ui() (which creates the inputs and output
# placeholders) and wire_keyword_search() (which fills them).
.keyword_search_ids <- function(id_prefix) {
  list(
    keyword      = paste0(id_prefix, "_keyword"),
    search_btn   = paste0(id_prefix, "_search_btn"),
    clear_btn    = paste0(id_prefix, "_clear_btn"),
    available    = paste0(id_prefix, "_available"),
    avail_ui     = paste0(id_prefix, "_availability_ui"),
    status_ui    = paste0(id_prefix, "_status_ui"),
    has_results  = paste0(id_prefix, "_has_results"),
    dt           = paste0(id_prefix, "_results_dt"),
    rows_sel     = paste0(id_prefix, "_results_dt_rows_selected"),
    instruction_ui = paste0(id_prefix, "_instruction_ui"),
    selection_ui = paste0(id_prefix, "_selection_ui"),
    action_ui    = paste0(id_prefix, "_action_ui"),
    open_single  = paste0(id_prefix, "_open_single_btn"),
    open_multi   = paste0(id_prefix, "_open_multi_btn"),
    dl_ui        = paste0(id_prefix, "_dl_ui"),
    download     = paste0(id_prefix, "_download")
  )
}

# UI for one "Keyword Search" sub-tab. Placed here (rather than inlined three
# times in ui.R) for the same reason as build_de_cutoffs_ui() /
# de_table_legend_html(): the markup is identical across omics. The *call*
# still lives in ui.R.
#
# `id_prefix` is the only input this builder needs: the panel contains no
# omics-specific literals. Every user-visible label that varies by omics
# (feature names, Explorer tab names, status and fallback wording) is rendered
# server-side by wire_keyword_search() into the uiOutput placeholders below.
#
# Two nested gates:
#   output.<prefix>_available   == 'yes'  -> the dataset can be searched at all
#   output.<prefix>_has_results == 'yes'  -> a search returned >= 1 row
keyword_search_panel_ui <- function(id_prefix) {
  ids <- .keyword_search_ids(id_prefix)

  shiny::fluidRow(
    shiny::column(
      width = 12,

      # Keyword Search styles, kept with the panel (rather than in ui.R) so the
      # feature stays self-contained and travels as one unit; both rules are
      # idempotent when several instances of this panel are on the page.
      #   mark.kws-hit  -- highlight for matched substrings in the results table
      #   .kws-results  -- breathing room between the status line
      #                    ("N matched ...") and DT's "Show N entries" selector.
      #                    Scoped to the results block, so the initial /
      #                    too-short / no-results states are untouched.
      shiny::tags$style(shiny::HTML(paste0(
        "mark.kws-hit{background-color:#FFF3A3;color:inherit;padding:0 1px;border-radius:2px;}",
        ".kws-results{margin-top:14px;}"
      ))),

      shiny::h4("Keyword Search", class = "plot-title full"),

      # Unavailable-state banner (also carries the "cannot open in the
      # Explorers" note when the data is searchable but not linkable).
      shiny::uiOutput(ids$avail_ui),

      shiny::conditionalPanel(
        condition = paste0("output.", ids$available, " == 'yes'"),

        # Label above the flex row so the box and buttons share one line.
        shiny::tags$label(
          "Search Annotation:", `for` = ids$keyword,
          style = paste("font-size:18px; font-weight:600; display:block;",
                        "margin-top:8px; margin-bottom:5px;")),
        shiny::div(
          style = "display:flex; gap:8px; align-items:center; flex-wrap:wrap;",
          htmltools::tagAppendAttributes(
            shiny::textInput(
              inputId     = ids$keyword,
              label       = NULL,
              value       = "",
              width       = "420px",
              placeholder = "e.g. a compound name, description, KEGG ID..."
            ),
            style = "margin-bottom: 0;"
          ),
          shiny::actionButton(ids$search_btn, "Search",
                              icon  = shiny::icon("magnifying-glass"),
                              class = "btn btn-primary btn-sm"),
          shiny::actionButton(ids$clear_btn, "Clear",
                              icon  = shiny::icon("trash"),
                              class = "btn btn-default btn-sm")
        ),
        # Status line is the ONLY message under the controls: the initial
        # "Enter a keyword..." prompt before any search, then the result count /
        # no-results / minimum-characters message afterwards. Deliberately no
        # static explanatory paragraph here.
        shiny::uiOutput(ids$status_ui),

        shiny::conditionalPanel(
          condition = paste0("output.", ids$has_results, " == 'yes'"),

          # §3: shared Explorer-handoff instruction ABOVE the table (same helper +
          # styling as the Enrichment gene tables). Rendered only when hits are
          # transferable to the Explorers.
          shiny::uiOutput(ids$instruction_ui),
          # .kws-results carries the top margin that separates the status line
          # from DT's page-length selector (see the style block above).
          shiny::div(class = "kws-results", DT::DTOutput(ids$dt)),

          # Selection hint (left) + the single context-appropriate action
          # button (right). Same flex idiom as the heatmap-zoom blocks.
          shiny::div(
            style = paste("display:flex; gap:14px; align-items:center;",
                          "flex-wrap:wrap; padding:8px 0; margin-top:4px;",
                          "border-top:1px solid #eee;"),
            shiny::div(style = "flex:1 1 auto;", shiny::uiOutput(ids$selection_ui)),
            shiny::div(style = "flex:0 0 auto;", shiny::uiOutput(ids$action_ui))
          ),

          shiny::uiOutput(ids$dl_ui),
          shiny::br()
        )
      )
    )
  )
}

# ===========================================================================
# Shared "selected rows -> Gene / Multi-Gene Explorer" hand-off (Stage 4F).
# Extracted VERBATIM from wire_keyword_search() so the Enrichment Member-Genes /
# Core-Genes tables reuse the IDENTICAL interaction (one action button that is
# 0 selected -> none, 1 -> "Open in <single>", 2+ -> "Open in <multi>"; the
# keep-openable skip toast; the Multi-Gene cap + warning; and the load-bearing
# priority = -100 deferred single-selectize write). Both wire_keyword_search()
# and wire_enrichment() call this, so the two features cannot drift apart.
#
#   action_id      : output id for the single action button (0->none/1/2+)
#   open_single_id : actionButton input id -> single Explorer
#   open_multi_id  : actionButton input id -> Multi Explorer
#   selected_ids   : reactive() -> stable selected IDs (character; maps the DT
#                    *_rows_selected to the data frame, NOT the visible page, so it
#                    is correct after sort/search/paging)
#   id_universe    : reactive() -> openable IDs (expr-matrix rownames)
#   has_table      : reactive() -> TRUE when the source table is present
#   transfer_ok    : reactive() -> TRUE when the hand-off is possible
#   enabled        : static logical; FALSE -> render nothing, register no observers
#   tabs_id / single_tab_label / single_selectize_id / multi_tab_label /
#   multi_controller : the same nav knobs wire_keyword_search uses.
wire_explorer_handoff <- function(session, input, output,
                                  action_id, open_single_id, open_multi_id,
                                  selected_ids, id_universe,
                                  has_table   = function() TRUE,
                                  transfer_ok = function() TRUE,
                                  enabled, tabs_id,
                                  single_tab_label, single_selectize_id,
                                  multi_tab_label, multi_controller,
                                  feature_label        = "gene",
                                  feature_label_plural = "genes") {
  # Exactly ONE navigation button is ever in the DOM:
  #   0 selected -> none, 1 -> single Explorer, 2+ -> Multi Explorer.
  output[[action_id]] <- shiny::renderUI({
    if (!isTRUE(has_table())) return(NULL)
    if (!isTRUE(enabled) || !isTRUE(transfer_ok())) return(NULL)
    n <- length(selected_ids())
    if (n == 0L) return(NULL)
    if (n == 1L) {
      shiny::actionButton(open_single_id, paste0("Open in ", single_tab_label),
                          icon  = shiny::icon("up-right-from-square"),
                          class = "btn btn-primary btn-sm")
    } else {
      shiny::actionButton(open_multi_id, paste0("Open in ", multi_tab_label),
                          icon  = shiny::icon("up-right-from-square"),
                          class = "btn btn-primary btn-sm")
    }
  })
  if (!isTRUE(enabled)) return(invisible(NULL))

  pending_single <- shiny::reactiveVal(NULL)

  # Drop IDs that are not in the Explorer ID universe; toast about any skipped.
  keep_openable <- function(sel) {
    universe <- id_universe()
    if (is.null(universe)) return(character(0))
    kept    <- sel[sel %in% universe]
    skipped <- length(sel) - length(kept)
    if (skipped > 0L) {
      shiny::showNotification(
        sprintf("%d selected %s %s not present in the expression data and %s skipped.",
                skipped, if (skipped == 1L) feature_label else feature_label_plural,
                if (skipped == 1L) "is" else "are",
                if (skipped == 1L) "was" else "were"),
        type = "warning", duration = 5)
    }
    kept
  }

  shiny::observeEvent(input[[open_single_id]], {
    if (!isTRUE(has_table())) return()
    sel <- selected_ids()
    if (length(sel) == 0L) {
      shiny::showNotification("Select a row first.",
                              type = "message", duration = 4); return()
    }
    if (length(sel) > 1L) {
      shiny::showNotification(
        sprintf("Select exactly one row to open the %s.", single_tab_label),
        type = "warning", duration = 4); return()
    }
    sel <- keep_openable(sel)
    if (length(sel) == 0L) return()
    # Switch tabs first; the selectize write lands in the post-nav observer below,
    # which fires AFTER the single Explorer's own choices resync.
    pending_single(sel[1])
    bslib::nav_select(id = tabs_id, selected = single_tab_label)
  }, ignoreInit = TRUE)

  # Phase 2: write the selectize AFTER the single Explorer's own choices resync,
  # which fires on the same tab-change flush and resets `selected` to character(0).
  # priority = -100 is load-bearing (see wire_keyword_search). Each hand-off guards
  # its OWN pending_single, so several of these observers can coexist on one
  # tabs_id (keyword search + ORA + GSEA) and only the triggered one writes.
  shiny::observeEvent(input[[tabs_id]], {
    if (!identical(input[[tabs_id]], single_tab_label)) return()
    val <- pending_single()
    if (is.null(val) || !nzchar(as.character(val))) return()
    shiny::updateSelectizeInput(session, single_selectize_id,
                                choices  = id_universe(),
                                selected = val,
                                server   = TRUE)
    pending_single(NULL)
  }, ignoreInit = TRUE, priority = -100)

  shiny::observeEvent(input[[open_multi_id]], {
    if (!isTRUE(has_table())) return()
    sel <- selected_ids()
    if (length(sel) == 0L) {
      shiny::showNotification("Select at least one row first.",
                              type = "message", duration = 4); return()
    }
    sel <- keep_openable(sel)
    if (length(sel) == 0L) return()
    max_n <- multi_controller$max_n
    if (!is.null(max_n) && length(sel) > max_n) {
      shiny::showNotification(
        sprintf("You selected %d %s; only the first %d will be opened.",
                length(sel), feature_label_plural, max_n),
        type = "warning", duration = 5)
      sel <- sel[seq_len(max_n)]
    }
    multi_controller$set(sel)
    bslib::nav_select(id = tabs_id, selected = multi_tab_label)
  }, ignoreInit = TRUE)

  invisible(NULL)
}

# Generic server wiring for one Keyword Search sub-tab.
#
#   env                 : omics env carrying feature_annot + expr_raw/expr_norm
#   id_prefix           : e.g. "metab_kws" (see .keyword_search_ids())
#   feature_col_name    : leading ID column label ("Gene"/"Protein"/"Metabolite")
#   tabs_id             : the omics navset_pill_list id ("metab_tabs")
#   single_tab_label    : nav_panel title of the single Explorer
#   single_selectize_id : that Explorer's selectize input id
#   multi_tab_label     : nav_panel title of the Multi Explorer
#   multi_controller    : the list() returned by wire_multi_gene_explorer()
#   min_chars           : keywords shorter than this are rejected (default 2)
#   warn_hits           : advisory toast above this many hits; results are NEVER
#                         capped or hidden
#
# The search index is built lazily on first use and cached for the session, so
# app start-up is unaffected.
wire_keyword_search <- function(session, input, output, env,
                                id_prefix,
                                feature_label        = "feature",
                                feature_label_plural = NULL,
                                feature_col_name     = "Feature",
                                tabs_id              = NULL,
                                single_tab_label     = NULL,
                                single_selectize_id  = NULL,
                                multi_tab_label      = NULL,
                                multi_controller     = NULL,
                                min_chars            = 2L,
                                warn_hits            = 500L) {
  ids    <- .keyword_search_ids(id_prefix)
  plural <- if (is.null(feature_label_plural)) paste0(feature_label, "s")
            else feature_label_plural

  nav_enabled <- !is.null(tabs_id) && nzchar(as.character(tabs_id)) &&
                 !is.null(single_tab_label) && !is.null(single_selectize_id) &&
                 !is.null(multi_tab_label)  && !is.null(multi_controller)

  # ---- lazy, session-cached status + index + ID universe --------------------
  ks_state <- shiny::reactive({
    keyword_search_status(env)
  })

  # Searchable = there is at least one annotation value to match against.
  # Transferable = hits can additionally be opened in the Explorers.
  is_searchable <- shiny::reactive({
    identical(ks_state()$status, "ok") ||
      identical(ks_state()$status, "no_id_universe") ||
      identical(ks_state()$status, "unlinked")
  })
  can_transfer <- shiny::reactive({
    nav_enabled && identical(ks_state()$status, "ok")
  })

  output[[ids$available]] <- shiny::renderText({
    if (isTRUE(is_searchable())) "yes" else "no"
  })
  shiny::outputOptions(output, ids$available, suspendWhenHidden = FALSE)

  output[[ids$avail_ui]] <- shiny::renderUI({
    st <- ks_state()$status
    # One wording for every "there is no annotation to search" state: the table
    # is absent, has zero annotation columns, or every value in it is blank.
    # The distinction is meaningless to the user, so do not re-split these.
    msg_no_annot <- paste("Keyword search is not available because no annotation",
                          "information is available for this dataset.")
    msg <- switch(
      st,
      "no_env"               = "Keyword search is not available for this dataset.",
      "no_annot"             = msg_no_annot,
      "no_columns"           = msg_no_annot,
      "no_searchable_values" = msg_no_annot,
      NULL
    )
    if (!is.null(msg)) {
      return(shiny::div(
        style = paste("padding:12px; margin-top:8px; border-radius:8px;",
                      "background:#EDF0F1; color:#2c3e50; font-size:14px;"),
        msg))
    }
    # Searchable, but hits cannot be handed to the Explorers.
    if (st %in% c("no_id_universe", "unlinked")) {
      return(shiny::div(
        style = paste("padding:10px; margin-top:8px; border-radius:8px;",
                      "background:#FCF3CF; color:#7D6608; font-size:13px;"),
        sprintf(paste("Annotation rows for this dataset cannot be linked to the",
                      "expression data, so results cannot be opened in the %s or",
                      "the %s."),
                if (is.null(single_tab_label)) "Explorer" else single_tab_label,
                if (is.null(multi_tab_label)) "Multi Explorer" else multi_tab_label)))
    }
    NULL
  })

  # ---- search ---------------------------------------------------------------
  # NULL = nothing searched yet. Otherwise list(kw, ids, state) where state is
  # "empty" / "short" / "none" / "hits".
  results_rv <- shiny::reactiveVal(NULL)

  shiny::observeEvent(input[[ids$search_btn]], {
    kw <- trimws(as.character(input[[ids$keyword]]))
    if (length(kw) == 0L || is.na(kw)) kw <- ""

    if (!nzchar(kw)) {
      results_rv(list(kw = "", ids = character(0), state = "empty")); return()
    }
    if (nchar(kw) < min_chars) {
      results_rv(list(kw = kw, ids = character(0), state = "short")); return()
    }

    st  <- ks_state()
    hit <- search_feature_annot_index(st$index, kw)

    # Advisory only -- the result table always shows every match.
    if (length(hit) > warn_hits) {
      shiny::showNotification(
        sprintf("Your search matched %d %s. Consider a more specific keyword.",
                length(hit), plural),
        type = "warning", duration = 6)
    }
    results_rv(list(kw = kw, ids = hit,
                    state = if (length(hit) > 0L) "hits" else "none"))
  }, ignoreInit = TRUE)

  shiny::observeEvent(input[[ids$clear_btn]], {
    shiny::updateTextInput(session, ids$keyword, value = "")
    results_rv(NULL)
  }, ignoreInit = TRUE)

  # ---- result table --------------------------------------------------------
  # Returns NULL (not a req() failure) when there is nothing to show, so the
  # action observers can distinguish "no table" from "no selection".
  results_df <- shiny::reactive({
    r <- results_rv()
    if (is.null(r) || length(r$ids) == 0L) return(NULL)
    fa <- ks_state()$fa
    if (is.null(fa)) return(NULL)
    tbl <- build_feature_annot_table(fa, r$ids,
                                     feature_col_name = feature_col_name)
    if (nrow(tbl) == 0L) return(NULL)
    # Same cleanup the Explorer Annotation sections apply: hide columns that are
    # entirely blank across the matched rows.
    drop_blank_annot_columns(tbl, keep_cols = feature_col_name)
  })

  output[[ids$has_results]] <- shiny::renderText({
    if (!is.null(results_df())) "yes" else "no"
  })
  shiny::outputOptions(output, ids$has_results, suspendWhenHidden = FALSE)

  output[[ids$status_ui]] <- shiny::renderUI({
    r <- results_rv()
    note <- function(txt, col = "#888") shiny::div(
      style = sprintf("color:%s; font-size:14px; margin-top:10px;", col), txt)

    if (is.null(r)) {
      return(note(sprintf("Enter a keyword to search the annotation table.")))
    }
    if (identical(r$state, "empty")) {
      return(note("Enter a keyword to search the annotation table."))
    }
    if (identical(r$state, "short")) {
      return(note(sprintf("Enter at least %d characters to search.", min_chars)))
    }
    if (identical(r$state, "none")) {
      return(note(sprintf('No %s matched "%s".', plural, r$kw), "#7D6608"))
    }
    n <- length(r$ids)
    shiny::div(
      style = "font-size:14px; margin-top:10px; color:#2c3e50; font-weight:600;",
      sprintf('%d %s matched "%s".', n,
              if (n == 1L) feature_label else plural, r$kw))
  })

  # DISPLAY-ONLY copy: same rows, same order, same count as results_df(), with
  # every cell HTML-escaped and each keyword occurrence wrapped in <mark>.
  # Recomputed only when a new search lands (results_rv changes) -- not while
  # typing, not on tab open, and not on sorting / paging / row selection, none
  # of which this reactive depends on.
  display_df <- shiny::reactive({
    df <- results_df()
    if (is.null(df)) return(NULL)
    r <- results_rv()
    build_keyword_highlight_df(df, if (is.null(r)) "" else r$kw,
                               id_cols = feature_col_name)
  })

  output[[ids$dt]] <- DT::renderDT({
    df <- display_df()
    shiny::req(!is.null(df))
    # searching = FALSE: no DT search field here (the panel has its own keyword
    # box). Sorting / paging / scrollX / row selection all stay enabled.
    # length_menu: page-length selector (10/20/50, default 10) -> dom "ltip".
    # escape = FALSE is safe ONLY because build_keyword_highlight_df() has
    # already escaped every cell and injected just the <mark> markup.
    render_feature_annot_dt(df, selection = "multiple", searching = FALSE,
                            escape = FALSE, length_menu = c(10, 20, 50))
  })

  # ---- selection ----------------------------------------------------------
  # DT's *_rows_selected indexes the data frame handed to renderDT, NOT the
  # sorted / filtered / paged view, so this mapping stays correct after any
  # combination of those.
  selected_ids <- shiny::reactive({
    df <- results_df()
    if (is.null(df)) return(character(0))
    sel <- input[[ids$rows_sel]]
    if (is.null(sel) || length(sel) == 0L) return(character(0))
    sel <- sel[!is.na(sel) & sel >= 1L & sel <= nrow(df)]
    if (length(sel) == 0L) return(character(0))
    out <- as.character(df[[feature_col_name]][sel])
    out <- out[!is.na(out) & nzchar(out)]
    unique(out)   # defensive: rownames are unique, so this is a no-op today
  })

  # §3: the shared one-line Explorer-handoff instruction (same helper + styling as
  # the Enrichment gene tables). Only when hits can actually be opened.
  output[[ids$instruction_ui]] <- shiny::renderUI({
    if (!isTRUE(can_transfer()) || is.null(results_df())) return(NULL)
    .explorer_instruction(single_tab_label, multi_tab_label, feature_label)
  })
  output[[ids$selection_ui]] <- shiny::renderUI({
    df <- results_df()
    shiny::req(!is.null(df))
    n <- length(selected_ids())
    if (n == 0L) {
      return(shiny::div(
        style = "color:#888; font-size:13px;",
        sprintf("Select one or more rows to open them in an Explorer.")))
    }
    shiny::div(
      style = "font-size:13px; color:#2c3e50; font-weight:600;",
      sprintf("Selected: %d of %d %s", n, nrow(df), plural))
  })

  # The selected-rows -> Explorer action button + navigation now live in the
  # shared wire_explorer_handoff() (same code the Enrichment tables use). Behavior
  # is unchanged: 0 selected -> none, 1 -> single Explorer, 2+ -> Multi Explorer;
  # keep-openable skip toast; Multi cap + warning; priority = -100 deferred write.
  wire_explorer_handoff(
    session, input, output,
    action_id      = ids$action_ui,
    open_single_id = ids$open_single,
    open_multi_id  = ids$open_multi,
    selected_ids   = selected_ids,
    id_universe    = function() ks_state()$id_universe,
    has_table      = function() !is.null(results_df()),
    transfer_ok    = function() isTRUE(can_transfer()),
    enabled        = nav_enabled,
    tabs_id        = tabs_id,
    single_tab_label    = single_tab_label,
    single_selectize_id = single_selectize_id,
    multi_tab_label     = multi_tab_label,
    multi_controller    = multi_controller,
    feature_label        = feature_label,
    feature_label_plural = plural)

  # ---- Excel export of exactly what is displayed --------------------------
  output[[ids$download]] <- shiny::downloadHandler(
    filename = function() {
      paste0(id_prefix, "_keyword_search_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      df <- results_df()
      writexl::write_xlsx(if (is.null(df)) data.frame() else df, path = file)
    }
  )

  output[[ids$dl_ui]] <- shiny::renderUI({
    shiny::req(!is.null(results_df()))
    shiny::div(
      style = "display:flex; justify-content:flex-end; margin-top:8px;",
      shiny::downloadButton(ids$download, "Download Results")
    )
  })

  # (navigation now wired by wire_explorer_handoff() above)

  invisible(list(current = function() selected_ids()))
}

# Wire Gene Explorer (expects NEW names in env: expr_raw, expr_norm, sample_meta, de_stats, contrasts, padj_cutoff, fc_cutoff)
# override_id_reactive: optional reactive() that supplies the gene ID directly (used in MultiOmics tab).
#   When provided, the selectize input and click-routing logic are bypassed entirely.
wire_gene_explorer <- function(session, input, output, env, tabs_id = "rna_tabs",
                               click_sources = c("rna_de_volcano", "rna_de_ma"),
                               id_prefix = "gene_explorer",
                               override_id_reactive = NULL,
                               feature_label = NULL,
                               feature_col_name = "Gene") {

  # Derive all input/output IDs from the prefix
  id_gene           <- paste0(id_prefix, "_gene")
  id_scale          <- paste0(id_prefix, "_scale")
  id_color_by       <- paste0(id_prefix, "_color_by")
  id_padj_cutoff    <- paste0(id_prefix, "_padj_cutoff")  # single significance cutoff input
  id_fc_cutoff      <- paste0(id_prefix, "_fc_cutoff")
  id_sig_mode       <- paste0(id_prefix, "_sig_mode")
  id_sel_title      <- paste0(id_prefix, "_selected_gene_title")
  id_plot_title     <- paste0(id_prefix, "_plot_title")
  id_counts_plot    <- paste0(id_prefix, "_counts_plot")
  id_data_table     <- paste0(id_prefix, "_data_table")
  id_de_table       <- paste0(id_prefix, "_de_table")
  id_cutoffs_ui     <- paste0(id_prefix, "_cutoffs_ui")
  id_de_legend_ui   <- paste0(id_prefix, "_de_legend_ui")
  id_dl_plot        <- paste0("download_", id_prefix, "_plot")
  id_dl_data        <- paste0("download_", id_prefix, "_data")
  id_view_btn       <- paste0(id_prefix, "_view_clicked")

  # Selectize input and click-routing are only used when there is no external override
  if (is.null(override_id_reactive)) {
    gene_choices_all <- NULL
    if (exists("expr_raw", envir = env, inherits = FALSE)) {
      gene_choices_all <- prepare_id_choices_from_matrix(get("expr_raw", envir = env, inherits = FALSE))
    } else if (exists("expr_norm", envir = env, inherits = FALSE)) {
      gene_choices_all <- prepare_id_choices_from_matrix(get("expr_norm", envir = env, inherits = FALSE))
    }

    ge_choices_ready <- shiny::reactiveVal(FALSE)
    ge_pending_gene  <- shiny::reactiveVal(NULL)
    ge_clicked_gene  <- shiny::reactiveVal(NULL)

    ge_sync_choices <- function() {
      if (is.null(gene_choices_all) || length(gene_choices_all) == 0) return(FALSE)

      cur  <- shiny::isolate(input[[id_gene]])
      keep <- if (!is.null(cur) && nzchar(cur) && cur %in% gene_choices_all) cur else character(0)

      shiny::updateSelectizeInput(
        session,
        inputId  = id_gene,
        choices  = gene_choices_all,
        server   = TRUE,
        selected = keep
      )
      ge_choices_ready(TRUE)
      TRUE
    }

    session$onFlushed(function() ge_sync_choices(), once = TRUE)

    shiny::observeEvent(input[[tabs_id]], {
      if (identical(input[[tabs_id]], "Gene Explorer")) ge_sync_choices()
    }, ignoreInit = TRUE)
  }

  shiny::observe({
    shiny::validate(
      shiny::need(exists("color", envir = env, inherits = FALSE), "color not available.")
    )

    color_col <- get("color", envir = env, inherits = FALSE)
    shape_col <- if (exists("shape", envir = env, inherits = FALSE)) get("shape", envir = env, inherits = FALSE) else NULL

    choices_vec <- unique(na.omit(c(color_col, shape_col)))

    shiny::updateSelectInput(
      session,
      id_scale,
      choices  = get_scale_choices_from_env(env),
      selected = "norm"
    )

    shiny::updateSelectInput(
      session,
      id_color_by,
      choices  = choices_vec,
      selected = color_col
    )
  })

  # Single source of truth for the selected gene ID.
  # When override_id_reactive is provided it is the true source; otherwise the selectize input is.
  selected_gene_id <- shiny::reactive({
    if (!is.null(override_id_reactive)) {
      val <- override_id_reactive()
      shiny::req(!is.null(val) && nzchar(trimws(as.character(val))))
      return(trimws(as.character(val)))
    }
    shiny::req(input[[id_gene]])
    trimws(as.character(input[[id_gene]]))
  })

  output[[id_sel_title]] <- shiny::renderUI({
    gene_id <- selected_gene_id()
    if (!is.null(feature_label)) {
      shiny::tags$div(
        style = "text-align:center; margin-top:6px; margin-bottom:2px;",
        shiny::tags$span(paste0(feature_label, ": "),
                         style = "color:#777; font-size:16px; font-weight:400;"),
        shiny::tags$span(gene_id,
                         style = "color:#18BC9C; font-size:24px; font-weight:800;")
      )
    } else {
      shiny::tags$div(gene_id, style = "color:#18BC9C;font-size:28px;font-weight:800;margin-top:6px;margin-bottom:2px;text-align:center;")
    }
  })

  output[[id_plot_title]] <- shiny::renderUI({
    selected_gene_id()  # req() guard
    shiny::h4("Expression Profile", class = "plot-title full")
  })

  ge_data <- shiny::reactive({
    shiny::validate(
      shiny::need(exists("expr_raw", envir = env, inherits = FALSE), "expr_raw not available."),
      shiny::need(exists("sample_meta", envir = env, inherits = FALSE), "sample_meta not available.")
    )

    gene_id   <- selected_gene_id()
    color_by  <- shiny::req(input[[id_color_by]])
    scale_sel <- shiny::req(input[[id_scale]])

    expr_raw    <- get("expr_raw", envir = env, inherits = FALSE)
    sample_meta <- get("sample_meta", envir = env, inherits = FALSE)

    expr_norm <- NULL
    if (scale_sel == "log2") {
      shiny::validate(shiny::need(exists("expr_norm", envir = env, inherits = FALSE), "expr_norm not available."))
      expr_norm <- get("expr_norm", envir = env, inherits = FALSE)
    } else if (exists("expr_norm", envir = env, inherits = FALSE)) {
      expr_norm <- get("expr_norm", envir = env, inherits = FALSE)
    }

    get_explorer_df_from_mats(
      feature_id  = gene_id,
      scale_sel   = scale_sel,
      color_by    = color_by,
      expr_raw    = expr_raw,
      expr_norm   = expr_norm,
      sample_meta = sample_meta
    )
  })

  ge_plot <- shiny::reactive({
    df <- ge_data()
    make_explorer_plot(
      df            = df,
      color_by      = shiny::req(input[[id_color_by]]),
      scale_sel     = shiny::req(input[[id_scale]]),
      ylab_log2     = get_norm_method_label_from_env(env),
      theme_variant = "gene"
    )
  })

  output[[id_counts_plot]] <- shiny::renderPlot({ ge_plot() })

  # Expression Data: wide-format table (one row for the selected feature,
  # one column per sample). Reuses the shared helpers so styling is
  # byte-identical to the Multi-Gene Explorer table.
  ge_table_df <- shiny::reactive({
    build_explorer_wide_df(ge_data(),
                           feature_id       = selected_gene_id(),
                           feature_col_name = feature_col_name)
  })

  output[[id_data_table]] <- DT::renderDT({
    render_explorer_wide_dt(ge_table_df(), feature_col_name = feature_col_name)
  })

  # Downloads. Plot keeps the existing PNG behavior. Data download is the
  # wide table (RAW unrounded numerics) — matches the Multi-Explorer and the
  # on-screen Expression Data table.
  wire_downloads(
    output            = output,
    id_plot_btn       = id_dl_plot,
    id_data_btn       = id_dl_data,
    get_plot          = ge_plot,
    get_df            = ge_table_df,
    prefix            = "GeneExplorer",
    get_selected_id   = selected_gene_id,
    get_scale_sel     = function() input[[id_scale]],
    norm_method_label = get_norm_method_label_from_env(env)
  )

  output[[id_cutoffs_ui]] <- shiny::renderUI({
    padj0 <- if (exists("padj_cutoff", envir = env, inherits = FALSE)) {
      get("padj_cutoff", envir = env, inherits = FALSE)
    } else 0.05

    fc0 <- if (exists("fc_cutoff", envir = env, inherits = FALSE)) {
      get("fc_cutoff", envir = env, inherits = FALSE)
    } else 1

    build_de_cutoffs_ui(
      id_sig_mode   = id_sig_mode,
      id_sig_cutoff = id_padj_cutoff,
      id_fc_cutoff  = id_fc_cutoff,
      sig0          = padj0,
      fc0           = fc0
    )
  })

  ge_de_stats <- shiny::reactive({
    gene_id <- selected_gene_id()
    shiny::validate(
      shiny::need(exists("de_stats", envir = env, inherits = FALSE), "de_stats not available."),
      shiny::need(exists("contrasts", envir = env, inherits = FALSE), "contrasts not available.")
    )
    get_gene_de_stats_for_table(
      de_stats  = get("de_stats", envir = env, inherits = FALSE),
      contrasts = get("contrasts", envir = env, inherits = FALSE),
      gene_id   = gene_id
    )
  })

  output[[id_de_table]] <- DT::renderDataTable({
    df <- ge_de_stats()
    shiny::validate(shiny::need(!is.null(df), "No DE statistics available for this gene."))

    padj0 <- if (exists("padj_cutoff", envir = env, inherits = FALSE)) get("padj_cutoff", envir = env, inherits = FALSE) else 0.05
    fc0   <- if (exists("fc_cutoff",   envir = env, inherits = FALSE)) get("fc_cutoff",   envir = env, inherits = FALSE) else 1

    render_de_stats_dt(
      df          = df,
      sig_thr     = input[[id_padj_cutoff]],
      fc_thr      = input[[id_fc_cutoff]],
      sig_default = padj0,
      fc_default  = fc0,
      sig_mode    = if (!is.null(input[[id_sig_mode]])) input[[id_sig_mode]] else "padj"
    )
  })

  output[[id_de_legend_ui]] <- shiny::renderUI({
    selected_gene_id()  # req() guard — ensures re-render when gene changes
    de_table_legend_html()
  })

  # Annotation section (canonical feature_annot from this env).
  wire_feature_annotation(
    output           = output,
    env              = env,
    id_prefix        = id_prefix,
    ids_reactive     = selected_gene_id,
    feature_col_name = feature_col_name,
    multiple         = FALSE
  )

  # Click routing is only meaningful when there is a real selectize input to navigate to
  if (is.null(override_id_reactive)) {
    wire_click_to_explorer(
      session            = session,
      input              = input,
      tabs_id            = tabs_id,
      selected_tab_label = "Gene Explorer",
      selectize_input_id = id_gene,
      choices_vec        = gene_choices_all,
      choices_ready_rv   = ge_choices_ready,
      pending_rv         = ge_pending_gene,
      clicked_rv         = ge_clicked_gene,
      view_btn_id        = id_view_btn,
      modal_title        = "Gene selected",
      click_sources      = click_sources
    )
  }
}

# Get proteomics aesthetics from env (uses env$color / env$shape)
prot_get_effects_from_env <- function(env) {
  color_col <- if (exists("color", envir = env, inherits = FALSE)) get("color", envir = env, inherits = FALSE) else NULL
  shape_col <- if (exists("shape", envir = env, inherits = FALSE)) get("shape", envir = env, inherits = FALSE) else NULL
  
  eff <- character(0)
  if (!is.null(shape_col) && is.character(shape_col) && length(shape_col) == 1 && nzchar(shape_col)) eff <- c(eff, shape_col)
  if (!is.null(color_col) && is.character(color_col) && length(color_col) == 1 && nzchar(color_col)) eff <- c(eff, color_col)
  unique(eff)
}

# Build Protein/Metabolite Explorer DE stats table from env-format.
# Supports two ID columns: "Protein" (proteomics) or "feature_id" (metabolomics).
# Supports two stat-column patterns:
#   linearFC.imputs.<contrast> (proteomics)  and  linearFC.<contrast> (metabolomics).
get_protein_de_stats_for_table <- function(stats_df, contrasts_data, protein_id) {
  if (is.null(stats_df) || is.null(contrasts_data)) return(NULL)
  if (!"Contrast_name" %in% colnames(contrasts_data)) return(NULL)

  # Detect ID column: prefer "Protein", fall back to "feature_id"
  id_col <- if ("Protein" %in% colnames(stats_df)) {
    "Protein"
  } else if ("feature_id" %in% colnames(stats_df)) {
    "feature_id"
  } else {
    return(NULL)
  }

  row <- stats_df[stats_df[[id_col]] == protein_id, , drop = FALSE]
  if (nrow(row) == 0) return(NULL)

  cnames <- contrasts_data$Contrast_name
  out_list <- list()

  for (cname in cnames) {
    # Try proteomics naming first (linearFC.imputs.<contrast>),
    # then fall back to plain naming (linearFC.<contrast>).
    col_fc   <- paste0("linearFC.imputs.", cname)
    col_p    <- paste0("pvalue.imputs.",  cname)
    col_padj <- paste0("padj.imputs.",    cname)

    if (!all(c(col_fc, col_p, col_padj) %in% colnames(row))) {
      col_fc   <- paste0("linearFC.", cname)
      col_p    <- paste0("pvalue.",   cname)
      col_padj <- paste0("padj.",     cname)
    }

    if (!all(c(col_fc, col_p, col_padj) %in% colnames(row))) next

    out_list[[cname]] <- data.frame(
      Contrast = cname,
      linearFC = as.numeric(row[[col_fc]]),
      pvalue   = as.numeric(row[[col_p]]),
      padj     = as.numeric(row[[col_padj]]),
      stringsAsFactors = FALSE
    )
  }

  if (length(out_list) == 0) return(NULL)
  do.call(rbind, out_list)
}

# Wire Protein Explorer (expects NEW names in env: expr_raw, expr_norm, sample_meta, de_stats, contrasts, padj_cutoff, fc_cutoff)
# override_id_reactive: optional reactive() that supplies the protein ID directly (used in MultiOmics tab).
#   When provided, the selectize input and click-routing logic are bypassed entirely.
wire_protein_explorer <- function(session, input, output, env, tabs_id = "prot_tabs",
                                  click_sources = c("prot_de_volcano", "prot_de_ma"),
                                  id_prefix = "prot_explorer",
                                  override_id_reactive = NULL,
                                  explorer_tab_label = "Protein Explorer",
                                  feature_label = NULL,
                                  feature_col_name = "Protein",
                                  # Optional named vector of scale-selector
                                  # choices passed to updateSelectInput.
                                  # Defaults to the Proteomics imputation
                                  # labels; the Metabolite Explorer overrides
                                  # this with the chosen_norm-based labels.
                                  scale_choices = c(
                                    "without imputation (log2 scale)" = "raw",
                                    "Imputed (log2 scale)"            = "imp1"
                                  ),
                                  # Force a fixed scale value instead of
                                  # reading input[[id_scale]]. When set
                                  # (e.g. "imp1"), the Normalization
                                  # selectInput is NOT updated and every
                                  # `input[[id_scale]]` read inside this
                                  # helper falls back to fixed_scale.
                                  # Used by Metabolomics to lock the data
                                  # path to the normalized matrix without
                                  # rendering a dropdown.
                                  fixed_scale = NULL) {

  # Derive all input/output IDs from the prefix
  id_protein        <- paste0(id_prefix, "_protein")
  id_scale          <- paste0(id_prefix, "_scale")
  id_color_by       <- paste0(id_prefix, "_color_by")
  id_padj_cutoff    <- paste0(id_prefix, "_padj_cutoff")  # single significance cutoff input
  id_fc_cutoff      <- paste0(id_prefix, "_fc_cutoff")
  id_sig_mode       <- paste0(id_prefix, "_sig_mode")
  id_sel_title      <- paste0(id_prefix, "_selected_title")
  id_plot_title     <- paste0(id_prefix, "_plot_title")
  id_counts_plot    <- paste0(id_prefix, "_counts_plot")
  id_data_table     <- paste0(id_prefix, "_data_table")
  id_de_table       <- paste0(id_prefix, "_de_table")
  id_cutoffs_ui     <- paste0(id_prefix, "_cutoffs_ui")
  id_de_legend_ui   <- paste0(id_prefix, "_de_legend_ui")
  id_dl_plot        <- paste0("download_", id_prefix, "_plot")
  id_dl_data        <- paste0("download_", id_prefix, "_data")
  id_view_btn       <- paste0(id_prefix, "_view_clicked")

  # Selectize input and click-routing are only used when there is no external override
  if (is.null(override_id_reactive)) {
    prot_choices_all <- NULL
    if (exists("expr_raw", envir = env, inherits = FALSE)) {
      prot_choices_all <- prepare_id_choices_from_matrix(get("expr_raw", envir = env, inherits = FALSE))
    } else if (exists("expr_norm", envir = env, inherits = FALSE)) {
      prot_choices_all <- prepare_id_choices_from_matrix(get("expr_norm", envir = env, inherits = FALSE))
    }

    prot_choices_ready   <- shiny::reactiveVal(FALSE)
    prot_pending_protein <- shiny::reactiveVal(NULL)
    prot_clicked_protein <- shiny::reactiveVal(NULL)

    prot_sync_choices <- function() {
      if (is.null(prot_choices_all) || length(prot_choices_all) == 0) return(FALSE)

      cur  <- shiny::isolate(input[[id_protein]])
      keep <- if (!is.null(cur) && nzchar(cur) && cur %in% prot_choices_all) cur else character(0)

      shiny::updateSelectizeInput(
        session,
        inputId  = id_protein,
        choices  = prot_choices_all,
        server   = TRUE,
        selected = keep
      )
      prot_choices_ready(TRUE)
      TRUE
    }

    session$onFlushed(function() prot_sync_choices(), once = TRUE)

    shiny::observeEvent(input[[tabs_id]], {
      if (identical(input[[tabs_id]], explorer_tab_label)) prot_sync_choices()
    }, ignoreInit = TRUE)
  }

  # Update UI inputs (scale + color-by) from env configuration
  shiny::observe({
    eff <- prot_get_effects_from_env(env)
    shiny::validate(shiny::need(length(eff) >= 1, "No 'color' found for Protein Explorer (expected env$color)."))

    # User-facing scale labels for the single Protein / Metabolite Explorer.
    # `scale_choices` is parameterized (see wire_protein_explorer signature):
    # Proteomics uses the default imputation labels; Metabolomics passes
    # chosen_norm-based labels. Values ("raw"/"imp1") are unchanged so the
    # downstream data path keeps using the same matrices.
    # When fixed_scale is set the widget is hidden in ui.R — skip the
    # update so we don't push values into a non-existent control.
    if (is.null(fixed_scale)) {
      shiny::updateSelectInput(
        session, id_scale,
        choices  = scale_choices,
        selected = "raw"
      )
    }
    shiny::updateSelectInput(session, id_color_by, choices = eff, selected = eff[1])
  })

  # Single source of truth for the selected protein ID.
  # When override_id_reactive is provided it is the true source; otherwise the selectize input is.
  selected_protein_id <- shiny::reactive({
    if (!is.null(override_id_reactive)) {
      val <- override_id_reactive()
      shiny::req(!is.null(val) && nzchar(trimws(as.character(val))))
      return(trimws(as.character(val)))
    }
    shiny::req(input[[id_protein]])
    trimws(as.character(input[[id_protein]]))
  })

  # Selected protein title
  output[[id_sel_title]] <- shiny::renderUI({
    pid <- selected_protein_id()
    if (!is.null(feature_label)) {
      shiny::tags$div(
        style = "text-align:center; margin-top:6px; margin-bottom:2px;",
        shiny::tags$span(paste0(feature_label, ": "),
                         style = "color:#777; font-size:16px; font-weight:400;"),
        shiny::tags$span(pid,
                         style = "color:#18BC9C; font-size:24px; font-weight:800;")
      )
    } else {
      shiny::tags$div(pid, style = "color:#18BC9C;font-size:28px;font-weight:800;margin-top:6px;margin-bottom:2px;text-align:center;")
    }
  })

  # Plot title
  output[[id_plot_title]] <- shiny::renderUI({
    selected_protein_id()  # req() guard
    shiny::h4("Expression Profile", class = "plot-title full")
  })

  # Build long-format data for the selected protein
  prot_data <- shiny::reactive({
    shiny::validate(
      shiny::need(exists("sample_meta", envir = env, inherits = FALSE), "sample_meta is not available."),
      shiny::need(
        exists("expr_raw", envir = env, inherits = FALSE) || exists("expr_norm", envir = env, inherits = FALSE),
        "Neither expr_raw nor expr_norm is available."
      )
    )

    pid       <- selected_protein_id()
    color_by  <- shiny::req(input[[id_color_by]])
    # Use fixed_scale when set (Metabolomics), otherwise read the visible
    # selectInput (Proteomics / MultiOmics protein).
    scale_sel <- if (!is.null(fixed_scale)) fixed_scale else shiny::req(input[[id_scale]])  # raw / imp1

    sample_meta <- get("sample_meta", envir = env, inherits = FALSE)

    mat <- if (identical(scale_sel, "raw")) {
      shiny::validate(shiny::need(exists("expr_raw", envir = env, inherits = FALSE), "expr_raw is not available."))
      get("expr_raw", envir = env, inherits = FALSE)
    } else {
      shiny::validate(shiny::need(exists("expr_norm", envir = env, inherits = FALSE), "expr_norm is not available."))
      get("expr_norm", envir = env, inherits = FALSE)
    }

    get_explorer_df_from_mat(
      feature_id  = pid,
      mat         = mat,
      color_by    = color_by,
      sample_meta = sample_meta
    )
  })

  # Build the protein expression plot
  prot_plot <- shiny::reactive({
    df <- prot_data()

    scale_now <- if (!is.null(fixed_scale)) fixed_scale else input[[id_scale]]
    ylab <- if (identical(scale_now, "imp1")) {
      "Log2 protein intensity (imputed)"
    } else {
      "Log2 protein intensity (raw)"
    }

    make_explorer_plot(
      df            = df,
      color_by      = shiny::req(input[[id_color_by]]),
      scale_sel     = "log2",
      ylab_log2     = ylab,
      theme_variant = "gene"
    )
  })

  output[[id_counts_plot]] <- shiny::renderPlot({ prot_plot() })

  # Expression Data: wide-format table (one row for the selected feature,
  # one column per sample). Same shared helpers as the Multi-Protein /
  # Multi-Metabolite Explorers.
  prot_table_df <- shiny::reactive({
    build_explorer_wide_df(prot_data(),
                           feature_id       = selected_protein_id(),
                           feature_col_name = feature_col_name)
  })

  output[[id_data_table]] <- DT::renderDT({
    render_explorer_wide_dt(prot_table_df(), feature_col_name = feature_col_name)
  })

  # Downloads (plot + data). Data is the wide RAW-numeric table — matches
  # what the on-screen Expression Data table displays.
  wire_downloads(
    output            = output,
    id_plot_btn       = id_dl_plot,
    id_data_btn       = id_dl_data,
    get_plot          = prot_plot,
    get_df            = prot_table_df,
    prefix            = "ProteinExplorer",
    get_selected_id   = selected_protein_id,
    get_scale_sel     = function() {
      if (!is.null(fixed_scale)) fixed_scale else input[[id_scale]]
    },
    norm_method_label = NULL,
    log2_levels       = c("raw", "imp1")
  )

  # Cutoffs UI (highlighting-mode toggle + single significance cutoff +
  # linearFC cutoff) built from env, so UI does not need global variables
  output[[id_cutoffs_ui]] <- shiny::renderUI({
    padj0 <- if (exists("padj_cutoff", envir = env, inherits = FALSE)) {
      get("padj_cutoff", envir = env, inherits = FALSE)
    } else 0.05

    fc0 <- if (exists("fc_cutoff", envir = env, inherits = FALSE)) {
      get("fc_cutoff", envir = env, inherits = FALSE)
    } else 1

    build_de_cutoffs_ui(
      id_sig_mode   = id_sig_mode,
      id_sig_cutoff = id_padj_cutoff,
      id_fc_cutoff  = id_fc_cutoff,
      sig0          = padj0,
      fc0           = fc0
    )
  })

  # DE stats for the selected protein
  prot_de_stats <- shiny::reactive({
    pid <- selected_protein_id()
    shiny::validate(
      shiny::need(exists("de_stats",  envir = env, inherits = FALSE), "de_stats is not available."),
      shiny::need(exists("contrasts", envir = env, inherits = FALSE), "contrasts is not available.")
    )

    get_protein_de_stats_for_table(
      stats_df       = get("de_stats",  envir = env, inherits = FALSE),
      contrasts_data = get("contrasts", envir = env, inherits = FALSE),
      protein_id     = pid
    )
  })

  output[[id_de_table]] <- DT::renderDataTable({
    df <- prot_de_stats()
    shiny::validate(shiny::need(!is.null(df), "No DE statistics available for this protein."))

    padj0 <- if (exists("padj_cutoff", envir = env, inherits = FALSE)) get("padj_cutoff", envir = env, inherits = FALSE) else 0.05
    fc0   <- if (exists("fc_cutoff",   envir = env, inherits = FALSE)) get("fc_cutoff",   envir = env, inherits = FALSE) else 1

    render_de_stats_dt(
      df          = df,
      sig_thr     = input[[id_padj_cutoff]],
      fc_thr      = input[[id_fc_cutoff]],
      sig_default = padj0,
      fc_default  = fc0,
      sig_mode    = if (!is.null(input[[id_sig_mode]])) input[[id_sig_mode]] else "padj"
    )
  })

  output[[id_de_legend_ui]] <- shiny::renderUI(de_table_legend_html())

  # Annotation section (canonical feature_annot from this env).
  wire_feature_annotation(
    output           = output,
    env              = env,
    id_prefix        = id_prefix,
    ids_reactive     = selected_protein_id,
    feature_col_name = feature_col_name,
    multiple         = FALSE
  )

  # Click routing is only meaningful when there is a real selectize input to navigate to
  if (is.null(override_id_reactive)) {
    wire_click_to_explorer(
      session            = session,
      input              = input,
      tabs_id            = tabs_id,
      selected_tab_label = explorer_tab_label,
      selectize_input_id = id_protein,
      choices_vec        = prot_choices_all,
      choices_ready_rv   = prot_choices_ready,
      pending_rv         = prot_pending_protein,
      clicked_rv         = prot_clicked_protein,
      view_btn_id        = id_view_btn,
      modal_title        = "Protein selected",
      click_sources      = click_sources
    )
  }
}









# ── MultiOmics helpers ─────────────────────────────────────────────────────────

# Resolve a gene-protein pair from the matching table given a user query.
# query        : gene ID or protein ID (with or without version suffix)
# matching_tbl : data.frame with columns "proteins" (no suffix) and "genes"
# prot_ids_full: character vector of full protein IDs from the proteomics data
# Returns list(gene = ..., protein_full = ...) or NULL if no match found.
resolve_multi_pair <- function(query, matching_tbl, prot_ids_full) {
  if (is.null(matching_tbl) || nrow(matching_tbl) == 0 || !nzchar(query)) return(NULL)

  query_stripped <- sub("\\.\\d+$", "", query)

  # Strip version suffixes from the matching table's protein column
  tbl_prot_stripped <- sub("\\.\\d+$", "", matching_tbl$proteins)

  # Try matching the query as a protein ID (with or without suffix)
  match_idx <- which(tbl_prot_stripped == query_stripped | matching_tbl$proteins == query)

  if (length(match_idx) == 0) {
    # Try matching as a gene ID
    match_idx <- which(matching_tbl$genes == query)
  }

  if (length(match_idx) == 0) return(NULL)

  row       <- matching_tbl[match_idx[1], , drop = FALSE]
  gene_id   <- as.character(row$genes)
  prot_stub <- sub("\\.\\d+$", "", as.character(row$proteins))  # stub from Excel

  # Find the full protein ID in the proteomics data that matches the stub
  prot_full_ids_stripped <- sub("\\.\\d+$", "", prot_ids_full)
  hit <- prot_ids_full[prot_full_ids_stripped == prot_stub]

  if (length(hit) == 0) return(NULL)

  list(gene = gene_id, protein_full = hit[1])
}

# Build a named choices vector for the MultiOmics search selectize input.
# Combines gene IDs and full protein IDs that appear in the matching table.
# matching_tbl : data.frame with columns "proteins" (no suffix) and "genes"
# prot_env     : environment containing expr_raw (protein expression matrix)
# Returns a named character vector (label -> value) for updateSelectizeInput.
build_multi_search_choices <- function(matching_tbl, prot_env) {
  if (is.null(matching_tbl) || nrow(matching_tbl) == 0) return(character(0))

  gene_ids <- unique(as.character(matching_tbl$genes))
  gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]

  prot_ids_full <- character(0)
  if (exists("expr_raw", envir = prot_env, inherits = FALSE)) {
    prot_ids_full <- rownames(get("expr_raw", envir = prot_env, inherits = FALSE))
  } else if (exists("expr_norm", envir = prot_env, inherits = FALSE)) {
    prot_ids_full <- rownames(get("expr_norm", envir = prot_env, inherits = FALSE))
  }

  if (length(prot_ids_full) > 0) {
    tbl_prot_stripped    <- sub("\\.\\d+$", "", as.character(matching_tbl$proteins))
    full_stripped        <- sub("\\.\\d+$", "", prot_ids_full)
    matched_full         <- prot_ids_full[full_stripped %in% tbl_prot_stripped]
    matched_full         <- unique(matched_full[!is.na(matched_full) & nzchar(matched_full)])
  } else {
    matched_full <- character(0)
  }

  choices <- c(
    setNames(gene_ids,    paste0("[Gene] ",    gene_ids)),
    setNames(matched_full, paste0("[Protein] ", matched_full))
  )
  choices
}

# ---------------------------------------------------------------------------
# MultiOmics Explorer — Gene–Protein correlation helpers
# ---------------------------------------------------------------------------
# Isolated to the MultiOmics Explorer: build the paired normalized expression
# of a gene and a protein across the samples shared by both omics, decide
# which correlation methods are computable, run cor.test safely, and draw the
# scatter. None of these touch the matching / Expression Profile / DE paths.

# Paired data frame of normalized gene vs protein expression across the
# samples shared by both matrices. Only complete (non-NA) pairs are kept.
# Returns NULL when either feature is absent, the matrices share no samples,
# or fewer than 2 complete pairs remain.
get_gene_protein_paired_df <- function(gene_id, protein_id,
                                       gene_mat, protein_mat) {
  if (is.null(gene_mat) || is.null(protein_mat)) return(NULL)
  if (is.null(rownames(gene_mat)) || is.null(rownames(protein_mat))) return(NULL)
  if (is.null(colnames(gene_mat)) || is.null(colnames(protein_mat))) return(NULL)
  if (!gene_id %in% rownames(gene_mat))       return(NULL)
  if (!protein_id %in% rownames(protein_mat)) return(NULL)

  gv <- as.numeric(gene_mat[gene_id, ]);        names(gv) <- colnames(gene_mat)
  pv <- as.numeric(protein_mat[protein_id, ]);  names(pv) <- colnames(protein_mat)

  common <- intersect(names(gv), names(pv))
  common <- common[!is.na(common) & nzchar(common)]
  if (length(common) == 0L) return(NULL)

  df <- data.frame(
    Sample    = common,
    gene      = gv[common],
    protein   = pv[common],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  df <- df[stats::complete.cases(df[, c("gene", "protein"), drop = FALSE]), ,
           drop = FALSE]
  if (nrow(df) < 2L) return(NULL)
  df
}

# Safely run cor.test for one method on the paired df. All errors and
# warnings (e.g. "cannot compute exact p-value with ties") are swallowed;
# returns the htest object, or NULL on any failure / too-few samples.
gene_protein_corr_test <- function(df, method) {
  if (is.null(df) || nrow(df) < 3L) return(NULL)
  tryCatch(
    suppressWarnings(stats::cor.test(df$gene, df$protein, method = method)),
    error = function(e) NULL
  )
}

# Detect which correlation methods can actually be computed for the pair.
# A method survives only when cor.test succeeds and returns a finite
# estimate; constant vectors, fewer than 3 paired samples, or other
# statistical limitations drop the method from the result.
gene_protein_corr_valid_methods <- function(df,
                                            methods = c("pearson", "spearman",
                                                        "kendall")) {
  if (is.null(df) || nrow(df) < 3L) return(character(0))
  if (stats::sd(df$gene) == 0 || stats::sd(df$protein) == 0) return(character(0))
  ok <- vapply(methods, function(m) {
    res <- gene_protein_corr_test(df, m)
    !is.null(res) && is.finite(unname(res$estimate))
  }, logical(1))
  methods[ok]
}

# Human-readable label for a correlation method key.
gene_protein_corr_method_label <- function(method) {
  unname(c(pearson = "Pearson", spearman = "Spearman",
           kendall = "Kendall")[method])
}

# Correlation statistics table for the MultiOmics section: one row per *valid*
# method (Pearson/Spearman/Kendall), with the coefficient and p-value. Reuses
# the same validation path as the rest of the section
# (gene_protein_corr_valid_methods + gene_protein_corr_test), so a method only
# appears when cor.test actually succeeds. Returns NULL when no method is
# valid. Columns are display-ready strings.
gene_protein_corr_stats_df <- function(df,
                                       methods = c("pearson", "spearman",
                                                   "kendall")) {
  valid <- gene_protein_corr_valid_methods(df, methods)
  if (length(valid) == 0L) return(NULL)
  rows <- lapply(valid, function(m) {
    res <- gene_protein_corr_test(df, m)
    if (is.null(res)) return(NULL)
    pval <- res$p.value
    # Clamp to [0, 1]: a perfectly monotonic pair can make cor.test return a
    # tiny negative p (floating-point artifact). Anything below machine
    # precision is shown as "< 2.2e-16" rather than a misleading exact figure.
    pval_disp <- if (is.null(pval) || is.na(pval)) {
      "not available"
    } else {
      pv <- min(1, max(0, pval))
      if (pv < .Machine$double.eps) "< 2.2e-16"
      else                          format_scientific_threshold(pv)
    }
    data.frame(
      Method                    = gene_protein_corr_method_label(m),
      `Correlation Coefficient` = formatC(unname(res$estimate),
                                          format = "f", digits = 3),
      `P-value`                 = pval_disp,
      check.names      = FALSE,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Interactive Plotly scatter of gene (x) vs protein (y) normalized expression
# for the MultiOmics correlation section. Each point is one biological sample;
# hover shows the sample name plus its gene and protein values. A single linear
# (lm) trend line is overlaid as a *visual aid only* — it is NOT specific to
# Pearson/Spearman/Kendall. Drawn only when >= 3 samples and gene is non-constant.
plotly_gene_protein_correlation <- function(df, gene_id, protein_id,
                                            add_trend = TRUE, source_id = NULL) {
  df$Sample <- as.character(df$Sample)

  # No trace names anywhere: both traces are unnamed and excluded from the
  # legend, so "Samples"/"Linear trend" never surface in the legend, hover, or
  # any Plotly control. The marker hovertemplate's <extra></extra> also strips
  # the trace-name box from tooltips; the trend line has hover disabled.
  p <- plotly::plot_ly(source = source_id)
  p <- plotly::add_trace(
    p,
    data       = df,
    x          = ~gene, y = ~protein,
    type       = "scatter", mode = "markers",
    showlegend = FALSE,
    marker = list(size = 10, color = "#18BC9C",
                  line = list(width = 1, color = "#0f7a66")),
    text   = ~paste0("Sample: ", Sample,
                     "<br>Gene expression: ",
                     formatC(gene, format = "f", digits = 3),
                     "<br>Protein expression: ",
                     formatC(protein, format = "f", digits = 3)),
    hovertemplate = "%{text}<extra></extra>"
  )

  # Simple linear trend (lm) — visual aid only, method-agnostic. Guarded so a
  # degenerate fit never errors the render.
  if (isTRUE(add_trend) && nrow(df) >= 3L && stats::sd(df$gene) > 0) {
    fit <- tryCatch(stats::lm(protein ~ gene, data = df), error = function(e) NULL)
    if (!is.null(fit)) {
      xr <- range(df$gene, na.rm = TRUE)
      xs <- seq(xr[1], xr[2], length.out = 50)
      ys <- tryCatch(
        as.numeric(stats::predict(fit, newdata = data.frame(gene = xs))),
        error = function(e) NULL)
      if (!is.null(ys) && all(is.finite(ys))) {
        p <- plotly::add_lines(
          p, x = xs, y = ys,
          showlegend = FALSE,
          line       = list(color = "#2C3E50", width = 2),
          hoverinfo  = "skip", inherit = FALSE)
      }
    }
  }

  p <- plotly::layout(
    p,
    xaxis  = list(title = paste0("Gene expression: ", gene_id),
                  zeroline = FALSE),
    yaxis  = list(title = paste0("Protein expression: ", protein_id),
                  zeroline = FALSE),
    showlegend = FALSE,
    autosize   = TRUE
  )

  # Exploration-only plot: strip the image-export / data-export modebar
  # actions (and the plotly logo) but keep hover, zoom, pan, and reset.
  plotly::config(
    p,
    displaylogo          = FALSE,
    modeBarButtonsToRemove = list("toImage", "sendDataToCloud")
  )
}


# ============================================================================
# ENRICHMENT tab (Stage 4) — OMICS-AGNOSTIC implementation
# ----------------------------------------------------------------------------
# This whole block is generic: it depends ONLY on the §14 enrichment payload
# contract (`env$enrichment`) and the standard omics-env fields it reuses
# (`env$de_stats`, `env$expr_norm`), NEVER on any RNA-seq-specific global. Every
# helper, resolver, table, plot, master-detail and download takes `env`,
# `enrichment`, `manifest`, or `u` as a parameter; every Shiny input/output id is
# namespaced by `id_prefix`. RNA-seq is simply the FIRST wired instance
# (`enrichment_panel_ui("rna_enrich")` + `wire_enrichment(env=rnaseq_env,
# id_prefix="rna_enrich", tabs_id="rna_tabs")`). A future omics whose RDS carries
# the SAME payload wires the identical UI with its own env + id_prefix + tabs_id
# — no logic is duplicated or rewritten. An omics whose env has no (or an
# unavailable) `enrichment` shows the friendly "not available" message via the
# same generic path (`enrichment_data_available_env()` → FALSE).
#
# The manifest (`enrichment$manifest`, columns analysis/database/group/item/
# evaluated/status/n_significant/has_simplify/storage_key) is the SOLE
# availability source. Internal state uses the manifest-aligned levels
# analysis→database→group→item (GSEA ranking maps to `group`, GSEA contrast to
# `item`); display labels vary by analysis. Pure helpers below are runtime-free
# and unit-tested directly.
# ============================================================================

.ENRICH_ANALYSIS_ORDER  <- c("ORA", "GSEA")
.ENRICH_DB_ORDER        <- c("GO_BP", "GO_MF", "GO_CC", "KEGG")
# §17 R2/S1: the contrast collections are UNIFIED under one group "contrasts"; the
# with-/without-direction distinction is a `subgroup` (up/down/up_down vs all), so
# the ORA hierarchy is  Database -> Group -> Item(contrast|cluster|pattern|compare)
# -> [Subgroup only for contrasts]. `contrasts_wo_direction` is retained in the
# maps ONLY for backward-compatibility with pre-redesign manifests.
.ENRICH_ORA_GROUP_ORDER <- c("contrasts", "all_DE", "partition", "binary_patterns")
# §23 H: user-facing order Up, Down, DE, Compare Up & Down. INTERNAL subgroup ids
# stay up/down/all/up_down (only labels + on-disk dir names change); "all" here is
# the subgroup id, NOT the unrelated Top-N sentinel.
.ENRICH_ORA_SUBGROUP_ORDER <- c("up", "down", "all", "up_down")
.ENRICH_RANKING_ORDER   <- c("pval_wo_direction", "pval_with_direction", "fc")
.ENRICH_ORA_GROUP_LABEL <- c(contrasts = "Contrasts",
                             contrasts_wo_direction = "Contrasts without direction",
                             all_DE = "All DE", partition = "Partition",
                             binary_patterns = "Binary patterns")
# §23 H: labels only — ids unchanged. all -> "DE"; up_down -> "Compare Up vs. Down"
# (user-facing rename; internal id stays `up_down`, dir stays `compare_up_down`).
.ENRICH_ORA_SUBGROUP_LABEL <- c(up = "Up", down = "Down",
                                all = "DE", up_down = "Compare Up vs. Down")
# §23 H: internal subgroup id -> on-disk directory name (up/down/de/compare_up_down).
.ENRICH_ORA_SUBGROUP_DIR <- c(up = "up", down = "down",
                              all = "de", up_down = "compare_up_down")
# §23 I: GSEA ranking labels — "P-value with direction" / "P-value without direction".
.ENRICH_RANKING_LABEL   <- c(fc = "FC",
                             pval_with_direction = "P-value with direction",
                             pval_wo_direction = "P-value without direction")
# manifest group value -> payload `enrichment$ora[[db]][[...]]` node name. For the
# unified "contrasts" group the node is subgroup-aware (see .enrich_ora_payload_group).
.ENRICH_ORA_PAYLOAD_GROUP <- c(contrasts = "contrasts_with_direction",
                               contrasts_wo_direction = "contrasts_without_direction",
                               all_DE = "all_DE", partition = "partition",
                               binary_patterns = "binary_patterns")
# Expression heatmaps (ORA member-gene, GSEA core-gene, GSEA all-gene) are NOT
# gene-count-limited (§16.43): they attempt to render every relevant gene, never
# truncating/sampling, and only render lazily when their view is opened. The old
# ENRICH_MAX_HEATMAP_GENES=500 heatmap guard was removed.
# The all-genes RIDGEPLOT is ALSO no longer gene-count-limited (§16.44): the old
# ENRICH_RIDGE_MAX_GENES=500 per-pathway exclusion was removed so the Shiny ridge
# includes the SAME top-N pathways as the pipeline (large GO root terms included).
# Y-axis term/pathway label size shared by ALL enrichment bar plots (ORA + GSEA
# + gene recurrence, §16.37/§16.39). Term names are the point of these plots, so
# they must be comfortably readable. Centralized (via .enrich_plot_theme) so the
# plots can never drift apart and the size is retunable in one place. Ridgeplots
# keep their own smaller scale (many stacked rows) and are NOT changed here.
ENRICH_PLOT_AXIS_TEXT_SIZE <- 13L
# Visualization-only display bounds (§16.39). None of these touch a biological
# value: they bound how much is DRAWN, never what is computed, tabled or
# downloaded.
ENRICH_BAR_TOP_N        <- 20L   # ORA/GSEA enrichment bar plot: top terms by FDR
ENRICH_BAR_ANNOT_SIZE   <- 4.4   # ggplot size of the bar-end annotation text
ENRICH_BAR_EXPAND       <- 0.42  # bar-axis headroom so annotations never clip
ENRICH_RECURRENCE_TOP_N <- 50L   # shared-genes CHART bound (the table shows all)
# §18 D18.3: shared-terms matrix default view = top 20 by adj p (Top 20/50/All
# control); also the threshold above which the control is shown. Was 40.
ENRICH_T2T_TOP_N        <- 20L   # shared-terms matrix: default term count
# UI refinement (post-S5): a bigger per-label legibility budget + larger minimums
# so the Shared Terms matrix never renders too small to read its axis labels.
ENRICH_T2T_PX_PER_TERM  <- 22L   # legibility budget for one axis label (was 14)
ENRICH_T2T_MIN_HEIGHT   <- 460L
ENRICH_T2T_MAX_HEIGHT   <- 6000L
ENRICH_T2T_MIN_WIDTH    <- 620L
ENRICH_T2T_MAX_WIDTH    <- 6000L
ENRICH_TERM_LABEL_CHARS <- 60L   # term-label truncation (was 45)
ENRICH_RIDGE_TOP_N      <- 20L   # GSEA ridgeplot: top pathways (pipeline show_category)
# GSEA ridgeplot X = fold-change (§16.42): the ridge distribution always uses the
# per-gene FC vector (log2FC) for the SAME contrast, even when a p-value ranking
# is selected — the ranking still governs the GSEA result (NES/FDR/membership),
# only the ridge X-values are remapped to FC. `any_contrast`/"All DE" has NO FC
# vector, so the ridge is shown as unavailable rather than falling back to ranks.
ENRICH_RIDGE_X_LABEL    <- "log2 Fold Change"   # verified scale (pipeline rank_by_fc = log2FC)
.ENRICH_RIDGE_FC_CAPTION <- paste(
  "Ridge distributions show fold-change values for genes in each pathway;",
  "GSEA significance is based on the selected ranking method.")
# One-line encoding legend shown ONCE BELOW each enrichment bar plot (§16.41),
# never repeated per bar. Same phrasing pattern for both analyses so they read
# as one family, while each names its OWN biological x metric — ORA is a gene
# count, GSEA is a signed NES, and those must never be conflated.
ENRICH_ORA_BAR_CAPTION <- paste(
  "Bar length = gene count  ·  Colour = -log10(FDR)  ·",
  "Label = genes in term / total term genes  ·  FE = fold enrichment")
ENRICH_GSEA_BAR_CAPTION <- paste(
  "Bar length = signed NES  ·  Colour = -log10(FDR)  ·",
  "Label = pathway set size  ·  LE = leading-edge genes")

#' Shared theme for the enrichment bar plots — minimal base + enlarged y-axis
#' term labels. Purely presentational.
#' @noRd
.enrich_plot_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(
      size = ENRICH_PLOT_AXIS_TEXT_SIZE),
      panel.grid.major.y = ggplot2::element_blank())
}

#' DISPLAY-ONLY -log10(FDR) (§16.39). The biological `p.adjust` is never
#' modified — this only produces a colour value. Zero / underflowed FDRs (fgsea
#' and clusterProfiler can emit exact 0) are floored at the smallest POSITIVE
#' adjusted p-value actually present, so the scale stays on the data's own
#' range instead of exploding to ~308. Non-finite input is treated as p = 1.
#' @noRd
.enrich_nlp <- function(padj) {
  p <- suppressWarnings(as.numeric(padj))
  p[!is.finite(p)] <- 1
  pos <- p[p > 0]
  eps <- if (length(pos)) max(min(pos), 1e-300) else 1e-300
  -log10(pmax(p, eps))
}

#' Upper colour limit for a -log10(FDR) fill: the 98th percentile, so one
#' hyper-significant term cannot flatten the whole scale. Used with
#' `oob = scales::squish` (values above the limit take the extreme colour).
#' @noRd
.enrich_nlp_limits <- function(nlp, q = 0.98) {
  nlp <- nlp[is.finite(nlp)]
  if (!length(nlp)) return(c(0, 1))
  hi <- suppressWarnings(stats::quantile(nlp, q, na.rm = TRUE, names = FALSE))
  if (!is.finite(hi) || hi <= 0) hi <- max(c(nlp, 1))
  c(0, hi)
}

#' Overlap-gene lists for a unit's ORA terms (split the `geneID` column).
#' Single source used by the recurrence chart and the shared-terms matrix.
#' @noRd
.enrich_gene_lists <- function(df) {
  lapply(strsplit(as.character(df$geneID), "/", fixed = TRUE),
         function(x) x[!is.na(x) & nzchar(x)])
}

#' Term size (genes annotated to the term in the background) for the bar
#' annotation: the BgRatio numerator ("29/9087" -> 29), falling back to the
#' `in_term` column. NA when neither is present — the annotation then omits it.
#' @noRd
.enrich_term_size <- function(df) {
  if (!is.null(df$BgRatio)) {
    n <- suppressWarnings(as.integer(sub("/.*$", "", as.character(df$BgRatio))))
    if (any(!is.na(n))) return(n)
  }
  if (!is.null(df$in_term)) return(suppressWarnings(as.integer(df$in_term)))
  rep(NA_integer_, nrow(df))
}

#' Full-precision fold enrichment, preferring clusterProfiler's unrounded
#' `FoldEnrichment` over the 1-dp `Fold_enrichment` mirror.
#' @noRd
.enrich_fold <- function(df) {
  v <- if (!is.null(df$FoldEnrichment)) df$FoldEnrichment else df$Fold_enrichment
  suppressWarnings(as.numeric(v))
}

#' Truncated, unique y-axis term labels.
#' @noRd
.enrich_term_labels <- function(x, chars = ENRICH_TERM_LABEL_CHARS) {
  make.unique(substr(as.character(x), 1, chars))
}

# Title of the unit-level plots panel in the flat content row (§16.38). Used to
# gate the selected-item detail, which belongs to the RESULT views only (§16.40).
.ENRICH_PLOTS_VIEW <- "Plots"

#' Is the Plots view the active content view? (§16.40) The selected term/pathway
#' detail renders under Full Results / GO Simplified only — under Plots the user
#' wants the unit-level plot area alone. NULL (before the client reports a value)
#' means the first panel, Full Results, is active -> not Plots.
#' @noRd
.enrich_is_plots_view <- function(v) identical(v, .ENRICH_PLOTS_VIEW)

.enr_lookup <- function(map, key) {
  if (is.null(key) || is.na(key)) return(key)
  v <- map[key]
  if (length(v) == 0 || is.na(v)) key else unname(v)
}

#' Is a usable enrichment payload present in `env`? (never errors)
#' @noRd
enrichment_data_available_env <- function(env) {
  if (!exists("enrichment", envir = env, inherits = FALSE)) return(FALSE)
  e <- tryCatch(get("enrichment", envir = env, inherits = FALSE),
                error = function(err) NULL)
  if (is.null(e) || !is.list(e)) return(FALSE)
  if (!isTRUE(e$available)) return(FALSE)
  m <- e$manifest
  if (!is.data.frame(m)) return(FALSE)
  req <- c("analysis", "database", "group", "item", "status", "n_significant")
  all(req %in% names(m))
}

# §17 S1: the selector levels, deepest last. `subgroup` (up/down/up_down/all) is
# the OPTIONAL 5th level, reached only for ORA contrasts; all other units are
# terminal at `item` (partition/binary) or `group` (all_DE).
.ENRICH_LEVEL_ORDER <- c("analysis", "database", "group", "item", "subgroup")

#' Initial (empty) selector state.
#' @noRd
enrichment_init_state <- function() {
  list(analysis = NULL, database = NULL, group = NULL, item = NULL,
       subgroup = NULL, selected_id = NULL)
}

#' Apply a selection at `level`; clears all downstream levels + selected_id.
#' @noRd
enrichment_apply_selection <- function(state, level, value) {
  order <- .ENRICH_LEVEL_ORDER
  idx <- match(level, order)
  if (is.na(idx)) return(state)
  state[[level]] <- value
  if (idx < length(order))
    for (lv in order[(idx + 1):length(order)]) state[[lv]] <- NULL
  state$selected_id <- NULL
  state
}

#' Set database, optionally restoring a pending (partial-nav) group/item/subgroup.
#' @noRd
enrichment_set_database <- function(state, database, pending = NULL) {
  new <- enrichment_apply_selection(state, "database", database)
  if (!is.null(pending)) {
    new$group <- pending$group; new$item <- pending$item
    if (!is.null(pending$subgroup)) new$subgroup <- pending$subgroup
  }
  new
}

#' Clear a level and everything downstream (breadcrumb step-back).
#' @noRd
enrichment_reset_level <- function(state, level) {
  order <- .ENRICH_LEVEL_ORDER
  idx <- match(level, order)
  if (is.na(idx)) return(state)
  for (lv in order[idx:length(order)]) state[[lv]] <- NULL
  state$selected_id <- NULL
  state
}

#' Build a (partial or full) state for controller `open()`.
#' GSEA callers may pass ranking=/contrast=; they map to group/item. `subgroup`
#' (up/down/up_down/all) is optional and only meaningful for an ORA contrast.
#' @noRd
enrichment_open_state <- function(analysis, database = NULL, group = NULL,
                                  item = NULL, ranking = NULL, contrast = NULL,
                                  subgroup = NULL) {
  g  <- if (!is.null(group)) group else ranking
  it <- if (!is.null(item))  item  else contrast
  s <- enrichment_init_state()
  s$analysis <- analysis; s$database <- database; s$group <- g; s$item <- it
  s$subgroup <- subgroup
  s
}

#' Classify a set of manifest rows into "enabled" / "greyed" / NA(=hidden).
#' Failed rows are excluded from existence AND significance (never "empty").
#' @noRd
.enrich_classify <- function(rows) {
  if (is.null(rows) || nrow(rows) == 0) return(NA_character_)
  nf <- rows[!is.na(rows$status) & rows$status != "failed", , drop = FALSE]
  if (nrow(nf) == 0) return(NA_character_)
  if (any(nf$status == "significant")) "enabled" else "greyed"
}

.enrich_option_label <- function(level, analysis, value) {
  if (level %in% c("analysis", "database")) return(value)
  if (level == "group") {
    if (identical(analysis, "GSEA")) return(.enr_lookup(.ENRICH_RANKING_LABEL, value))
    return(.enr_lookup(.ENRICH_ORA_GROUP_LABEL, value))
  }
  if (level == "subgroup") return(.enr_lookup(.ENRICH_ORA_SUBGROUP_LABEL, value))
  # item level: the aggregate-contrast unit. GSEA shows "Any Contrast" (§D29 §1);
  # ORA's all_DE keeps "All DE". Internal id stays `any_contrast` in both cases.
  if (identical(value, "any_contrast"))
    return(if (identical(analysis, "GSEA")) "Any Contrast" else "All DE")
  value
}

#' First non-NA `group_label` among a set of manifest rows (the §17 R2 display
#' string, e.g. "Cluster 1" / "Compare Clusters" / "Up + Down"); NA if none.
#' @noRd
.enrich_group_label_of <- function(rows) {
  if (is.null(rows) || !("group_label" %in% names(rows)) || nrow(rows) == 0)
    return(NA_character_)
  gl <- as.character(rows$group_label); gl <- gl[!is.na(gl) & nzchar(gl)]
  if (length(gl)) gl[1] else NA_character_
}

#' Three-state options for one selector level given the current upstream state.
#' Returns data.frame(value, label, status) in canonical order; hidden options
#' (no non-failed rows in scope) are omitted entirely. §17 S1 adds the `subgroup`
#' level (ORA contrasts: up/down/up_down/all) and prefers the manifest
#' `group_label` for ORA item (cluster/pattern/compare) + subgroup display.
#' @noRd
resolve_enrichment_options <- function(manifest, level, state) {
  empty <- data.frame(value = character(0), label = character(0),
                      status = character(0), stringsAsFactors = FALSE)
  if (is.null(manifest) || !is.data.frame(manifest) || nrow(manifest) == 0) return(empty)
  if (level == "subgroup" && !("subgroup" %in% names(manifest))) return(empty)
  m <- manifest
  if (level %in% c("database", "group", "item", "subgroup") && !is.null(state$analysis))
    m <- m[m$analysis == state$analysis, , drop = FALSE]
  if (level %in% c("group", "item", "subgroup") && !is.null(state$database))
    m <- m[m$database == state$database, , drop = FALSE]
  if (level %in% c("item", "subgroup") && !is.null(state$group))
    m <- m[m$group == state$group, , drop = FALSE]
  if (level == "subgroup" && !is.null(state$item))
    m <- m[m$item == state$item, , drop = FALSE]
  if (nrow(m) == 0) return(empty)
  cand <- switch(level, analysis = "analysis", database = "database",
                 group = "group", item = "item", subgroup = "subgroup")
  present <- unique(as.character(m[[cand]]))
  present <- present[!is.na(present)]            # subgroup NA (non-contrast rows) never an option
  universe <- switch(level,
    analysis = intersect(.ENRICH_ANALYSIS_ORDER, present),
    database = c(intersect(.ENRICH_DB_ORDER, present), setdiff(present, .ENRICH_DB_ORDER)),
    group    = {
      ord <- if (identical(state$analysis, "GSEA")) .ENRICH_RANKING_ORDER else .ENRICH_ORA_GROUP_ORDER
      c(intersect(ord, present), setdiff(present, ord))
    },
    subgroup = c(intersect(.ENRICH_ORA_SUBGROUP_ORDER, present),
                 setdiff(present, .ENRICH_ORA_SUBGROUP_ORDER)),
    item     = present)
  # §23 H: SUBGROUP labels come from `.ENRICH_ORA_SUBGROUP_LABEL` (id -> new label,
  # e.g. all -> "DE", up_down -> "Compare Up & Down") via `.enrich_option_label`,
  # NOT the manifest's baked `group_label` — so the rename applies to existing
  # payloads with no regeneration. Only ORA partition/binary ITEM labels (data-
  # derived cluster/pattern names) still read the manifest group_label.
  use_group_label <- identical(state$analysis, "ORA") &&
    level == "item" && !is.null(state$group) &&
    state$group %in% c("partition", "binary_patterns")
  rows <- lapply(universe, function(v) {
    sub <- m[as.character(m[[cand]]) == v, , drop = FALSE]
    stt <- .enrich_classify(sub)
    if (is.na(stt)) return(NULL)
    lbl <- if (use_group_label) {
      gl <- .enrich_group_label_of(sub)
      if (!is.na(gl)) gl else .enrich_option_label(level, state$analysis, v)
    } else .enrich_option_label(level, state$analysis, v)
    data.frame(value = v, label = lbl, status = stt, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  # §5/§13: Compare Contrasts — a synthetic cross-contrast ITEM option appended
  # AFTER the individual contrasts, shown when >=2 contrasts have a significant
  # direction-agnostic DE (subgroup "all") result. Derived from the manifest alone
  # (no payload read, no new manifest row required for existing payloads); skipped
  # if the payload already supplies a `compare_contrasts` item. Generic: it keys on
  # the generic `contrasts` group + `all` subgroup, never on hard-coded contrast
  # names, so the resolver stays omics-agnostic.
  if (identical(level, "item") && identical(state$analysis, "ORA") &&
      identical(.enrich_scalar(state$group), "contrasts") &&
      !("compare_contrasts" %in% universe) && "subgroup" %in% names(m)) {
    de <- m[!is.na(m$subgroup) & as.character(m$subgroup) == "all" &
            !is.na(m$status) & m$status == "significant", , drop = FALSE]
    if (length(unique(as.character(de$item))) >= 2L)
      rows[[length(rows) + 1L]] <- data.frame(
        value = "compare_contrasts", label = "Compare Contrasts",
        status = "enabled", stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) return(empty)
  do.call(rbind, rows)
}

#' Count of "failed" units within the currently-fixed scope (for the note).
#' @noRd
enrichment_scope_failed_count <- function(manifest, state) {
  if (is.null(manifest) || !is.data.frame(manifest) || nrow(manifest) == 0) return(0L)
  m <- manifest
  if (!is.null(state$analysis)) m <- m[m$analysis == state$analysis, , drop = FALSE]
  if (!is.null(state$database)) m <- m[m$database == state$database, , drop = FALSE]
  if (!is.null(state$group))    m <- m[m$group    == state$group,    , drop = FALSE]
  sum(!is.na(m$status) & m$status == "failed")
}

#' Resolved unit (or NULL). Terminal levels (§17 S1):
#'   - ORA all_DE      -> terminal at group (item forced "any_contrast");
#'   - ORA contrasts   -> terminal at SUBGROUP (up/down/up_down/all) when the
#'                        manifest carries a `subgroup` column (the redesign);
#'   - everything else -> terminal at item (partition/binary cluster/pattern/compare,
#'                        GSEA contrast). Pre-redesign manifests (no `subgroup`
#'                        column) keep resolving ORA contrasts at item.
#' Carries `subgroup`/`group_type`/`group_label` for the downstream resolver + UI.
#' @noRd
enrichment_resolved_unit <- function(manifest, state) {
  if (is.null(manifest) || !is.data.frame(manifest)) return(NULL)
  a <- state$analysis; db <- state$database; g <- state$group
  if (is.null(a) || is.null(db) || is.null(g)) return(NULL)
  it <- if (identical(a, "ORA") && identical(g, "all_DE")) "any_contrast" else state$item
  if (is.null(it)) return(NULL)
  # §5: Compare Contrasts — a synthetic TERMINAL cross-contrast compare unit (no
  # subgroup, no per-item manifest row). group_type = "contrast" so the dot matrix
  # + Download Table label the groups with the real contrast names.
  if (identical(a, "ORA") && identical(g, "contrasts") &&
      identical(.enrich_scalar(it), "compare_contrasts")) {
    return(list(analysis = a, database = db, group = g, item = "compare_contrasts",
                subgroup = NA_character_, group_type = "contrast",
                group_label = NA_character_, n_significant = NA_integer_,
                status = "significant", storage_key = NA_character_,
                has_simplify = FALSE))
  }
  has_subgroup <- identical(a, "ORA") && identical(g, "contrasts") &&
                  ("subgroup" %in% names(manifest))
  sg <- if (has_subgroup) state$subgroup else NA_character_
  if (has_subgroup && is.null(sg)) return(NULL)   # contrast chosen, subgroup pending
  rows <- manifest[manifest$analysis == a & manifest$database == db &
                   manifest$group == g & manifest$item == it, , drop = FALSE]
  if (has_subgroup)
    rows <- rows[as.character(rows$subgroup) == as.character(sg), , drop = FALSE]
  if (nrow(rows) == 0) return(NULL)
  row <- rows[1, ]
  fld <- function(nm) if (nm %in% names(row)) as.character(row[[nm]]) else NA_character_
  list(analysis = a, database = db, group = g, item = it,
       subgroup   = if (has_subgroup) as.character(sg) else NA_character_,
       group_type = fld("group_type"), group_label = fld("group_label"),
       n_significant = row$n_significant, status = row$status,
       storage_key = if ("storage_key" %in% names(row)) row$storage_key else NA_character_,
       has_simplify = if ("has_simplify" %in% names(row)) isTRUE(row$has_simplify) else FALSE)
}

# §17 R2/R4 resolution contract, MIRRORED on the Shiny side (a pure filter-on-read
# of the compact whole-table payload nodes — NO biology is recomputed in Shiny;
# the same up/down/up_down/all + compare semantics the pipeline validated).

#' Payload `enrichment$ora[[db]]` node name for an ORA (group, subgroup). For the
#' unified "contrasts" group it is subgroup-aware: subgroup "all" -> the DISTINCT
#' `contrasts_without_direction` node; up/down/up_down -> `contrasts_with_direction`.
#' @noRd
# Normalize a possibly-NULL / length-0 / NA level value to a scalar (NA when absent).
.enrich_scalar <- function(x) {
  x <- as.character(x)
  if (length(x) == 0) NA_character_ else x[1]
}
.enrich_ora_payload_group <- function(group, subgroup) {
  if (identical(group, "contrasts"))
    return(if (identical(.enrich_scalar(subgroup), "all")) "contrasts_without_direction"
           else "contrasts_with_direction")
  pg <- .ENRICH_ORA_PAYLOAD_GROUP[group]
  if (is.na(pg)) group else unname(pg)
}

#' `Cluster` value to keep for an ORA selection, or NA = the WHOLE table
#' (UP+DOWN / All / Compare / all_DE). Fully determined by (group, subgroup, item).
#' Robust to a NULL / length-0 subgroup (a manually built or pre-redesign unit) —
#' which resolves to the whole table, the pre-subgroup behavior.
#' @noRd
.enrich_ora_cluster_filter <- function(group, subgroup, item) {
  sg <- .enrich_scalar(subgroup)
  if (identical(group, "contrasts"))
    return(if (!is.na(sg) && sg %in% c("up", "down")) sg else NA_character_)
  if (group %in% c("partition", "binary_patterns")) {
    it <- .enrich_scalar(item)
    return(if (identical(it, "compare")) NA_character_ else it)
  }
  NA_character_
}

#' The single biologically-defined contrast a resolved unit is tied to, or NA
#' when it is NOT tied to one contrast. Used to restrict the selected-term /
#' core-gene DE-stat table to just that contrast's columns (pass-3, item 1). A
#' contrast exists for:
#'   - an ORA `contrasts` unit (item = the contrast, for ANY direction subgroup
#'     up/down/up_down/all — the direction is a slice of the SAME contrast);
#'   - a GSEA unit whose item is a real contrast.
#' It does NOT exist (returns NA) for ORA partition / binary_patterns / all_DE,
#' GSEA `any_contrast`, or any `compare` unit — those aggregate DE sets or
#' rankings with no single contrast, so the table keeps EVERY contrast. Never
#' invents a contrast: the value returned is only ever the unit's own `item`.
#' @noRd
.enrich_unit_contrast <- function(u) {
  if (is.null(u)) return(NA_character_)
  it <- .enrich_scalar(u$item)
  if (is.na(it) || it %in% c("any_contrast", "compare")) return(NA_character_)
  if (identical(u$analysis, "GSEA")) return(it)
  if (identical(u$analysis, "ORA") &&
      identical(.enrich_scalar(u$group), "contrasts")) return(it)
  NA_character_
}

#' Resolve the ORA result data.frame for a resolved ORA unit (Stage 4B, §17 S1).
#' UP/DOWN filter the ONE `contrasts_with_direction[[contrast]]` table on `Cluster`;
#' UP+DOWN returns it whole; All returns the DISTINCT `contrasts_without_direction[[
#' contrast]]` table (never a union); a single cluster/pattern filters the one
#' partition/binary table; Compare returns it whole; all_DE returns its one table.
#' Kept VERBATIM (all columns, incl. geneID/Cluster/… for master-detail). NULL if
#' not an ORA unit or no table/rows exist.
#' @noRd
enrichment_ora_table <- function(enrichment, u) {
  if (is.null(enrichment) || is.null(u) || !identical(u$analysis, "ORA")) return(NULL)
  ora <- enrichment$ora
  if (is.null(ora)) return(NULL)
  db_node <- ora[[u$database]]
  if (is.null(db_node)) return(NULL)
  # §5/§6: Compare Contrasts — combine every contrast's direction-agnostic DE
  # (`contrasts_without_direction`) into ONE table, with Cluster = contrast name
  # (the compare group). No recompute, no Up/Down rows, no unions — just the
  # stored DE tables stacked. Only contrasts with actual significant rows appear.
  if (identical(u$group, "contrasts") &&
      identical(.enrich_scalar(u$item), "compare_contrasts")) {
    node <- db_node[["contrasts_without_direction"]]
    if (is.null(node) || !length(node)) return(NULL)
    parts <- list()
    for (ct in names(node)) {
      d <- node[[ct]]
      if (is.data.frame(d) && nrow(d) > 0) { d$Cluster <- ct; parts[[length(parts) + 1L]] <- d }
    }
    if (!length(parts)) return(NULL)
    return(do.call(rbind, parts))
  }
  pg   <- .enrich_ora_payload_group(u$group, u$subgroup)
  node <- db_node[[pg]]
  if (is.null(node)) return(NULL)
  # contrast-keyed nodes hold one table per contrast; the others hold one table.
  keyed_by_contrast <- pg %in% c("contrasts_with_direction", "contrasts_without_direction")
  df <- if (keyed_by_contrast) node[[u$item]] else node
  if (is.null(df) || !is.data.frame(df)) return(NULL)
  filt <- .enrich_ora_cluster_filter(u$group, u$subgroup, u$item)
  if (!is.na(filt) && "Cluster" %in% names(df))
    df <- df[as.character(df$Cluster) == filt, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df
}

#' Resolve the GSEA results table for a resolved GSEA unit, verbatim (NULL for
#' non-GSEA / missing). Layout is simple: `enrichment$gsea[[db]][[ranking]][[
#' contrast]]` is the df directly (group = ranking, item = contrast; no
#' direction split, no Cluster filter). The SAME object feeds the DT and the
#' selection/download so selected-row -> biological ID resolution is stable.
#' @noRd
enrichment_gsea_table <- function(enrichment, u) {
  if (is.null(enrichment) || is.null(u) || !identical(u$analysis, "GSEA")) return(NULL)
  g <- enrichment$gsea
  if (is.null(g)) return(NULL)
  db_node <- g[[u$database]]
  if (is.null(db_node)) return(NULL)
  rk_node <- db_node[[u$group]]                    # group = ranking
  if (is.null(rk_node)) return(NULL)
  df <- rk_node[[u$item]]                           # item = contrast
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
  df
}

#' Resolve the *simplified GO* ORA table for a resolved unit (GO only; NULL for
#' KEGG or when no simplify node exists). Same node-mapping/keying/Cluster-filter
#' logic as enrichment_ora_table(), applied to `enrichment$ora[[db]]$simplify`.
#' @noRd
enrichment_ora_simplify_table <- function(enrichment, u) {
  if (is.null(enrichment) || is.null(u) || !identical(u$analysis, "ORA")) return(NULL)
  db_node <- enrichment$ora[[u$database]]
  sm <- if (is.null(db_node)) NULL else db_node$simplify
  if (is.null(sm)) return(NULL)
  pg   <- .enrich_ora_payload_group(u$group, u$subgroup)
  node <- sm[[pg]]
  if (is.null(node)) return(NULL)
  keyed_by_contrast <- pg %in% c("contrasts_with_direction", "contrasts_without_direction")
  df <- if (keyed_by_contrast) node[[u$item]] else node
  if (is.null(df) || !is.data.frame(df)) return(NULL)
  filt <- .enrich_ora_cluster_filter(u$group, u$subgroup, u$item)
  if (!is.na(filt) && "Cluster" %in% names(df))
    df <- df[as.character(df$Cluster) == filt, , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df
}

#' The EXACT row of an ORA table for a selection (§17 S3/S4). In a COMPARE table
#' the same term `ID` appears once per group, so `id` ALONE is ambiguous — the
#' selection is keyed by (ID + Group), where `group` is the row's internal
#' `Cluster` value. When `group` is NULL/NA (single-group unit) the first ID match
#' is used (the pre-compare behavior). Returns a 1-row df or NULL.
#' @noRd
.enrich_ora_row <- function(df, id, group = NULL) {
  if (is.null(df) || is.null(id)) return(NULL)
  hit <- as.character(df$ID) == as.character(id)
  g <- .enrich_scalar(group)
  if (!is.na(g) && "Cluster" %in% names(df))
    hit <- hit & (as.character(df$Cluster) == g)
  row <- df[hit, , drop = FALSE]
  if (nrow(row) == 0) NULL else row[1, , drop = FALSE]
}

#' Member (overlap) gene IDs for one selected ORA term (split `geneID`). §17 S3/S4:
#' `group` disambiguates a compare table (a term significant in >1 group) so the
#' member genes are the SELECTED group's, never a union across groups.
#' @noRd
enrichment_ora_member_genes <- function(ora_df, id, group = NULL) {
  if (is.null(ora_df) || is.null(id)) return(character(0))
  row <- .enrich_ora_row(ora_df, id, group)
  if (is.null(row) || is.null(row$geneID)) return(character(0))
  g <- unlist(strsplit(paste(row$geneID, collapse = "/"), "/", fixed = TRUE))
  unique(g[!is.na(g) & nzchar(g)])
}

#' The contrasts data.frame (`Contrast_name` column) used to build the per-gene
#' DE-stat rows. Prefers the payload's canonical `env$contrasts` (so the contrast
#' set + ORDER match the Multi Explorer exactly); falls back to the contrast
#' names discovered in `de_stats` (`linearFC.<c>` suffixes) when it is absent.
#' @noRd
.enrich_contrasts_frame <- function(env) {
  if (exists("contrasts", envir = env, inherits = FALSE)) {
    ct <- get("contrasts", envir = env, inherits = FALSE)
    if (is.data.frame(ct) && "Contrast_name" %in% names(ct)) return(ct)
  }
  ds <- if (exists("de_stats", envir = env, inherits = FALSE))
    get("de_stats", envir = env, inherits = FALSE) else NULL
  nm <- if (is.data.frame(ds))
    sub("^linearFC\\.", "", grep("^linearFC\\.", names(ds), value = TRUE)) else
    character(0)
  data.frame(Contrast_name = nm, stringsAsFactors = FALSE)
}

#' Member-gene / core-gene DE-statistics table, built to the SAME long shape and
#' via the SAME per-gene helper (`get_gene_de_stats_for_table`) as the Multi
#' Explorer "Differential Expression Statistics" table, so the two present
#' identically (columns `Gene`, `Contrast`, `linearFC`, `pvalue`, `padj`; the
#' render helper drops `Contrast` when a single contrast is shown). Genes absent
#' from de_stats are KEPT as NA rows (complete gene set, never truncated).
#'
#' `contrast` (pass-3, item 1): a single-contrast unit passes its contrast name
#' to restrict to that ONE contrast; NA keeps ALL contrasts (partition / binary /
#' all_DE / any_contrast / compare). Defensive: a name with no matching contrast
#' falls back to all contrasts rather than an empty table. `de_stats` is never
#' modified. NULL when de_stats is unavailable.
#' @noRd
ora_member_de_table <- function(env, member_ids, contrast = NA_character_) {
  if (length(member_ids) == 0 ||
      !exists("de_stats", envir = env, inherits = FALSE)) return(NULL)
  ds <- get("de_stats", envir = env, inherits = FALSE)
  if (is.null(ds) || !is.data.frame(ds) || !"FeatureID" %in% names(ds)) return(NULL)
  cf <- .enrich_contrasts_frame(env)
  ct <- .enrich_scalar(contrast)
  if (!is.na(ct) && ct %in% cf$Contrast_name)      # restrict to the one contrast
    cf <- cf[cf$Contrast_name == ct, , drop = FALSE]
  if (nrow(cf) == 0) return(NULL)
  empty_one <- data.frame(Contrast = cf$Contrast_name, linearFC = NA_real_,
                          pvalue = NA_real_, padj = NA_real_,
                          stringsAsFactors = FALSE)
  parts <- lapply(member_ids, function(g) {
    one <- tryCatch(get_gene_de_stats_for_table(ds, cf, g),
                    error = function(e) NULL)
    # absent gene (or no usable contrast rows) -> NA rows, so no gene is dropped
    if (is.null(one) || nrow(one) == 0) one <- empty_one
    cbind(Gene = g, one, stringsAsFactors = FALSE)
  })
  m <- do.call(rbind, parts)
  rownames(m) <- NULL
  m
}

#' Rows/page for the paginated enrichment DE tables (§21). Display-only.
ENRICH_DE_PAGE_LEN <- 20L

#' Display-shaped DE-stat df shared by the render + the download so both match.
#' SINGLE contrast (one unique `Contrast`): drop it -> `Gene, linearFC, pvalue,
#' padj` (the Multi Explorer single-contrast layout). MULTIPLE contrasts (§21):
#' PIVOT to a WIDE table — ONE row per Gene, a `linearFC.<c> / pvalue.<c> /
#' padj.<c>` block per contrast (contrast order preserved). Genes with missing
#' statistics stay as NA; no gene or contrast is dropped, never gene×contrast rows.
#' @noRd
.enrich_de_display_df <- function(m) {
  if (is.null(m) || !"Contrast" %in% names(m)) return(m)
  cts <- unique(as.character(m$Contrast))
  if (length(cts) <= 1L) { m$Contrast <- NULL; return(m) }
  genes <- unique(as.character(m$Gene))            # first-appearance (member) order
  out <- data.frame(Gene = genes, stringsAsFactors = FALSE, check.names = FALSE)
  for (ct in cts) {
    sub <- m[as.character(m$Contrast) == ct, , drop = FALSE]
    idx <- match(genes, as.character(sub$Gene))     # NA where a gene is absent
    out[[paste0("linearFC.", ct)]] <- sub$linearFC[idx]
    out[[paste0("pvalue.",   ct)]] <- sub$pvalue[idx]
    out[[paste0("padj.",     ct)]] <- sub$padj[idx]
  }
  out
}

#' TRUE when a display df is the multi-contrast WIDE layout (per-contrast stat
#' columns like `linearFC.<c>`) rather than the single-contrast narrow layout.
#' @noRd
.enrich_de_is_wide <- function(d) {
  !is.null(d) && any(grepl("^linearFC\\.", names(d)))
}

#' TRUE when a member/core DE table (long builder output) resolves to the SINGLE-
#' contrast coloured layout — i.e. it is the only table that carries the Up/Down
#' legend (§21). NULL/empty -> FALSE.
#' @noRd
.enrich_de_single <- function(m) {
  !is.null(m) && nrow(m) > 0 && !.enrich_de_is_wide(.enrich_de_display_df(m))
}

#' Render the MULTI-contrast WIDE DE table (§21): one row per gene, a per-contrast
#' `linearFC/pvalue/padj` block for every contrast. NO Up/Down row colouring — a
#' gene may move in opposite directions across contrasts, so a single row colour
#' would be ambiguous. Same number formatting (`format_scientific_threshold`) and
#' cell-border styling as the single-contrast table, paginated (20/page) with
#' search + sort. `de_stats` is untouched; genes/contrasts are never dropped.
#' @noRd
render_enrichment_wide_de_dt <- function(d, page_length = ENRICH_DE_PAGE_LEN,
                                         selection = "none") {
  disp <- d
  num_cols <- grep("^(linearFC|pvalue|padj)\\.", names(disp), value = TRUE)
  for (cn in num_cols) disp[[cn]] <- format_scientific_threshold(disp[[cn]])
  DT::datatable(
    disp, rownames = FALSE, selection = selection, class = "cell-border",
    options = list(dom = "ftip", paging = TRUE, pageLength = page_length,
                   searching = TRUE, ordering = TRUE, order = list(),
                   info = TRUE, scrollX = TRUE, stripeClasses = c("", "")))
}

#' Render the enrichment member/core-gene DE table.
#'  * SINGLE contrast -> the EXACT Multi Explorer presentation via
#'    `render_de_stats_dt()` (columns, numeric / adjusted-p formatting, header
#'    tooltips, Up/Down row colouring) using the dataset's own cutoffs, now
#'    PAGINATED (20/page) + search + sort (§21).
#'  * MULTIPLE contrasts -> the WIDE table (`render_enrichment_wide_de_dt`), same
#'    formatting/paging but no row colouring.
#' @noRd
render_enrichment_member_dt <- function(env, m, selection = "none") {
  d <- .enrich_de_display_df(m)
  if (.enrich_de_is_wide(d))
    return(render_enrichment_wide_de_dt(d, selection = selection))
  padj0 <- if (exists("padj_cutoff", envir = env, inherits = FALSE))
    get("padj_cutoff", envir = env, inherits = FALSE) else 0.05
  fc0 <- if (exists("fc_cutoff", envir = env, inherits = FALSE))
    get("fc_cutoff", envir = env, inherits = FALSE) else 1
  render_de_stats_dt(d, sig_thr = padj0, fc_thr = fc0,
                     sig_default = padj0, fc_default = fc0, sig_mode = "padj",
                     paginate = TRUE, page_length = ENRICH_DE_PAGE_LEN,
                     searching = TRUE, ordering = TRUE, selection = selection,
                     colour_rows = FALSE)   # §1: no Up/Down row colouring in Enrichment
}

#' ORA dotplot (top-N terms by adjusted p-value). Static ggplot.
#' @noRd
#' ORA Enrichment Plot (§23 B) — horizontal bars. X = gene Count, fill =
#' display-only -log10(FDR), bar annotation "<count>/<term size> · FE <fold>".
#' SELECTION is the top-N terms by adjusted p-value (most significant); the
#' SELECTED terms are then DISPLAYED ordered by Count DESCENDING (highest Count on
#' top) — selection and display order are intentionally different. Nothing here
#' recomputes biology: every value is read from the ORA table verbatim.
#' @noRd
plot_ora_enrichment_bar <- function(df, max_terms = ENRICH_BAR_TOP_N) {
  # (1) SELECT top-N by FDR ascending
  d <- df[order(df$p.adjust), , drop = FALSE]
  d <- utils::head(d, max_terms)
  d$.count <- suppressWarnings(as.numeric(d$Count))
  # (2) DISPLAY order = Count descending (tie-break FDR asc, then Description)
  d <- d[order(-d$.count, suppressWarnings(as.numeric(d$p.adjust)),
               as.character(d$Description)), , drop = FALSE]
  d$.nlp   <- .enrich_nlp(d$p.adjust)
  tsz <- .enrich_term_size(d); fe <- .enrich_fold(d)
  d$.ann <- ifelse(
    is.na(tsz),
    ifelse(is.na(fe), format(d$.count),
           sprintf("%s · FE %.1f", format(d$.count), fe)),
    ifelse(is.na(fe), sprintf("%s/%s", format(d$.count), tsz),
           sprintf("%s/%s · FE %.1f", format(d$.count), tsz, fe)))
  d$.lab <- .enrich_term_labels(d$Description)
  d$.lab <- factor(d$.lab, levels = rev(d$.lab))   # most significant on top
  lim <- .enrich_nlp_limits(d$.nlp)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$.count, y = .data$.lab,
                                  fill = .data$.nlp)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_text(ggplot2::aes(label = .data$.ann), hjust = -0.06,
                       size = ENRICH_BAR_ANNOT_SIZE, colour = "grey25") +
    ggplot2::scale_fill_viridis_c(name = "-log10(adj p)", option = "D",
                                  limits = lim, oob = scales::squish) +
    # generous right-hand headroom so the enlarged "<count>/<size> · FE x.x"
    # labels are never clipped, even for 4-digit counts (§16.40)
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, ENRICH_BAR_EXPAND))) +
    ggplot2::labs(x = "Gene count", y = NULL) +
    .enrich_plot_theme()
}

# ============================================================================
# §17 S3 — COMPARE PLOTS (interactive counterpart of the pipeline R4B artifacts).
# UP+DOWN -> diverging bars (D9); Compare Clusters/Patterns -> dot matrix (D10).
# Shiny Top-N control (D11): Top 20 / 50 / All, default 20 — TOP-N PER GROUP by
# adj p, then UNION across groups (a term appears once on Y, present in the groups
# that contain it). Viz-only — never truncates the table/payload/downloads.
# ============================================================================

ENRICH_COMPARE_TOPN_CHOICES <- c("Top 20" = "20", "Top 50" = "50", "All" = "all")

#' Map the Top-N radio value ("20"/"50"/"all"/NULL) to a numeric per-group cap
#' (Inf = All; default 20 before the control has reported a value). §17 S3 (D11).
#' @noRd
.enrich_topn_val <- function(x) {
  if (is.null(x)) return(ENRICH_BAR_TOP_N)
  if (identical(as.character(x), "all")) return(Inf)
  n <- suppressWarnings(as.numeric(x)); if (is.na(n)) ENRICH_BAR_TOP_N else n
}

#' Kept term IDs = union over groups of each group's top-N terms by adj p-value.
#' `top_n` NULL / non-finite / "all" -> every term. §17 S3 (D11).
#' @noRd
.enrich_compare_topn_ids <- function(df, top_n) {
  ids_all <- unique(as.character(df$ID))
  n <- suppressWarnings(as.numeric(top_n))
  if (length(n) == 0 || is.na(n) || !is.finite(n)) return(ids_all)
  parts <- split(seq_len(nrow(df)), as.character(df$Cluster))
  keep <- unlist(lapply(parts, function(ix) {
    sub <- df[ix, , drop = FALSE]
    o <- order(suppressWarnings(as.numeric(sub$p.adjust)))
    as.character(sub$ID[o])[seq_len(min(n, length(o)))]
  }), use.names = FALSE)
  unique(keep)
}
#' Compare df filtered to the Top-N union of terms (all group rows for kept IDs).
#' @noRd
.enrich_compare_topn_df <- function(df, top_n) {
  keep <- .enrich_compare_topn_ids(df, top_n)
  df[as.character(df$ID) %in% keep, , drop = FALSE]
}
#' Which compare plot a unit uses. §23 G: ALWAYS the compare dot matrix — Compare
#' Up & Down now uses the same dot matrix as Compare Clusters/Patterns (the
#' diverging-bars figure is retired as the canonical up/down comparison).
#' @noRd
.enrich_compare_plot_kind <- function(u) "dotmatrix"

# UI refinement (post-S5): compare plots size DYNAMICALLY to the number of terms
# shown so every label stays readable — the plot grows taller than the page and
# scrolls, especially at Top-N = All (matching the pipeline's full static figure).
#' Distinct terms shown after the Top-N union filter.
#' @noRd
.enrich_compare_n_terms <- function(df, top_n) length(.enrich_compare_topn_ids(df, top_n))
#' Pixel height for the diverging bar plot (one readable row band per term).
#' @noRd
.enrich_compare_bar_height <- function(n_terms)
  as.integer(max(420, 30 * max(n_terms, 1L) + 150))
#' Pixel height / width for the dot matrix (per term / per group).
#' @noRd
.enrich_dotmatrix_height <- function(n_terms)
  as.integer(max(420, 26 * max(n_terms, 1L) + 170))
# §13: wider per-group horizontal spacing so group/cluster labels are clearly
# readable and dots are not crowded — width scales with the number of groups, so a
# many-cluster comparison gets more horizontal room than a two-group Up/Down one.
.enrich_dotmatrix_width <- function(n_groups)
  as.integer(max(640, 150 * max(n_groups, 1L) + 420))

# §18 D18.1: the single-group ORA/GSEA Enrichment Plot gains the SAME Top 20/50/
# All viewing control as the compare plots (ENRICH_COMPARE_TOPN_CHOICES, default
# "20"), and — like them — sizes dynamically to the number of terms shown so the
# figure stays readable at Top-N = All (grows/scrolls vertically). Viz-only: the
# results table, payload and downloads are never truncated by this control.
#' Number of terms/pathways actually shown for the current Top-N selection.
#' @noRd
.enrich_bar_shown_n <- function(df, top_n)
  as.integer(min(nrow(df), .enrich_topn_val(top_n)))
#' Pixel height for the single-group enrichment bar (one readable band per term).
#' @noRd
.enrich_bar_height <- function(n_terms)
  as.integer(max(360, 26 * max(n_terms, 1L) + 150))

# §4/§5: explicit, unit-scoped memory of the inner content tab for each ORA result
# set (Full Results vs GO Simplified). Pure helpers so the behavior is unit-testable
# independently of the reactive session; wire_enrichment() drives them from observers.
#' Fresh memory — both result sets start at the first content tab ("Results").
#' @noRd
.enrich_inner_tab_init <- function(default = "Results Table")
  list(full = default, simplified = default)
#' Normalize a result-set key to "full" / "simplified".
#' @noRd
.enrich_inner_tab_key <- function(active_set)
  if (identical(active_set, "simplified")) "simplified" else "full"
#' Remember `tab` as the inner tab last used in the active result set. A NULL/empty
#' tab leaves the memory unchanged (never erases a remembered choice).
#' @noRd
.enrich_inner_tab_remember <- function(mem, active_set, tab) {
  if (is.null(mem)) mem <- .enrich_inner_tab_init()
  if (is.null(tab) || !nzchar(tab)) return(mem)
  mem[[.enrich_inner_tab_key(active_set)]] <- tab
  mem
}
#' The inner tab to restore for `active_set` (first visit / unknown -> default).
#' @noRd
.enrich_inner_tab_restore <- function(mem, active_set, default = "Results Table") {
  if (is.null(mem)) return(default)
  t <- mem[[.enrich_inner_tab_key(active_set)]]
  if (is.null(t) || !nzchar(t)) default else t
}

# §9/§10: enrichment download filenames carry the unit context and (for Top-N
# figures) the selected variant (top20/top50/all) so Top 20/50/All downloads never
# overwrite each other — and never a download date. Pure helpers (no reactive deps)
# so the naming rule is unit-testable; wire_enrichment()'s .dl_name delegates here.
#' Build an enrichment download filename. `variant` is appended only when non-empty.
#' @noRd
.enrich_dl_filename <- function(kind, ctx, ext, variant = NULL) {
  var <- if (is.null(variant) || !nzchar(variant)) "" else paste0("_", variant)
  sprintf("enrichment_%s_%s%s.%s", kind, ctx, var, ext)
}
#' Top-N variant tag for a filename from a raw Top-N input value
#' ("20"->top20, "50"->top50, "all"/"Inf"->all; NULL/NA -> the on-screen default).
#' @noRd
.enrich_dl_topn_tag <- function(x) {
  if (is.null(x) || (length(x) == 1 && is.na(x))) return("top20")
  xs <- as.character(x)
  if (identical(xs, "all") || identical(xs, "Inf")) "all" else paste0("top", xs)
}

# §D25.2: the app's canonical loading indicator is `shinycssloaders::withSpinner`
# (type 4; heatmaps use the app's #2c3e50 colour, e.g. ui.R de_heatmap). Reuse it —
# never a bespoke spinner — on every enrichment plot/heatmap output so the module
# matches the rest of the site while a figure renders. Centralized here so there is
# ONE wrapper (no nested/duplicate spinners) and it degrades to the bare output when
# the package is unavailable (e.g. headless tests). `heatmap=TRUE` selects the
# heatmap colour used elsewhere in the app.
#' @noRd
.enrich_spin <- function(output, heatmap = FALSE) {
  if (!requireNamespace("shinycssloaders", quietly = TRUE)) return(output)
  if (isTRUE(heatmap))
    shinycssloaders::withSpinner(output, type = 4, color = "#2c3e50")
  else
    shinycssloaders::withSpinner(output, type = 4)
}

#' §17 S4 (D13): per-gene group label for a compare unit's GROUPED heatmap, from
#' the term `geneID` membership + `Cluster` (self-contained; the compare groups
#' are gene-disjoint by construction — up vs down DE genes, one cluster/pattern per
#' gene). Returns a named character vector gene -> group label (first group wins on
#' the rare overlap). `group_type` gives the display labels (Up/Down, Cluster N…).
#' @noRd
.enrich_compare_gene_groups <- function(df, group_type = NULL) {
  if (is.null(df) || is.null(df$geneID) || is.null(df$Cluster)) return(character(0))
  gl    <- .enrich_row_group_labels(df$Cluster, group_type)
  glist <- .enrich_gene_lists(df)
  out <- character(0)
  for (i in seq_along(glist)) {
    new <- setdiff(glist[[i]], names(out))
    if (length(new)) out[new] <- gl[i]
  }
  out
}

#' §17 S3 / D9: UP+DOWN diverging horizontal bar plot. Y=term, RIGHT=UP / LEFT=DOWN,
#' bar length=Gene Count, fill=-log10(FDR) (a magnitude gradient, NOT a red/blue
#' direction encoding), per-bar annotation "count/term size · FE". ggplot.
#' @noRd
plot_ora_compare_diverging <- function(df, top_n = ENRICH_BAR_TOP_N) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
  d <- .enrich_compare_topn_df(df, top_n)
  d <- d[as.character(d$Cluster) %in% c("up", "down"), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d$.count  <- suppressWarnings(as.numeric(d$Count))
  d$.signed <- ifelse(as.character(d$Cluster) == "up", d$.count, -d$.count)
  d$.nlp    <- .enrich_nlp(d$p.adjust)
  tsz <- .enrich_term_size(d); fe <- .enrich_fold(d)
  d$.ann <- ifelse(is.na(tsz), format(d$.count), sprintf("%s/%s", format(d$.count), tsz))
  d$.ann <- ifelse(is.na(fe), d$.ann, sprintf("%s · FE %.1f", d$.ann, fe))
  d$.lab <- .enrich_term_labels(d$Description)
  ord <- tapply(d$.nlp, d$.lab, max)                       # strongest term on top
  d$.lab <- factor(d$.lab, levels = names(sort(ord)))
  lim <- .enrich_nlp_limits(d$.nlp)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$.signed, y = .data$.lab, fill = .data$.nlp)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
    ggplot2::geom_text(ggplot2::aes(label = .data$.ann,
                                    hjust = ifelse(.data$.signed >= 0, -0.06, 1.06)),
                       size = ENRICH_BAR_ANNOT_SIZE, colour = "grey25") +
    ggplot2::scale_x_continuous(labels = function(x) abs(x),
                                expand = ggplot2::expansion(mult = ENRICH_BAR_EXPAND)) +
    ggplot2::scale_fill_viridis_c(name = "-log10(adj p)", option = "D",
                                  limits = lim, oob = scales::squish) +
    ggplot2::labs(x = "Gene count   (DOWN <-  |  -> UP)", y = NULL) +
    .enrich_plot_theme()
}

#' §17 S3 / D10: comparison dot matrix (interactive). Y=term, X=group, dot size =
#' Gene Count, dot colour = -log10(FDR). Hover: Term, Group, Gene Count, GeneRatio,
#' FDR, Fold Enrichment. Native plotly. `group_type` maps Cluster -> display group.
#' @noRd
plot_ora_compare_dotmatrix <- function(df, top_n = ENRICH_BAR_TOP_N,
                                       group_type = NULL, source_id = NULL) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
  d <- .enrich_compare_topn_df(df, top_n)
  if (nrow(d) == 0) return(NULL)
  d$.grp   <- .enrich_row_group_labels(d$Cluster, group_type)
  # §23 F: the Y-axis ROW IDENTITY is the biological Term ID — a term shared across
  # groups is ONE row with multiple points, and two different IDs sharing a
  # Description never collapse. Descriptions remain the display tick labels
  # (make.unique only on the ORDERED per-ID labels, purely to keep ticks distinct).
  d$.id    <- as.character(if (!is.null(d$ID)) d$ID else d$Description)
  d$.count <- suppressWarnings(as.numeric(d$Count))
  d$.nlp   <- .enrich_nlp(d$p.adjust)
  fe <- .enrich_fold(d)
  gr <- if (!is.null(d$GeneRatio)) as.character(d$GeneRatio) else rep(NA_character_, nrow(d))
  ord   <- tapply(d$.nlp, d$.id, max)
  yl_id <- names(sort(ord))                                # strongest term at top
  id2lab <- tapply(as.character(d$Description), d$.id, function(x) x[1])
  yl  <- make.unique(unname(id2lab[yl_id]))               # display labels (Descriptions)
  xg  <- sort(unique(d$.grp))
  yi  <- match(d$.id, yl_id); xi <- match(d$.grp, xg)
  # Dot size = Gene Count, mapped into a readable pixel range. Keep the mapping so
  # the size legend below uses the SAME transform for its reference markers.
  crng <- range(d$.count, na.rm = TRUE); srng <- c(7, 26)
  to_size <- function(v) if (crng[1] == crng[2]) rep(mean(srng), length(v))
                         else scales::rescale(v, from = crng, to = srng)
  sz  <- to_size(d$.count)
  hover <- paste0("Term: ", d$Description, "<br>Group: ", d$.grp,
                  "<br>Gene count: ", d$.count, "<br>Gene ratio: ", gr,
                  "<br>FDR: ", signif(d$p.adjust, 3),
                  "<br>Fold enrichment: ", signif(fe, 2))
  n_terms <- length(yl); n_grp <- length(xg)
  p <- plotly::plot_ly(                                    # dynamic size on the widget
    x = xi, y = yi, type = "scatter", mode = "markers", source = source_id,
    height = .enrich_dotmatrix_height(n_terms), width = .enrich_dotmatrix_width(n_grp),
    marker = list(size = sz, color = d$.nlp, colorscale = "Viridis",
                  showscale = TRUE,
                  colorbar = list(title = "-log10(FDR)", len = 0.55, y = 0.78),
                  line = list(width = 0.5, color = "grey40")),
    text = hover, hoverinfo = "text", showlegend = FALSE)
  # Item 5: a Gene Count SIZE legend — reference markers at representative counts,
  # drawn off the data range (clipped by the fixed axes) so only the legend shows.
  brks <- unique(round(stats::quantile(d$.count, c(0.1, 0.5, 0.95),
                                       names = FALSE, na.rm = TRUE)))
  brks <- brks[brks > 0]
  for (b in brks)
    p <- plotly::add_trace(p, x = -99, y = -99, type = "scatter", mode = "markers",
                           marker = list(size = to_size(b), color = "grey55",
                                         line = list(width = 0.5, color = "grey40")),
                           name = paste0(b, " genes"), showlegend = TRUE,
                           hoverinfo = "skip", inherit = FALSE)
  plotly::layout(p,
    legend = list(title = list(text = "Gene count"), x = 1.02, y = 0.30,
                  itemsizing = "trace", bgcolor = "rgba(255,255,255,0.6)"),
    xaxis = list(tickmode = "array", tickvals = seq_along(xg), ticktext = xg,
                 title = "", zeroline = FALSE, tickfont = list(size = 13),
                 tickangle = if (n_grp > 6) 45 else 0,
                 range = c(0.5, length(xg) + 0.5)),
    yaxis = list(tickmode = "array", tickvals = seq_along(yl), ticktext = yl,
                 title = "", automargin = TRUE, tickfont = list(size = 13),
                 range = c(0.5, length(yl) + 0.5)),
    margin = list(l = 10, r = 40))
}

#' §17 (post-S5 item 11): STATIC ggplot dot matrix for the FULL-figure download —
#' every term, dot size = Gene Count (with a size legend), colour = -log10(FDR).
#' @noRd
plot_ora_compare_dotmatrix_gg <- function(df, group_type = NULL, top_n = Inf) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
  d <- .enrich_compare_topn_df(df, top_n)
  if (nrow(d) == 0) return(NULL)
  d$.grp   <- .enrich_row_group_labels(d$Cluster, group_type)
  # §23 F: row identity = Term ID; Description is the display label.
  d$.id    <- as.character(if (!is.null(d$ID)) d$ID else d$Description)
  d$.count <- suppressWarnings(as.numeric(d$Count))
  d$.nlp   <- .enrich_nlp(d$p.adjust)
  ord    <- tapply(d$.nlp, d$.id, max)
  yl_id  <- names(sort(ord))
  id2lab <- tapply(as.character(d$Description), d$.id, function(x) x[1])
  d$.id  <- factor(d$.id, levels = yl_id)
  d$.grp <- factor(d$.grp, levels = sort(unique(d$.grp)))
  lim <- .enrich_nlp_limits(d$.nlp)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$.grp, y = .data$.id,
                                  size = .data$.count, colour = .data$.nlp)) +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(labels = function(ids) unname(id2lab[ids])) +
    ggplot2::scale_colour_viridis_c(name = "-log10(FDR)", limits = lim, oob = scales::squish) +
    ggplot2::scale_size_continuous(name = "Gene count", range = c(1.5, 9)) +
    ggplot2::labs(x = NULL, y = NULL) + .enrich_plot_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Full-figure (ALL significant terms) compare plot for download — diverging bars
#' for UP+DOWN, else the static dot matrix. Independent of the on-screen Top-N.
#' Returns list(plot, width, height) sized to fit every term, or NULL.
#' @noRd
.enrich_compare_download_figure <- function(df, u, top_n = Inf) {
  if (is.null(df) || is.null(u)) return(NULL)
  # §11: the downloaded figure MATCHES the on-screen Top-N selection (Top 20/50/All)
  # — same term union as the displayed dot matrix, not always "all terms".
  d <- .enrich_compare_topn_df(df, top_n)
  if (nrow(d) == 0) return(NULL)
  n_terms <- length(unique(as.character(d$ID)))
  if (identical(.enrich_compare_plot_kind(u), "diverging")) {
    p <- plot_ora_compare_diverging(d, top_n = Inf)
    h <- max(5, 0.28 * n_terms + 2); w <- 11
  } else {
    p <- plot_ora_compare_dotmatrix_gg(d, u$group_type, top_n = Inf)
    n_grp <- length(unique(as.character(d$Cluster)))
    # §13: wider per-group download width so group labels stay readable.
    h <- max(5, 0.26 * n_terms + 2); w <- max(8, 1.7 * n_grp + 5.5)
  }
  if (is.null(p)) NULL else list(plot = p, width = w, height = h)
}

#' Build shared-gene matrices for a unit's terms: gene x term membership (0/1)
#' and term x term Jaccard overlap, from the `geneID` column.
#' @noRd
#' Shared-gene structure for one ORA unit. VECTORIZED (§16.39): the pairwise
#' intersections come from `crossprod()` on the binary gene x term membership
#' matrix (|A n B| = t(G) %*% G, |A u B| = |A| + |B| - |A n B|), replacing the
#' former O(n^2) R double loop over union()/intersect().
#'
#' §18 D18.3 — the shared-terms overlap metric is standardized on Sørensen–Dice
#' `100 · 2|A∩B|/(|A|+|B|)` (percent 0–100), matching the pipeline term×term CSV
#' EXACTLY. `dice` is what the display now colours by, what clustering orders by,
#' and what both downloads export. `term_term` (Jaccard) is retained only as a
#' legacy field for the regression test; nothing user-facing reads it any more.
#' `shared` (raw count) is kept for the hover tooltip. `ord_sig` is the
#' significance order used for the top-N default.
#' @noRd
build_ora_shared_matrices <- function(df) {
  labs <- .enrich_term_labels(df$Description)
  gl <- .enrich_gene_lists(df)
  all_genes <- sort(unique(unlist(gl)))
  g2t <- matrix(0L, nrow = length(all_genes), ncol = length(labs),
                dimnames = list(all_genes, labs))
  for (j in seq_along(gl)) if (length(gl[[j]])) g2t[gl[[j]], j] <- 1L
  inter <- crossprod(g2t)                       # |A n B| for every term pair
  sz    <- diag(inter)                          # |A| (overlap genes per term)
  uni   <- outer(sz, sz, "+") - inter           # |A u B|
  jac   <- ifelse(uni == 0, 0, inter / uni)     # legacy Jaccard (regression only)
  szsum <- outer(sz, sz, "+")                   # |A| + |B|
  dice  <- ifelse(szsum == 0, 0, 100 * 2 * inter / szsum)   # Sørensen–Dice %
  dimnames(inter) <- dimnames(jac) <- dimnames(dice) <- list(labs, labs)
  # Top-N default is by FDR; a table without p.adjust keeps its own row order
  # (the matrix itself never depended on p.adjust and still must not).
  ord <- if (is.null(df$p.adjust)) seq_len(nrow(df)) else order(df$p.adjust)
  list(gene_term = g2t, term_term = jac, dice = dice, shared = inter,
       term_sizes = sz, n_genes = length(all_genes), n_terms = length(labs),
       ord_sig = ord, full_labels = as.character(df$Description))
}

#' Deterministic cluster order for the shared-terms matrix: average-linkage
#' hclust on the Sørensen–Dice DISTANCE (1 - Dice fraction), so redundant GO
#' terms fall into visible blocks (§18 D18.3 — clustering is by Dice, matching
#' the display metric). Accepts the Dice matrix on EITHER the 0–1 fraction or the
#' 0–100 percent scale (percent is detected and rescaled) so the distance is a
#' proper 0–1 dissimilarity. Identity order for < 3 terms or on any failure. The
#' underlying values are never altered — this only permutes rows/columns.
#' @noRd
.enrich_term_order <- function(sim) {
  n <- ncol(sim)
  if (is.null(n) || n < 3) return(seq_len(max(n, 0)))
  frac <- if (length(sim) && max(sim, na.rm = TRUE) > 1) sim / 100 else sim
  h <- tryCatch(stats::hclust(stats::as.dist(1 - frac), method = "average"),
                error = function(e) NULL)
  if (is.null(h)) seq_len(n) else h$order
}

#' Legibility-driven height for the shared-terms matrix: ~14 px per term plus a
#' margin allowance, clamped to sane bounds (§16.39).
#' @noRd
.enrich_t2t_height <- function(n_terms) {
  h <- as.integer(round(n_terms * ENRICH_T2T_PX_PER_TERM)) + 170L
  max(ENRICH_T2T_MIN_HEIGHT, min(ENRICH_T2T_MAX_HEIGHT, h))
}

#' Matching width rule (§16.40). The matrix is symmetric, so the x axis needs the
#' same per-term budget as the y axis or the (rotated) column labels crush
#' together exactly when "Show all terms" is on. The larger constant offset
#' covers the long left-hand y labels. When the result exceeds the page the
#' panel scrolls horizontally rather than squeezing the plot.
#' @noRd
.enrich_t2t_width <- function(n_terms) {
  w <- as.integer(round(n_terms * ENRICH_T2T_PX_PER_TERM)) + 320L
  max(ENRICH_T2T_MIN_WIDTH, min(ENRICH_T2T_MAX_WIDTH, w))
}

#' Which terms the shared-terms matrix shows (§18 D18.3): the top-N terms by
#' adjusted p-value (significance-based selection), or ALL of them when `top_n`
#' is Inf/NA/NULL or the unit is already <= N. Returned in the ORIGINAL matrix
#' order; clustering is applied afterwards. Top-N is a VIEWING choice only — the
#' full-table view and the table download are never restricted by it.
#' @noRd
.enrich_t2t_keep <- function(s, top_n = ENRICH_T2T_TOP_N) {
  n <- s$n_terms
  k <- suppressWarnings(as.numeric(top_n))
  if (length(k) == 0 || is.na(k) || !is.finite(k) || n <= k) return(seq_len(n))
  sort(s$ord_sig[seq_len(k)])
}

#' Shared Terms matrix (§18 D18.3): term x term overlap heatmap coloured by the
#' Sørensen–Dice percent `100 · 2|A∩B|/(|A|+|B|)` — the SAME metric, scale (0–100)
#' and colouring intent as the pipeline term×term heatmap. Terms are cluster-
#' ordered by Dice, height scales with the number of terms, and every displayed
#' term is explicitly ticked (no Plotly auto-thinning). Hover reports the exact
#' Dice percent AND the raw shared-gene count. `top_n` (Inf = all) selects the
#' most significant terms to show — a VIEWING choice only.
#' @noRd
plot_ora_term_overlap <- function(s, top_n = ENRICH_T2T_TOP_N) {
  keep <- .enrich_t2t_keep(s, top_n)
  D <- s$dice[keep, keep, drop = FALSE]           # Sørensen–Dice percent (0–100)
  M <- s$shared[keep, keep, drop = FALSE]         # raw shared-gene count (hover)
  ord <- .enrich_term_order(D)
  D <- D[ord, ord, drop = FALSE]; M <- M[ord, ord, drop = FALSE]
  labs <- rownames(D)
  n <- length(labs)
  h <- .enrich_t2t_height(n); w <- .enrich_t2t_width(n)
  # Colour is the Dice percent itself (no transform); hover carries both the
  # Dice percent and the exact integer shared-gene count. plotly drops
  # `customdata` on heatmap traces, so the text matrix + hoverinfo="text" is the
  # reliable route.
  txt <- matrix(sprintf("%s<br>%s<br>Sørensen–Dice: %.1f%%<br>shared genes: %s",
                        matrix(rep(labs, times = n), n, n),   # row (y) term
                        matrix(rep(labs, each  = n), n, n),   # column (x) term
                        D, M), n, n)
  ax <- function(angle = 0) list(
    tickmode = "array", tickvals = labs, ticktext = labs,
    tickangle = angle, tickfont = list(size = 10), automargin = TRUE)
  plotly::plot_ly(
    z = D, x = labs, y = labs, type = "heatmap",
    text = txt, hoverinfo = "text", colorscale = "Viridis", zmin = 0, zmax = 100,
    height = h, width = w,
    colorbar = list(title = list(text = "Sørensen–\nDice (%)"))) |>
    plotly::layout(xaxis = ax(-40), yaxis = ax(0),
                   margin = list(l = 20, b = 20)) |>
    plotly::config(displaylogo = FALSE,
                   modeBarButtonsToRemove = list("sendDataToCloud"))
}

#' DT for the Shared Terms FULL TABLE view (§18 D18.4) — the COMPLETE term × term
#' Sørensen–Dice percent matrix (all terms, never truncated by the Top-N view),
#' scrollable and searchable. Terms are cluster-ordered by Dice so the on-screen
#' table and the heatmap read in the same order; the first column is the term
#' name. Values are numeric percents rounded for display (client sort stays
#' numeric). The download serves the same complete matrix.
#' @noRd
render_term_overlap_dt <- function(s) {
  D <- s$dice
  ord <- .enrich_term_order(D)
  D <- D[ord, ord, drop = FALSE]
  d <- cbind(Term = rownames(D),
             as.data.frame(round(D, 1), check.names = FALSE))
  DT::datatable(
    d, rownames = FALSE, selection = "none",
    class = "cell-border stripe hover compact",
    options = list(pageLength = 20, dom = "ftip", scrollX = TRUE))
}

#' STATIC ggplot version of the Shared Terms Dice heatmap, for the figure
#' download (§18 D18.3). Same metric (Sørensen–Dice %), same Dice cluster order,
#' same 0–100 viridis scale and Top-N selection as the on-screen plotly view, so
#' the downloaded figure MATCHES the current selection. Returns NULL for < 2
#' terms. Sized (inches) to the number of terms shown.
#' @noRd
plot_ora_term_overlap_static <- function(s, top_n = ENRICH_T2T_TOP_N) {
  if (is.null(s) || is.null(s$dice) || s$n_terms < 2) return(NULL)
  keep <- .enrich_t2t_keep(s, top_n)
  D <- s$dice[keep, keep, drop = FALSE]
  ord <- .enrich_term_order(D)
  D <- D[ord, ord, drop = FALSE]
  labs <- rownames(D)
  long <- data.frame(
    row = factor(rep(labs, times = ncol(D)), levels = rev(labs)),
    col = factor(rep(labs, each = nrow(D)),  levels = labs),
    dice = as.vector(D))
  ggplot2::ggplot(long, ggplot2::aes(.data$col, .data$row, fill = .data$dice)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "Sørensen–Dice (%)", option = "D",
                                  limits = c(0, 100)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 7),
      panel.grid = ggplot2::element_blank())
}

#' Inches size for the static Dice heatmap download (square, grows with N).
#' @noRd
.enrich_t2t_static_in <- function(n) max(5, min(40, 0.28 * max(n, 1L) + 2.5))

#' Shared Genes (§16.39) — replaces the unreadable genes x terms binary matrix.
#' Answers "which genes recur across the most enriched terms?": one row per
#' gene, X = number of enriched terms containing it, fill = the BEST (minimum)
#' adjusted p-value among those terms, as a display-only -log10. Bounded to the
#' top `top_n` genes, so it always renders regardless of unit size — no
#' gene-count guard. Ordering is fully deterministic: recurrence desc, then best
#' FDR asc, then gene id asc (stable tie-break). All underlying data is intact;
#' top-N is a visualization choice only.
#' @noRd
build_ora_gene_recurrence <- function(df, top_n = ENRICH_RECURRENCE_TOP_N) {
  if (is.null(df) || is.null(df$geneID) || nrow(df) == 0) return(NULL)
  gl <- .enrich_gene_lists(df)
  genes <- unlist(gl, use.names = FALSE)
  if (!length(genes)) return(NULL)
  padj <- suppressWarnings(as.numeric(df$p.adjust))
  per_gene_padj <- rep(padj, lengths(gl))
  # Per gene, the enriched TERMS containing it (§18 D18.2). Term index repeated
  # per member gene, so we can gather each gene's term set. Deterministic order:
  # by the term's adjusted p ascending, then term label — the same "most
  # significant first" ordering as best_padj, so the column reads predictably.
  term_lab <- .enrich_term_labels(df$Description)
  per_gene_term_ix <- rep(seq_len(nrow(df)), lengths(gl))
  ord_terms <- order(padj, term_lab)
  term_rank <- integer(nrow(df)); term_rank[ord_terms] <- seq_len(nrow(df))
  gene_terms <- tapply(per_gene_term_ix, genes, function(ix) {
    ix <- unique(ix)
    paste(term_lab[ix[order(term_rank[ix])]], collapse = "; ")
  })
  n_terms <- table(genes)
  best <- tapply(per_gene_padj, genes, min, na.rm = TRUE)
  g <- names(n_terms)
  out <- data.frame(gene = g, n_terms = as.integer(n_terms[g]),
                    best_padj = as.numeric(best[g]),
                    terms = as.character(gene_terms[g]),
                    stringsAsFactors = FALSE)
  out <- out[order(-out$n_terms, out$best_padj, out$gene), , drop = FALSE]
  rownames(out) <- NULL
  # top_n = NULL/Inf/NA -> the COMPLETE table (used by the Full Table view,
  # which must never truncate); the chart passes the top-N bound.
  if (is.null(top_n) || !is.finite(top_n)) out else utils::head(out, top_n)
}

#' DT for the complete Shared Genes table (§16.40) — EVERY gene of the current
#' ORA unit, never truncated. Same house style as the other enrichment tables;
#' counts and FDR stay numeric so client-side sorting is correct. Data order is
#' the deterministic recurrence order, which is also the default display sort.
#' @noRd
render_gene_recurrence_dt <- function(rec) {
  # §18 D18.2: the readable "Terms" column lists the enriched terms containing
  # each gene (deterministic order), so no membership information is lost.
  d <- data.frame(Gene = rec$gene, `Enriched terms` = rec$n_terms,
                  `Best FDR` = rec$best_padj,
                  Terms = if (is.null(rec$terms)) "" else rec$terms,
                  check.names = FALSE, stringsAsFactors = FALSE)
  DT::datatable(
    d, rownames = FALSE, selection = "none",
    class = "cell-border stripe hover compact",
    options = list(pageLength = 20, dom = "ftip", scrollX = TRUE,
                   order = list(list(1, "desc"), list(2, "asc")),
                   columnDefs = list(list(targets = 3, className = "dt-left")))) |>
    DT::formatSignif("Best FDR", digits = 3)
}

#' Bar chart for build_ora_gene_recurrence().
#' @noRd
plot_ora_gene_recurrence <- function(rec) {
  rec$.nlp <- .enrich_nlp(rec$best_padj)
  rec$.lab <- factor(rec$gene, levels = rev(rec$gene))
  lim <- .enrich_nlp_limits(rec$.nlp)
  ggplot2::ggplot(rec, ggplot2::aes(x = .data$n_terms, y = .data$.lab,
                                    fill = .data$.nlp)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::scale_fill_viridis_c(name = "-log10(best adj p)", option = "D",
                                  limits = lim, oob = scales::squish) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = "Enriched terms containing the gene", y = NULL) +
    .enrich_plot_theme()
}

# ============================================================================
# §23 C/D — Shared Genes = Term × Gene binary MEMBERSHIP MATRIX (replaces the
# recurrence bar as the displayed plot). Rows = enriched terms, columns = genes,
# cells = 0/1 membership. Complete matrix always preserved (table download);
# display Top-N selects GENES (columns) by the approved ranking (recurrence desc,
# best-FDR asc, gene id) reused from build_ora_gene_recurrence().
# ============================================================================

ENRICH_SG_PX_PER_GENE <- 15L    # per-column legibility budget
ENRICH_SG_PX_PER_TERM <- 20L    # per-row legibility budget

#' Ranked gene IDs for the Shared Genes matrix columns (§23 C2): recurrence desc,
#' best-FDR asc, gene id — exactly build_ora_gene_recurrence()'s order.
#' @noRd
.enrich_sg_gene_order <- function(df) {
  rec <- build_ora_gene_recurrence(df, top_n = NULL)
  if (is.null(rec)) character(0) else as.character(rec$gene)
}
#' Gene IDs kept for the current Top-N view (Inf/NA/NULL or <= N -> all).
#' @noRd
.enrich_sg_keep_genes <- function(df, top_n = ENRICH_BAR_TOP_N) {
  g <- .enrich_sg_gene_order(df)
  n <- suppressWarnings(as.numeric(top_n))
  if (length(n) == 0 || is.na(n) || !is.finite(n) || length(g) <= n) g else g[seq_len(n)]
}
#' §3: gene–gene JACCARD distance over a binary Term × Gene membership matrix `m`
#' (rows = Terms, columns = Genes). Jaccard = |Terms(A) ∩ Terms(B)| /
#' |Terms(A) ∪ Terms(B)| (shared ABSENCES do not count), distance = 1 − Jaccard.
#' Vectorized via `crossprod` (same pattern as build_ora_shared_matrices). Two
#' genes with no memberships at all (empty ∪) are treated as identical (distance 0).
#' @noRd
.enrich_gene_jaccard_dist <- function(m) {
  x <- (m > 0) * 1L                              # Terms x Genes, 0/1
  inter <- crossprod(x)                          # Genes x Genes: |A ∩ B|
  cs    <- colSums(x)                            # |A| per gene
  uni   <- outer(cs, cs, "+") - inter            # |A ∪ B|
  jac   <- ifelse(uni == 0, 1, inter / uni)      # both-empty -> similarity 1
  stats::as.dist(1 - jac)
}
#' §3/§7: DISPLAY column order of the genes in a binary Term × Gene membership
#' matrix `m` by similarity of each gene's Term-membership vector — genes that
#' belong to the SAME/SIMILAR sets of enriched Terms sit next to each other. Jaccard
#' distance (§3) + hierarchical clustering with COMPLETE linkage (§4: the user's
#' stated fallback and the linkage the pipeline's old pheatmap gene clustering used).
#' Deterministic (§10): identical-membership genes form one exact-distance-0 group
#' that is clustered as ONE representative and expanded back IN THE INPUT ORDER (the
#' recurrence/best-FDR/Gene-ID selection order), so ties never reorder between runs.
#' Uses ONLY the Terms present in `m` (the active result view, §7). Changes column
#' order only — never the membership values; never drops a gene. Degenerate cases
#' (< 3 genes, or all genes share one membership pattern) keep the input order.
#' @noRd
.enrich_membership_gene_order <- function(m) {
  genes <- colnames(m)
  if (is.null(genes) || length(genes) < 3L) return(genes)   # 0/1/2 genes: unchanged
  keys <- apply(m, 2, function(col) paste0(as.integer(col > 0), collapse = ""))
  uniq <- unique(keys)
  if (length(uniq) < 2L) return(genes)                      # all identical: keep order
  rep_ix <- vapply(uniq, function(k) which(keys == k)[1], integer(1))
  M <- m[, rep_ix, drop = FALSE]                            # one column per unique pattern
  D <- .enrich_gene_jaccard_dist(M)
  ord <- tryCatch(stats::hclust(D, method = "complete")$order,
                  error = function(e) seq_along(uniq))
  # expand each clustered pattern back to its genes, IN INPUT (selection) ORDER
  unlist(lapply(uniq[ord], function(k) genes[keys == k]), use.names = FALSE)
}
#' Complete Term × Gene binary matrix (rows = terms, cols = genes) — the source of
#' truth for the table download; never truncated. Rows labelled by (unique)
#' Description, columns are gene IDs in the approved ranking order.
#' @noRd
build_ora_gene_term_matrix <- function(df) {
  s <- build_ora_shared_matrices(df)                 # gene_term is genes x terms
  m <- t(s$gene_term)                                # -> terms x genes
  gord <- .enrich_sg_gene_order(df)
  gord <- gord[gord %in% colnames(m)]
  if (length(gord)) m <- m[, gord, drop = FALSE]     # rank the gene columns
  m
}
.enrich_sg_height <- function(n_terms)
  as.integer(max(320, ENRICH_SG_PX_PER_TERM * max(n_terms, 1L) + 150))
.enrich_sg_width <- function(n_genes)
  as.integer(max(420, ENRICH_SG_PX_PER_GENE * max(n_genes, 1L) + 260))

#' Shared Genes Term × Gene binary heatmap (§23 C/D) — plotly. Rows = terms, cols
#' = the Top-N ranked genes. Cell = membership (0/1); hover names the term + gene.
#' Dynamic width/height so the matrix stays as readable as practical.
#' @noRd
plot_ora_shared_genes_matrix <- function(df, top_n = ENRICH_BAR_TOP_N,
                                         source_id = NULL) {
  full <- build_ora_gene_term_matrix(df)
  if (is.null(full) || nrow(full) < 1 || ncol(full) < 1) return(NULL)
  keep <- .enrich_sg_keep_genes(df, top_n)      # SELECTION: recurrence/FDR/id (unchanged)
  keep <- keep[keep %in% colnames(full)]
  if (!length(keep)) return(NULL)
  m <- full[, keep, drop = FALSE]               # Terms x selected Genes (membership)
  m <- m[, .enrich_membership_gene_order(m), drop = FALSE]  # §3: Term-membership (Jaccard) order
  terms <- rownames(m); genes <- colnames(m)
  nT <- length(terms); nG <- length(genes)
  txt <- matrix(sprintf("Term: %s<br>Gene: %s<br>%s",
                        matrix(rep(terms, times = nG), nT, nG),
                        matrix(rep(genes, each = nT), nT, nG),
                        ifelse(m == 1L, "member", "—")), nT, nG)
  ax_x <- list(tickmode = "array", tickvals = genes, ticktext = genes,
               tickangle = -60, tickfont = list(size = 9), automargin = TRUE)
  ax_y <- list(tickmode = "array", tickvals = terms, ticktext = terms,
               tickfont = list(size = 10), automargin = TRUE)
  plotly::plot_ly(
    z = m, x = genes, y = terms, type = "heatmap",
    text = txt, hoverinfo = "text", source = source_id,
    colorscale = list(c(0, "#f3f4f6"), c(1, "#3b6fb5")),
    showscale = FALSE, zmin = 0, zmax = 1,
    xgap = 1, ygap = 1,
    height = .enrich_sg_height(nT), width = .enrich_sg_width(nG)) |>
    plotly::layout(xaxis = ax_x, yaxis = ax_y, margin = list(l = 20, b = 20)) |>
    plotly::config(displaylogo = FALSE,
                   modeBarButtonsToRemove = list("sendDataToCloud"))
}

#' STATIC ggplot Shared Genes matrix for the figure download (§23 C/D) — same
#' rows/cols/selection as the on-screen plotly view.
#' @noRd
plot_ora_shared_genes_matrix_static <- function(df, top_n = ENRICH_BAR_TOP_N) {
  full <- build_ora_gene_term_matrix(df)
  if (is.null(full) || nrow(full) < 1 || ncol(full) < 1) return(NULL)
  keep <- .enrich_sg_keep_genes(df, top_n); keep <- keep[keep %in% colnames(full)]
  if (!length(keep)) return(NULL)
  m <- full[, keep, drop = FALSE]
  m <- m[, .enrich_membership_gene_order(m), drop = FALSE]  # §3: same Jaccard order as on-screen
  long <- data.frame(
    term = factor(rep(rownames(m), times = ncol(m)), levels = rev(rownames(m))),
    gene = factor(rep(colnames(m), each = nrow(m)), levels = colnames(m)),
    member = as.vector(m))
  ggplot2::ggplot(long, ggplot2::aes(.data$gene, .data$term,
                                     fill = factor(.data$member))) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = c("0" = "#f3f4f6", "1" = "#3b6fb5"),
                               guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 7),
                   axis.text.y = ggplot2::element_text(size = 7),
                   panel.grid = ggplot2::element_blank())
}
#' Inches size for the static Shared Genes matrix download.
#' @noRd
.enrich_sg_static_in <- function(n) max(4, min(50, 0.16 * max(n, 1L) + 2))

#' DT of the COMPLETE Term × Gene matrix for the table download (§23 D) — all
#' genes, never restricted by the Top-N view. First column = Term.
#' @noRd
shared_genes_matrix_df <- function(df) {
  m <- build_ora_gene_term_matrix(df)
  if (is.null(m)) return(data.frame())
  cbind(Term = rownames(m), as.data.frame(m, check.names = FALSE))
}

.enrich_group_display <- function(analysis, group) {
  if (is.null(group)) return(NULL)
  if (identical(analysis, "GSEA")) return(.enr_lookup(.ENRICH_RANKING_LABEL, group))
  .enr_lookup(.ENRICH_ORA_GROUP_LABEL, group)
}

.enrich_item_row_label <- function(analysis, group) {
  if (identical(analysis, "GSEA")) return("Contrast:")
  gg <- if (is.null(group)) "" else group
  switch(gg,
    contrasts = "Contrast:", contrasts_wo_direction = "Contrast:",
    partition = "Cluster:",  binary_patterns = "Pattern:", "Item:")
}

#' Breadcrumb (tag) with clickable step-back segments; "(choose X)" for a
#' pending middle level (partial navigation). `any_contrast` shows "Any Contrast"
#' for GSEA and "All DE" for ORA (§D29 §1).
#' @noRd
# ---- Shared presentation helpers (§16.38) ---------------------------------
# Purely cosmetic label/markup builders used by BOTH the ORA and GSEA views, so
# the two branches cannot drift apart visually. Omics-agnostic: everything is
# derived from the resolved unit / manifest fields, never from an omics global.

# Friendly, human-facing database names. Unknown keys pass through unchanged, so
# a future database needs no code change here.
.ENRICH_DB_LABEL <- c(GO_BP = "GO Biological Process",
                      GO_MF = "GO Molecular Function",
                      GO_CC = "GO Cellular Component",
                      KEGG  = "KEGG")

#' Friendly database name for the resolved-unit metadata line.
#' @noRd
.enrich_db_label <- function(db) .enr_lookup(.ENRICH_DB_LABEL, db)

#' Human-facing context for a RESOLVED unit — the "<context>" in
#' "Enrichment for <context>". Group-driven (not omics- or RNA-specific):
#' partition -> "Cluster <x>", binary_patterns -> "Pattern <x>", the terminal
#' ORA all_DE / any_contrast -> "All DE", the terminal GSEA any_contrast ->
#' "Any Contrast" (§D29 §1), anything else (contrasts, GSEA rankings)
#' -> the item verbatim, which is already the human contrast name.
#' @noRd
.enrich_context_label <- function(u) {
  if (is.null(u)) return(NULL)
  it <- if (is.null(u$item)) NA_character_ else as.character(u$item)
  # §D29 §1: the GSEA aggregate unit reads "Any Contrast"; ORA's all_DE keeps "All DE".
  if (identical(u$analysis, "GSEA") && identical(it, "any_contrast")) return("Any Contrast")
  if (identical(u$group, "all_DE") || identical(it, "any_contrast")) return("All DE")
  gl <- if (!is.null(u$group_label) && !is.na(u$group_label) && nzchar(u$group_label))
          as.character(u$group_label) else NA_character_
  # §17 S1: contrasts append the direction subgroup (e.g. "G.vs.Sy — Up + Down").
  if (identical(u$group, "contrasts")) {
    sg <- .enrich_scalar(u$subgroup)
    sg_lbl <- if (!is.na(gl)) gl
              else if (!is.na(sg)) .enr_lookup(.ENRICH_ORA_SUBGROUP_LABEL, sg)
              else NA_character_
    if (!is.na(it) && !is.na(sg_lbl)) return(paste0(it, " — ", sg_lbl))
    return(it)
  }
  if (u$group %in% c("partition", "binary_patterns")) {
    if (!is.na(gl)) return(gl)                                  # "Cluster 1" / "Compare Clusters"
    if (identical(it, "compare"))
      return(if (identical(u$group, "partition")) "Compare Clusters" else "Compare Patterns")
    return(paste(if (identical(u$group, "partition")) "Cluster" else "Pattern", it))
  }
  if (is.na(it)) .enrich_group_display(u$analysis, u$group) else it
}

#' Result noun for the metadata line: GSEA always reports pathways; ORA reports
#' GO *terms* and KEGG *pathways*. Singular when n == 1.
#' @noRd
.enrich_result_noun <- function(analysis, database, n = 2L) {
  pathway_like <- identical(analysis, "GSEA") || identical(database, "KEGG")
  base <- if (pathway_like) "pathway" else "term"
  if (!is.na(n) && n == 1L) base else paste0(base, "s")
}

#' Main section title — the same visual language as the panel's "Enrichment"
#' heading (used for "Enrichment for <context>" and for the selected
#' term/pathway name).
#' @noRd
.enrich_title <- function(text) shiny::h4(text, class = "plot-title full")

#' Compact muted metadata line rendered directly under a .enrich_title(), and
#' CENTERED beneath it (§16.39) — applies to the resolved-unit line and to both
#' the ORA and GSEA selected-item lines, since all three go through here.
#' @noRd
.enrich_meta_line <- function(text) {
  shiny::div(text, class = "text-muted",
             style = "margin:-2px 0 12px; font-size:0.92em; text-align:center;")
}

#' Encoding legend rendered directly BELOW a plot (§16.41). One generic helper
#' used by both the ORA and GSEA enrichment bar plots, so their captions cannot
#' drift apart in typography or alignment; it takes the text as a parameter and
#' contains no analysis- or omics-specific logic, so any future omics on the
#' generic enrichment view gets the identical treatment.
#' @noRd
.enrich_plot_caption <- function(text) {
  shiny::helpText(text, style = paste(
    "text-align:center; margin:8px auto 0; max-width:900px;",
    "font-size:0.88em; line-height:1.35;"))
}

#' §4: the configured enrichment FDR cutoff from the payload (`enr$config$
#' pvalue_cutoff`), or NULL when absent/degenerate. Both ORA and GSEA store ONLY
#' terms significant at this cutoff (ORA: p.adjust < cutoff; GSEA: p.adjust <=
#' cutoff AND pvalue <= cutoff), so the Enrichment Plot — including "All" — never
#' shows non-significant tested terms.
#' @noRd
.enrich_sig_cutoff <- function(enr) {
  if (is.null(enr) || is.null(enr$config)) return(NULL)
  cu <- suppressWarnings(as.numeric(enr$config$pvalue_cutoff))
  if (length(cu) != 1L || is.na(cu) || cu <= 0) NULL else cu
}
#' §4: a short, unobtrusive note above the Enrichment Plot stating that only
#' significant results are shown (so Top 20 / Top 50 / All are not misread as "all
#' tested terms"). Renders the ACTUAL configured FDR cutoff when known, else
#' accurate threshold-independent wording. `≤` (<=) is accurate for both ORA
#' (strictly < cutoff) and GSEA (<= cutoff).
#' @noRd
.enrich_sig_note <- function(cutoff = NULL) {
  txt <- if (is.null(cutoff)) "Showing significant enrichment results only."
         else sprintf("Showing significant enrichment results (FDR ≤ %s).",
                      format(cutoff))
  shiny::div(txt, class = "text-muted",
             style = "font-size:12.5px; margin:2px 0 6px;")
}

#' Secondary content heading — clearly a section, but visually subordinate to
#' .enrich_title() (used for Member Genes / Expression Heatmap / Core Genes /
#' the GSEA heatmap sections).
#' @noRd
.enrich_section_h <- function(text, top = "16px") {
  shiny::tags$h5(text, style = sprintf(
    "margin:%s 0 6px; font-weight:600; font-size:1.02em;", top))
}

#' Bullet-joined metadata pieces ("A • B • C"); NULL/empty pieces are dropped.
#' @noRd
.enrich_meta_join <- function(...) {
  p <- unlist(list(...)); p <- p[!is.na(p) & nzchar(p)]
  paste(p, collapse = " • ")
}

#' "<n> <noun>" with big.mark and singular/plural handling.
#' @noRd
.enrich_count_phrase <- function(n, noun) {
  sprintf("%s %s%s", format(n, big.mark = ","), noun, if (n == 1L) "" else "s")
}

enrichment_breadcrumb_html <- function(state, change_channel = NULL) {
  if (is.null(state$analysis)) return(NULL)
  order <- c("analysis", "database", "group", "item")
  disp <- list(
    analysis = state$analysis,
    database = state$database,
    group    = .enrich_group_display(state$analysis, state$group),
    item     = if (is.null(state$item)) NULL             # §D29 §1: GSEA -> "Any Contrast"
               else if (identical(state$item, "any_contrast"))
                 (if (identical(state$analysis, "GSEA")) "Any Contrast" else "All DE")
               else state$item)
  set_flags <- vapply(order, function(l) !is.null(state[[l]]), logical(1))
  deepest <- if (any(set_flags)) max(which(set_flags)) else 0L
  segs <- list()
  for (i in seq_len(deepest)) {
    l <- order[i]
    if (!is.null(state[[l]])) {
      d <- disp[[l]]
      seg <- if (!is.null(change_channel)) {
        shiny::tags$a(d, href = "#", style = "text-decoration:none;",
          onclick = sprintf("Shiny.setInputValue('%s', %s, {priority:'event'}); return false;",
                            change_channel, jsonlite::toJSON(l, auto_unbox = TRUE)))
      } else d
      segs[[length(segs) + 1]] <- seg
    } else {
      segs[[length(segs) + 1]] <- shiny::tags$span(sprintf("(choose %s)", l),
                                                   style = "color:#a06000;")
    }
  }
  inner <- list()
  for (i in seq_along(segs)) {
    if (i > 1) inner[[length(inner) + 1]] <- shiny::HTML(" &rsaquo; ")
    inner[[length(inner) + 1]] <- segs[[i]]
  }
  shiny::div(style = "font-weight:500; margin-top:4px;",
             do.call(shiny::tagList, inner))
}

#' Render one selector level as an application-style nav bar (Bootstrap-3, to
#' match the flatly theme + bslib navsets used elsewhere). Enabled options are
#' clickable tabs/pills that fire `Shiny.setInputValue` (NOT a bslib navset — so
#' the initial state has NO active tab, preserving the no-auto-selection rule);
#' the currently-selected value gets `.active`; greyed (evaluated-but-empty)
#' options render as BS3 `.disabled` tabs (greyed, no suffix text — the disabled
#' state alone conveys "no significant results", §16.42); hidden options are
#' omitted. `nav_class` = "nav-tabs" (primary) or "nav-pills".
#' @noRd
.enrich_nav_row <- function(row_label, channel, opts, selected = NULL,
                            nav_class = "nav-pills", badges = NULL) {
  if (is.null(opts) || nrow(opts) == 0) return(NULL)
  # §2: optional small "Recommended"-style badge next to a specific option, keyed
  # by its VALUE. Looks like an annotation, not part of the label; reuses the
  # site's green accent (#00A651).
  mk_label <- function(v, lab) {
    bdg <- if (!is.null(badges) && !is.null(badges[[v]])) badges[[v]] else NULL
    if (is.null(bdg)) return(lab)
    shiny::tagList(lab, shiny::tags$span(bdg, class = "badge",
      style = paste("margin-left:6px; font-size:0.72em; font-weight:600;",
                    "color:#fff; background:#00A651; border-radius:10px;",
                    "padding:1px 7px; vertical-align:middle;")))
  }
  items <- lapply(seq_len(nrow(opts)), function(i) {
    v <- opts$value[i]; lab <- opts$label[i]; stt <- opts$status[i]
    if (identical(stt, "greyed")) {
      # Evaluated-but-empty: keep the option VISIBLE and disabled/greyed; the
      # disabled state alone conveys "no significant results" (no suffix text).
      # Failed/not-run options are hidden upstream, so this stays distinct.
      shiny::tags$li(class = "disabled", role = "presentation",
        shiny::tags$a(href = "#", onclick = "return false;", `aria-disabled` = "true",
          shiny::span(mk_label(v, lab))))
    } else {
      is_sel <- !is.null(selected) && identical(as.character(selected), v)
      shiny::tags$li(class = if (is_sel) "active" else NULL, role = "presentation",
        shiny::tags$a(href = "#",
          onclick = sprintf("Shiny.setInputValue('%s', %s, {priority:'event'}); return false;",
                            channel, jsonlite::toJSON(v, auto_unbox = TRUE)),
          mk_label(v, lab)))
    }
  })
  nav <- shiny::tags$ul(class = paste("nav", nav_class), role = "tablist",
                        style = "margin-bottom:0;", do.call(shiny::tagList, items))
  if (is.null(row_label)) return(nav)
  shiny::div(style = "margin:8px 0;",
    shiny::tags$div(row_label, class = "text-muted",
      style = "font-size:0.82em; font-weight:600; text-transform:uppercase; letter-spacing:.03em; margin-bottom:3px;"),
    nav)
}

#' Render the ORA results DT for a resolved ORA unit (Stage 4B).
#'
#' Visible columns (renamed): Term (Description) · ID · Gene ratio (GeneRatio) ·
#' Fold enrichment (Fold_enrichment) · Adj. p-value (p.adjust) · Count. ALL other
#' columns (geneID, Cluster, pvalue, qvalue, RichFactor, zScore, BgRatio, padj,
#' pathway, …) are RETAINED in the widget data but hidden via `columnDefs`, so
#' `input$<id>_rows_selected` maps to the full, original-order df and later
#' master-detail can read geneID etc. No truncation. `selection = "single"`.
#' Default display sort by adjusted p-value ascending (a client-side DataTables
#' order — the widget DATA stays in original order, which is what makes the
#' selected-row -> biological-ID resolution stable across sort/search/paging).
#' @noRd
# §17 S2 (D8): a COMPARE unit (UP+DOWN, or a partition/binary "compare" item) is
# resolved to the WHOLE multi-group table; is_compare identifies it so the result
# table preserves — never collapses — the originating group.
.enrich_is_compare_unit <- function(u) {
  if (is.null(u)) return(FALSE)
  identical(.enrich_scalar(u$subgroup), "up_down") ||
    identical(.enrich_scalar(u$item), "compare") ||
    identical(.enrich_scalar(u$item), "compare_contrasts")   # §5: cross-contrast DE compare
}

#' Per-row display Group labels for a compare table, from `Cluster` + `group_type`
#' (D8): direction up/down -> Up/Down; partition Cluster N -> "Cluster N"; pattern
#' -> "Pattern <x>". `Cluster` stays the internal key; this is UI-only.
#' @noRd
.enrich_row_group_labels <- function(cluster, group_type) {
  cl <- as.character(cluster)
  gt <- .enrich_scalar(group_type)
  if (identical(gt, "direction")) {
    lbl <- unname(.ENRICH_ORA_SUBGROUP_LABEL[cl])          # up->Up, down->Down
    lbl[is.na(lbl)] <- cl[is.na(lbl)]
    return(lbl)
  }
  if (identical(gt, "partition")) return(paste("Cluster", cl))
  if (identical(gt, "pattern"))   return(paste("Pattern", cl))
  cl
}

#' Group labels to show for a compare unit's result df, or NULL when the unit is
#' single-group (no Group column) or has no `Cluster` column.
#' @noRd
.enrich_compare_group_labels <- function(df, u) {
  if (is.null(df) || !.enrich_is_compare_unit(u) || !("Cluster" %in% names(df)))
    return(NULL)
  .enrich_row_group_labels(df$Cluster, u$group_type)
}

#' §12: the group-identifying column name for a compare unit's Download Table,
#' from `group_type`: Up/Down -> "Direction", partition -> "Cluster", pattern ->
#' "Pattern", contrast -> "Contrast" (fallback "Group").
#' @noRd
.enrich_compare_group_colname <- function(group_type) {
  gt <- .enrich_scalar(group_type); if (is.na(gt)) gt <- ""
  switch(gt, direction = "Direction", partition = "Cluster",
         pattern = "Pattern", contrast = "Contrast", "Group")
}
#' §12: COMPLETE compare download table — a leading group column (readable
#' Up/Down · Cluster · Pattern · Contrast label) followed by the per-term stats,
#' every row of the compare df (never Top-N filtered). Friendly column names.
#' @noRd
.enrich_compare_table_df <- function(df, u) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(data.frame())
  grp  <- .enrich_row_group_labels(df$Cluster, u$group_type)
  keep <- intersect(c("Description", "ID", "GeneRatio", "Fold_enrichment",
                      "p.adjust", "Count", "geneID"), names(df))
  out <- data.frame(.group = grp, df[, keep, drop = FALSE],
                    check.names = FALSE, stringsAsFactors = FALSE)
  ren <- c(.group = .enrich_compare_group_colname(u$group_type),
           Description = "Term", p.adjust = "FDR",
           Fold_enrichment = "Fold enrichment", GeneRatio = "Gene ratio",
           Count = "Gene count", geneID = "Genes")
  for (k in names(ren)) if (k %in% names(out)) names(out)[names(out) == k] <- ren[[k]]
  out
}

#' Render the ORA results DT for a resolved ORA unit (Stage 4B / §17 S2).
#' `group_labels` (per-row, length == nrow(df)) adds a leading, visible **Group**
#' column for COMPARE tables (D8) — preserving Up/Down or Cluster/Pattern per term.
#' Rows stay in original order so `input$<id>_rows_selected` still maps to df.
#' @noRd
render_ora_results_dt <- function(df, group_labels = NULL) {
  # §2/§5: the ORA Results Table shows EXACTLY these columns, in this order. This is
  # a DISPLAY mapping only — the underlying data-frame columns (and their values /
  # meanings) are UNCHANGED; `Count` is still the in-cluster-AND-in-term gene count,
  # just relabelled. All other source columns stay in the widget data but hidden, so
  # `*_rows_selected` still maps to the full df and Term-ID selection stays stable.
  visible <- c("Description", "ID", "Count", "in_cluster", "in_term", "in_genome",
               "FoldEnrichment", "pvalue", "p.adjust")
  real_vis <- visible[visible %in% names(df)]
  hidden <- setdiff(names(df), real_vis)
  disp <- df[, c(real_vis, hidden), drop = FALSE]             # visible first; rows UNCHANGED
  add_group <- !is.null(group_labels) && length(group_labels) == nrow(disp)
  # A compare "Group" column (used only by compare tables, which are now dot-matrix
  # only — kept for safety) sits third, after Term and ID. `Cluster` stays hidden.
  vis_cols <- real_vis
  if (add_group) {
    disp$Group <- as.character(group_labels)
    pos <- match("ID", real_vis); if (is.na(pos)) pos <- length(real_vis)
    vis_cols <- append(real_vis, "Group", after = pos)       # Term, ID, Group, …
    disp <- disp[, c(vis_cols, hidden), drop = FALSE]         # reorder Group into place
  }
  pretty <- c(Group = "Group", Description = "Term", ID = "ID",
              Count = "Count (in cluster in term)", in_cluster = "In cluster",
              in_term = "In term", in_genome = "In genome",
              FoldEnrichment = "Fold Enrichment", pvalue = "P-value",
              p.adjust = "FDR")
  disp_names <- unname(vapply(names(disp), function(c) {
    p <- pretty[c]; if (is.na(p)) c else unname(p) }, character(1)))
  n_vis <- length(vis_cols)
  # a compare table sorts by Term initially; a single-group table keeps FDR ascending.
  sort_col <- if (add_group) match("Description", vis_cols) else match("p.adjust", vis_cols)
  order_opt <- if (!is.na(sort_col)) list(list(sort_col - 1L, "asc")) else list()
  col_defs <- list()
  if (n_vis < ncol(disp))
    col_defs <- list(list(targets = seq.int(n_vis, ncol(disp) - 1L),  # 0-based hidden cols
                          visible = FALSE, searchable = FALSE))
  dt <- DT::datatable(
    disp, rownames = FALSE, colnames = disp_names, selection = "single",
    class = "cell-border stripe hover compact",
    options = list(
      pageLength = 20,
      lengthMenu = list(c(10, 20, 50, 100), c("10", "20", "50", "100")),
      dom = "ftip", scrollX = TRUE, order = order_opt, columnDefs = col_defs
    )
  )
  if ("p.adjust" %in% vis_cols)       dt <- DT::formatSignif(dt, "p.adjust", digits = 3)
  if ("pvalue" %in% vis_cols)         dt <- DT::formatSignif(dt, "pvalue", digits = 3)
  if ("FoldEnrichment" %in% vis_cols) dt <- DT::formatRound(dt, "FoldEnrichment", digits = 2)
  dt
}

#' Render the GSEA results DT for a resolved GSEA unit (Stage 4D).
#'
#' §D29 §2: visible columns, in order — Pathway (Description) · ID · NES · Set
#' size (setSize) · Leading edge (leading_edge — the clusterProfiler
#' "tags/list/signal" summary string, kept as-is, NOT a numeric count) · P-value
#' (pvalue) · FDR (p.adjust). This is a DISPLAY mapping only — the underlying
#' data-frame columns and their values/meanings are UNCHANGED. ALL other columns
#' (enrichmentScore, qvalue, rank, contrast, ranking_method, padj, pathway, …)
#' are RETAINED but hidden, so `input$<id>_rows_selected` maps to the full
#' original-order df and 4E can read rank/enrichmentScore etc. NES/setSize/
#' pvalue/p.adjust stay NUMERIC (correct client-side numeric sort). No truncation;
#' `selection = "single"`. Default display sort = p.adjust (FDR) asc, secondary
#' abs(NES) desc (§14.5) via a client-side order over a hidden `.absNES` key —
#' the widget DATA stays in original row order, which is what makes the
#' selected-row -> biological-ID resolution stable across sort/search/paging.
#' @noRd
render_gsea_results_dt <- function(df) {
  visible <- c("Description", "ID", "NES", "setSize", "leading_edge", "pvalue", "p.adjust")
  present_vis <- visible[visible %in% names(df)]
  hidden <- setdiff(names(df), present_vis)
  disp <- df[, c(present_vis, hidden), drop = FALSE]           # visible first; rows UNCHANGED
  disp$.absNES <- if ("NES" %in% names(disp)) abs(disp$NES) else NA_real_  # hidden 2nd-sort key
  pretty <- c(Description = "Pathway", ID = "ID", NES = "NES", setSize = "Set size",
              leading_edge = "Leading edge", pvalue = "P-value", p.adjust = "FDR")
  disp_names <- unname(vapply(names(disp), function(c) {
    p <- pretty[c]; if (is.na(p)) c else unname(p) }, character(1)))
  n_vis <- length(present_vis)
  padj_pos   <- match("p.adjust", present_vis)                 # 0-based among visible
  absnes_pos <- match(".absNES", names(disp))                  # 1-based over all cols
  order_opt <- list()
  if (!is.na(padj_pos))   order_opt <- c(order_opt, list(list(padj_pos - 1L, "asc")))
  if (!is.na(absnes_pos)) order_opt <- c(order_opt, list(list(absnes_pos - 1L, "desc")))
  col_defs <- list()
  if (n_vis < ncol(disp))
    col_defs <- list(list(targets = seq.int(n_vis, ncol(disp) - 1L),  # 0-based hidden cols
                          visible = FALSE, searchable = FALSE))
  dt <- DT::datatable(
    disp, rownames = FALSE, colnames = disp_names, selection = "single",
    class = "cell-border stripe hover compact",
    options = list(
      pageLength = 20,
      lengthMenu = list(c(10, 20, 50, 100), c("10", "20", "50", "100")),
      dom = "ftip", scrollX = TRUE, order = order_opt, columnDefs = col_defs
    )
  )
  if ("p.adjust" %in% present_vis) dt <- DT::formatSignif(dt, "p.adjust", digits = 3)
  if ("pvalue" %in% present_vis)   dt <- DT::formatSignif(dt, "pvalue", digits = 3)
  if ("NES" %in% present_vis)      dt <- DT::formatRound(dt, "NES", digits = 3)
  dt
}

# ==== GSEA master-detail + plot reconstruction (Stage 4E) ====================
# All reconstruction uses the compact Stage-3 payload components verbatim; Shiny
# never recomputes the biological leading edge or membership. Integer indices
# resolve through `gene_index`.

#' The enrichment gene universe (`gene_index`), or NULL.
#' @noRd
enrichment_gene_index <- function(enrichment) {
  gi <- enrichment$gene_index
  if (is.null(gi) || !length(gi)) NULL else as.character(gi)
}

#' Raw integer indices (into gene_index) for one pathway's leading edge, from
#' `gsea_leading_edge[[db]][[ranking]][[contrast]][[pathway]]`. integer(0) if absent.
#' @noRd
.gsea_le_idx <- function(enrichment, u, pathway_id) {
  if (is.null(u) || !identical(u$analysis, "GSEA")) return(integer(0))
  node <- enrichment$gsea_leading_edge[[u$database]][[u$group]][[u$item]][[as.character(pathway_id)]]
  if (is.null(node)) integer(0) else as.integer(node)
}

#' Raw integer indices for one pathway's full membership, from
#' `pathway_membership[[db]][[pathway]]` (per-db; shared across ranking/contrast).
#' @noRd
.gsea_membership_idx <- function(enrichment, u, pathway_id) {
  if (is.null(u) || !identical(u$analysis, "GSEA")) return(integer(0))
  node <- enrichment$pathway_membership[[u$database]][[as.character(pathway_id)]]
  if (is.null(node)) integer(0) else as.integer(node)
}

#' Leading-edge / membership gene IDs for one pathway (decoded via gene_index).
#' @noRd
enrichment_gsea_leading_edge_ids <- function(enrichment, u, pathway_id) {
  gi <- enrichment_gene_index(enrichment); if (is.null(gi)) return(character(0))
  idx <- .gsea_le_idx(enrichment, u, pathway_id)
  idx <- idx[idx >= 1L & idx <= length(gi)]
  unique(gi[idx])
}
#' @noRd
enrichment_gsea_membership_ids <- function(enrichment, u, pathway_id) {
  gi <- enrichment_gene_index(enrichment); if (is.null(gi)) return(character(0))
  idx <- .gsea_membership_idx(enrichment, u, pathway_id)
  idx <- idx[idx >= 1L & idx <= length(gi)]
  unique(gi[idx])
}

#' §23 J: ALL genes annotated to a term, for ANY analysis — decode
#' `pathway_membership[[db]][[term_id]]` (integer indices) via `gene_index`. This
#' is the same compact membership the GSEA All-Gene heatmap uses; ORA reuses it so
#' a term selected through Full Results or GO Simplified resolves to the SAME
#' complete pathway membership. Returns character(0) when the term has no stored
#' membership (older payloads / unsupported db).
#' @noRd
enrichment_term_membership_ids <- function(enrichment, db, term_id) {
  gi <- enrichment_gene_index(enrichment); if (is.null(gi)) return(character(0))
  node <- enrichment$pathway_membership[[db]][[as.character(term_id)]]
  if (is.null(node)) return(character(0))
  idx <- as.integer(node); idx <- idx[idx >= 1L & idx <= length(gi)]
  unique(gi[idx])
}

#' The ranked-gene vector for a resolved GSEA unit: `list(idx, score)` from
#' `gsea_rankings[[ranking]][[contrast]]` (group=ranking, item=contrast; db-
#' independent). NULL if absent/malformed.
#' @noRd
enrichment_gsea_ranking <- function(enrichment, u) {
  if (is.null(u) || !identical(u$analysis, "GSEA")) return(NULL)
  rk <- enrichment$gsea_rankings[[u$group]][[u$item]]
  if (is.null(rk) || is.null(rk$idx) || is.null(rk$score) ||
      length(rk$idx) != length(rk$score) || !length(rk$idx)) return(NULL)
  list(idx = as.integer(rk$idx), score = as.numeric(rk$score))
}

#' The FOLD-CHANGE ranking vector for the SAME contrast as a resolved GSEA unit
#' (`gsea_rankings$fc[[contrast]]`, log2FC; §16.42). Used ONLY to remap the
#' ridgeplot X-values — never the GSEA result, running score or NES. NULL when no
#' FC vector exists for the contrast (e.g. `any_contrast`/"All DE"), which the
#' caller must surface as unavailable rather than falling back to p-value ranks.
#' @noRd
enrichment_fc_ranking <- function(enrichment, u) {
  if (is.null(u) || !identical(u$analysis, "GSEA")) return(NULL)
  rk <- enrichment$gsea_rankings[["fc"]][[u$item]]
  if (is.null(rk) || is.null(rk$idx) || is.null(rk$score) ||
      length(rk$idx) != length(rk$score) || !length(rk$idx)) return(NULL)
  list(idx = as.integer(rk$idx), score = as.numeric(rk$score))
}

#' Reconstruct the weighted (p=1) GSEA running-enrichment-score walk over the
#' ranked list, hitting a pathway's membership indices. Reproduces clusterProfiler's
#' enrichmentScore (peak deviation) and `rank` (peak position). NULL if unusable.
#' @noRd
prepare_gsea_running_score <- function(rank_idx, rank_score, member_idx) {
  if (length(rank_idx) == 0 || length(member_idx) == 0) return(NULL)
  hits <- rank_idx %in% member_idx
  Nh <- sum(hits); N <- length(rank_idx); Nmiss <- N - Nh
  if (Nh == 0 || Nmiss == 0) return(NULL)
  Nr <- sum(abs(rank_score[hits])); if (Nr == 0) return(NULL)
  run <- cumsum(ifelse(hits, abs(rank_score) / Nr, -1 / Nmiss))
  peak <- which.max(abs(run))
  # `score` (the ranked-list metric per position) is carried for the gseaplot2-
  # style metric panel + colour strip; it does not affect the ES walk.
  list(curve = data.frame(position = seq_len(N), running_es = run, hit = hits,
                          score = rank_score),
       peak_position = peak, es = run[peak], n_hits = Nh)
}

#' Ranked-list colour-strip rectangles for the gseaplot2-style Enrichment Plot,
#' replicating `enrichplot::gseaplot2` EXACTLY (§16.43). enrichplot colours the
#' strip by binning `rev(cumsum(position))` — the cumulative pathway-hit count,
#' counted from the bottom of the ranked list — into 9 equal breakpoints, then
#' maps those bins onto a fixed 10-colour blue→red brewer palette. Because
#' `cumsum(hits)` is monotone in rank, the bins run red (top of the list) → light
#' → blue (bottom) left-to-right regardless of the metric distribution; band
#' widths reflect where the pathway's genes concentrate. `pos` = the 0/1 hit
#' indicator over the ranked list (from the reconstructed running score). Returns
#' data.frame(xmin, xmax, col) or NULL. Uses only the hit indicator — fully
#' reproducible from the compact payload; no enrichplot call at runtime.
#' @noRd
.gsea_running_strip <- function(pos) {
  total <- sum(pos)
  if (total < 1) return(NULL)
  v   <- seq(1, total, length.out = 9)
  inv <- findInterval(rev(cumsum(pos)), v)
  if (min(inv) == 0) inv <- inv + 1
  # c(rev(brewer.pal(5,"Blues")), brewer.pal(5,"Reds")) — hard-coded to avoid an
  # RColorBrewer dependency; identical hex values.
  col <- c("#08519C", "#3182BD", "#6BAED6", "#BDD7E7", "#EFF3FF",
           "#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15")
  xmin <- which(!duplicated(inv))
  xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
  data.frame(xmin = xmin, xmax = xmax, col = col[unique(inv)], stringsAsFactors = FALSE)
}

#' Selected-pathway "Enrichment Plot" — a gseaplot2-style reconstruction of the
#' running enrichment score, restyled (§16.42) to resemble the pipeline's
#' `enrichplot::gseaplot2` output WITHOUT enrichplot and WITHOUT changing any
#' calculation (ES/NES/leading edge/ranking are all untouched). Three coordinated
#' panels sharing the rank x-axis, composed with patchwork:
#'   (A) running enrichment-score curve (green) + dashed line at the ES peak;
#'   (B) pathway-gene hit marks (black ticks) over the enrichplot red→blue
#'       ranked-list colour strip (`.gsea_running_strip`, exact algorithm);
#'   (C) the ranked-list metric (grey area) with a zero line.
#' `prep$curve` must carry `score` (the ranked metric per position). Falls back to
#' the single ES-curve panel if patchwork is unavailable.
#' @noRd
plot_gsea_running_score <- function(prep, title = "", nes = NA_real_, padj = NA_real_) {
  d <- prep$curve; peak <- prep$peak_position
  hitpos <- d$position[d$hit]
  # prep$n_hits = pathway genes present in the ranked list (the GSEA set size),
  # NOT the leading-edge/core count — label it accurately (§16.43, item 6).
  sub <- paste0(if (!is.na(nes)) sprintf("NES = %.2f", nes) else "",
                if (!is.na(padj)) sprintf("    FDR = %s", format_scientific_threshold(padj)) else "",
                sprintf("    %d genes in set", prep$n_hits))
  x_scale <- ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.01)))

  # (A) running enrichment score
  pA <- ggplot2::ggplot(d, ggplot2::aes(.data$position, .data$running_es)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = peak, linetype = "dashed",
                        color = "grey55", linewidth = 0.4) +
    ggplot2::geom_line(color = "#00A651", linewidth = 0.9) + x_scale +
    ggplot2::labs(x = NULL, y = "Running\nEnrichment Score", title = title, subtitle = sub) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   panel.grid.minor = ggplot2::element_blank())
  if (!requireNamespace("patchwork", quietly = TRUE)) return(pA)

  # (B) black hit ticks over the enrichplot ranked-list colour strip. The strip is
  # the EXACT enrichplot algorithm (§16.43): red (top of ranked list) → light →
  # blue (bottom), band widths from the cumulative hit distribution. This replaces
  # the previous metric-value gradient, which mapped most genes (|score|≈0 for
  # p-value rankings) to the pale midpoint and lost the red→blue strip.
  strip <- .gsea_running_strip(as.integer(d$hit))
  pB <- ggplot2::ggplot()
  if (!is.null(strip))
    pB <- pB + ggplot2::geom_rect(
      data = strip, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                 ymin = 0, ymax = 0.24, fill = .data$col)) +
      ggplot2::scale_fill_identity()
  # §23 L: taller hit ticks — they now occupy 0.3→1.0 of the (also taller) middle
  # panel (was 0.5→1.0), above a slimmer colour strip. Hit POSITIONS are unchanged.
  pB <- pB +
    ggplot2::geom_segment(data = data.frame(position = hitpos),
                          ggplot2::aes(x = .data$position, xend = .data$position,
                                       y = 0.3, yend = 1),
                          color = "black", linewidth = 0.2) +
    x_scale + ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    ggplot2::labs(x = NULL, y = NULL) + ggplot2::theme_void()

  # (C) ranked-list metric
  pC <- ggplot2::ggplot(d, ggplot2::aes(.data$position, .data$score)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    ggplot2::geom_area(fill = "grey70") + x_scale +
    ggplot2::labs(x = "Rank in Ordered Dataset", y = "Ranked List\nMetric") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  # §23 L: middle hit-tick band enlarged (0.5 -> 1.1), visually comparable to the
  # pipeline gseaplot2 rel_heights c(1.5, 1.1, 1). NES / running-ES / ranking /
  # hit positions / ranked metric unchanged — layout only.
  patchwork::wrap_plots(pA, pB, pC, ncol = 1, heights = c(3, 1.1, 1.4))
}

#' GSEA dotplot (top-N by p.adjust, secondary abs(NES) desc). Static ggplot.
#' @noRd
#' GSEA Enrichment Plot (§16.39) — same bar-plot language as ORA, replacing the
#' old dotplot. X = SIGNED NES (never abs(): positive and negative enrichment
#' must stay distinguishable) with a reference line at 0; fill = display-only
#' -log10(FDR) on its OWN limits (GSEA's FDR range is far narrower than ORA's,
#' so the two plots must not share a scale); annotation "<set size> genes - LE
#' <n>". `le_counts` is supplied by the caller from the compact payload's stored
#' leading edge — the biological leading edge is never recomputed here.
#' @noRd
plot_gsea_enrichment_bar <- function(df, le_counts = NULL,
                                     max_terms = ENRICH_BAR_TOP_N) {
  # (1) SELECT top-N by FDR ascending (tie-break |NES| desc) — selection unchanged.
  d <- df[order(df$p.adjust, -abs(df$NES)), , drop = FALSE]
  d <- utils::head(d, max_terms)
  # (2) DISPLAY order = SIGNED NES descending (§D25.1) — the bar X axis is signed
  # NES, so ordering the SELECTED pathways by signed NES gives the most natural,
  # interpretable arrangement (highest +NES at top → lowest −NES at bottom).
  # Selection stays FDR-based; NOT setSize/gene count. Tie-break FDR asc, Description.
  .nesd <- suppressWarnings(as.numeric(d$NES))
  d <- d[order(-.nesd, suppressWarnings(as.numeric(d$p.adjust)),
               as.character(d$Description)), , drop = FALSE]
  d$.nes <- suppressWarnings(as.numeric(d$NES))
  d$.nlp <- .enrich_nlp(d$p.adjust)
  le <- if (is.null(le_counts)) rep(NA_integer_, nrow(d)) else
          suppressWarnings(as.integer(le_counts[as.character(d$ID)]))
  d$.ann <- ifelse(is.na(le), sprintf("%s genes", d$setSize),
                   sprintf("%s genes · LE %d", d$setSize, le))
  d$.lab <- .enrich_term_labels(d$Description)
  d$.lab <- factor(d$.lab, levels = rev(d$.lab))   # highest signed NES on top
  lim <- .enrich_nlp_limits(d$.nlp)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$.nes, y = .data$.lab,
                                  fill = .data$.nlp)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$.ann,
                                    hjust = ifelse(.data$.nes >= 0, -0.06, 1.06)),
                       size = ENRICH_BAR_ANNOT_SIZE, colour = "grey25") +
    ggplot2::scale_fill_viridis_c(name = "-log10(adj p)", option = "D",
                                  limits = lim, oob = scales::squish) +
    # same headroom rule as ORA (§16.41) — on BOTH sides here, since bars run
    # left and right of zero and each carries an annotation at its far end
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(
      mult = c(ENRICH_BAR_EXPAND, ENRICH_BAR_EXPAND))) +
    ggplot2::labs(x = "NES", y = NULL) +
    .enrich_plot_theme()
}

#' Build the long df for a unit-level GSEA ridgeplot over the top-N pathways
#' (top-N by p.adjust, secondary abs(NES) desc — the pathway-selection logic is
#' unchanged). NO gene-count limit any more (§16.44): every selected pathway is
#' included, so the Shiny ridge covers the SAME pathways as the pipeline (large GO
#' root terms included). Genes are never truncated or sampled.
#'
#' Gene population per pathway = exactly the pipeline's (`plot_gsea_ridgeplot` /
#' legacy `ridgeplot_edited`): "all_genes" uses `pathway_membership ∩ the ranked
#' GSEA universe` = the GSEA geneSet (setSize genes), NOT the broader raw
#' membership; "leading_edge" uses the stored leading edge (already ⊆ the ranked
#' list). X-VALUES are FOLD CHANGE for the same contrast (the legacy FC hack,
#' §16.42/§16.44) — the GSEA result/NES/FDR/membership are untouched. The current
#' *pipeline* writes the RANKING metric (`x@geneList`) to its data_all_genes.csv,
#' so for a p-value ranking the Shiny X-values (FC) match the LEGACY intent, not
#' the current pipeline CSV; for the `fc` ranking both are FC and identical (see
#' §16.44 for the audit). No FC vector (e.g. any_contrast) -> fc_available = FALSE
#' and NO fallback to ranks. Returns list(df, n_paths, fc_available).
#' @noRd
build_gsea_ridge_df <- function(enrichment, u, table_df,
                                mode = c("leading_edge", "all_genes"),
                                max_paths = ENRICH_RIDGE_TOP_N) {
  mode <- match.arg(mode)
  fc  <- enrichment_fc_ranking(enrichment, u)                   # FC (log2FC) -> ridge X-axis
  sel <- enrichment_gsea_ranking(enrichment, u)                 # selected ranking -> ranked universe
  gi  <- enrichment_gene_index(enrichment)
  if (is.null(gi) || is.null(table_df) || nrow(table_df) == 0)
    return(list(df = NULL, n_paths = 0L, fc_available = !is.null(fc)))
  if (is.null(fc))
    return(list(df = NULL, n_paths = 0L, fc_available = FALSE))
  sc <- rep(NA_real_, length(gi)); sc[fc$idx] <- fc$score       # log2FC by gene_index position
  ranked <- if (!is.null(sel)) sel$idx else fc$idx              # the ranked GSEA list (geneSet universe)
  d <- table_df[order(table_df$p.adjust, -abs(table_df$NES)), , drop = FALSE]
  d <- utils::head(d, max_paths)
  rows <- list()
  for (i in seq_len(nrow(d))) {
    pid <- d$ID[i]
    if (mode == "leading_edge") {
      idx <- .gsea_le_idx(enrichment, u, pid)                   # leading edge (already ⊆ ranked)
    } else {
      idx <- .gsea_membership_idx(enrichment, u, pid)           # all pathway genes
      idx <- idx[idx %in% ranked]                               # ∩ ranked list = geneSet (setSize) — NO cap
    }
    if (length(idx) == 0) next
    vals <- sc[idx]; vals <- vals[!is.na(vals)]
    if (length(vals) < 2) next
    rows[[length(rows) + 1L]] <- data.frame(
      pathway = substr(as.character(d$Description[i]), 1, ENRICH_TERM_LABEL_CHARS),
      value = vals, neglp = -log10(d$p.adjust[i]),
      nes = suppressWarnings(as.numeric(d$NES[i])), stringsAsFactors = FALSE)
  }
  if (length(rows) == 0)
    return(list(df = NULL, n_paths = 0L, fc_available = TRUE))
  out <- do.call(rbind, rows)
  # Pipeline ordering: NES ascending — most-negative NES at the BOTTOM of the
  # ridge, most-positive at the top (matches plot_gsea_ridgeplot's factor levels).
  paths <- unique(out$pathway)
  nesv  <- vapply(paths, function(p) out$nes[out$pathway == p][1], numeric(1))
  out$pathway <- factor(out$pathway, levels = paths[order(nesv, decreasing = FALSE)])
  list(df = out, n_paths = nlevels(out$pathway), fc_available = TRUE)
}

#' Heatmap floor message (§16.43): heatmaps are NOT gene-count-limited any more —
#' the only remaining note is the correctness floor (a heatmap needs >= 2 genes to
#' cluster). NULL when >= 2 genes, or when `n` is NULL (no selection). `unaffected`
#' is retained for call-site compatibility (no longer used).
#' @noRd
.enrich_hm_guard_msg <- function(n, unaffected = "gene table") {
  if (is.null(n) || n >= 2) return(NULL)
  shiny::div("At least 2 genes are needed to draw an expression heatmap.",
             class = "text-muted", style = "padding:4px 0;")
}

#' §23 J3: how many of `ids` are present in the expression matrix (env$expr_norm
#' rownames). Used to report how many pathway genes the All-Gene Heatmap can draw.
#' @noRd
.enrich_expr_present_count <- function(env, ids) {
  if (!exists("expr_norm", envir = env, inherits = FALSE)) return(NA_integer_)
  rn <- rownames(get("expr_norm", envir = env, inherits = FALSE))
  if (is.null(rn)) return(NA_integer_)
  sum(unique(as.character(ids)) %in% rn)
}

#' §23 K: dynamic pixel height for a ridgeplot of `n` pathways, so each term has
#' comfortable vertical separation (used by the on-screen render and the download).
#' @noRd
.enrich_ridge_height <- function(n) as.integer(max(360, 62 * max(as.integer(n), 1L) + 120))
#' Distinct pathways in a ridge df (0 when NULL/empty).
#' @noRd
.enrich_ridge_n <- function(ridge_df)
  if (is.null(ridge_df) || is.null(ridge_df$pathway)) 0L else
    length(unique(as.character(ridge_df$pathway)))

#' Real `ggridges` ridgeplot from a build_gsea_ridge_df() long df, styled to match
#' the pipeline `plot_gsea_ridgeplot` (§16.42): default `geom_density_ridges()`
#' geometry, blue→red `-log10(p.adjust)` fill, NES-ascending order (baked into the
#' `pathway` factor levels by the builder), `theme_minimal()`, no y-axis title. X =
#' fold change (log2FC). Optional concise `caption` when a p-value ranking is shown.
#' @noRd
plot_gsea_ridge <- function(ridge_df, xlab = ENRICH_RIDGE_X_LABEL, caption = NULL) {
  # UI refinement (post-S5): enlarge the displayed text so pathway names, axis
  # titles, tick labels and the legend are comfortably readable. Presentation
  # only — the ridge data/geometry/fill are unchanged.
  # §23 K: more vertical separation per pathway — reduce ridge overlap (scale) and
  # add y-axis expansion; canvas height is set dynamically by the caller. Data /
  # FC / selection / ranking unchanged (layout only).
  ggplot2::ggplot(ridge_df,
                  ggplot2::aes(x = .data$value, y = .data$pathway, fill = .data$neglp)) +
    ggridges::geom_density_ridges(scale = 0.9, rel_min_height = 0.01) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.4, 1.2))) +
    ggplot2::scale_fill_gradient(low = "blue", high = "red", name = "-log10(p.adjust)") +
    ggplot2::labs(x = xlab, y = NULL, caption = caption) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = 14, colour = "#2c3e50"),  # pathway labels
      axis.text.x  = ggplot2::element_text(size = 13),
      axis.title.x = ggplot2::element_text(size = 15),
      legend.title = ggplot2::element_text(size = 13),
      legend.text  = ggplot2::element_text(size = 12),
      plot.caption = ggplot2::element_text(size = 12, hjust = 0, colour = "grey45"))
}

#' Enrichment panel markup (called from ui.R). OMICS-AGNOSTIC: `id_prefix`
#' namespaces every input/output id, so each omics instance gets an isolated UI
#' (RNA uses "rna_enrich"; a future omics passes e.g. "prot_enrich"). Required —
#' no default — so a second instance can never silently collide with RNA's ids.
#' @noRd
enrichment_panel_ui <- function(id_prefix) {
  id <- function(s) paste0(id_prefix, "_", s)
  bslib::nav_panel(
    "Enrichment",
    shiny::div(
      style = "padding:14px 16px; max-width:1100px;",
      # Section title in the app's standard section-heading style.
      shiny::h4("Enrichment", class = "plot-title full"),
      shiny::uiOutput(id("availability_ui")),
      shiny::conditionalPanel(
        condition = sprintf("output.%s == 'yes'", id("available")),
        # PRIMARY / HIERARCHICAL SELECTION — all tab-style (§16.36):
        # ORA/GSEA, then Database, Analysis group/Ranking, Cluster/Pattern/Contrast.
        shiny::uiOutput(id("analysis_ui")),
        shiny::div(
          style = "margin-top:12px;",
          shiny::uiOutput(id("database_ui")),
          shiny::uiOutput(id("group_ui")),
          shiny::uiOutput(id("item_ui")),
          shiny::uiOutput(id("subgroup_ui"))   # §17 S1: ORA contrasts only (Up/Down/Up+Down/All)
        ),
        # (The "Selection:" breadcrumb was removed — the hierarchy is clear from
        # the tab controls; selection/reset state logic is unchanged. Only the
        # non-redundant failed-units note (§14.3) is kept, as a standalone notice.)
        shiny::uiOutput(id("notice_ui")),
        # resolved unit: count + results/plots, set off by a rule
        shiny::conditionalPanel(
          condition = sprintf("output.%s == 'yes'", id("resolved")),
          shiny::div(
            style = "margin-top:14px; padding-top:12px; border-top:1px solid #e5e5e5;",
            shiny::uiOutput(id("count_ui")),
            shiny::uiOutput(id("results_tabs_ui"))
          )
        )
      )
    )
  )
}

#' Wire an Enrichment tab (Stage 4). OMICS-AGNOSTIC: works for ANY omics `env`
#' carrying the §14 payload contract — reads `env$enrichment` (+ reuses
#' `env$de_stats`/`env$expr_norm` for detail tables/heatmaps) via the `env`
#' parameter only; never touches an RNA-seq global. `id_prefix` (required)
#' namespaces every input/output id, and `tabs_id` (required) is the parent
#' tabset used for cross-tab `open()` navigation — both are per-instance so
#' multiple omics can be wired side by side. `mge_controller`/`single_selectize_id`
#' are the (optional) Explorer hand-off knobs, wired later per omics (4F). RNA-seq
#' is the current instance (`env=rnaseq_env`, `id_prefix="rna_enrich"`,
#' `tabs_id="rna_tabs"`). Returns the cross-tab controller
#' `list(open, reset, state, max_n)`.
#' @noRd
wire_enrichment <- function(session, input, output, env, id_prefix, tabs_id,
                            mge_controller = NULL, single_selectize_id = NULL,
                            single_tab_label = "Gene Explorer",
                            multi_tab_label  = "Multi-Gene Explorer") {
  id <- function(s) paste0(id_prefix, "_", s)
  st <- shiny::reactiveVal(enrichment_init_state())
  pending_target <- shiny::reactiveVal(NULL)
  # Selected term/pathway ID is held SEPARATELY from `st` so that a row-click
  # (which only changes sel_id) does NOT invalidate resolved_r() / the results
  # navset / the DT — keeping the table + its selection stable while the inline
  # detail updates. Cleared whenever the resolved unit changes (below).
  sel_id <- shiny::reactiveVal(NULL)
  # §17 S3/S4: the selected row's internal group key (Cluster) — disambiguates a
  # compare table where a term ID appears in more than one group. NA/NULL for
  # single-group units. Cleared with sel_id whenever the resolved unit changes.
  sel_group <- shiny::reactiveVal(NULL)

  enr_r      <- shiny::reactive(
    if (enrichment_data_available_env(env))
      tryCatch(get("enrichment", envir = env, inherits = FALSE), error = function(e) NULL)
    else NULL)
  manifest_r <- shiny::reactive({ e <- enr_r(); if (is.null(e)) NULL else e$manifest })
  avail_r    <- shiny::reactive(!is.null(enr_r()))
  resolved_r <- shiny::reactive(enrichment_resolved_unit(manifest_r(), st()))
  # ORA results df for the resolved unit (NULL unless a resolved ORA unit). Kept
  # verbatim; the SAME object feeds the DT and the selection/download, so the
  # selected row -> biological ID resolution is stable.
  ora_table_r <- shiny::reactive({
    u <- resolved_r()
    if (is.null(u) || !identical(u$analysis, "ORA")) return(NULL)
    tryCatch(enrichment_ora_table(enr_r(), u), error = function(e) NULL)
  })
  # GSEA results df for the resolved unit (NULL unless a resolved GSEA unit).
  # Same verbatim/stable-object contract as ora_table_r (4D).
  gsea_table_r <- shiny::reactive({
    u <- resolved_r()
    if (is.null(u) || !identical(u$analysis, "GSEA")) return(NULL)
    tryCatch(enrichment_gsea_table(enr_r(), u), error = function(e) NULL)
  })
  # ---- GSEA master-detail inputs (4E). All key off sel_id → single-pathway,
  # computed only on selection (lazy); never loop over pathways. ---------------
  gsea_row_r <- shiny::reactive({
    df <- gsea_table_r(); id <- sel_id()
    if (is.null(df) || is.null(id)) return(NULL)
    r <- df[as.character(df$ID) == as.character(id), , drop = FALSE]
    if (nrow(r) == 0) NULL else r[1, , drop = FALSE]
  })
  gsea_core_ids_r <- shiny::reactive({
    u <- resolved_r(); id <- sel_id()
    if (is.null(u) || is.null(id)) return(character(0))
    enrichment_gsea_leading_edge_ids(enr_r(), u, id)
  })
  gsea_member_ids_r <- shiny::reactive({
    u <- resolved_r(); id <- sel_id()
    if (is.null(u) || is.null(id)) return(character(0))
    enrichment_gsea_membership_ids(enr_r(), u, id)
  })
  gsea_running_r <- shiny::reactive({
    u <- resolved_r(); id <- sel_id()
    if (is.null(u) || is.null(id)) return(NULL)
    rk <- enrichment_gsea_ranking(enr_r(), u); if (is.null(rk)) return(NULL)
    midx <- .gsea_membership_idx(enr_r(), u, id); if (length(midx) == 0) return(NULL)
    tryCatch(prepare_gsea_running_score(rk$idx, rk$score, midx), error = function(e) NULL)
  })
  # ---- ORA master-detail + unit-plot inputs (4C). Detail keys off sel_id, so
  # it updates on a row-click without rebuilding the table/navset. -----------
  # §17 S3/S4: the ORA selection is keyed by (ID + Group) so a compare row picks
  # the EXACT group's term, not the first ID match (a term can be significant in
  # more than one group). sel_group() carries the selected row's internal Cluster.
  ora_row_r <- shiny::reactive({
    .enrich_ora_row(ora_view_table_r(), sel_id(), sel_group())
  })
  ora_member_ids_r <- shiny::reactive(
    enrichment_ora_member_genes(ora_view_table_r(), sel_id(), sel_group()))
  # §23 J: ALL genes annotated to the selected ORA term (complete pathway
  # membership), decoded from the compact payload; same for Full or GO Simplified.
  ora_all_gene_ids_r <- shiny::reactive({
    u <- resolved_r(); id <- sel_id()
    if (is.null(u) || !identical(u$analysis, "ORA") || is.null(id)) return(character(0))
    enrichment_term_membership_ids(enr_r(), u$database, id)
  })
  # §21: shared DE-stat tables (long builder + contrast scope). Used by the DT
  # render, the Up/Down legend gate (single-contrast only) and the download, so
  # all three stay consistent and the table is built once per selection.
  ora_member_de_r <- shiny::reactive(
    ora_member_de_table(env, ora_member_ids_r(), .enrich_unit_contrast(resolved_r())))
  gsea_core_de_r <- shiny::reactive(
    ora_member_de_table(env, gsea_core_ids_r(), .enrich_unit_contrast(resolved_r())))

  # ---- 4F: Enrichment -> Explorer hand-off (selection-driven) ----------------
  # Reuse the SAME mechanism as Keyword Search (wire_explorer_handoff). The gene
  # hand-off is enabled only when the RNA instance supplies the Explorer knobs
  # (prot/metab omit them -> disabled, exactly as for Keyword Search).
  .handoff_enabled <- !is.null(mge_controller) &&
    !is.null(single_selectize_id) && nzchar(as.character(single_selectize_id)) &&
    !is.null(tabs_id) && nzchar(as.character(tabs_id)) &&
    !is.null(single_tab_label) && !is.null(multi_tab_label)
  # Openable-gene universe: the SAME expr_raw -> expr_norm rownames the Explorers
  # and Keyword Search use (prepare_id_choices_from_matrix). Session-stable.
  enrich_id_universe <- shiny::reactive({
    ids <- NULL
    if (exists("expr_raw", envir = env, inherits = FALSE))
      ids <- prepare_id_choices_from_matrix(get("expr_raw", envir = env, inherits = FALSE))
    if ((is.null(ids) || length(ids) == 0L) &&
        exists("expr_norm", envir = env, inherits = FALSE))
      ids <- prepare_id_choices_from_matrix(get("expr_norm", envir = env, inherits = FALSE))
    ids
  })
  # Stable selected Gene IDs from a detail DE table. `*_rows_selected` indexes the
  # data frame passed to renderDT (the .enrich_de_display_df output, whose first
  # column is `Gene`), NOT the sorted/searched/paged view — so the biological IDs
  # stay correct after any of those. Bounds-checked against the CURRENT table, so a
  # stale index after a Term/Pathway change can never hand off the wrong gene.
  .selected_detail_genes <- function(m, sel_input) {
    d <- .enrich_de_display_df(m)
    if (is.null(d) || is.null(d$Gene) || nrow(d) == 0L) return(character(0))
    sel <- sel_input
    sel <- sel[!is.na(sel) & sel >= 1L & sel <= nrow(d)]
    if (length(sel) == 0L) return(character(0))
    out <- as.character(d$Gene[sel]); out <- out[!is.na(out) & nzchar(out)]
    unique(out)
  }
  ora_selected_genes_r <- shiny::reactive(
    .selected_detail_genes(ora_member_de_r(), input[[id("detail_member_dt_rows_selected")]]))
  gsea_selected_genes_r <- shiny::reactive(
    .selected_detail_genes(gsea_core_de_r(), input[[id("gsea_core_dt_rows_selected")]]))

  # ORA Member Genes -> Explorer hand-off (reuses the Keyword Search mechanism).
  wire_explorer_handoff(
    session, input, output,
    action_id      = id("ora_explorer_action"),
    open_single_id = id("ora_open_single_btn"),
    open_multi_id  = id("ora_open_multi_btn"),
    selected_ids   = ora_selected_genes_r,
    id_universe    = enrich_id_universe,
    has_table      = function() !is.null(ora_row_r()) &&
      { m <- ora_member_de_r(); !is.null(m) && nrow(m) > 0 },
    transfer_ok    = function() length(enrich_id_universe()) > 0,
    enabled        = .handoff_enabled, tabs_id = tabs_id,
    single_tab_label = single_tab_label, single_selectize_id = single_selectize_id,
    multi_tab_label  = multi_tab_label,  multi_controller    = mge_controller)
  # GSEA Core Genes -> Explorer hand-off (Core Genes ONLY; All Genes is heatmap-only).
  wire_explorer_handoff(
    session, input, output,
    action_id      = id("gsea_explorer_action"),
    open_single_id = id("gsea_open_single_btn"),
    open_multi_id  = id("gsea_open_multi_btn"),
    selected_ids   = gsea_selected_genes_r,
    id_universe    = enrich_id_universe,
    has_table      = function() !is.null(gsea_row_r()) &&
      { m <- gsea_core_de_r(); !is.null(m) && nrow(m) > 0 },
    transfer_ok    = function() length(enrich_id_universe()) > 0,
    enabled        = .handoff_enabled, tabs_id = tabs_id,
    single_tab_label = single_tab_label, single_selectize_id = single_selectize_id,
    multi_tab_label  = multi_tab_label,  multi_controller    = mge_controller)
  # Simplified-GO slice for the resolved unit (§16.37: a RESULT table, so it
  # lives under Results, not Plots). NULL unless the manifest says the unit has
  # a simplify node AND the payload actually yields non-empty rows — that pair
  # is what gates the "Simplified results" view, so no empty placeholder is ever
  # shown. KEGG has no simplify node and has_simplify is FALSE there, so it can
  # never reach this. Unit-level and sub-millisecond (§15.6), keyed on
  # resolved_r() only — it does NOT depend on sel_id, so a row-click never
  # rebuilds it, and the payload/simplify computation itself is untouched.
  ora_simplify_r <- shiny::reactive({
    u <- resolved_r()
    if (is.null(u) || !identical(u$analysis, "ORA") || !isTRUE(u$has_simplify))
      return(NULL)
    sm <- tryCatch(enrichment_ora_simplify_table(enr_r(), u), error = function(e) NULL)
    if (is.null(sm) || !is.data.frame(sm) || nrow(sm) == 0) NULL else sm
  })
  # §23 A: the ACTIVE ORA result view — Full Results or GO Simplified — drives the
  # Results table, Enrichment Plot, Shared Genes, Shared Terms and the selected-term
  # detail. Defaults to Full; falls back to Full when no simplified set exists.
  ora_view_table_r <- shiny::reactive({
    full <- ora_table_r()
    if (identical(rv_resultset(), "simplified")) {
      sm <- ora_simplify_r()
      if (!is.null(sm)) return(sm)
    }
    full
  })
  # Shared-gene matrices are unit-level; built from the ACTIVE result view (§23 A).
  ora_shared_r <- shiny::reactive({
    df <- ora_view_table_r()
    if (is.null(df) || is.null(df$geneID)) return(NULL)
    tryCatch(build_ora_shared_matrices(df), error = function(e) NULL)
  })
  # Gene-recurrence table for the Shared Genes chart (§16.39). Unit-level and
  # bounded to the top-N genes; like ora_shared_r it is only pulled when the
  # corresponding Plots sub-tab renders.
  # COMPLETE recurrence table for the unit (every gene, no bound) — the single
  # source for both the Full Table view and the chart, which just takes its
  # head(). One computation per unit; still only pulled when the sub-tab renders.
  ora_recurrence_all_r <- shiny::reactive({
    df <- ora_view_table_r()
    if (is.null(df) || is.null(df$geneID)) return(NULL)
    tryCatch(build_ora_gene_recurrence(df, top_n = NULL), error = function(e) NULL)
  })
  ora_recurrence_r <- shiny::reactive({
    all <- ora_recurrence_all_r()
    if (is.null(all)) NULL else utils::head(all, ENRICH_RECURRENCE_TOP_N)
  })

  # ---- flags driving the conditionalPanels (must compute while hidden) ----
  output[[id("available")]] <- shiny::renderText(if (isTRUE(avail_r())) "yes" else "no")
  shiny::outputOptions(output, id("available"), suspendWhenHidden = FALSE)
  output[[id("resolved")]]  <- shiny::renderText(if (!is.null(resolved_r())) "yes" else "no")
  shiny::outputOptions(output, id("resolved"), suspendWhenHidden = FALSE)

  # ---- availability banner / selectors / breadcrumb / count ----
  output[[id("availability_ui")]] <- shiny::renderUI({
    if (isTRUE(avail_r())) return(NULL)
    shiny::div("Enrichment results are not available for this dataset.",
               style = "color:#888; padding:8px 0;")
  })
  output[[id("analysis_ui")]] <- shiny::renderUI(
    .enrich_nav_row(NULL, id("analysis"),
                    resolve_enrichment_options(manifest_r(), "analysis", st()),
                    st()$analysis, nav_class = "nav-tabs"))
  # Hierarchical selectors are all TAB-STYLE (§16.36) — same nav-tabs language as
  # the ORA/GSEA bar. No auto-selection / greyed / hidden behavior is changed.
  output[[id("database_ui")]] <- shiny::renderUI({
    s <- st(); if (is.null(s$analysis)) return(NULL)
    .enrich_nav_row("Database", id("database"),
                    resolve_enrichment_options(manifest_r(), "database", s),
                    s$database, nav_class = "nav-tabs")
  })
  output[[id("group_ui")]] <- shiny::renderUI({
    s <- st(); if (is.null(s$database)) return(NULL)
    lbl <- if (identical(s$analysis, "GSEA")) "Ranking" else "Analysis group"
    # §2: mark the recommended GSEA ranking ("P-value without direction") with a
    # small Recommended badge; the label and internal id (`pval_wo_direction`) are
    # unchanged. ORA groups get no badge.
    badges <- if (identical(s$analysis, "GSEA"))
      stats::setNames(list("Recommended"), "pval_wo_direction") else NULL
    .enrich_nav_row(lbl, id("group"),
                    resolve_enrichment_options(manifest_r(), "group", s),
                    s$group, nav_class = "nav-tabs", badges = badges)
  })
  output[[id("item_ui")]] <- shiny::renderUI({
    s <- st(); if (is.null(s$group)) return(NULL)
    if (identical(s$analysis, "ORA") && identical(s$group, "all_DE")) return(NULL)  # terminal
    lbl <- sub(":$", "", .enrich_item_row_label(s$analysis, s$group))
    .enrich_nav_row(lbl, id("item"),
                    resolve_enrichment_options(manifest_r(), "item", s),
                    s$item, nav_class = "nav-tabs")
  })
  # §17 S1: the direction subgroup (Up / Down / Up + Down / All) — shown ONLY for
  # an ORA contrast once a contrast (item) is chosen; contrasts are terminal here.
  output[[id("subgroup_ui")]] <- shiny::renderUI({
    s <- st()
    if (!(identical(s$analysis, "ORA") && identical(s$group, "contrasts"))) return(NULL)
    if (is.null(s$item)) return(NULL)
    # §5/§13: Compare Contrasts is a TERMINAL cross-contrast view (not a direction
    # subgroup of one contrast) — it has no Up/Down/DE/Compare-Up-vs-Down selector.
    if (identical(.enrich_scalar(s$item), "compare_contrasts")) return(NULL)
    opts <- resolve_enrichment_options(manifest_r(), "subgroup", s)
    if (nrow(opts) == 0) return(NULL)   # pre-redesign payload: no subgroup level
    .enrich_nav_row("Direction", id("subgroup"), opts, s$subgroup, nav_class = "nav-tabs")
  })
  # Non-redundant failed-units note only (the "Selection:" breadcrumb was removed;
  # enrichment_breadcrumb_html + reset_level state logic remain intact).
  output[[id("notice_ui")]] <- shiny::renderUI({
    s <- st(); if (is.null(s$analysis)) return(NULL)
    nf <- enrichment_scope_failed_count(manifest_r(), s)
    if (nf == 0) return(NULL)
    shiny::div(
      shiny::icon("triangle-exclamation"),
      sprintf(" %d enrichment unit%s could not be evaluated because of a technical error.",
              nf, if (nf == 1) "" else "s"),
      style = "color:#a06000; font-size:0.88em; margin-top:10px;")
  })
  # Main context title + polished metadata line for the resolved unit (§16.38).
  # "Enrichment for <context>" in the same heading style as the panel title,
  # then "<friendly database> • <n> significant <terms|pathways>". Counts come
  # straight from the manifest — no counting logic changed.
  output[[id("count_ui")]] <- shiny::renderUI({
    u <- resolved_r(); if (is.null(u)) return(NULL)
    n <- u$n_significant
    shiny::tagList(
      .enrich_title(sprintf("Enrichment for %s", .enrich_context_label(u))),
      .enrich_meta_line(.enrich_meta_join(
        .enrich_db_label(u$database),
        sprintf("%s significant %s", format(n, big.mark = ","),
                .enrich_result_noun(u$analysis, u$database, n)))))
  })
  # §4/§5: result-set navigation state. `rv_resultset` is the active ORA result set
  # (full | simplified); `.rv_inner` remembers the inner content tab INDEPENDENTLY
  # for each set, so Full Results and GO Simplified never share an inner-tab input.
  # Both are unit-scoped: a change of biological unit (contrast/cluster/pattern/
  # database/analysis unit — all captured by resolved_r()) resets them. The pure
  # helpers `.enrich_inner_tab_*` carry the logic (unit-tested separately).
  rv_resultset <- shiny::reactiveVal("full")
  .rv_inner     <- shiny::reactiveVal(.enrich_inner_tab_init())
  # A click on the Full Results / GO Simplified tab selector sets the active set.
  shiny::observeEvent(input[[id("result_view")]], {
    v <- input[[id("result_view")]]
    if (!is.null(v) && !identical(v, rv_resultset())) rv_resultset(v)
  }, ignoreInit = TRUE)
  # §5: changing the biological unit resets the active set (Full) and both remembered
  # inner tabs (Results) so no state leaks between contrasts/clusters/patterns/dbs.
  shiny::observeEvent(resolved_r(), {
    rv_resultset("full")
    .rv_inner(.enrich_inner_tab_init())
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  # Remember the inner content tab for whichever result set is active as the user
  # navigates Results / Enrichment Plot / Shared Genes / Shared Terms.
  shiny::observeEvent(input[[id("plots_view")]], {
    .rv_inner(.enrich_inner_tab_remember(.rv_inner(), rv_resultset(),
                                         input[[id("plots_view")]]))
  }, ignoreInit = TRUE)
  # §4: switching result set restores THAT set's last inner tab (first visit ->
  # "Results"); the other set keeps its own remembered tab untouched.
  shiny::observeEvent(rv_resultset(), {
    want <- .enrich_inner_tab_restore(.rv_inner(), rv_resultset())
    bslib::nav_select(id("plots_view"), selected = want, session = session)
  }, ignoreInit = TRUE)

  # Results/Plots navset. Rebuilt ONLY when the resolved unit changes (depends
  # on resolved_r, not sel_id) so the table + its selection persist across a
  # row-click. Plot outputs live under the Plots tab and stay suspended (lazy)
  # until that tab — and their sub-tab — is opened.
  # Right-aligned, default-gray download button BELOW the table (app convention),
  # label without "Excel" (still writes .xlsx). Shared by ORA + GSEA.
  dl_row <- shiny::div(
    style = "display:flex; justify-content:flex-end; margin-top:8px;",
    shiny::downloadButton(paste0("download_", id_prefix, "_results"), "Download results"))
  # The flattened content row (§16.38) makes GO Simplified a sibling panel with
  # its own visible table, so it gets its own matching download button — the
  # download always serves the table the user is looking at. Same style/label
  # convention and the same .xlsx writer; the full-results button is unchanged.
  dl_simplify_row <- shiny::div(
    style = "display:flex; justify-content:flex-end; margin-top:8px;",
    shiny::downloadButton(paste0("download_", id_prefix, "_simplify"), "Download results"))

  output[[id("results_tabs_ui")]] <- shiny::renderUI({
    u <- resolved_r(); if (is.null(u)) return(NULL)
    # Consistent right-aligned download button (same style everywhere, item 10).
    dl_btn <- function(sfx, label = "Download figure") shiny::div(
      style = "display:flex; justify-content:flex-end; margin-top:8px;",
      shiny::downloadButton(paste0("download_", id_prefix, "_", sfx), label))
    # §2/§22: every unit type now uses ONE flat tab-style content navset (navset_tab,
    # id = plots_view) whose first tab is "Results" — the same visual language as the
    # hierarchical selectors above. The selected-term detail is a separate output
    # below it, shown only on the "Results" tab.
    if (identical(u$analysis, "GSEA")) {
      # GSEA content navset. Rebuilt only on unit change; the inline master-detail
      # (gsea_detail_ui) is a separate output keyed on sel_id.
      gplot_tabs <- list(
        bslib::nav_panel("Enrichment Plot",
          # §18 D18.1: same Top 20/50/All viewing control as ORA (default 20),
          # dynamic height; download matches the on-screen selection.
          shiny::radioButtons(id("gsea_enrich_topn"), "Show pathways:",
            choices = ENRICH_COMPARE_TOPN_CHOICES, selected = "20", inline = TRUE),
          # §4: the GSEA plot shows ONLY significant pathways (audit: p.adjust <=
          # cutoff AND pvalue <= cutoff), so a note names that so "All" is not read
          # as "all tested pathways".
          .enrich_sig_note(.enrich_sig_cutoff(enr_r())),
          # same treatment as ORA (§16.41): legend BELOW the plot, one line,
          # naming GSEA's own metrics (signed NES / leading edge)
          shiny::div(style = "overflow-x:auto; width:100%;",
            .enrich_spin(shiny::plotOutput(id("gsea_enrich_bar"), width = "100%", height = "auto"))),
          .enrich_plot_caption(ENRICH_GSEA_BAR_CAPTION),
          dl_btn("gsea_bar")),
        bslib::nav_panel("Leading-edge ridgeplot",
          shiny::helpText("Ranking-metric distribution of each top pathway's leading-edge genes."),
          # §18/§22: same full-width plot wrapper + below-plot "Download figure" row
          # as every other enrichment figure (height/spacing unchanged — §23 K keeps
          # the approved ridge spacing).
          shiny::div(style = "overflow-x:auto; width:100%;",
            .enrich_spin(shiny::plotOutput(id("gsea_le_ridge"), width = "100%", height = "480px"))),
          dl_btn("gsea_le_ridge")),
        bslib::nav_panel("All-genes ridgeplot",
          shiny::helpText("Ranking-metric distribution of each top pathway's full membership (large pathways omitted)."),
          shiny::uiOutput(id("gsea_all_ridge_msg")),
          shiny::div(style = "overflow-x:auto; width:100%;",
            .enrich_spin(shiny::plotOutput(id("gsea_all_ridge"), width = "100%", height = "480px"))),
          dl_btn("gsea_all_ridge")))
      # §2/§22: flat content tabs — "Results" is the first tab, with no outer
      # "Full Results | Plots" wrapper (GSEA has no simplification, so Full Results
      # is implicit and the user proceeds straight to the content-tab level). Same
      # single-level navset_tab language as single-group ORA. The master-detail sits
      # BELOW the navset so the selected pathway reads as its own page section.
      gsea_tabs <- c(list(
        bslib::nav_panel("Results Table", DT::DTOutput(id("gsea_results_dt")), dl_row)),
        gplot_tabs)
      gsea_navset <- do.call(bslib::navset_tab,
                             c(list(id = id("plots_view")), gsea_tabs))
      return(shiny::tagList(gsea_navset, shiny::uiOutput(id("gsea_detail_ui"))))
    }
    # §4/§15/§16: a COMPARE unit shows ONLY the comparison dot matrix + Download
    # Figure + Download Table — NO Results Table, NO Grouped Heatmap, NO
    # selected-term detail, no inner content navset. Applies to Compare Up vs. Down
    # / Compare Clusters / Compare Patterns / Compare Contrasts. (Ordinary non-compare
    # units are unchanged — they still have their Results Table + detail + heatmaps.)
    if (.enrich_is_compare_unit(u)) {
      topn_ctrl <- shiny::radioButtons(
        id("compare_topn"), "Show terms:", choices = ENRICH_COMPARE_TOPN_CHOICES,
        selected = "20", inline = TRUE)                        # D11 default = Top 20
      cmp_out <- shiny::div(style = "overflow-x:auto; width:100%;",
        .enrich_spin(plotly::plotlyOutput(id("compare_plot_ly"), width = "auto", height = "auto")))
      cmp_caption <- "Dot size = gene count  ·  colour = -log10(FDR)  ·  one column per group; hover for details"
      # §11 Download Figure (MATCHES the on-screen Top-N) + §12 Download Table (the
      # COMPLETE compare data incl. the group column), both below the dot matrix.
      cmp_dl <- shiny::div(style = "display:flex; gap:8px; justify-content:flex-end; margin-top:8px;",
        shiny::downloadButton(paste0("download_", id_prefix, "_compare_plot"), "Download figure"),
        shiny::downloadButton(paste0("download_", id_prefix, "_compare_table"), "Download table"))
      return(shiny::div(topn_ctrl, cmp_out, .enrich_plot_caption(cmp_caption), cmp_dl))
    }
    # §23 A/C/D/E — single-group ORA content. When GO simplification is available a
    # Result-set selector (Full Results / GO Simplified) drives EVERY view (Results,
    # Enrichment Plot, Shared Genes, Shared Terms, detail). Shared Genes = Term ×
    # Gene matrix (§23 C/D); Shared Terms = Dice heatmap (§23 E). Both: Top 20/50/All
    # + Download Figure + Download Table (complete), NO visible full-table tab.
    enrich_topn_ctrl <- shiny::radioButtons(
      id("enrich_topn"), "Show terms:", choices = ENRICH_COMPARE_TOPN_CHOICES,
      selected = "20", inline = TRUE)                        # default = Top 20
    dl_pair <- function(fig_sfx, tbl_sfx, tbl_label = "Download table")
      shiny::div(style = "display:flex; gap:8px; justify-content:flex-end; margin-top:8px;",
        shiny::downloadButton(paste0("download_", id_prefix, "_", fig_sfx), "Download figure"),
        shiny::downloadButton(paste0("download_", id_prefix, "_", tbl_sfx), tbl_label))
    content_tabs <- list(
      bslib::nav_panel("Results Table", DT::DTOutput(id("results_dt")), dl_row),
      bslib::nav_panel("Enrichment Plot",
        enrich_topn_ctrl,
        # §4: the plot shows ONLY significant terms (audit: ORA table is p.adjust <
        # cutoff), so a note names that so "All" is not read as "all tested terms".
        .enrich_sig_note(.enrich_sig_cutoff(enr_r())),
        shiny::div(style = "overflow-x:auto; width:100%;",
          .enrich_spin(shiny::plotOutput(id("plot_enrich_bar"), width = "100%", height = "auto"))),
        .enrich_plot_caption(ENRICH_ORA_BAR_CAPTION),
        dl_btn("enrich_bar")),
      # §5: Gene-Term Matrix (was "Shared Genes") = Term × Gene binary matrix; Top-N
      # selects GENES; Download Figure (matches view) + Download Table (complete).
      bslib::nav_panel("Gene-Term Matrix",
        shiny::uiOutput(id("sg_controls")),
        shiny::div(style = "overflow-x:auto; width:100%;",
          .enrich_spin(plotly::plotlyOutput(id("shared_genes_matrix"), width = "auto", height = "auto"),
                       heatmap = TRUE)),
        dl_pair("shared_genes_fig", "shared_genes")),
      # §5: Term-Term Matrix (was "Shared Terms") = Dice heatmap; Top-N; Download
      # Figure (matches view) + Download Table (complete Dice matrix).
      bslib::nav_panel("Term-Term Matrix",
        shiny::uiOutput(id("t2t_controls")),
        shiny::div(style = "overflow-x:auto; width:100%;",
          .enrich_spin(plotly::plotlyOutput(id("shared_terms"), width = "auto", height = "auto"),
                       heatmap = TRUE)),
        dl_pair("shared_terms_fig", "shared_terms")))
    content_navset <- do.call(bslib::navset_tab,
                              c(list(id = id("plots_view")), content_tabs))
    # §1/§3: the Full Results / GO Simplified selector is a real tab-style
    # navigation level (same nav-tabs language as every other enrichment selector),
    # rendered by result_view_ui and present ONLY when a non-empty simplified term
    # set exists (so KEGG / non-GO units show no selector). The explanatory
    # GO-Simplified sentence is removed. §4/§5: each result set keeps its OWN inner
    # content tab (see the result-set state observers above).
    rv_ui <- if (is.null(ora_simplify_r())) NULL
             else shiny::uiOutput(id("result_view_ui"))
    shiny::tagList(rv_ui, content_navset, shiny::uiOutput(id("detail_ui")))
  })

  # §1: tab-style Full Results / GO Simplified selector (result-set navigation
  # level) — only when the unit actually has a non-empty simplified set. Uses the
  # shared tab language (.enrich_nav_row / nav-tabs) and the server-controlled
  # active state (rv_resultset) so switching sets never couples their inner tabs.
  output[[id("result_view_ui")]] <- shiny::renderUI({
    if (is.null(ora_simplify_r())) return(NULL)
    opts <- data.frame(
      value  = c("full", "simplified"),
      label  = c("Full Results", "GO Simplified"),
      status = c("ok", "ok"), stringsAsFactors = FALSE)
    shiny::div(style = "margin-bottom:10px;",
      .enrich_nav_row(NULL, id("result_view"), opts,
                      selected = rv_resultset(), nav_class = "nav-tabs"))
  })

  # ---- ORA master-detail (inline, below the table) — lazy on sel_id (4C) ----
  output[[id("detail_ui")]] <- shiny::renderUI({
    # The detail belongs to the RESULT table and shows only on the "Results" content
    # tab (§2/§22: all unit types now use a single flat content navset `plots_view`
    # whose first tab is "Results"). sel_id is untouched, so returning to Results
    # restores the selection.
    pv <- input[[id("plots_view")]]
    if (!is.null(pv) && !identical(pv, "Results Table")) return(NULL)
    row <- ora_row_r(); if (is.null(row)) return(NULL)   # nothing until a term is selected
    n <- length(ora_member_ids_r())
    # Selected-term section (§16.38): the term NAME as a prominent title in the
    # same language as the main heading, one compact "<ID> • <n> genes" line,
    # then subordinate content headings. Shares .enrich_title/.enrich_meta_line/
    # .enrich_section_h with the GSEA detail so both read identically.
    # Same secondary-detail navigation as GSEA (§16.42): a subtle pill navset with
    # Member Genes (default, renders first) and Expression Heatmap (lazy — only
    # when its pill is opened, so the table and heatmap are never shown at once).
    detail_navset <- bslib::navset_pill(
      id = id("ora_detail_view"),
      bslib::nav_panel("Member Genes",
        # §2/§3: no Up/Down colour legend; one shared instruction ABOVE the table.
        if (.handoff_enabled) .explorer_instruction() else NULL,
        DT::DTOutput(id("detail_member_dt")),
        # 4F.1: select gene rows -> open them in the Gene / Multi-Gene Explorer
        # (same selection-driven action as Keyword Search). The dynamic action
        # button sits bottom-right, immediately ABOVE the Download control; it is
        # absent when 0 genes are selected (no hint text).
        shiny::div(style = "display:flex; justify-content:flex-end; margin-top:8px;",
          shiny::uiOutput(id("ora_explorer_action"))),
        shiny::div(style = "display:flex; justify-content:flex-end; margin-top:8px;",
          shiny::downloadButton(paste0("download_", id_prefix, "_member_genes"),
                                "Download results"))),
      # §23 J: DE-Gene Heatmap = the term's DE member genes (drove the ORA term).
      bslib::nav_panel("DE-Gene Heatmap",
        .enrich_section_h("DE-Gene Heatmap", top = "10px"),
        shiny::uiOutput(id("detail_hm_msg")),
        .enrich_spin(plotly::plotlyOutput(id("detail_member_hm"), height = "auto"), heatmap = TRUE),
        shiny::div(style = "display:flex; justify-content:flex-end; margin-top:8px;",
          shiny::downloadButton(paste0("download_", id_prefix, "_member_hm"),
                                "Download data"))),
      # §23 J: All-Gene Heatmap = ALL genes annotated to the term (DE + non-DE)
      # that are present in the expression matrix. Lazy; reuses the same heatmap.
      bslib::nav_panel("All-Gene Heatmap",
        .enrich_section_h("All-Gene Heatmap", top = "10px"),
        shiny::uiOutput(id("detail_all_hm_msg")),
        .enrich_spin(plotly::plotlyOutput(id("detail_all_hm"), height = "auto"), heatmap = TRUE),
        shiny::div(style = "display:flex; justify-content:flex-end; margin-top:8px;",
          shiny::downloadButton(paste0("download_", id_prefix, "_all_hm"),
                                "Download data"))))
    n_all <- length(ora_all_gene_ids_r())
    shiny::div(
      style = "margin-top:18px; padding-top:14px; border-top:1px solid #e5e5e5;",
      .enrich_title(as.character(row$Description)),
      # §23 J: distinguish DE genes from all pathway genes in the metadata line.
      .enrich_meta_line(.enrich_meta_join(
        as.character(row$ID), .enrich_count_phrase(n, "DE gene"),
        if (n_all > 0) .enrich_count_phrase(n_all, "pathway gene") else NULL)),
      detail_navset)
  })
  output[[id("detail_member_dt")]] <- DT::renderDT({
    shiny::req(!is.null(ora_row_r()))
    # pass-3: single-contrast ORA unit -> that contrast only; else all contrasts.
    m <- ora_member_de_r()
    shiny::validate(shiny::need(!is.null(m) && nrow(m) > 0,
                                "No member genes available for this term."))
    # §20/§21: single contrast -> Multi Explorer coloured table (paginated);
    # multiple contrasts -> wide one-row-per-gene table (no colouring).
    # 4F: selectable rows drive the Explorer hand-off (selection is read by the
    # stable Gene column, so it survives sort/search/paging).
    render_enrichment_member_dt(env, m, selection = "multiple")
  }, server = FALSE)
  # §2: the Up/Down colour legend is removed from the Enrichment gene tables
  # (rows are no longer colour-coded).
  # No gene-count limit (§16.43): the heatmap attempts to render ALL member genes;
  # the only note is the >= 2-gene correctness floor. Still lazy (renders only when
  # the Expression Heatmap pill is opened).
  output[[id("detail_hm_msg")]] <- shiny::renderUI({
    shiny::req(!is.null(ora_row_r()))
    .enrich_hm_guard_msg(length(ora_member_ids_r()))
  })
  output[[id("detail_member_hm")]] <- plotly::renderPlotly({
    shiny::req(!is.null(ora_row_r()))
    ids <- ora_member_ids_r()
    shiny::req(length(ids) >= 2)
    prepped <- tryCatch(prepare_selected_features_heatmap_inputs(env, ids),
                        error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(prepped),
      "Not enough of these genes are present in the expression matrix to draw a heatmap."))
    plot_selected_features_heatmap_plotly(prepped)
  })
  # §23 J3: All-Gene Heatmap — ALL genes annotated to the term (DE + non-DE) that
  # are present in the expression matrix. Reports total membership, number mapped,
  # and any omitted for lack of expression. Lazy; reuses the same heatmap builder.
  output[[id("detail_all_hm_msg")]] <- shiny::renderUI({
    shiny::req(!is.null(ora_row_r()))
    all_ids <- ora_all_gene_ids_r()
    if (length(all_ids) == 0)
      return(shiny::div("Complete pathway membership is not available for this term.",
                        class = "text-muted", style = "padding:4px 0;"))
    present <- tryCatch(.enrich_expr_present_count(env, all_ids), error = function(e) NA_integer_)
    # §16: the routine "<n> pathway genes · <n> in the expression matrix" summary
    # sentence is removed; only genuine guard messages (too few genes to draw a
    # heatmap) remain.
    .enrich_hm_guard_msg(present)
  })
  output[[id("detail_all_hm")]] <- plotly::renderPlotly({
    shiny::req(!is.null(ora_row_r()))
    ids <- ora_all_gene_ids_r()
    shiny::validate(shiny::need(length(ids) >= 2,
      "Complete pathway membership is not available for this term."))
    prepped <- tryCatch(prepare_selected_features_heatmap_inputs(env, ids),
                        error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(prepped),
      "Not enough pathway genes are present in the expression matrix to draw a heatmap."))
    plot_selected_features_heatmap_plotly(prepped)
  })

  # ---- ORA unit-level plots (lazy: only render when the Plots sub-tab shows) -
  # §18 D18.1: Top 20/50/All viewing control (default 20); dynamic height so the
  # figure stays readable at All. The download regenerates at the SAME selection.
  output[[id("plot_enrich_bar")]] <- shiny::renderPlot({
    df <- ora_view_table_r(); shiny::req(!is.null(df), nrow(df) >= 1)
    plot_ora_enrichment_bar(df, max_terms = .enrich_topn_val(input[[id("enrich_topn")]]))
  }, height = function() {
    df <- ora_view_table_r()
    if (is.null(df) || nrow(df) < 1) return(400L)
    .enrich_bar_height(.enrich_bar_shown_n(df, input[[id("enrich_topn")]]))
  })
  # §17 S3: COMPARE plots (lazy — only render when the Comparison sub-tab is open).
  # UP+DOWN -> diverging bars (ggplot); Compare Clusters/Patterns -> dot matrix
  # (plotly). Top-N per group by adj p, union across groups (D11).
  output[[id("compare_plot")]] <- shiny::renderPlot({          # UP+DOWN diverging
    u <- resolved_r(); df <- ora_table_r()
    shiny::req(.enrich_is_compare_unit(u), !is.null(df))
    p <- plot_ora_compare_diverging(df, .enrich_topn_val(input[[id("compare_topn")]]))
    shiny::validate(shiny::need(!is.null(p), "No comparison to plot for this unit."))
    p
  }, height = function() {                                     # item 3: dynamic height
    u <- resolved_r(); df <- ora_table_r()
    if (is.null(df) || !.enrich_is_compare_unit(u)) return(440L)
    .enrich_compare_bar_height(.enrich_compare_n_terms(
      df, .enrich_topn_val(input[[id("compare_topn")]])))
  })
  output[[id("compare_plot_ly")]] <- plotly::renderPlotly({    # Compare dot matrix
    u <- resolved_r(); df <- ora_table_r()
    shiny::req(.enrich_is_compare_unit(u), !is.null(df))
    p <- plot_ora_compare_dotmatrix(df, .enrich_topn_val(input[[id("compare_topn")]]),
                                    group_type = u$group_type, source_id = id("compare_dot"))
    shiny::req(!is.null(p)); p
  })
  # §17 S4: GROUPED compare heatmap (lazy). Genes in this comparison's terms,
  # row-grouped by their group (Up/Down or cluster/pattern), from expr_norm.
  output[[id("compare_heatmap")]] <- plotly::renderPlotly({
    u <- resolved_r(); df <- ora_table_r()
    shiny::req(.enrich_is_compare_unit(u), !is.null(df))
    gg <- .enrich_compare_gene_groups(df, u$group_type)
    prepped <- prepare_grouped_heatmap_inputs(env, gg)
    shiny::validate(shiny::need(
      !is.null(prepped), "Not enough genes with expression data for a grouped heatmap."))
    plot_grouped_heatmap_plotly(prepped, source_id = id("compare_hm"))
  })
  # §23 C/D: Shared Genes = Term × Gene binary MATRIX. Top-N selects GENES
  # (columns) by recurrence desc / best-FDR asc / gene id. Dynamic size + scroll.
  output[[id("sg_controls")]] <- shiny::renderUI({
    df <- ora_view_table_r()
    ng <- if (is.null(df)) 0L else length(.enrich_sg_gene_order(df))
    if (ng <= ENRICH_BAR_TOP_N) return(NULL)
    shiny::div(
      style = "display:flex; align-items:center; gap:10px; margin-bottom:6px;",
      shiny::radioButtons(id("sg_topn"), "Show genes:",
        choices = c("Top 20" = "20", "Top 50" = "50", "All" = "all"),
        selected = "20", inline = TRUE),
      shiny::span(sprintf("%s genes total; top-N by recurrence then best FDR.",
                          format(ng, big.mark = ",")),
                  class = "text-muted", style = "font-size:0.9em;"))
  })
  .sg_topn <- function() {
    v <- input[[id("sg_topn")]]
    if (is.null(v)) return(ENRICH_BAR_TOP_N)
    if (identical(as.character(v), "all")) return(Inf)
    n <- suppressWarnings(as.numeric(v)); if (is.na(n)) ENRICH_BAR_TOP_N else n
  }
  output[[id("shared_genes_matrix")]] <- plotly::renderPlotly({
    df <- ora_view_table_r()
    shiny::validate(shiny::need(!is.null(df) && !is.null(df$geneID) && nrow(df) >= 1,
                                "No overlap genes available for this unit."))
    p <- plot_ora_shared_genes_matrix(df, top_n = .sg_topn(),
                                      source_id = id("sg_matrix"))
    shiny::validate(shiny::need(!is.null(p),
      "Not enough gene membership to draw a Term × Gene matrix."))
    p
  })
  # §18 D18.3: Top 20/50/All viewing control (default 20, significance-based) —
  # rendered only inside the Shared Terms Plot sub-tab, so it stays lazy; resets
  # to the default whenever the unit changes. Only shown when the unit has more
  # than the default number of terms (otherwise all terms already fit).
  output[[id("t2t_controls")]] <- shiny::renderUI({
    s <- ora_shared_r(); if (is.null(s) || s$n_terms <= ENRICH_T2T_TOP_N) return(NULL)
    shiny::div(
      style = "display:flex; align-items:center; gap:10px; margin-bottom:6px;",
      shiny::radioButtons(id("t2t_topn"), "Show terms:",
        choices = c("Top 20" = "20", "Top 50" = "50", "All" = "all"),
        selected = "20", inline = TRUE),
      shiny::span(sprintf("%s terms total; top-N selected by adjusted p-value.",
                          format(s$n_terms, big.mark = ",")),
                  class = "text-muted", style = "font-size:0.9em;"))
  })
  # Map the radio to a numeric cap (default 20 before the control reports; Inf =
  # all). Reused by the render AND the figure download so they stay in lockstep.
  .t2t_topn <- function() {
    v <- input[[id("t2t_topn")]]
    if (is.null(v)) return(ENRICH_T2T_TOP_N)
    if (identical(as.character(v), "all")) return(Inf)
    n <- suppressWarnings(as.numeric(v))
    if (is.na(n)) ENRICH_T2T_TOP_N else n
  }
  output[[id("shared_terms")]] <- plotly::renderPlotly({
    s <- ora_shared_r(); shiny::req(!is.null(s))
    shiny::validate(shiny::need(s$n_terms >= 2,
      "At least 2 terms are needed for the shared-terms view."))
    plot_ora_term_overlap(s, top_n = .t2t_topn())
  })
  # §18 D18.4: Full Table view — the COMPLETE Dice matrix (all terms), never
  # restricted by the Top-N view. Scrollable/searchable; download serves the same.
  output[[id("shared_terms_dt")]] <- DT::renderDT({
    s <- ora_shared_r(); shiny::req(!is.null(s))
    shiny::validate(shiny::need(s$n_terms >= 2,
      "At least 2 terms are needed for the shared-terms view."))
    render_term_overlap_dt(s)
  }, server = FALSE)
  # ---- Simplified GO results table (§16.37: under Results, not Plots) -------
  # Same renderer as the main ORA table, so columns/sort/search/paging and the
  # hidden-but-retained columns behave identically. Only rendered when the
  # "Simplified results" sub-view exists (ora_simplify_r() non-NULL).
  output[[id("simplify_dt")]] <- DT::renderDT({
    sm <- ora_simplify_r(); shiny::req(!is.null(sm))
    render_ora_results_dt(sm, .enrich_compare_group_labels(sm, resolved_r()))
  }, server = FALSE)

  # Stable biological-ID selection from the simplified table (§16.14), identical
  # contract to the main ORA table: DT reports the index into the df handed to
  # renderDT, so sm$ID[sel] is the biological ID — correct after sort/search/
  # paging. Simplified GO terms are a SUBSET of the unit's ORA rows, so the
  # shared master-detail (ora_row_r/ora_member_ids_r, which resolve by ID
  # against ora_table_r()) renders the same term detail either table selects it
  # from. An ID somehow absent from the full table simply yields no detail.
  shiny::observeEvent(input[[id("simplify_dt_rows_selected")]], {
    sel <- input[[id("simplify_dt_rows_selected")]]
    sm  <- ora_simplify_r()
    ok  <- !is.null(sm) && length(sel) == 1L && sel >= 1L && sel <= nrow(sm)
    sel_id(if (ok) as.character(sm$ID[sel]) else NULL)
    sel_group(if (ok && "Cluster" %in% names(sm)) as.character(sm$Cluster[sel]) else NULL)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  # ---- ORA results table + stable-ID row selection + Excel download (4B) ----
  # §23 A: renders the ACTIVE result view (Full Results or GO Simplified).
  output[[id("results_dt")]] <- DT::renderDT({
    df <- ora_view_table_r(); shiny::req(!is.null(df))
    render_ora_results_dt(df, .enrich_compare_group_labels(df, resolved_r()))
  }, server = FALSE)

  # DT reports the selected index into the ORIGINAL df passed to renderDT
  # (sort/search/page-aware), so df$ID[sel] is the stable biological ID. It is
  # stored in sel_id (NOT st) so it does not invalidate the table/navset.
  shiny::observeEvent(input[[id("results_dt_rows_selected")]], {
    sel <- input[[id("results_dt_rows_selected")]]
    df  <- ora_view_table_r()
    ok  <- !is.null(df) && length(sel) == 1L && sel >= 1L && sel <= nrow(df)
    sel_id(if (ok) as.character(df$ID[sel]) else NULL)
    # §17 S3/S4: capture the selected row's group (Cluster) for the composite key.
    sel_group(if (ok && "Cluster" %in% names(df)) as.character(df$Cluster[sel]) else NULL)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  # ---- GSEA results table + stable-ID row selection (4D) ----
  output[[id("gsea_results_dt")]] <- DT::renderDT({
    df <- gsea_table_r(); shiny::req(!is.null(df))
    render_gsea_results_dt(df)
  }, server = FALSE)

  # Same stable-ID contract as ORA: DT index -> gsea_table_r()$ID[sel]. Only one
  # of the two tables is rendered at a time, so only its rows_selected fires.
  shiny::observeEvent(input[[id("gsea_results_dt_rows_selected")]], {
    sel <- input[[id("gsea_results_dt_rows_selected")]]
    df  <- gsea_table_r()
    sel_id(if (!is.null(df) && length(sel) == 1L &&
               sel >= 1L && sel <= nrow(df)) as.character(df$ID[sel]) else NULL)
    sel_group(NULL)                                   # GSEA has no group dimension
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  # Any change to the resolved unit (a selector change) clears the term
  # selection so no stale detail is shown. A row-click does NOT change st /
  # resolved_r, so this never fires on selection itself.
  shiny::observeEvent(resolved_r(), { sel_id(NULL); sel_group(NULL) },
                      ignoreNULL = FALSE, ignoreInit = TRUE)

  # ---- GSEA master-detail (inline, below the table) — lazy on sel_id (4E) ----
  output[[id("gsea_detail_ui")]] <- shiny::renderUI({
    # Same rule as ORA (§2/§22): no selected-pathway detail off the "Results" tab —
    # the running score, core-gene table and heatmaps belong to the plot views. The
    # GSEA content is now a single flat navset `plots_view` (first tab "Results").
    # sel_id is preserved, so returning to Results restores the detail.
    pv <- input[[id("plots_view")]]
    if (!is.null(pv) && !identical(pv, "Results Table")) return(NULL)
    row <- gsea_row_r(); if (is.null(row)) return(NULL)   # nothing until a pathway is selected
    nc <- length(gsea_core_ids_r()); nm <- length(gsea_member_ids_r())
    # Detail subviews stay a LAZY nested navset: stacking all four would render
    # the running-score plot, the core-gene table AND both heatmaps at once for
    # every selected pathway, which §16.23 forbids. Each panel carries the
    # §16.38 secondary content heading, and both heatmaps are headed immediately
    # before the plot as required. Subtle pill-style (content switching).
    # Order (§16.42): Core Genes (default) → Core-Gene Heatmap → All-Gene Heatmap
    # → Enrichment Plot. Core Genes renders first; heatmaps and the running-score
    # Enrichment Plot each render only when their pill is opened (lazy, §16.23).
    gdl <- function(sfx, label) shiny::div(
      style = "display:flex; justify-content:flex-end; margin-top:8px;",
      shiny::downloadButton(paste0("download_", id_prefix, "_", sfx), label))
    detail_navset <- bslib::navset_pill(
      id = id("detail_view"),
      bslib::nav_panel("Core Genes",
        # §2/§3: no Up/Down colour legend; one shared instruction ABOVE the table.
        if (.handoff_enabled) .explorer_instruction() else NULL,
        DT::DTOutput(id("gsea_core_dt")),
        # 4F.1: select Core-Gene rows -> Gene / Multi-Gene Explorer (same as KWS).
        # GSEA hands off ONLY from Core Genes; All Genes stays heatmap-only. The
        # dynamic action button sits bottom-right, immediately ABOVE the Download
        # control; absent when 0 genes are selected (no hint text).
        shiny::div(style = "display:flex; justify-content:flex-end; margin-top:8px;",
          shiny::uiOutput(id("gsea_explorer_action"))),
        gdl("gsea_core_genes", "Download results")),
      bslib::nav_panel("Core-Gene Heatmap",
        .enrich_section_h("Core-Gene Heatmap", top = "10px"),
        shiny::uiOutput(id("gsea_core_hm_msg")),
        .enrich_spin(plotly::plotlyOutput(id("gsea_core_hm"), height = "auto"), heatmap = TRUE),
        gdl("gsea_core_hm", "Download data")),
      bslib::nav_panel("All-Gene Heatmap",
        .enrich_section_h("All-Gene Heatmap", top = "10px"),
        shiny::uiOutput(id("gsea_all_hm_msg")),
        .enrich_spin(plotly::plotlyOutput(id("gsea_all_hm"), height = "auto"), heatmap = TRUE),
        gdl("gsea_all_hm", "Download data")),
      bslib::nav_panel("Enrichment Plot",
        .enrich_spin(shiny::plotOutput(id("gsea_running_plot"), height = "460px")),
        gdl("gsea_running", "Download figure")))
    # Same title/metadata language as the ORA selected-term section (§16.38):
    # pathway NAME as the prominent title, then "<ID> • <n> pathway genes •
    # <n> core genes". Counts are the existing membership/leading-edge lengths.
    shiny::div(
      style = "margin-top:18px; padding-top:14px; border-top:1px solid #e5e5e5;",
      .enrich_title(as.character(row$Description)),
      .enrich_meta_line(.enrich_meta_join(
        as.character(row$ID),
        .enrich_count_phrase(nm, "pathway gene"),
        .enrich_count_phrase(nc, "core gene"))),
      detail_navset)
  })

  # Running-score / gseaplot2-style plot (single selected pathway).
  output[[id("gsea_running_plot")]] <- shiny::renderPlot({
    row <- gsea_row_r(); shiny::req(!is.null(row))
    prep <- gsea_running_r()
    shiny::validate(shiny::need(!is.null(prep),
      "A running-score plot is not available for this pathway (ranking or membership missing)."))
    plot_gsea_running_score(prep, title = as.character(row$Description),
                            nes = suppressWarnings(as.numeric(row$NES)),
                            padj = suppressWarnings(as.numeric(row$p.adjust)))
  })

  # Core (leading-edge) gene table — reuses the shared de_stats-join renderer.
  # pass-3 item 1: a specific-contrast GSEA unit shows ONLY that contrast's DE
  # stats; any_contrast keeps all contrasts (NA -> all).
  output[[id("gsea_core_dt")]] <- DT::renderDT({
    shiny::req(!is.null(gsea_row_r()))
    m <- gsea_core_de_r()
    shiny::validate(shiny::need(!is.null(m) && nrow(m) > 0,
                                "No leading-edge genes available for this pathway."))
    # §20/§21: single contrast -> Multi Explorer coloured table (paginated);
    # multiple contrasts (any_contrast) -> wide one-row-per-gene table (no colour).
    # 4F: selectable Core-Gene rows drive the Explorer hand-off (stable Gene col).
    render_enrichment_member_dt(env, m, selection = "multiple")
  }, server = FALSE)
  # §2: the Up/Down colour legend is removed from the Enrichment gene tables.

  # Core-gene heatmap (leading-edge genes) — no gene-count limit (§16.43), lazy.
  output[[id("gsea_core_hm_msg")]] <- shiny::renderUI(
    .enrich_hm_guard_msg(if (is.null(gsea_row_r())) NULL else length(gsea_core_ids_r())))
  output[[id("gsea_core_hm")]] <- plotly::renderPlotly({
    shiny::req(!is.null(gsea_row_r()))
    ids <- gsea_core_ids_r()
    shiny::req(length(ids) >= 2)
    prepped <- tryCatch(prepare_selected_features_heatmap_inputs(env, ids), error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(prepped),
      "Not enough of these genes are present in the expression matrix to draw a heatmap."))
    plot_selected_features_heatmap_plotly(prepped)
  })

  # All-gene heatmap (full membership) — no gene-count limit (§16.43): renders ALL
  # pathway genes (large GO terms included), lazily.
  output[[id("gsea_all_hm_msg")]] <- shiny::renderUI(
    .enrich_hm_guard_msg(if (is.null(gsea_row_r())) NULL else length(gsea_member_ids_r())))
  output[[id("gsea_all_hm")]] <- plotly::renderPlotly({
    shiny::req(!is.null(gsea_row_r()))
    ids <- gsea_member_ids_r()
    shiny::req(length(ids) >= 2)
    prepped <- tryCatch(prepare_selected_features_heatmap_inputs(env, ids), error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(prepped),
      "Not enough of these genes are present in the expression matrix to draw a heatmap."))
    plot_selected_features_heatmap_plotly(prepped)
  })

  # ---- GSEA unit-level plots (lazy: only render when the Plots sub-tab shows) -
  # §18 D18.1: Top 20/50/All viewing control (default 20); dynamic height. Leading-
  # edge counts for the bar annotation come from the stored compact payload
  # (decode only), and ONLY for the pathways actually drawn (the current
  # selection) — no recomputation of the biological leading edge. The download
  # regenerates at the SAME selection.
  gsea_bar_le <- function(u, df, n) {
    top <- utils::head(df[order(df$p.adjust, -abs(df$NES)), , drop = FALSE], n)
    vapply(as.character(top$ID), function(pid)
      length(enrichment_gsea_leading_edge_ids(enr_r(), u, pid)), integer(1))
  }
  output[[id("gsea_enrich_bar")]] <- shiny::renderPlot({
    u <- resolved_r(); df <- gsea_table_r()
    shiny::req(!is.null(u), !is.null(df), nrow(df) >= 1)
    n <- .enrich_bar_shown_n(df, input[[id("gsea_enrich_topn")]])
    plot_gsea_enrichment_bar(df, le_counts = gsea_bar_le(u, df, n), max_terms = n)
  }, height = function() {
    df <- gsea_table_r()
    if (is.null(df) || nrow(df) < 1) return(400L)
    .enrich_bar_height(.enrich_bar_shown_n(df, input[[id("gsea_enrich_topn")]]))
  })
  # The ridge X-axis is FOLD CHANGE for the contrast; caption clarifies this when
  # a p-value ranking is selected (the ranking still governs the GSEA result).
  .ridge_caption <- function(u) if (identical(u$group, "fc")) NULL else .ENRICH_RIDGE_FC_CAPTION
  # Shared "no FC vector" message (e.g. any_contrast / "Any Contrast") — never falls back.
  .ridge_no_fc <- paste("Fold-change ridge distributions are unavailable for this unit:",
    "no per-gene fold-change vector exists for \"Any Contrast\" / any_contrast.",
    "Select a specific contrast to view fold-change ridges.")
  # §23 K: dynamic height so each pathway has comfortable vertical separation.
  .le_ridge_r <- shiny::reactive({
    u <- resolved_r(); df <- gsea_table_r()
    if (is.null(u) || is.null(df)) return(NULL)
    build_gsea_ridge_df(enr_r(), u, df, mode = "leading_edge")
  })
  output[[id("gsea_le_ridge")]] <- shiny::renderPlot({
    u <- resolved_r(); shiny::req(!is.null(u))
    rd <- .le_ridge_r()
    shiny::validate(shiny::need(!is.null(rd) && isTRUE(rd$fc_available), .ridge_no_fc))
    shiny::validate(shiny::need(!is.null(rd$df),
      "No leading-edge fold-change values available to plot for this unit."))
    plot_gsea_ridge(rd$df, caption = .ridge_caption(u))
  }, height = function() .enrich_ridge_height(.enrich_ridge_n(
      tryCatch(.le_ridge_r()$df, error = function(e) NULL))))
  gsea_all_ridge_r <- shiny::reactive({
    u <- resolved_r(); df <- gsea_table_r()
    if (is.null(u) || is.null(df)) return(NULL)
    build_gsea_ridge_df(enr_r(), u, df, mode = "all_genes")
  })
  # No pathways are excluded any more (§16.44) — the note output stays as a slot
  # for layout stability but renders nothing.
  output[[id("gsea_all_ridge_msg")]] <- shiny::renderUI(NULL)
  output[[id("gsea_all_ridge")]] <- shiny::renderPlot({
    rd <- gsea_all_ridge_r(); u <- resolved_r()
    shiny::validate(shiny::need(!is.null(rd) && isTRUE(rd$fc_available), .ridge_no_fc))
    shiny::validate(shiny::need(!is.null(rd$df),
      "No all-genes fold-change values available to plot for this unit."))
    plot_gsea_ridge(rd$df, caption = .ridge_caption(u))
  }, height = function() .enrich_ridge_height(.enrich_ridge_n(
      tryCatch(gsea_all_ridge_r()$df, error = function(e) NULL))))

  # One download handler for whichever analysis is resolved (full table, no
  # truncation). ORA -> ora_table_r; GSEA -> gsea_table_r.
  output[[paste0("download_", id_prefix, "_results")]] <- shiny::downloadHandler(
    filename = function() {   # §10/§22: no date (consistent with all enrichment downloads)
      u <- resolved_r()
      sprintf("enrichment_%s_%s_%s.xlsx",
              if (is.null(u)) "unit" else u$analysis,
              if (is.null(u)) "unit" else u$database,
              if (is.null(u)) "unit" else paste0(u$group, "_", u$item))
    },
    content = function(file) {
      u  <- resolved_r()
      df <- if (!is.null(u) && identical(u$analysis, "GSEA")) gsea_table_r() else ora_table_r()
      writexl::write_xlsx(if (is.null(df)) data.frame() else df, path = file)
    }
  )

  # Matching handler for the GO Simplified panel's own button (§16.38). Same
  # writer/format as the full-results download, serving the simplify slice.
  output[[paste0("download_", id_prefix, "_simplify")]] <- shiny::downloadHandler(
    filename = function() {   # §10/§22: no date
      u <- resolved_r()
      sprintf("enrichment_%s_%s_%s_simplified.xlsx",
              if (is.null(u)) "unit" else u$analysis,
              if (is.null(u)) "unit" else u$database,
              if (is.null(u)) "unit" else paste0(u$group, "_", u$item))
    },
    content = function(file) {
      sm <- ora_simplify_r()
      writexl::write_xlsx(if (is.null(sm)) data.frame() else sm, path = file)
    }
  )

  # Full Shared Genes table download (§16.40) — same style/format as the other
  # enrichment downloads, serving every gene of the unit.
  output[[paste0("download_", id_prefix, "_shared_genes")]] <- shiny::downloadHandler(
    filename = function() {   # §10/§22: no date
      u <- resolved_r()
      sprintf("enrichment_%s_%s_%s_gene_term_matrix.xlsx",
              if (is.null(u)) "unit" else u$analysis,
              if (is.null(u)) "unit" else u$database,
              if (is.null(u)) "unit" else paste0(u$group, "_", u$item))
    },
    content = function(file) {
      # §23 D: the Download Table is the COMPLETE Term × Gene binary matrix of the
      # ACTIVE result view — every gene, never restricted by the Top-N plot view.
      df <- ora_view_table_r()
      out <- if (is.null(df)) data.frame() else shared_genes_matrix_df(df)
      writexl::write_xlsx(out, path = file)
    }
  )

  # ---- post-S5 download coverage --------------------------------------------
  # §9/§10: filenames carry the unit context and (for Top-N figures) the selected
  # variant (top20/top50/all) so Top 20/50/All downloads never overwrite each other
  # — and NO download date. Delegates to the pure .enrich_dl_* helpers (tested).
  .dl_name <- function(kind, ext, variant = NULL) {
    u <- resolved_r()
    ctx <- if (is.null(u)) "unit" else paste0(u$database, "_", u$group, "_", u$item,
      if (!is.na(.enrich_scalar(u$subgroup))) paste0("_", u$subgroup) else "")
    .enrich_dl_filename(kind, ctx, ext, variant)
  }
  .dl_topn_tag <- .enrich_dl_topn_tag
  .ggsave_pdf <- function(file, p, w = 11, h = 8) {
    if (is.null(p)) { grDevices::pdf(file, width = 6, height = 4); grDevices::dev.off() }
    else ggplot2::ggsave(file, p, device = "pdf", width = w, height = h, limitsize = FALSE)
  }
  .hm_xlsx <- function(file, ids) {
    pr <- if (length(ids) >= 2)
            tryCatch(prepare_selected_features_heatmap_inputs(env, ids), error = function(e) NULL) else NULL
    out <- if (is.null(pr)) data.frame() else cbind(Gene = pr$rn, as.data.frame(pr$mat, check.names = FALSE))
    writexl::write_xlsx(out, path = file)
  }
  # Item 11: the compare figure ALWAYS contains ALL significant terms as a full,
  # readable static figure — independent of the on-screen Top-N (viz-only).
  output[[paste0("download_", id_prefix, "_compare_plot")]] <- shiny::downloadHandler(
    filename = function() .dl_name("compare", "pdf", .dl_topn_tag(input[[id("compare_topn")]])),
    content = function(file) {
      u <- resolved_r(); df <- ora_table_r()
      tn  <- .enrich_topn_val(input[[id("compare_topn")]])          # §11: match display
      fig <- if (is.null(df)) NULL else .enrich_compare_download_figure(df, u, top_n = tn)
      if (is.null(fig)) { grDevices::pdf(file, width = 6, height = 4); grDevices::dev.off(); return() }
      ggplot2::ggsave(file, fig$plot, device = "pdf",
                      width = fig$width, height = fig$height, limitsize = FALSE)
    })
  # §12: Download Table — the COMPLETE compare data (every group's every
  # significant term, NOT the on-screen Top-N subset), with the group-identifying
  # column (Up/Down · Cluster · Pattern · Contrast). Excel, same convention as the
  # other enrichment table downloads.
  output[[paste0("download_", id_prefix, "_compare_table")]] <- shiny::downloadHandler(
    filename = function() .dl_name("compare", "xlsx"),
    content = function(file) {
      u <- resolved_r(); df <- ora_table_r()
      out <- if (is.null(df)) data.frame() else .enrich_compare_table_df(df, u)
      writexl::write_xlsx(out, path = file)
    })
  # Single-group ORA Enrichment Plot — figure MATCHES the on-screen Top-N
  # selection (§18 D18.1): Top 20 -> top-20 figure, All -> all terms.
  output[[paste0("download_", id_prefix, "_enrich_bar")]] <- shiny::downloadHandler(
    filename = function() .dl_name("enrichment_plot", "pdf", .dl_topn_tag(input[[id("enrich_topn")]])),
    content = function(file) {
      df <- ora_view_table_r()                                 # §23 A: active view
      n  <- if (is.null(df)) 0L else .enrich_bar_shown_n(df, input[[id("enrich_topn")]])
      p  <- if (is.null(df) || nrow(df) < 1) NULL else plot_ora_enrichment_bar(df, max_terms = n)
      if (is.null(p)) { grDevices::pdf(file, width = 6, height = 4); grDevices::dev.off(); return() }
      ggplot2::ggsave(file, p, device = "pdf", width = 11,
                      height = max(5, 0.28 * n + 2), limitsize = FALSE)
    })
  # §23 C/D: Shared Genes FIGURE — Term × Gene matrix MATCHING the on-screen Top-N
  # gene selection. The complete matrix is the separate Download Table.
  output[[paste0("download_", id_prefix, "_shared_genes_fig")]] <- shiny::downloadHandler(
    filename = function() .dl_name("gene_term_matrix", "pdf", .dl_topn_tag(input[[id("sg_topn")]])),
    content = function(file) {
      df <- ora_view_table_r()
      p <- if (is.null(df)) NULL else
        plot_ora_shared_genes_matrix_static(df, .sg_topn())
      if (is.null(p)) { grDevices::pdf(file, width = 6, height = 4); grDevices::dev.off(); return() }
      n_g <- length(.enrich_sg_keep_genes(df, .sg_topn())); n_t <- nrow(df)
      ggplot2::ggsave(file, p, device = "pdf",
                      width = .enrich_sg_static_in(n_g),
                      height = .enrich_sg_static_in(n_t), limitsize = FALSE)
    })
  # Shared Terms FIGURE — Dice heatmap MATCHING the on-screen Top-N selection
  # (§18 D18.3). Static ggplot mirror of the plotly view.
  output[[paste0("download_", id_prefix, "_shared_terms_fig")]] <- shiny::downloadHandler(
    filename = function() .dl_name("term_term_matrix", "pdf", .dl_topn_tag(input[[id("t2t_topn")]])),
    content = function(file) {
      s <- ora_shared_r()
      p <- if (is.null(s)) NULL else plot_ora_term_overlap_static(s, .t2t_topn())
      if (is.null(p)) { grDevices::pdf(file, width = 6, height = 4); grDevices::dev.off(); return() }
      n <- length(.enrich_t2t_keep(s, .t2t_topn()))
      din <- .enrich_t2t_static_in(n)
      ggplot2::ggsave(file, p, device = "pdf", width = din, height = din,
                      limitsize = FALSE)
    })
  # Shared Terms TABLE — the COMPLETE term×term Sørensen–Dice matrix (all terms,
  # never restricted by the Top-N view; §18 D18.3/D18.4), Dice cluster-ordered.
  output[[paste0("download_", id_prefix, "_shared_terms")]] <- shiny::downloadHandler(
    filename = function() .dl_name("term_term_matrix", "xlsx"),
    content = function(file) {
      s <- ora_shared_r()
      m <- if (is.null(s) || is.null(s$dice)) NULL else s$dice
      if (!is.null(m)) { ord <- .enrich_term_order(m); m <- m[ord, ord, drop = FALSE] }
      out <- if (is.null(m)) data.frame() else
        cbind(Term = rownames(m), as.data.frame(round(m, 4), check.names = FALSE))
      writexl::write_xlsx(out, path = file)
    })
  # Selected-term member genes (DE table). Download follows the on-screen rule:
  # single-contrast unit -> only that contrast's stats (pass-3 item 1).
  output[[paste0("download_", id_prefix, "_member_genes")]] <- shiny::downloadHandler(
    filename = function() .dl_name("member_genes", "xlsx"),
    content = function(file) {
      m <- ora_member_de_r()
      # §20/§21: same structure as the on-screen table — single-contrast narrow
      # (Gene|linearFC|pvalue|padj) or multi-contrast WIDE (per-contrast blocks).
      # Complete gene rows; pagination is display-only.
      writexl::write_xlsx(if (is.null(m)) data.frame() else .enrich_de_display_df(m),
                          path = file)
    })
  # --- GSEA downloads (item 10): unit plots + selected-pathway detail ---------
  # §18 D18.1: the GSEA enrichment figure MATCHES the on-screen Top-N selection.
  output[[paste0("download_", id_prefix, "_gsea_bar")]] <- shiny::downloadHandler(
    filename = function() .dl_name("gsea_plot", "pdf", .dl_topn_tag(input[[id("gsea_enrich_topn")]])),
    content = function(file) {
      u <- resolved_r(); df <- gsea_table_r()
      n <- if (is.null(df)) 0L else .enrich_bar_shown_n(df, input[[id("gsea_enrich_topn")]])
      p <- if (is.null(df) || nrow(df) < 1) NULL else {
        top <- utils::head(df[order(df$p.adjust, -abs(df$NES)), , drop = FALSE], n)
        le  <- vapply(as.character(top$ID), function(pid)
                 length(enrichment_gsea_leading_edge_ids(enr_r(), u, pid)), integer(1))
        plot_gsea_enrichment_bar(df, le_counts = le, max_terms = n) }
      .ggsave_pdf(file, p, h = max(5, 0.28 * n + 2))
    })
  output[[paste0("download_", id_prefix, "_gsea_le_ridge")]] <- shiny::downloadHandler(
    filename = function() .dl_name("gsea_leading_edge_ridge", "pdf"),
    content = function(file) {
      u <- resolved_r(); df <- gsea_table_r()
      rd <- if (is.null(u) || is.null(df)) NULL else build_gsea_ridge_df(enr_r(), u, df, mode = "leading_edge")
      p  <- if (is.null(rd) || !isTRUE(rd$fc_available) || is.null(rd$df)) NULL else plot_gsea_ridge(rd$df)
      # §23 K: taller canvas so the downloaded ridge has the same spacing on screen
      .ggsave_pdf(file, p, h = max(5, 0.62 * .enrich_ridge_n(rd$df) + 1.5))
    })
  output[[paste0("download_", id_prefix, "_gsea_all_ridge")]] <- shiny::downloadHandler(
    filename = function() .dl_name("gsea_all_genes_ridge", "pdf"),
    content = function(file) {
      rd <- gsea_all_ridge_r()
      p  <- if (is.null(rd) || !isTRUE(rd$fc_available) || is.null(rd$df)) NULL else plot_gsea_ridge(rd$df)
      .ggsave_pdf(file, p, h = max(5, 0.62 * .enrich_ridge_n(rd$df) + 1.5))
    })
  output[[paste0("download_", id_prefix, "_gsea_running")]] <- shiny::downloadHandler(
    filename = function() .dl_name("gsea_running_score", "pdf"),
    content = function(file) {
      row <- gsea_row_r(); prep <- gsea_running_r()
      p <- if (is.null(row) || is.null(prep)) NULL else
        plot_gsea_running_score(prep, title = as.character(row$Description),
          nes = suppressWarnings(as.numeric(row$NES)),
          padj = suppressWarnings(as.numeric(row$p.adjust)))
      .ggsave_pdf(file, p, w = 10, h = 7)
    })
  output[[paste0("download_", id_prefix, "_gsea_core_genes")]] <- shiny::downloadHandler(
    filename = function() .dl_name("gsea_core_genes", "xlsx"),
    content = function(file) {
      m <- gsea_core_de_r()
      # §20/§21: same structure as the on-screen table (narrow single-contrast or
      # multi-contrast WIDE); complete gene rows; pagination is display-only.
      writexl::write_xlsx(if (is.null(m)) data.frame() else .enrich_de_display_df(m),
                          path = file)
    })
  output[[paste0("download_", id_prefix, "_gsea_core_hm")]] <- shiny::downloadHandler(
    filename = function() .dl_name("gsea_core_heatmap", "xlsx"),
    content = function(file) .hm_xlsx(file, gsea_core_ids_r()))
  output[[paste0("download_", id_prefix, "_gsea_all_hm")]] <- shiny::downloadHandler(
    filename = function() .dl_name("gsea_all_heatmap", "xlsx"),
    content = function(file) .hm_xlsx(file, gsea_member_ids_r()))
  output[[paste0("download_", id_prefix, "_member_hm")]] <- shiny::downloadHandler(
    filename = function() .dl_name("de_gene_heatmap", "xlsx"),
    content = function(file) .hm_xlsx(file, ora_member_ids_r()))
  # §23 J3: All-Gene Heatmap data (complete pathway membership ∩ expression).
  output[[paste0("download_", id_prefix, "_all_hm")]] <- shiny::downloadHandler(
    filename = function() .dl_name("all_gene_heatmap", "xlsx"),
    content = function(file) .hm_xlsx(file, ora_all_gene_ids_r()))

  # ---- selector observers (progressive reveal + downstream reset) ----
  shiny::observeEvent(input[[id("analysis")]], {
    pending_target(NULL)
    st(enrichment_apply_selection(st(), "analysis", input[[id("analysis")]]))
  }, ignoreInit = TRUE)
  shiny::observeEvent(input[[id("database")]], {
    pend <- pending_target()
    st(enrichment_set_database(st(), input[[id("database")]], pend))
    if (!is.null(pend)) pending_target(NULL)
  }, ignoreInit = TRUE)
  shiny::observeEvent(input[[id("group")]],
    st(enrichment_apply_selection(st(), "group", input[[id("group")]])), ignoreInit = TRUE)
  shiny::observeEvent(input[[id("item")]],
    st(enrichment_apply_selection(st(), "item", input[[id("item")]])), ignoreInit = TRUE)
  shiny::observeEvent(input[[id("subgroup")]],
    st(enrichment_apply_selection(st(), "subgroup", input[[id("subgroup")]])), ignoreInit = TRUE)
  shiny::observeEvent(input[[id("change")]], {
    pending_target(NULL)
    st(enrichment_reset_level(st(), input[[id("change")]]))
  }, ignoreInit = TRUE)

  # ---- cross-tab controller (partial navigation) ----
  open_fn <- function(analysis, database = NULL, group = NULL, item = NULL,
                      ranking = NULL, contrast = NULL) {
    s <- enrichment_open_state(analysis, database, group, item, ranking, contrast)
    st(s)
    if (is.null(database) && !is.null(s$group))
      pending_target(list(group = s$group, item = s$item))
    else pending_target(NULL)
    tryCatch(bslib::nav_select(id = tabs_id, selected = "Enrichment"),
             error = function(e) NULL)
    invisible(NULL)
  }
  reset_fn <- function() {
    pending_target(NULL); sel_id(NULL); sel_group(NULL)
    st(enrichment_init_state()); invisible(NULL)
  }

  list(open = open_fn, reset = reset_fn,
       # read-only accessor; selected_id is merged from its own reactiveVal
       state = function() { s <- st(); s$selected_id <- sel_id(); s },
       # §4/§5: read-only accessors for the result-set navigation state (the active
       # ORA result set + the per-set remembered inner content tab) — for tests.
       resultset = function() rv_resultset(),
       inner_tabs = function() .rv_inner(),
       max_n = if (!is.null(mge_controller)) mge_controller$max_n else 100L)
}

