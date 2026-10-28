# Define server logic (started by server() at the end of this file, after the Appstream check)
app_server <- function(input, output, session) {

  # Pre-register plotly click events for all source IDs that will be used with

  # event_data("plotly_click", source = ...).  The plots themselves call
  # event_register(p, "plotly_click") on their plotly objects, but that
  # registration only reaches the session after the widget renders in the
  # browser.  Because wire_click_to_explorer creates observeEvent(event_data())
  # at server start, plotly's onFlushed check fires before the plots render and
  # emits a warning for each unregistered source.  Populating the session-level
  # registry here prevents those startup warnings.
  plotly_click_sources <- c(
    "rna_de_volcano",  "rna_de_ma",
    "clust_hc",        "clust_pc",        "clust_bp",
    "prot_de_volcano", "prot_de_ma",
    "prot_clust_hc",   "prot_clust_pc",   "prot_clust_bp",
    "metab_de_volcano","metab_de_ma",
    "metab_clust_hc",  "metab_clust_pc",  "metab_clust_bp"
  )
  session$userData$plotlyShinyEventIDs <- unique(c(
    session$userData$plotlyShinyEventIDs,
    paste("plotly_click", plotly_click_sources, sep = "-")
  ))

  # Sources that also support box/lasso selection (Multi-* Feature flows).
  plotly_selected_sources <- c(
    "rna_de_volcano",   "rna_de_ma",
    "prot_de_volcano",  "prot_de_ma",
    "metab_de_volcano", "metab_de_ma"
  )
  session$userData$plotlyShinyEventIDs <- unique(c(
    session$userData$plotlyShinyEventIDs,
    paste("plotly_selected", plotly_selected_sources, sep = "-")
  ))

  # ---------- tab panel 1 - Home ----------

  # Researcher-specific Home content (home_info_ui() in global.R). Inserted here, inside
  # app_server(), so it is only sent after a successful Appstream check (or locally).
  # app_server() runs once per session, so this inserts exactly once. immediate = TRUE
  # sends it right away instead of with the session's first (slow) batch of outputs.
  insertUI("#home_info_ui", where = "afterBegin", ui = home_info_ui(), immediate = TRUE)

  # ---------- tab panel 2 - Experiment design ----------
  
  # Load the experiment design Excel (sample annotations table)
  excel_path <- "data/experiment_design.xlsx"
  experiment_design <- if (file.exists(excel_path)) read_excel(excel_path) else NULL
  
  # Render the Experiment design table with YES checkmarks
  output$col_table <- renderDT({
    validate(
      need(file.exists(excel_path), "Experiment design table not available"),
      need(!is.null(experiment_design), "Experiment design table not available")
    )
    
    df <- experiment_design
    assay_cols <- c("RNAseq", "Proteomics", "Metabolomics", "Lipidomics")
    
    # If expected assay columns are missing, show the raw table
    if (!all(assay_cols %in% colnames(df))) {
      return(datatable(df, options = list(dom = "t"), rownames = FALSE))
    }
    
    # Replace YES with a green checkmark
    df[assay_cols] <- lapply(df[assay_cols], function(x) {
      ifelse(
        x == "YES",
        "<span style='color: green; font-size: 24px;'>&#10004;</span>",
        ""
      )
    })
    
    datatable(
      df,
      options  = list(dom = "t", pageLength = nrow(df)),
      rownames = FALSE,
      escape   = FALSE
    )
  })
  
  
  # ---------- tab panel 3 - Data & Files ----------
  
  # Read the documentation manifest (Excel / CSV) as a reactive object
  doc_manifest <- reactive({
    read_Data_and_Files_manifest(DOC_MANIFEST)
  })
  
  # Build the Data & Files page UI dynamically from the manifest
  output$doc_page <- renderUI({
    df <- doc_manifest()
    
    validate(
      need(!is.null(df) && nrow(df) > 0, "No Data & files are available right now.")
    )
    
    # 1) Files with NO category -> show without a category header
    df_no_cat <- df[is.na(df$category), , drop = FALSE]
    
    block_no_cat <- NULL
    if (nrow(df_no_cat) > 0) {
      rows <- lapply(seq_len(nrow(df_no_cat)), function(i) {
        item <- df_no_cat[i, , drop = FALSE]
        out_id <- paste0("doc_dl_", item$.id)
        doc_file_row_ui(item, output_id_download = out_id)
      })
      
      block_no_cat <- tags$div(
        tags$div(rows),
        tags$div(style = "height: 10px;")
      )
    }
    
    # 2) Categories in order of FIRST appearance in the Excel
    cats <- df$category[!is.na(df$category)]
    cats <- cats[!duplicated(cats)]
    
    # 3) For each category, show ALL its files (even if scattered in Excel)
    blocks_cat <- lapply(cats, function(cat) {
      sub <- df[!is.na(df$category) & df$category == cat, , drop = FALSE]
      if (nrow(sub) == 0) return(NULL)
      
      rows <- lapply(seq_len(nrow(sub)), function(i) {
        item <- sub[i, , drop = FALSE]
        out_id <- paste0("doc_dl_", item$.id)
        doc_file_row_ui(item, output_id_download = out_id)
      })
      
      tags$div(
        style = "margin-bottom: 18px;",
        tags$hr(style = "margin-top: 6px; margin-bottom: 12px;"),
        h4(strong(tags$u(cat))),
        tags$div(rows)
      )
    })
    
    tagList(block_no_cat, do.call(tagList, blocks_cat))
  })
  
  # Create one downloadHandler per manifest row (suitable for a small number of files)
  observe({
    df <- doc_manifest()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    for (i in seq_len(nrow(df))) {
      local({
        item <- df[i, , drop = FALSE]
        out_id <- paste0("doc_dl_", item$.id)
        src <- item$.full_path[1]
        
        # Download a documentation file from the server filesystem
        output[[out_id]] <- downloadHandler(
          filename = function() {
            if (is.na(src) || !nzchar(src)) return("file")
            basename(src)
          },
          content = function(file) {
            validate(need(!is.na(src) && file.exists(src), "File not found on server."))
            file.copy(src, file, overwrite = TRUE)
          }
        )
      })
    }
  })
  

  
  
  # ---------- tab panel 4 - RNAseq ----------
  
  ## ---- Subtab: Exploratory Analyses ----
  ## ---- Subtab: Data structure ----
  
  # Render 3D PCA iframe
  output$pca_3d_ui <- plotly::renderPlotly({
    shiny::validate(
      shiny::need(exists("pca_3d", envir = rnaseq_env, inherits = FALSE),
                  "pca_3d not available.")
    )
    get("pca_3d", envir = rnaseq_env, inherits = FALSE)
  }) |> shiny::bindCache("rna_pca_3d")

  output$pca_pc1_pc2_plot <- plotly::renderPlotly(
    make_pca_plotly(rnaseq_env, 1, 2)
  ) |> shiny::bindCache("rna_pca_12")

  output$pca_pc1_pc3_plot <- plotly::renderPlotly(
    make_pca_plotly(rnaseq_env, 1, 3)
  ) |> shiny::bindCache("rna_pca_13")

  output$sample_corr_heatmap <- renderPlotly(
    make_sample_corr_plotly(rnaseq_env)
  ) |> shiny::bindCache("rna_corr_heatmap")
  
  
  ## ---- Subtab: DE Genes ----
  
  observe({
    update_de_contrast_choices_from_env(session, rnaseq_env, "de_contrast")
  })
  
  output$volcano_plot <- plotly::renderPlotly({
    req(input$de_contrast)
    plot_volcano_plotly_from_env(rnaseq_env, input$de_contrast, backend = "deseq2", source_id = "rna_de_volcano")
  })
  

  output$ma_plot <- plotly::renderPlotly({
    req(input$de_contrast)
    plot_ma_plotly_from_env(rnaseq_env, input$de_contrast, backend = "deseq2", source_id = "rna_de_ma")
  })
  
  

  
  ## ---- Subtab: Stats results ----

  wire_stats_results(
    output   = output,
    env      = rnaseq_env,
    data_dir = "data/RNAseq",
    id_prefix = "RNAseq",
    
    obj_de_stats = "de_summary",
    obj_de_mat   = "de_expr_norm",
    obj_all_mats = list(norm = "expr_norm"),
    
    out_de_table_id         = "DE_genes_stats",
    out_results_buttons_id  = "results_buttons",
    btn_results_all_id      = "dl_results_all",
    btn_results_de_id       = "dl_results_de",
    dl_all_mat_btn_ids      = list(norm = "dl_norm_all_txt"),
    dl_de_mat_btn_id        = "dl_norm_de_txt",
    # Excel exports now consume the preformatted, Excel-ready objects
    # `all_final_xlsx` / `de_final_xlsx` produced upstream. write.xlsx
    # handles either a single data frame or a named list of sheets, so
    # the same handler works whether the upstream object is one table or
    # multiple.
    obj_results_all = "all_final_xlsx",
    obj_results_de  = "de_final_xlsx",
    file_label_all = "genes",
    file_label_de  = "genes"
  )

 
  
  ## ---- Subtab: Clustering ----
  #  static figures
    output$hc_static_fig_rna <- shiny::renderPlot({
      shiny::validate(shiny::need(
        clust_static_fig_available_env(rnaseq_env, "clust_heatmap_hier_fig"),
        "Clustering results are not available for this dataset."
      ))
      tryCatch({
        grid::grid.newpage()
        grid::grid.draw(get("clust_heatmap_hier_fig", envir = rnaseq_env, inherits = FALSE))
      }, error = function(e) {
        shiny::validate(shiny::need(FALSE, "Clustering figure could not be rendered."))
      })
    })

  output$pc_static_fig_rna <- shiny::renderPlot({
    shiny::validate(shiny::need(
      clust_static_fig_available_env(rnaseq_env, "clust_heatmap_partition_fig"),
      "Clustering results are not available for this dataset."
    ))
    tryCatch({
      grid::grid.newpage()
      grid::grid.draw(get("clust_heatmap_partition_fig", envir = rnaseq_env, inherits = FALSE))
    }, error = function(e) {
      shiny::validate(shiny::need(FALSE, "Clustering figure could not be rendered."))
    })
  })
  
  
  rna_clustering_ctrl <- wire_clustering(
    session = session,
    input   = input,
    output  = output,
    env     = rnaseq_env,
    data_dir = "data/RNAseq",
    id = list(
      heatmap_hc            = "de_heatmap",
      pc_title_ui           = "pc_title_ui",
      pc_cluster_pick       = "pc_cluster_pick",
      pc_cluster_title_ui   = "pc_cluster_title_ui",
      pc_cluster_note_ui    = "pc_cluster_note_ui",
      pc_heatmap_cluster    = "pc_heatmap_cluster",
      pc_cluster_subgraph   = "pc_cluster_subgraph",
      pc_results_title_ui   = "pc_results_title_ui",
      dl_results_de_ui      = "dl_results_de_filtered_ui",
      dl_results_de_btn     = "dl_results_de_filtered",
      bp_pattern_ui         = "bp_pattern_ui",
      bp_static_title_ui    = "bp_static_title_ui",
      bp_interactive_title_ui = "bp_interactive_title_ui",
      bp_static             = "bp_static",
      bp_interactive        = "bp_interactive",
      bp_results_title_ui   = "bp_results_title_ui",
      dl_results_de_bp_ui   = "dl_results_de_bp_filtered_ui",
      dl_results_de_bp_btn  = "dl_results_de_bp_filtered"
    ),
    source = list(
      hc = "clust_hc",
      pc = "clust_pc",
      bp = "clust_bp"
    ),
    label = list(
      feature = "genes",
      results_de_obj = "de_final_table"
    )
  )

  ## ---- Subtab: Enrichment (ORA/GSEA — Stage 4) ----
  ## Wired below (after the Gene/Multi-Gene Explorer wirings) via
  ## wire_enrichment(); see the call after the Multi-Gene Explorer block.


  ## ---- Subtab: Gene Explorer ----
  wire_gene_explorer(
    session = session,
    input   = input,
    output  = output,
    env     = rnaseq_env,
    tabs_id = "rna_tabs",
    # plotly sources that should route into Gene Explorer:
    click_sources = c("rna_de_volcano", "rna_de_ma", "clust_hc", "clust_pc", "clust_bp")
  )

  # "Clear selection" button next to the Gene Explorer selectize. Resetting
  # the input to character(0) flips the conditionalPanel condition
  # ("input.gene_explorer_gene != '' && != null") to false, so the panels
  # naturally return to the empty state.
  #
  # IMPORTANT: pass `choices` alongside `selected = character(0)` and
  # `server = TRUE` — same pattern wire_multi_gene_explorer uses after Add.
  # Calling updateSelectizeInput with server=TRUE WITHOUT choices leaves
  # the server-side data source unregistered, after which the selectize's
  # dropdown can no longer open ("can't open after clear" bug).
  shiny::observeEvent(input$gene_explorer_clear_btn, {
    ch <- if (exists("expr_raw", envir = rnaseq_env, inherits = FALSE)) {
      prepare_id_choices_from_matrix(get("expr_raw", envir = rnaseq_env, inherits = FALSE))
    } else NULL
    shiny::updateSelectizeInput(session, "gene_explorer_gene",
                                choices  = ch,
                                selected = character(0),
                                server   = TRUE)
  }, ignoreInit = TRUE)

  ## ---- Subtab: Multi-Gene Explorer (RNAseq only) ----
  mge_controller <- wire_multi_gene_explorer(
    session             = session,
    input               = input,
    output              = output,
    env                 = rnaseq_env,
    id_prefix           = "mge",
    max_n               = 100,
    # Clickable feature names in the Expression Data table navigate to the
    # single Gene Explorer using the same modal UX as volcano/MA clicks.
    single_tab_label    = "Gene Explorer",
    tabs_id             = "rna_tabs",
    single_selectize_id = "gene_explorer_gene"
  )

  ## ---- Subtab: Enrichment (ORA/GSEA — Stage 4A) ----
  # Registered after the Gene/Multi-Gene Explorer wirings so any future deferred
  # (priority=-100) selectize hand-off lands after the Explorer's own resync.
  enrich_controller <- wire_enrichment(
    session             = session,
    input               = input,
    output              = output,
    env                 = rnaseq_env,
    id_prefix           = "rna_enrich",
    tabs_id             = "rna_tabs",
    mge_controller      = mge_controller,
    single_selectize_id = "gene_explorer_gene",
    # 4F: selection-driven Enrichment -> Explorer hand-off (same tabs as Keyword
    # Search). Gene Explorer for 1 selected gene, Multi-Gene Explorer for 2+.
    single_tab_label    = "Gene Explorer",
    multi_tab_label     = "Multi-Gene Explorer"
  )

  ## ---- Proteomics & Metabolomics Enrichment (SAME generic code) ----
  # Identical omics-agnostic wiring, only env/id_prefix/tabs_id differ — no
  # per-omics enrichment implementation, no ORA/GSEA duplication. These envs
  # currently carry NO `enrichment` payload, so the generic panel shows only the
  # friendly "not available" message (enrichment_data_available_env() -> FALSE).
  # When the same §14 payload is added to their RDS upstream, this exact code
  # activates the full ORA/GSEA interface with no further change. The Explorer /
  # clustering hand-off knobs (mge_controller, single_selectize_id) are omics-
  # specific and deferred to 4F — omitted here (the view works unavailable and,
  # once a payload exists, resolved WITHOUT an Explorer controller).
  prot_enrich_controller <- wire_enrichment(
    session   = session,
    input     = input,
    output    = output,
    env       = proteomics_env,
    id_prefix = "prot_enrich",
    tabs_id   = "prot_tabs"
  )
  metab_enrich_controller <- wire_enrichment(
    session   = session,
    input     = input,
    output    = output,
    env       = metabolomics_env,
    id_prefix = "metab_enrich",
    tabs_id   = "metab_tabs"
  )

  ## ---- Subtab: Keyword Search (RNAseq) ----
  # Registered after wire_gene_explorer + wire_multi_gene_explorer above so the
  # post-nav selectize write lands after the Explorer's own choices resync.
  #
  # NOTE: the current RNAseq payload ships feature_annot with 4928 rows and ZERO
  # annotation columns, so this panel renders its unavailable-state banner
  # ("...contains no annotation columns, so there is nothing to search"). That is
  # expected: per the agreed design, Keyword Search never falls back to
  # searching technical feature IDs. It becomes functional as soon as a payload
  # carries real annotation columns -- no code change required.
  wire_keyword_search(
    session              = session,
    input                = input,
    output               = output,
    env                  = rnaseq_env,
    id_prefix            = "rna_kws",
    feature_label        = "gene",
    feature_label_plural = "genes",
    feature_col_name     = "Gene",
    tabs_id              = "rna_tabs",
    single_tab_label     = "Gene Explorer",
    single_selectize_id  = "gene_explorer_gene",
    multi_tab_label      = "Multi-Gene Explorer",
    multi_controller     = mge_controller
  )

  # Box/lasso selection on RNAseq volcano + MA plots:
  #  - 1 gene  -> existing single-gene Gene Explorer prompt (different button id
  #               than the click flow, but same UX: modal -> "View in Gene Explorer")
  #  - 2..20 genes -> new Multi-Gene Explorer prompt
  wire_multi_select_to_explorer(
    session              = session,
    input                = input,
    output               = output,
    tabs_id              = "rna_tabs",
    single_tab_label     = "Gene Explorer",
    multi_tab_label      = "Multi-Gene Explorer",
    single_selectize_id  = "gene_explorer_gene",
    single_choices_vec   = if (exists("expr_raw", envir = rnaseq_env, inherits = FALSE))
                             prepare_id_choices_from_matrix(
                               get("expr_raw", envir = rnaseq_env, inherits = FALSE)
                             ) else NULL,
    mge_controller       = mge_controller,
    select_sources       = c("rna_de_volcano", "rna_de_ma")
  )

  # ============================================================================
  # HEATMAP-ZOOM MULTI-SELECT (revertible) — HC + Partition + Binary Pattern
  # When the user zooms a clustering heatmap on the y-axis, capture the visible
  # gene range and let them open those genes in the Multi-Gene Explorer.
  #
  # To revert: remove this block, the matching UI blocks in ui.R, the small
  # additions to wire_clustering() that expose pc_gene_order/bp_gene_order,
  # and the helper `wire_heatmap_zoom_to_multi_gene_explorer()` in
  # functions_shiny.R.
  # ============================================================================

  # HC: gene order is constant — derived from rnaseq_env objects only.
  hc_zoom_gene_order <- shiny::reactive({
    if (!exists("de_expr_norm",       envir = rnaseq_env, inherits = FALSE)) return(character(0))
    if (!exists("clust_heatmap_hier", envir = rnaseq_env, inherits = FALSE)) return(character(0))
    prepped <- tryCatch(
      prepare_de_heatmap_inputs(
        de_expr_norm       = get("de_expr_norm",       envir = rnaseq_env, inherits = FALSE),
        clust_heatmap_hier = get("clust_heatmap_hier", envir = rnaseq_env, inherits = FALSE)
      ),
      error = function(e) NULL
    )
    if (is.null(prepped) || is.null(prepped$rn)) return(character(0))
    as.character(prepped$rn)
  })

  # HC for Proteomics — derived from proteomics_env objects only.
  prot_hc_zoom_gene_order <- shiny::reactive({
    if (!exists("de_expr_norm",       envir = proteomics_env, inherits = FALSE)) return(character(0))
    if (!exists("clust_heatmap_hier", envir = proteomics_env, inherits = FALSE)) return(character(0))
    prepped <- tryCatch(
      prepare_de_heatmap_inputs(
        de_expr_norm       = get("de_expr_norm",       envir = proteomics_env, inherits = FALSE),
        clust_heatmap_hier = get("clust_heatmap_hier", envir = proteomics_env, inherits = FALSE)
      ),
      error = function(e) NULL
    )
    if (is.null(prepped) || is.null(prepped$rn)) return(character(0))
    as.character(prepped$rn)
  })

  # HC for Metabolomics — derived from metabolomics_env objects only.
  metab_hc_zoom_gene_order <- shiny::reactive({
    if (!exists("de_expr_norm",       envir = metabolomics_env, inherits = FALSE)) return(character(0))
    if (!exists("clust_heatmap_hier", envir = metabolomics_env, inherits = FALSE)) return(character(0))
    prepped <- tryCatch(
      prepare_de_heatmap_inputs(
        de_expr_norm       = get("de_expr_norm",       envir = metabolomics_env, inherits = FALSE),
        clust_heatmap_hier = get("clust_heatmap_hier", envir = metabolomics_env, inherits = FALSE)
      ),
      error = function(e) NULL
    )
    if (is.null(prepped) || is.null(prepped$rn)) return(character(0))
    as.character(prepped$rn)
  })

  # Pre-register plotly_relayout for all clustering sources (RNA + prot + metab)
  # to suppress the "event-not-registered" warning at startup.
  session$userData$plotlyShinyEventIDs <- unique(c(
    session$userData$plotlyShinyEventIDs,
    "plotly_relayout-clust_hc",
    "plotly_relayout-clust_pc",
    "plotly_relayout-clust_bp",
    "plotly_relayout-prot_clust_hc",
    "plotly_relayout-prot_clust_pc",
    "plotly_relayout-prot_clust_bp",
    "plotly_relayout-metab_clust_hc",
    "plotly_relayout-metab_clust_pc",
    "plotly_relayout-metab_clust_bp"
  ))

  wire_heatmap_zoom_to_multi_gene_explorer(
    session             = session,
    input               = input,
    output              = output,
    source_id           = "clust_hc",
    gene_order_reactive = hc_zoom_gene_order,
    mge_controller      = mge_controller,
    tabs_id             = "rna_tabs",
    multi_tab_label     = "Multi-Gene Explorer",
    id_prefix           = "hc_zoom",
    empty_msg           = "Zoom into the heatmap to detect visible genes."
  )

  # Partition: gene order depends on the user's cluster pick — exposed via
  # the wire_clustering controller's pc_gene_order reactive.
  wire_heatmap_zoom_to_multi_gene_explorer(
    session             = session,
    input               = input,
    output              = output,
    source_id           = "clust_pc",
    gene_order_reactive = rna_clustering_ctrl$pc_gene_order,
    mge_controller      = mge_controller,
    tabs_id             = "rna_tabs",
    multi_tab_label     = "Multi-Gene Explorer",
    id_prefix           = "pc_zoom",
    empty_msg           = "Zoom into the partition heatmap to detect visible genes."
  )

  # Binary Pattern: gene order depends on the user's pattern pick — exposed
  # via the wire_clustering controller's bp_gene_order reactive.
  wire_heatmap_zoom_to_multi_gene_explorer(
    session             = session,
    input               = input,
    output              = output,
    source_id           = "clust_bp",
    gene_order_reactive = rna_clustering_ctrl$bp_gene_order,
    mge_controller      = mge_controller,
    tabs_id             = "rna_tabs",
    multi_tab_label     = "Multi-Gene Explorer",
    id_prefix           = "bp_zoom",
    empty_msg           = "Zoom into the binary pattern heatmap to detect visible genes."
  )
  # ============================================================================
  # END HEATMAP-ZOOM MULTI-SELECT (RNA)
  # ----------------------------------------------------------------------------
  # NOTE: Proteomics and Metabolomics clustering zoom invocations live near
  # their respective Multi-*-Explorer controllers (mpe_controller /
  # mme_controller) further down, because controllers must exist before they
  # are passed to wire_heatmap_zoom_to_multi_gene_explorer().
  # ============================================================================

  # Navigate from Gene Explorer to Multi-omics Explorer.
  # Navigation always happens; if the matching table is unavailable, the
  # MultiOmics Explorer page will display its "information not available"
  # message instead of preselecting a query.
  shiny::observeEvent(input$gene_explorer_multi_btn, {
    gene_id <- shiny::isolate(input$gene_explorer_gene)
    shiny::req(!is.null(gene_id) && nzchar(trimws(as.character(gene_id))))
    if (!is.null(multi_matching_tbl)) {
      shiny::updateSelectizeInput(session, "multi_search_query",
                                  choices  = build_multi_search_choices(multi_matching_tbl, proteomics_env),
                                  selected = as.character(gene_id),
                                  server   = TRUE)
    }
    shiny::updateNavbarPage(session, "main_navbar", selected = "MultiOmics")
  }, ignoreInit = TRUE)



  # ---------- tab panel 5 - Proteomics ----------
  
  ## ---- Subtab: Exploratory Analyses ----
  ## ---- Subtab: Data distributions ----  
  
  
  # QC plot(s): observed vs imputed distributions
  output$imp_hist_samp_plot <- shiny::renderPlot({
    shiny::validate(
      shiny::need(exists("imp_hist_samp", envir = proteomics_env, inherits = FALSE),
                  "Proteomics data not available (imp_hist_samp missing).")
    )
    
    get("imp_hist_samp", envir = proteomics_env, inherits = FALSE)
  }) |> shiny::bindCache("qc_plots_static")
  
  
  # Normalized expression boxplots: controls (only if shape exists)
  output$prot_norm_boxplot_controls <- shiny::renderUI({
    shiny::validate(
      shiny::need(exists("expr_long", envir = proteomics_env, inherits = FALSE),
                  "Proteomics data not available (expr_long missing)."),
      shiny::need(exists("color", envir = proteomics_env, inherits = FALSE),
                  "Proteomics configuration missing (color).")
    )
    
    expr_long <- get("expr_long", envir = proteomics_env, inherits = FALSE)
    color_var <- normalize_opt_var(get("color", envir = proteomics_env, inherits = FALSE))
    
    shiny::validate(
      shiny::need(!is.null(color_var) && color_var %in% colnames(expr_long),
                  paste0("expr_long is missing required color column: ", color_var))
    )
    
    # shape is optional
    shape_var <- NULL
    if (exists("shape", envir = proteomics_env, inherits = FALSE)) {
      shape_var <- normalize_opt_var(get("shape", envir = proteomics_env, inherits = FALSE))
    }
    
    build_norm_boxplot_controls(expr_long, color_var, shape_var)
  })


  # Normalized expression boxplots: static ggplot (recomputes on color-by change)
  output$prot_norm_boxplot <- shiny::renderPlot({
    shiny::validate(
      shiny::need(exists("expr_long", envir = proteomics_env, inherits = FALSE),
                  "Proteomics data not available (expr_long missing)."),
      shiny::need(exists("color", envir = proteomics_env, inherits = FALSE),
                  "Proteomics configuration missing (color).")
    )

    expr_long <- get("expr_long", envir = proteomics_env, inherits = FALSE)
    color_var <- normalize_opt_var(get("color", envir = proteomics_env, inherits = FALSE))

    shiny::validate(
      shiny::need(all(c("SampleID", "value") %in% colnames(expr_long)),
                  "expr_long must contain SampleID and value columns."),
      shiny::need(!is.null(color_var) && color_var %in% colnames(expr_long),
                  paste0("expr_long is missing required color column: ", color_var))
    )

    # Default mapping is color; allow override from selector
    map_var <- color_var
    if (!is.null(input$prot_norm_boxplot_color_by) &&
        input$prot_norm_boxplot_color_by %in% colnames(expr_long)) {
      map_var <- input$prot_norm_boxplot_color_by
    }

    build_norm_boxplot_static(expr_long, map_var)
  }, res = 110)
  
  
  # Density plot: precompute once (fast render later)
  prot_dens_df_rv <- shiny::reactiveVal(NULL)
  
  shiny::observeEvent(TRUE, {
    if (!exists("expr_long", envir = proteomics_env, inherits = FALSE)) return()
    expr_long <- get("expr_long", envir = proteomics_env, inherits = FALSE)
    
    dens_df <- compute_density_df(expr_long, sample_col = "SampleID", value_col = "value")
    prot_dens_df_rv(dens_df)
  }, once = TRUE)
  
  
  output$prot_density_plot <- plotly::renderPlotly({
    dens_df <- prot_dens_df_rv()
    
    shiny::validate(
      shiny::need(!is.null(dens_df) && nrow(dens_df) > 0,
                  "Not enough data to compute densities.")
    )
    
    build_density_plotly(dens_df)
  }) |> shiny::bindCache("prot_density_static")
  

  ## ---- Subtab: Data structure ----
  
  # Render 3D PCA iframe 
  # Proteomics
  output$prot_pca_3d_ui <- plotly::renderPlotly({
    shiny::validate(
      shiny::need(exists("pca_3d", envir = proteomics_env, inherits = FALSE),
                  "pca_3d not available.")
    )
    get("pca_3d", envir = proteomics_env, inherits = FALSE)
  }) |> shiny::bindCache("prot_pca_3d")

  output$prot_pca_pc1_pc2_plot <- plotly::renderPlotly(
    make_pca_plotly(proteomics_env, 1, 2)
  ) |> shiny::bindCache("prot_pca_12")

  output$prot_pca_pc1_pc3_plot <- plotly::renderPlotly(
    make_pca_plotly(proteomics_env, 1, 3)
  ) |> shiny::bindCache("prot_pca_13")

  output$prot_sample_corr_heatmap <- renderPlotly(
    make_sample_corr_plotly(proteomics_env)
  ) |> shiny::bindCache("prot_corr_heatmap")
  
  
  ## ---- Subtab: DE Protein ----

  observe({
    update_de_contrast_choices_from_env(session, proteomics_env, "prot_de_contrast")
  })
  
  output$prot_volcano_plot <- plotly::renderPlotly({
    req(input$prot_de_contrast)
    plot_volcano_plotly_from_env(proteomics_env, input$prot_de_contrast, backend = "limma", source_id = "prot_de_volcano")
  })
  
  output$prot_ma_plot <- plotly::renderPlotly({
    req(input$prot_de_contrast)
    plot_ma_plotly_from_env(proteomics_env, input$prot_de_contrast, backend = "limma", source_id = "prot_de_ma")
  })
  
  
  ## ---- Subtab: Stats results ----
  wire_stats_results(
    output = output,
    env    = proteomics_env,
    data_dir = "data/Proteomics",
    id_prefix = "Proteomics",
    
    obj_de_stats = "de_summary",
    obj_de_mat   = "de_expr_norm",
    
    # IMPORTANT: names here ("raw","imputed") become part of the downloaded filename
    obj_all_mats = list(raw = "expr_raw", imputed = "expr_norm"),
    
    out_de_table_id         = "prot_DE_genes_stats",
    out_results_buttons_id  = "prot_results_buttons",
    btn_results_all_id      = "prot_dl_results_all",
    btn_results_de_id       = "prot_dl_results_de",
    
    # IMPORTANT: keys must match obj_all_mats keys exactly
    dl_all_mat_btn_ids      = list(
      raw     = "prot_dl_log2_raw_all_txt",
      imputed = "prot_dl_log2_imp1_all_txt"
    ),
    
    dl_de_mat_btn_id        = "prot_dl_norm_de_txt",
    
    label_results_all = "Result tables - All proteins",
    label_results_de  = "Result tables - DE proteins",
    
    file_label_all = "proteins",
    file_label_de  = "proteins",
    # Excel exports now consume the preformatted, Excel-ready objects.
    obj_results_all = "all_final_xlsx",
    obj_results_de  = "de_final_xlsx"
  )
  
  
  ## ---- Subtab: Clustering ----
  
  # static figures
  output$hc_static_fig_prot <- shiny::renderPlot({
    shiny::validate(shiny::need(
      clust_static_fig_available_env(proteomics_env, "clust_heatmap_hier_fig"),
      "Clustering results are not available for this dataset."
    ))
    tryCatch({
      grid::grid.newpage()
      grid::grid.draw(get("clust_heatmap_hier_fig", envir = proteomics_env, inherits = FALSE))
    }, error = function(e) {
      shiny::validate(shiny::need(FALSE, "Clustering figure could not be rendered."))
    })
  })

  output$pc_static_fig_prot <- shiny::renderPlot({
    shiny::validate(shiny::need(
      clust_static_fig_available_env(proteomics_env, "clust_heatmap_partition_fig"),
      "Clustering results are not available for this dataset."
    ))
    tryCatch({
      grid::grid.newpage()
      grid::grid.draw(get("clust_heatmap_partition_fig", envir = proteomics_env, inherits = FALSE))
    }, error = function(e) {
      shiny::validate(shiny::need(FALSE, "Clustering figure could not be rendered."))
    })
  })
  
  
  prot_clustering_ctrl <- wire_clustering(
    session = session,
    input   = input,
    output  = output,
    env     = proteomics_env,
    data_dir = "data/Proteomics",
    id = list(
      heatmap_hc            = "prot_de_heatmap",
      pc_title_ui           = "prot_pc_title_ui",
      pc_cluster_pick       = "prot_pc_cluster_pick",
      pc_cluster_title_ui   = "prot_pc_cluster_title_ui",
      pc_cluster_note_ui    = "prot_pc_cluster_note_ui",
      pc_heatmap_cluster    = "prot_pc_heatmap_cluster",
      pc_cluster_subgraph   = "prot_pc_cluster_subgraph",
      pc_results_title_ui   = "prot_pc_results_title_ui",
      dl_results_de_ui      = "prot_dl_results_de_filtered_ui",
      dl_results_de_btn     = "prot_dl_results_de_filtered",
      bp_pattern_ui         = "prot_bp_pattern_ui",
      bp_static_title_ui    = "prot_bp_static_title_ui",
      bp_interactive_title_ui = "prot_bp_interactive_title_ui",
      bp_static             = "prot_bp_static",
      bp_interactive        = "prot_bp_interactive",
      bp_results_title_ui   = "prot_bp_results_title_ui",
      dl_results_de_bp_ui   = "prot_dl_results_de_bp_filtered_ui",
      dl_results_de_bp_btn  = "prot_dl_results_de_bp_filtered"
    ),
    source = list(
      hc = "prot_clust_hc",
      pc = "prot_clust_pc",
      bp = "prot_clust_bp"
    ),
    label = list(
      feature = "proteins",
      results_de_obj = "de_final_table"
    )
  )
  
  ## ---- Subtab: Protein Explorer ----
  wire_protein_explorer(
    session = session,
    input   = input,
    output  = output,
    env     = proteomics_env,
    tabs_id = "prot_tabs",
    # plotly sources that should route into Protein Explorer:
    click_sources = c("prot_de_volcano", "prot_de_ma","prot_clust_hc", "prot_clust_pc", "prot_clust_bp")
  )

  # "Clear selection" button next to the Protein Explorer selectize.
  # See gene_explorer_clear_btn comment for why choices must be re-passed.
  shiny::observeEvent(input$prot_explorer_clear_btn, {
    ch <- if (exists("expr_raw", envir = proteomics_env, inherits = FALSE)) {
      prepare_id_choices_from_matrix(get("expr_raw", envir = proteomics_env, inherits = FALSE))
    } else NULL
    shiny::updateSelectizeInput(session, "prot_explorer_protein",
                                choices  = ch,
                                selected = character(0),
                                server   = TRUE)
  }, ignoreInit = TRUE)

  # Navigate from Protein Explorer to Multi-omics Explorer.
  # Navigation always happens; if the matching table is unavailable, the
  # MultiOmics Explorer page will display its "information not available"
  # message instead of preselecting a query.
  shiny::observeEvent(input$prot_explorer_multi_btn, {
    prot_id <- shiny::isolate(input$prot_explorer_protein)
    shiny::req(!is.null(prot_id) && nzchar(trimws(as.character(prot_id))))
    if (!is.null(multi_matching_tbl)) {
      shiny::updateSelectizeInput(session, "multi_search_query",
                                  choices  = build_multi_search_choices(multi_matching_tbl, proteomics_env),
                                  selected = as.character(prot_id),
                                  server   = TRUE)
    }
    shiny::updateNavbarPage(session, "main_navbar", selected = "MultiOmics")
  }, ignoreInit = TRUE)

  ## ---- Subtab: Multi-Protein Explorer (Proteomics) ----
  mpe_controller <- wire_multi_gene_explorer(
    session              = session,
    input                = input,
    output               = output,
    env                  = proteomics_env,
    id_prefix            = "mpe",
    max_n                = 100,
    feature_label        = "protein",
    feature_label_plural = "proteins",
    feature_col_name     = "Protein",
    expression_data_mode = "single_matrix",
    de_stats_fn          = get_protein_de_stats_for_table,
    scale_choices_fn     = function(env) c("without imputation (log2 scale)" = "raw",
                                            "Imputed (log2 scale)"            = "imp1"),
    default_scale        = "raw",
    ylab_log2_fn         = function(env, scale_sel) {
      if (identical(scale_sel, "imp1")) "Log2 protein intensity (imputed)"
      else                              "Log2 protein intensity (raw)"
    },
    theme_variant        = "protein",
    # Clickable feature names in the Expression Data table -> single
    # Protein Explorer, mirrors the volcano/MA single-click modal UX.
    single_tab_label     = "Protein Explorer",
    tabs_id              = "prot_tabs",
    single_selectize_id  = "prot_explorer_protein"
  )

  ## ---- Subtab: Keyword Search (Proteomics) ----
  # Registered after wire_protein_explorer + wire_multi_gene_explorer above so
  # the post-nav selectize write lands after the Explorer's own choices resync.
  wire_keyword_search(
    session              = session,
    input                = input,
    output               = output,
    env                  = proteomics_env,
    id_prefix            = "prot_kws",
    feature_label        = "protein",
    feature_label_plural = "proteins",
    feature_col_name     = "Protein",
    tabs_id              = "prot_tabs",
    single_tab_label     = "Protein Explorer",
    single_selectize_id  = "prot_explorer_protein",
    multi_tab_label      = "Multi-Protein Explorer",
    multi_controller     = mpe_controller
  )

  # Box/lasso multi-selection on Proteomics volcano + MA
  wire_multi_select_to_explorer(
    session                 = session,
    input                   = input,
    output                  = output,
    tabs_id                 = "prot_tabs",
    single_tab_label        = "Protein Explorer",
    multi_tab_label         = "Multi-Protein Explorer",
    single_selectize_id     = "prot_explorer_protein",
    single_choices_vec      = if (exists("expr_raw", envir = proteomics_env, inherits = FALSE))
                                prepare_id_choices_from_matrix(
                                  get("expr_raw", envir = proteomics_env, inherits = FALSE)
                                ) else NULL,
    mge_controller          = mpe_controller,
    select_sources          = c("prot_de_volcano", "prot_de_ma"),
    id_prefix               = "mpe_select",
    feature_label_singular  = "protein",
    feature_label_plural    = "proteins"
  )

  # ============================================================================
  # HEATMAP-ZOOM MULTI-SELECT (Proteomics) -> Multi-Protein Explorer
  # ============================================================================
  wire_heatmap_zoom_to_multi_gene_explorer(
    session             = session,
    input               = input,
    output              = output,
    source_id           = "prot_clust_hc",
    gene_order_reactive = prot_hc_zoom_gene_order,
    mge_controller      = mpe_controller,
    tabs_id             = "prot_tabs",
    multi_tab_label     = "Multi-Protein Explorer",
    id_prefix           = "prot_hc_zoom",
    empty_msg           = "Zoom into the heatmap to detect visible proteins."
  )
  wire_heatmap_zoom_to_multi_gene_explorer(
    session             = session,
    input               = input,
    output              = output,
    source_id           = "prot_clust_pc",
    gene_order_reactive = prot_clustering_ctrl$pc_gene_order,
    mge_controller      = mpe_controller,
    tabs_id             = "prot_tabs",
    multi_tab_label     = "Multi-Protein Explorer",
    id_prefix           = "prot_pc_zoom",
    empty_msg           = "Zoom into the partition heatmap to detect visible proteins."
  )
  wire_heatmap_zoom_to_multi_gene_explorer(
    session             = session,
    input               = input,
    output              = output,
    source_id           = "prot_clust_bp",
    gene_order_reactive = prot_clustering_ctrl$bp_gene_order,
    mge_controller      = mpe_controller,
    tabs_id             = "prot_tabs",
    multi_tab_label     = "Multi-Protein Explorer",
    id_prefix           = "prot_bp_zoom",
    empty_msg           = "Zoom into the binary pattern heatmap to detect visible proteins."
  )
  # ============================================================================
  # END HEATMAP-ZOOM (Proteomics)
  # ============================================================================



  # ---------- tab panel 6 - Metabolomics ----------

  ## ---- Subtab: Exploratory Analyses > Data distributions ----

  output$metab_norm_boxplot_controls <- shiny::renderUI({
    shiny::validate(
      shiny::need(exists("expr_long", envir = metabolomics_env, inherits = FALSE),
                  "Metabolomics data not available (expr_long missing)."),
      shiny::need(exists("color", envir = metabolomics_env, inherits = FALSE),
                  "Metabolomics configuration missing (color).")
    )

    expr_long <- get("expr_long", envir = metabolomics_env, inherits = FALSE)
    color_var <- normalize_opt_var(get("color", envir = metabolomics_env, inherits = FALSE))

    shiny::validate(
      shiny::need(!is.null(color_var) && color_var %in% colnames(expr_long),
                  paste0("expr_long is missing required color column: ", color_var))
    )

    shape_var <- NULL
    if (exists("shape", envir = metabolomics_env, inherits = FALSE)) {
      shape_var <- normalize_opt_var(get("shape", envir = metabolomics_env, inherits = FALSE))
    }

    build_norm_boxplot_controls(expr_long, color_var, shape_var,
                                color_by_input_id = "metab_norm_boxplot_color_by")
  })

  output$metab_norm_boxplot <- shiny::renderPlot({
    shiny::validate(
      shiny::need(exists("expr_long", envir = metabolomics_env, inherits = FALSE),
                  "Metabolomics data not available (expr_long missing)."),
      shiny::need(exists("color", envir = metabolomics_env, inherits = FALSE),
                  "Metabolomics configuration missing (color).")
    )

    expr_long <- get("expr_long", envir = metabolomics_env, inherits = FALSE)
    color_var <- normalize_opt_var(get("color", envir = metabolomics_env, inherits = FALSE))

    shiny::validate(
      shiny::need(all(c("SampleID", "value") %in% colnames(expr_long)),
                  "expr_long must contain SampleID and value columns."),
      shiny::need(!is.null(color_var) && color_var %in% colnames(expr_long),
                  paste0("expr_long is missing required color column: ", color_var))
    )

    map_var <- color_var
    if (!is.null(input$metab_norm_boxplot_color_by) &&
        input$metab_norm_boxplot_color_by %in% colnames(expr_long)) {
      map_var <- input$metab_norm_boxplot_color_by
    }

    build_norm_boxplot_static(expr_long, map_var)
  }, res = 110)

  metab_dens_df_rv <- shiny::reactiveVal(NULL)

  shiny::observeEvent(TRUE, {
    if (!exists("expr_long", envir = metabolomics_env, inherits = FALSE)) return()
    expr_long <- get("expr_long", envir = metabolomics_env, inherits = FALSE)
    dens_df <- compute_density_df(expr_long, sample_col = "SampleID", value_col = "value")
    metab_dens_df_rv(dens_df)
  }, once = TRUE)

  output$metab_density_plot <- plotly::renderPlotly({
    dens_df <- metab_dens_df_rv()
    shiny::validate(
      shiny::need(!is.null(dens_df) && nrow(dens_df) > 0,
                  "Not enough data to compute densities.")
    )
    build_density_plotly(dens_df)
  }) |> shiny::bindCache("metab_density_static")


  ## ---- Subtab: Exploratory Analyses > Data structure ----

  output$metab_pca_3d_ui <- plotly::renderPlotly({
    shiny::validate(
      shiny::need(exists("pca_3d", envir = metabolomics_env, inherits = FALSE),
                  "pca_3d not available.")
    )
    get("pca_3d", envir = metabolomics_env, inherits = FALSE)
  }) |> shiny::bindCache("metab_pca_3d")

  output$metab_pca_pc1_pc2_plot <- plotly::renderPlotly(
    make_pca_plotly(metabolomics_env, 1, 2)
  ) |> shiny::bindCache("metab_pca_12")

  output$metab_pca_pc1_pc3_plot <- plotly::renderPlotly(
    make_pca_plotly(metabolomics_env, 1, 3)
  ) |> shiny::bindCache("metab_pca_13")

  output$metab_sample_corr_heatmap <- renderPlotly(
    make_sample_corr_plotly(metabolomics_env)
  ) |> shiny::bindCache("metab_corr_heatmap")


  ## ---- Subtab: Exploratory Analyses > DE Metabolite ----

  observe({
    update_de_contrast_choices_from_env(session, metabolomics_env, "metab_de_contrast")
  })

  output$metab_volcano_plot <- plotly::renderPlotly({
    req(input$metab_de_contrast)
    plot_volcano_plotly_from_env(metabolomics_env, input$metab_de_contrast, backend = "limma", source_id = "metab_de_volcano")
  }) |> shiny::bindCache(input$metab_de_contrast, "metab_volcano")

  output$metab_ma_plot <- plotly::renderPlotly({
    req(input$metab_de_contrast)
    plot_ma_plotly_from_env(metabolomics_env, input$metab_de_contrast, backend = "limma", source_id = "metab_de_ma")
  }) |> shiny::bindCache(input$metab_de_contrast, "metab_ma")


  ## ---- Subtab: Stats results ----

  wire_stats_results(
    output   = output,
    env      = metabolomics_env,
    data_dir = "data/Metabolomics",
    id_prefix = "Metabolomics",

    obj_de_stats = "de_summary",
    obj_de_mat   = "de_expr_norm",

    obj_all_mats = list(raw = "expr_raw", norm = "expr_norm"),

    out_de_table_id         = "metab_DE_stats",
    out_results_buttons_id  = "metab_results_buttons",
    btn_results_all_id      = "metab_dl_results_all",
    btn_results_de_id       = "metab_dl_results_de",

    dl_all_mat_btn_ids = list(
      raw  = "metab_dl_log2_raw_all_txt",
      norm = "metab_dl_log2_norm_all_txt"
    ),

    dl_de_mat_btn_id = "metab_dl_norm_de_txt",

    label_results_all = "Result tables - All metabolites",
    label_results_de  = "Result tables - DE metabolites",

    file_label_all  = "metabolites",
    file_label_de   = "metabolites",
    # Excel exports now consume the preformatted, Excel-ready objects.
    obj_results_all = "all_final_xlsx",
    obj_results_de  = "de_final_xlsx"
  )


  ## ---- Subtab: Clustering ----

  output$hc_static_fig_metab <- shiny::renderPlot({
    shiny::validate(shiny::need(
      clust_static_fig_available_env(metabolomics_env, "clust_heatmap_hier_fig"),
      "Clustering results are not available for this dataset."
    ))
    tryCatch({
      grid::grid.newpage()
      grid::grid.draw(get("clust_heatmap_hier_fig", envir = metabolomics_env, inherits = FALSE))
    }, error = function(e) {
      shiny::validate(shiny::need(FALSE, "Clustering figure could not be rendered."))
    })
  }) |> shiny::bindCache("metab_hc_static")

  output$pc_static_fig_metab <- shiny::renderPlot({
    shiny::validate(shiny::need(
      clust_static_fig_available_env(metabolomics_env, "clust_heatmap_partition_fig"),
      "Clustering results are not available for this dataset."
    ))
    tryCatch({
      grid::grid.newpage()
      grid::grid.draw(get("clust_heatmap_partition_fig", envir = metabolomics_env, inherits = FALSE))
    }, error = function(e) {
      shiny::validate(shiny::need(FALSE, "Clustering figure could not be rendered."))
    })
  }) |> shiny::bindCache("metab_pc_static")

  metab_clustering_ctrl <- wire_clustering(
    session  = session,
    input    = input,
    output   = output,
    env      = metabolomics_env,
    data_dir = "data/Metabolomics",
    id = list(
      heatmap_hc              = "metab_de_heatmap",
      pc_title_ui             = "metab_pc_title_ui",
      pc_cluster_pick         = "metab_pc_cluster_pick",
      pc_cluster_title_ui     = "metab_pc_cluster_title_ui",
      pc_cluster_note_ui      = "metab_pc_cluster_note_ui",
      pc_heatmap_cluster      = "metab_pc_heatmap_cluster",
      pc_cluster_subgraph     = "metab_pc_cluster_subgraph",
      pc_results_title_ui     = "metab_pc_results_title_ui",
      dl_results_de_ui        = "metab_dl_results_de_filtered_ui",
      dl_results_de_btn       = "metab_dl_results_de_filtered",
      bp_pattern_ui           = "metab_bp_pattern_ui",
      bp_static_title_ui      = "metab_bp_static_title_ui",
      bp_interactive_title_ui = "metab_bp_interactive_title_ui",
      bp_static               = "metab_bp_static",
      bp_interactive          = "metab_bp_interactive",
      bp_results_title_ui     = "metab_bp_results_title_ui",
      dl_results_de_bp_ui     = "metab_dl_results_de_bp_filtered_ui",
      dl_results_de_bp_btn    = "metab_dl_results_de_bp_filtered"
    ),
    source = list(
      hc = "metab_clust_hc",
      pc = "metab_clust_pc",
      bp = "metab_clust_bp"
    ),
    label = list(
      feature        = "metabolites",
      results_de_obj = "de_final_table"
    )
  )


  ## ---- Subtab: Metabolite Explorer ----

  wire_protein_explorer(
    session            = session,
    input              = input,
    output             = output,
    env                = metabolomics_env,
    tabs_id            = "metab_tabs",
    click_sources      = c("metab_de_volcano", "metab_de_ma",
                           "metab_clust_hc", "metab_clust_pc", "metab_clust_bp"),
    id_prefix          = "metab_explorer",
    explorer_tab_label = "Metabolite Explorer",
    feature_col_name   = "Metabolite",
    # Metabolomics no longer exposes the raw/non-normalized option to the
    # user. fixed_scale tells wire_protein_explorer to skip rendering the
    # selectInput AND to use "imp1" everywhere it would have read
    # input[[id_scale]]. The data path is therefore locked to the
    # normalized matrix (expr_norm) without needing any hidden UI widget.
    fixed_scale        = "imp1"
  )

  # Static "Normalization Method:" value shown under the title in the
  # Metabolite Explorer (replaces the former dropdown). Falls back via
  # get_chosen_norm_label_from_env when chosen_norm is absent.
  output$metab_explorer_norm_text <- shiny::renderUI({
    shiny::tags$div(
      get_chosen_norm_label_from_env(metabolomics_env),
      style = "color:#333; font-size:14px; padding-top:2px;"
    )
  })

  # "Clear selection" button next to the Metabolite Explorer selectize.
  # See gene_explorer_clear_btn comment for why choices must be re-passed.
  shiny::observeEvent(input$metab_explorer_clear_btn, {
    ch <- if (exists("expr_raw", envir = metabolomics_env, inherits = FALSE)) {
      prepare_id_choices_from_matrix(get("expr_raw", envir = metabolomics_env, inherits = FALSE))
    } else NULL
    shiny::updateSelectizeInput(session, "metab_explorer_protein",
                                choices  = ch,
                                selected = character(0),
                                server   = TRUE)
  }, ignoreInit = TRUE)

  ## ---- Subtab: Multi-Metabolite Explorer (Metabolomics) ----
  mme_controller <- wire_multi_gene_explorer(
    session              = session,
    input                = input,
    output               = output,
    env                  = metabolomics_env,
    id_prefix            = "mme",
    max_n                = 100,
    feature_label        = "metabolite",
    feature_label_plural = "metabolites",
    feature_col_name     = "Metabolite",
    expression_data_mode = "single_matrix",  # one matrix per choice (expr_raw=raw, expr_norm=imp1)
    de_stats_fn          = get_protein_de_stats_for_table,
    # Metabolomics: no scale dropdown is rendered in ui.R. fixed_scale
    # below locks the data path to the normalized matrix; scale_choices_fn
    # / default_scale are unused when fixed_scale is set (the helper skips
    # the updateSelectInput and the input read), but we keep harmless
    # 1-element defaults for symmetry.
    scale_choices_fn     = function(env) setNames(
      "imp1",
      get_chosen_norm_label_from_env(env)
    ),
    default_scale        = "imp1",
    fixed_scale          = "imp1",
    ylab_log2_fn         = function(env, scale_sel) {
      if (identical(scale_sel, "imp1")) "Log2 metabolite intensity (imputed)"
      else                              "Log2 metabolite intensity (raw)"
    },
    theme_variant        = "protein",
    # Clickable feature names in the Expression Data table -> single
    # Metabolite Explorer, mirrors the volcano/MA single-click modal UX.
    single_tab_label     = "Metabolite Explorer",
    tabs_id              = "metab_tabs",
    single_selectize_id  = "metab_explorer_protein"
  )

  ## ---- Subtab: Keyword Search (Metabolomics) ----
  # Registered AFTER wire_protein_explorer + wire_multi_gene_explorer above:
  # the multi controller is needed, and Shiny delivers queued
  # updateSelectizeInput messages in observer-registration order, so the
  # post-nav selectize write must be registered after the Explorer's own resync.
  wire_keyword_search(
    session              = session,
    input                = input,
    output               = output,
    env                  = metabolomics_env,
    id_prefix            = "metab_kws",
    feature_label        = "metabolite",
    feature_label_plural = "metabolites",
    feature_col_name     = "Metabolite",
    tabs_id              = "metab_tabs",
    single_tab_label     = "Metabolite Explorer",
    single_selectize_id  = "metab_explorer_protein",
    multi_tab_label      = "Multi-Metabolite Explorer",
    multi_controller     = mme_controller
  )

  # Static "Normalization Method:" value shown under the title in the
  # Multi-Metabolite Explorer (replaces the former mme_scale dropdown).
  output$mme_norm_text <- shiny::renderUI({
    shiny::tags$div(
      get_chosen_norm_label_from_env(metabolomics_env),
      style = "color:#333; font-size:14px; padding-top:2px;"
    )
  })

  # Box/lasso multi-selection on Metabolomics volcano + MA
  wire_multi_select_to_explorer(
    session                 = session,
    input                   = input,
    output                  = output,
    tabs_id                 = "metab_tabs",
    single_tab_label        = "Metabolite Explorer",
    multi_tab_label         = "Multi-Metabolite Explorer",
    single_selectize_id     = "metab_explorer_protein",
    single_choices_vec      = if (exists("expr_raw", envir = metabolomics_env, inherits = FALSE))
                                prepare_id_choices_from_matrix(
                                  get("expr_raw", envir = metabolomics_env, inherits = FALSE)
                                ) else NULL,
    mge_controller          = mme_controller,
    select_sources          = c("metab_de_volcano", "metab_de_ma"),
    id_prefix               = "mme_select",
    feature_label_singular  = "metabolite",
    feature_label_plural    = "metabolites"
  )

  # ============================================================================
  # HEATMAP-ZOOM MULTI-SELECT (Metabolomics) -> Multi-Metabolite Explorer
  # ============================================================================
  wire_heatmap_zoom_to_multi_gene_explorer(
    session             = session,
    input               = input,
    output              = output,
    source_id           = "metab_clust_hc",
    gene_order_reactive = metab_hc_zoom_gene_order,
    mge_controller      = mme_controller,
    tabs_id             = "metab_tabs",
    multi_tab_label     = "Multi-Metabolite Explorer",
    id_prefix           = "metab_hc_zoom",
    empty_msg           = "Zoom into the heatmap to detect visible metabolites."
  )
  wire_heatmap_zoom_to_multi_gene_explorer(
    session             = session,
    input               = input,
    output              = output,
    source_id           = "metab_clust_pc",
    gene_order_reactive = metab_clustering_ctrl$pc_gene_order,
    mge_controller      = mme_controller,
    tabs_id             = "metab_tabs",
    multi_tab_label     = "Multi-Metabolite Explorer",
    id_prefix           = "metab_pc_zoom",
    empty_msg           = "Zoom into the partition heatmap to detect visible metabolites."
  )
  wire_heatmap_zoom_to_multi_gene_explorer(
    session             = session,
    input               = input,
    output              = output,
    source_id           = "metab_clust_bp",
    gene_order_reactive = metab_clustering_ctrl$bp_gene_order,
    mge_controller      = mme_controller,
    tabs_id             = "metab_tabs",
    multi_tab_label     = "Multi-Metabolite Explorer",
    id_prefix           = "metab_bp_zoom",
    empty_msg           = "Zoom into the binary pattern heatmap to detect visible metabolites."
  )
  # ============================================================================
  # END HEATMAP-ZOOM (Metabolomics)
  # ============================================================================


  # ---------- tab panel 7 - Lipidomics ----------
  # (No server logic yet — placeholder in UI)
  
  # ---------- tab panel 8 - MultiOmics ----------

  # Show an info message when multiomics matching data is not usable.
  output$multi_unavailable_msg <- shiny::renderUI({
    if (!is.null(multi_matching_tbl)) return(NULL)
    tags$div(
      class = "alert alert-info",
      style = "margin-top:10px; margin-bottom:10px; text-align:center; font-size:15px;",
      "No relevant multiomics matching information is available for the currently loaded datasets."
    )
  })

  # Populate search choices from the gene-protein matching table.
  # Explicitly pass selected = character(0) to prevent server-side selectize
  # from auto-selecting the first choice when choices are first loaded.
  observe({
    req(!is.null(multi_matching_tbl))
    choices <- build_multi_search_choices(multi_matching_tbl, proteomics_env)
    updateSelectizeInput(session, "multi_search_query",
                         choices  = choices,
                         selected = character(0),
                         server   = TRUE)
  })

  # "Clear selection" button next to the MultiOmics search input. Clearing
  # multi_search_query makes multi_pair() return NULL and flips the
  # conditionalPanel ("input.multi_search_query != null && != ''") to
  # false, returning the tab to its empty state.
  # See gene_explorer_clear_btn comment for why choices must be re-passed.
  shiny::observeEvent(input$multi_search_clear_btn, {
    ch <- if (!is.null(multi_matching_tbl)) {
      build_multi_search_choices(multi_matching_tbl, proteomics_env)
    } else NULL
    shiny::updateSelectizeInput(session, "multi_search_query",
                                choices  = ch,
                                selected = character(0),
                                server   = TRUE)
  }, ignoreInit = TRUE)

  # Resolve the gene-protein pair from the user's search query.
  # Returns NULL (not req()) so callers can distinguish "no query" from "no match".
  multi_pair <- reactive({
    query <- input$multi_search_query
    if (is.null(query) || !nzchar(trimws(query))) return(NULL)
    if (is.null(multi_matching_tbl)) return(NULL)
    prot_ids <- NULL
    if (exists("expr_raw", envir = proteomics_env, inherits = FALSE)) {
      prot_ids <- rownames(
        get("expr_raw", envir = proteomics_env, inherits = FALSE))
    } else if (exists("expr_norm", envir = proteomics_env, inherits = FALSE)) {
      prot_ids <- rownames(
        get("expr_norm", envir = proteomics_env, inherits = FALSE))
    }
    resolve_multi_pair(query, multi_matching_tbl, prot_ids)
  })

  # Show resolved gene <-> protein pair as a label line, or a warning when
  # the query is non-empty but no matching pair was found.
  output$multi_pair_info <- renderUI({
    query <- input$multi_search_query
    if (is.null(query) || !nzchar(trimws(query))) return(NULL)
    pair <- multi_pair()
    if (is.null(pair)) {
      return(tags$div(
        class = "alert alert-warning",
        style = "margin-top:10px; margin-bottom:6px;",
        "No matching gene-protein pair found for the selected ID."
      ))
    }
    tags$div(
      style = "text-align:center; margin-top:12px; margin-bottom:2px;",
      tags$span("Gene: ", style = "color:#777; font-size:15px; font-weight:400;"),
      tags$span(pair$gene, style = "color:#18BC9C; font-size:18px; font-weight:700;"),
      tags$span("  \u2194  ", style = "margin:0 10px; font-weight:bold; font-size:16px; color:#555;"),
      tags$span("Protein: ", style = "color:#777; font-size:15px; font-weight:400;"),
      tags$span(pair$protein_full, style = "color:#18BC9C; font-size:18px; font-weight:700;")
    )
  })

  # Empty-state messages for MultiOmics Gene/Protein panels when data is unavailable.
  output$multi_ge_no_data_msg <- shiny::renderUI({
    pair <- multi_pair()
    shiny::req(!is.null(pair))
    gene_id <- pair$gene
    expr <- if (exists("expr_raw",  envir = rnaseq_env, inherits = FALSE)) {
      get("expr_raw",  envir = rnaseq_env, inherits = FALSE)
    } else if (exists("expr_norm", envir = rnaseq_env, inherits = FALSE)) {
      get("expr_norm", envir = rnaseq_env, inherits = FALSE)
    } else NULL
    if (is.null(expr) || !gene_id %in% rownames(expr)) {
      tags$div(class = "alert alert-warning", style = "margin-top:8px; margin-bottom:4px;",
               paste0("No expression data available for gene: ", gene_id))
    } else {
      NULL
    }
  })

  output$multi_pe_no_data_msg <- shiny::renderUI({
    pair <- multi_pair()
    shiny::req(!is.null(pair))
    prot_id <- pair$protein_full
    expr <- if (exists("expr_raw",  envir = proteomics_env, inherits = FALSE)) {
      get("expr_raw",  envir = proteomics_env, inherits = FALSE)
    } else if (exists("expr_norm", envir = proteomics_env, inherits = FALSE)) {
      get("expr_norm", envir = proteomics_env, inherits = FALSE)
    } else NULL
    if (is.null(expr) || !prot_id %in% rownames(expr)) {
      tags$div(class = "alert alert-warning", style = "margin-top:8px; margin-bottom:4px;",
               paste0("No expression data available for protein: ", prot_id))
    } else {
      NULL
    }
  })

  # Gene/Protein ID labels above the DE statistics tables
  output$multi_ge_de_label <- shiny::renderUI({
    pair <- multi_pair()
    shiny::req(!is.null(pair))
    tags$span(pair$gene, style = "color:#18BC9C; font-size:16px; font-weight:700;")
  })
  output$multi_pe_de_label <- shiny::renderUI({
    pair <- multi_pair()
    shiny::req(!is.null(pair))
    tags$span(pair$protein_full, style = "color:#18BC9C; font-size:16px; font-weight:700;")
  })

  # Wire Gene Explorer for MultiOmics.
  # multi_pair() is the single source of truth — no secondary selectize input is used.
  wire_gene_explorer(
    session              = session,
    input                = input,
    output               = output,
    env                  = rnaseq_env,
    tabs_id              = "multi_tabs",
    click_sources        = c(),
    id_prefix            = "multi_ge",
    feature_label        = "Gene",
    override_id_reactive = reactive({
      pair <- multi_pair()
      req(!is.null(pair))
      pair$gene
    })
  )

  # Wire Protein Explorer for MultiOmics.
  # multi_pair() is the single source of truth — no secondary selectize input is used.
  wire_protein_explorer(
    session              = session,
    input                = input,
    output               = output,
    env                  = proteomics_env,
    tabs_id              = "multi_tabs",
    click_sources        = c(),
    id_prefix            = "multi_pe",
    feature_label        = "Protein",
    override_id_reactive = reactive({
      pair <- multi_pair()
      req(!is.null(pair))
      pair$protein_full
    })
  )

  # --- Shared "Color by" for MultiOmics ---
  # Populate the shared selector from the intersection of both envs' color/shape columns.
  observe({
    ge_color <- if (exists("color", envir = rnaseq_env, inherits = FALSE))
      get("color", envir = rnaseq_env, inherits = FALSE) else NULL
    ge_shape <- if (exists("shape", envir = rnaseq_env, inherits = FALSE))
      get("shape", envir = rnaseq_env, inherits = FALSE) else NULL

    pe_eff <- prot_get_effects_from_env(proteomics_env)

    # Build the union of grouping columns present in both envs
    ge_choices <- unique(na.omit(c(ge_color, ge_shape)))
    common     <- intersect(ge_choices, pe_eff)
    choices    <- if (length(common) > 0) common else ge_choices

    shiny::updateSelectInput(session, "multi_shared_color_by",
                             choices  = choices,
                             selected = choices[1])
  })

  # Sync: whenever the shared selector changes, push the value to both per-panel selectors.
  shiny::observeEvent(input$multi_shared_color_by, {
    val <- input$multi_shared_color_by
    shiny::updateSelectInput(session, "multi_ge_color_by", selected = val)
    shiny::updateSelectInput(session, "multi_pe_color_by", selected = val)
  }, ignoreNULL = TRUE)

  # --- Combined plot download ---
  output$download_multi_combined_plot <- shiny::downloadHandler(
    filename = function() {
      pair <- multi_pair()
      gene_id <- if (!is.null(pair)) sanitize_filename(pair$gene)       else "gene"
      prot_id <- if (!is.null(pair)) sanitize_filename(pair$protein_full) else "protein"
      paste0("MultiOmics_", gene_id, "_", prot_id, ".png")
    },
    content = function(file) {
      # Build gene plot
      ge_p <- tryCatch({
        gene_id <- multi_pair()$gene
        color_by  <- shiny::req(input$multi_shared_color_by)
        scale_sel <- shiny::req(input$multi_ge_scale)

        expr_raw    <- get("expr_raw",    envir = rnaseq_env, inherits = FALSE)
        sample_meta <- get("sample_meta", envir = rnaseq_env, inherits = FALSE)
        expr_norm   <- if (exists("expr_norm", envir = rnaseq_env, inherits = FALSE))
          get("expr_norm", envir = rnaseq_env, inherits = FALSE) else NULL

        df <- get_explorer_df_from_mats(
          feature_id = gene_id, scale_sel = scale_sel, color_by = color_by,
          expr_raw = expr_raw, expr_norm = expr_norm, sample_meta = sample_meta)

        make_explorer_plot(df, color_by = color_by, scale_sel = scale_sel,
                           ylab_log2 = get_norm_method_label_from_env(rnaseq_env),
                           theme_variant = "gene") +
          ggplot2::ggtitle(paste("Gene:", gene_id))
      }, error = function(e) NULL)

      # Build protein plot
      pe_p <- tryCatch({
        prot_id   <- multi_pair()$protein_full
        color_by  <- shiny::req(input$multi_shared_color_by)
        scale_sel <- shiny::req(input$multi_pe_scale)

        sample_meta <- get("sample_meta", envir = proteomics_env, inherits = FALSE)
        mat <- if (identical(scale_sel, "raw")) {
          get("expr_raw", envir = proteomics_env, inherits = FALSE)
        } else {
          get("expr_norm", envir = proteomics_env, inherits = FALSE)
        }

        df <- get_explorer_df_from_mat(
          feature_id = prot_id, mat = mat, color_by = color_by,
          sample_meta = sample_meta)

        ylab <- if (identical(scale_sel, "imp1")) "Log2 protein intensity (imputed)" else "Log2 protein intensity (raw)"
        make_explorer_plot(df, color_by = color_by, scale_sel = "log2",
                           ylab_log2 = ylab, theme_variant = "gene") +
          ggplot2::ggtitle(paste("Protein:", prot_id))
      }, error = function(e) NULL)

      # Combine side-by-side using gridExtra (grid is already loaded in global.R)
      combined <- if (!is.null(ge_p) && !is.null(pe_p)) {
        gridExtra::arrangeGrob(ge_p, pe_p, ncol = 2)
      } else if (!is.null(ge_p)) {
        ge_p
      } else if (!is.null(pe_p)) {
        pe_p
      } else {
        ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data available")
      }
      ggplot2::ggsave(filename = file, plot = combined, width = 16, height = 5, dpi = 300)
    }
  )

  # ---- Combined RNA + Protein Expression Data ------------------------------
  # Builds one wide table with an `Omics` column identifying RNA vs Protein.
  # Sample columns are the union of both omics; cells are NA where a sample
  # doesn't exist in the corresponding omics.
  #
  # We replicate the same data-fetch path the unused single tables use
  # (get_explorer_df_from_mat[s] -> build_explorer_wide_df) so styling and
  # numeric values match what the single Gene/Protein Explorers show.
  multi_combined_table_df <- reactive({
    pair <- multi_pair()
    req(!is.null(pair))

    color_by <- req(input$multi_shared_color_by)

    # RNA row(s)
    ge_wide <- tryCatch({
      gene_id   <- pair$gene
      scale_sel <- req(input$multi_ge_scale)
      expr_raw    <- get("expr_raw",    envir = rnaseq_env, inherits = FALSE)
      sample_meta <- get("sample_meta", envir = rnaseq_env, inherits = FALSE)
      expr_norm   <- if (exists("expr_norm", envir = rnaseq_env, inherits = FALSE))
        get("expr_norm", envir = rnaseq_env, inherits = FALSE) else NULL
      long_df <- get_explorer_df_from_mats(
        feature_id = gene_id, scale_sel = scale_sel, color_by = color_by,
        expr_raw = expr_raw, expr_norm = expr_norm, sample_meta = sample_meta)
      build_explorer_wide_df(long_df,
                             feature_id       = gene_id,
                             feature_col_name = "Feature")
    }, error = function(e) NULL)

    # Protein row(s)
    pe_wide <- tryCatch({
      prot_id   <- pair$protein_full
      scale_sel <- req(input$multi_pe_scale)
      sample_meta <- get("sample_meta", envir = proteomics_env, inherits = FALSE)
      mat <- if (identical(scale_sel, "raw")) {
        get("expr_raw", envir = proteomics_env, inherits = FALSE)
      } else {
        get("expr_norm", envir = proteomics_env, inherits = FALSE)
      }
      long_df <- get_explorer_df_from_mat(
        feature_id  = prot_id, mat = mat, color_by = color_by,
        sample_meta = sample_meta)
      build_explorer_wide_df(long_df,
                             feature_id       = prot_id,
                             feature_col_name = "Feature")
    }, error = function(e) NULL)

    rows <- list()
    if (!is.null(ge_wide) && nrow(ge_wide) > 0L) {
      rows[["RNA"]] <- cbind(Omics = "RNA",
                             ge_wide,
                             stringsAsFactors = FALSE)
    }
    if (!is.null(pe_wide) && nrow(pe_wide) > 0L) {
      rows[["Protein"]] <- cbind(Omics = "Protein",
                                 pe_wide,
                                 stringsAsFactors = FALSE)
    }
    if (length(rows) == 0L) return(NULL)

    # Align columns across both omics: Omics + Feature first, then the union
    # of sample columns. Missing samples in one omics become NA.
    other <- unique(unlist(lapply(rows, function(r)
      setdiff(colnames(r), c("Omics", "Feature")))))
    all_cols <- c("Omics", "Feature", other)

    combined <- do.call(rbind, lapply(rows, function(r) {
      miss <- setdiff(all_cols, colnames(r))
      for (cc in miss) r[[cc]] <- NA_real_
      r[, all_cols, drop = FALSE]
    }))
    rownames(combined) <- NULL
    combined
  })

  output$multi_combined_data_table <- DT::renderDT({
    # Reuse the shared renderer: setdiff treats c("Omics","Feature") as the
    # non-numeric leading columns, so dynamic decimals + DT options match
    # every other Expression Data table in the app.
    render_explorer_wide_dt(multi_combined_table_df(),
                            feature_col_name = c("Omics", "Feature"))
  })

  output$download_multi_combined_data <- downloadHandler(
    filename = function() paste0("multiomics_expression_data_", Sys.Date(), ".csv"),
    content  = function(file) {
      df <- multi_combined_table_df()
      utils::write.csv(if (is.null(df)) data.frame() else df,
                       file = file, row.names = FALSE)
    }
  )

  # ---- Gene–Protein Correlation (MultiOmics Explorer) ----------------------
  # Isolated to this tab. Uses the normalized matrices (rnaseq_env$expr_norm /
  # proteomics_env$expr_norm) for the currently selected pair, paired across
  # the samples shared by both omics. Methods are detected dynamically and
  # every cor.test is guarded — a failed calculation never crashes the app.

  # Paired normalized gene vs protein expression (complete pairs only).
  multi_corr_df <- reactive({
    pair <- multi_pair()
    req(!is.null(pair))
    gene_mat <- if (exists("expr_norm", envir = rnaseq_env, inherits = FALSE))
      get("expr_norm", envir = rnaseq_env, inherits = FALSE) else NULL
    prot_mat <- if (exists("expr_norm", envir = proteomics_env, inherits = FALSE))
      get("expr_norm", envir = proteomics_env, inherits = FALSE) else NULL
    get_gene_protein_paired_df(pair$gene, pair$protein_full, gene_mat, prot_mat)
  })

  # Correlation methods that can actually be computed for the current pair.
  multi_corr_methods <- reactive({
    gene_protein_corr_valid_methods(multi_corr_df())
  })

  # conditionalPanel flag: "yes" when at least one method is computable.
  output$multi_corr_available <- shiny::reactive({
    if (length(multi_corr_methods()) > 0) "yes" else "no"
  })
  shiny::outputOptions(output, "multi_corr_available", suspendWhenHidden = FALSE)

  # Interactive Plotly scatter (X = gene, Y = protein). The plot is
  # method-agnostic (no dropdown) — switching method never changed the points,
  # only the reported statistics, which now live in the table below.
  output$multi_corr_plot <- plotly::renderPlotly({
    df <- multi_corr_df()
    req(!is.null(df))
    req(length(multi_corr_methods()) > 0)
    pair <- multi_pair()
    plotly_gene_protein_correlation(df, pair$gene, pair$protein_full)
  })

  # Statistics table: one row per valid method (coefficient / p-value /
  # paired-sample count). Built once per pair via the shared validation path.
  output$multi_corr_stats_table <- DT::renderDT({
    df <- multi_corr_df()
    req(!is.null(df))
    stats_df <- gene_protein_corr_stats_df(df)
    req(!is.null(stats_df))
    DT::datatable(
      stats_df,
      rownames = FALSE,
      class    = "compact stripe hover",
      options  = list(dom = "t", paging = FALSE, searching = FALSE,
                      ordering = FALSE, info = FALSE,
                      columnDefs = list(list(className = "dt-center",
                                             targets = "_all")))
    )
  })

  # Friendly message shown when no method can be calculated.
  output$multi_corr_msg <- shiny::renderUI({
    pair <- multi_pair()
    req(!is.null(pair))
    if (length(multi_corr_methods()) > 0) return(NULL)
    tags$div(
      class = "alert alert-info",
      style = "margin-top:8px;",
      paste0("Correlation cannot be calculated for the selected gene–",
             "protein pair because there are not enough valid paired samples ",
             "or the expression values are constant.")
    )
  })

}


# ---------- Appstream authentication ----------
# On Shiny Server, the application (app_server above) starts only after the Appstream check.
server <- function(input, output, session) {

  # Local development: no Appstream check
  if (!APPSTREAM_REQUIRED) {
    app_server(input, output, session)
    return(invisible())
  }

  # Refuse only this session (instead of shiny::stopApp(), which stops the app for all users)
  deny <- function(reason) {
    message("Appstream: access denied to ", APPSTREAM_APP_PATH, " (", reason, ")")
    session$sendCustomMessage("appstream-denied", list(login_url = APPSTREAM_LOGIN_URL))
  }

  shiny::observeEvent(input$cookies, {

    cookies <- input$cookies

    required_cookies <- c(
      "Appstream_user",
      "Appstream_hash",
      "Appstream_md5"
    )

    if (!all(required_cookies %in% names(cookies))) {
      return(deny("missing Appstream cookies"))
    }

    user  <- as.character(cookies$Appstream_user)
    token <- as.character(cookies$Appstream_hash)
    md5   <- as.character(cookies$Appstream_md5)

    if (length(user) != 1 || !grepl("^[A-Za-z0-9._-]+$", user) || user %in% c(".", "..")) {
      return(deny("invalid user name"))
    }

    if (is.na(APPSTREAM_APP_PATH)) {
      return(deny("application path not detected"))
    }

    user_path <- file.path(APPSTREAM_USERS_DIR, user)
    yaml_path <- file.path(user_path, "Order.yaml")

    if (!file.exists(yaml_path)) {
      return(deny(paste("no Order.yaml for", user)))
    }

    yaml_lines <- base::readLines(
      yaml_path,
      warn = FALSE
    )

    app_url <- paste0("/", APPSTREAM_APP_PATH, "/")

    if (!any(grepl(app_url, yaml_lines, fixed = TRUE))) {
      return(deny(paste(user, "has no permission")))
    }

    token_file <- NA_character_

    user_files <- list.files(
      path = user_path,
      full.names = TRUE
    )
    user_files <- user_files[basename(user_files) != "Order.yaml"]

    for (file_path in user_files) {
      file_md5 <- as.character(
        tools::md5sum(file_path)
      )

      if (identical(file_md5, md5)) {
        token_file <- basename(file_path)
        break
      }
    }

    if (is.na(token_file)) {
      return(deny(paste("no matching token for", user)))
    }

    if (!isTRUE(tryCatch(bcrypt::checkpw(token_file, token), error = function(e) FALSE))) {
      return(deny(paste("invalid token for", user)))
    }

    for (name in required_cookies) {
      session$sendCustomMessage("cookie-remove", list(name = name))
    }

    message("Appstream: ", user, " authenticated for ", APPSTREAM_APP_PATH)

    # The unchanged bdika application starts here
    app_server(input, output, session)

  }, ignoreNULL = TRUE, ignoreInit = FALSE, once = TRUE)
}
