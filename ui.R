
ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(
    tags$link(
      rel = "stylesheet",
      href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.0/css/all.min.css"
    ),
    tags$style(HTML(
      ".plot-title{
     background:#EDF0F1;
     color:#2c3e50;
     padding:8px 12px;
     border-radius:10px;
     text-align:center;
     margin:10px 0 15px;
     font-weight:600
   }
   .plot-title.full{display:block;width:100%}

/* shrink the numericInput boxes themselves */
input#gene_explorer_padj_cutoff,
input#gene_explorer_fc_cutoff {
  height: 26px !important;
  padding: 2px 4px !important;
  font-size: 12px !important;
}

/* linearFC cutoff label — match the significance-dropdown label style
   (same font size / weight / color / spacing) across all explorer tabs */
label[for$='_fc_cutoff'] {
  font-size: 13px !important;
  font-weight: 600 !important;
  color: #2c3e50 !important;
  margin-bottom: 2px !important;
}

/* DE cutoffs: make the sig-mode <select> look like a compact label with a ▼ */
.sig-mode-dropdown .form-group { margin-bottom: 2px !important; }
.sig-mode-dropdown select.form-control {
  border: none !important;
  background: transparent !important;
  padding: 0 2px !important;
  height: auto !important;
  min-height: 20px !important;
  font-size: 13px !important;
  font-weight: 600 !important;
  color: #2c3e50 !important;
  box-shadow: none !important;
  cursor: pointer !important;
  width: auto !important;
}
.sig-mode-dropdown select.form-control:focus {
  outline: none !important;
  box-shadow: none !important;
}

/* Multi-Gene Explorer: removable gene tags */
.mge-tag {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  background: #18BC9C;
  color: white;
  border-radius: 12px;
  padding: 3px 10px;
  font-size: 12px;
  font-weight: 600;
  user-select: none;
}
.mge-tag .mge-tag-x {
  cursor: pointer;
  opacity: 0.85;
  font-weight: 700;
  font-size: 14px;
  line-height: 1;
}
.mge-tag .mge-tag-x:hover {
  opacity: 1;
  color: #ffe;
}"
    )),

    # Hover tooltip for the significance-cutoff mode dropdown.
    # Uses the browser's native `title` attribute — lightweight, no extra libs.
    # The tooltip text reflects the OTHER mode the user would switch TO.
    tags$script(HTML("
      (function() {
        function sigTipFor(val) {
          return val === 'pvalue'
            ? 'Switch to adjusted p-value'
            : 'Switch to p-value';
        }
        function applyTip(sel) {
          var tip = sigTipFor(sel.value);
          sel.setAttribute('title', tip);
          var wrap = sel.closest('.sig-mode-dropdown');
          if (wrap) wrap.setAttribute('title', tip);
        }
        function init() {
          document.querySelectorAll('.sig-mode-dropdown select').forEach(function(sel) {
            if (sel.dataset.sigTipBound === '1') return;
            sel.dataset.sigTipBound = '1';
            applyTip(sel);
            sel.addEventListener('change', function() { applyTip(sel); });
          });
        }
        // Run on initial load and whenever Shiny renders new outputs
        document.addEventListener('DOMContentLoaded', init);
        if (window.$) { $(document).on('shiny:value shiny:inputchanged', init); }
        // Fallback: poll briefly in case of late render
        var tries = 0; var t = setInterval(function() {
          init(); if (++tries > 20) clearInterval(t);
        }, 300);
      })();
    ")),

    # Authentication cookies
    tags$script(HTML("
  (function() {

    function readCookies() {
      var cookies = {};

      document.cookie.split(';').forEach(function(cookie) {
        var parts = cookie.split('=');
        var name = parts.shift().trim();

        if (!name) return;

        var value = parts.join('=');

        try {
          value = decodeURIComponent(value);
        } catch (e) {
          // Keep the original value if decoding fails.
        }

        cookies[name] = value;
      });

      return cookies;
    }

    $(document).on('shiny:connected', function() {
      Shiny.setInputValue(
        'cookies',
        readCookies(),
        {priority: 'event'}
      );
    });

    Shiny.addCustomMessageHandler('cookie-remove', function(message) {
      document.cookie =
        message.name +
        '=; Max-Age=0; expires=Thu, 01 Jan 1970 00:00:00 GMT; path=/';
    });

    // Failed Appstream check: replace this browser's page with a login message
    Shiny.addCustomMessageHandler('appstream-denied', function(message) {
      document.body.innerHTML =
        '<div style=\"max-width:600px; margin:80px auto; text-align:center;\">' +
        '<h3>Access denied</h3>' +
        '<p>Please log in through <a href=\"' + message.login_url + '\">Appstream</a>' +
        ' and open this application from there.</p></div>';
    });

  })();
"))
  ),
  
  navbarPage(
    "",
    id = "main_navbar",

    # ---------- tab panel 1 - Home ----------
    tabPanel(
      "Home",
      fluidPage(
        # Hero image with overlay text
        tags$div(
          style = "width: 100%; height: 250px; position: relative; overflow: hidden; margin-bottom: 30px;",
          tags$img(
            src   = "home_picture4.png",
            style = "width: 100%; height: 100%; object-fit: cover; display: block;"
          ),
          tags$div(
            "Multiomics Explorer",
            style = paste(
              "position: absolute; top: 50%; left: 65%;",
              "transform: translate(-50%, -50%);",
              "color: #2C3E50; font-size: 5vw; font-weight: bold;",
              "text-shadow: 2px 2px 6px black; text-align: center;",
              "white-space: nowrap;"
            )
          )
        ),
        
        # Researcher-specific content (researcher, project, description, team):
        # inserted here by app_server() only after a successful Appstream check.
        # See home_info_ui() in global.R.
        tags$div(id = "home_info_ui")
      )
    ),
    
    # ---------- tab panel 2 - Experiment design ----------

    tabPanel(
      "Experiment design",
      # Experiment design table (sample annotations)
      h4("Experiment design", class = "plot-title full"),
      fluidRow(
        column(
          width = 8,
          DTOutput("col_table")
        )
      )
    ),
    
    # ---------- tab panel 3 - Data & Files ----------
    
    tabPanel(
      "Data & Files",
      div(
        h4("Data & Files", class = "plot-title full"),
        tags$hr(),
        uiOutput("doc_page")
      )
    ),
    

    # ---------- tab panel 4 - RNAseq ----------

    # Shown only when data/shiny_payload_rnaseq.rds exists (SHOW_TAB in global.R)
    if (SHOW_TAB$rnaseq) tabPanel(
      "RNAseq",
      titlePanel("RNA-seq Analysis"),
      bslib::navset_pill_list(
        id     = "rna_tabs",
        widths = c(2, 9),
        
        
        ## ---- Subtab: Exploratory Analyses ----
        
        bslib::nav_panel(
          "Exploratory Analyses",
          bslib::navset_tab(
            id = "explore_subtab",
            
            ## ---- Subtab: Data structure ----
            bslib::nav_panel(
              "Data structure",
              
              # 3D PCA - iframe output
              h4("3D PCA", class = "plot-title full"),
              shinycssloaders::withSpinner(plotly::plotlyOutput("pca_3d_ui", height = "600px"), type = 4),
              br(), hr(),
              
              # 2D PCA - PC1 vs PC2
              h4("2D PCA", class = "plot-title full"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  h4("PC1 vs PC2"),
                  shinycssloaders::withSpinner(
                    plotlyOutput("pca_pc1_pc2_plot", width = "100%"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),
              br(),
              
              # 2D PCA - PC1 vs PC3
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  h4("PC1 vs PC3"),
                  shinycssloaders::withSpinner(
                    plotlyOutput("pca_pc1_pc3_plot", width = "100%"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),
              
              br(), hr(),
              
              # Sample distance heatmap (interactive)
              h4("Sample Distance Heatmap", class = "plot-title full"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotlyOutput("sample_corr_heatmap", height = "700px"),
                    type = 4
                  )
                ),
                column(width = 1)
              )
            ),
            
            ## ---- Subtab: DE Genes ----

            bslib::nav_panel(
              "DE Genes",
              fluidRow(
                column(
                  width = 3,
                  tags$h4("Choose contrast:"),
                  selectInput("de_contrast", label = NULL, choices = NULL)
                ),
                column(
                  width = 9,
                  # Volcano plot (DE analysis)
                  h4("Volcano plot", class = "plot-title full"),
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("volcano_plot", height = "520px"),
                    type = 4
                  ),
                  tags$hr(),
                  # MA plot (DE analysis)
                  h4("MA plot", class = "plot-title full"),
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("ma_plot", height = "520px"),
                    type = 4
                  )
                )
              )
            )
          )
        ),
        
        ## ---- Subtab: Stats results ----

        bslib::nav_panel(
          "Stats results",
          # Summary DE statistics table
          h4("Summary DE statistics", class = "plot-title full"),
          tableOutput("DE_genes_stats"),
          tags$hr(),
          # Download result tables (Excel)
          h4("Download result tables", class = "plot-title full"),
          uiOutput("results_buttons"),
          tags$hr(),
          # Download normalized count matrices
          h4("Download normalized counts matrix", class = "plot-title full"),
          div(
            style = "margin-bottom: 8px;",
            downloadButton(
              "dl_norm_all_txt",
              "Normalized counts matrix - All genes",
              class = "btn btn-success btn-sm",
              style = "margin:4px;"
            )
          ),
          div(
            downloadButton(
              "dl_norm_de_txt",
              "Normalized counts matrix - DE genes",
              class = "btn btn-success btn-sm",
              style = "margin:4px;"
            )
          )
        ),
        
        ## ---- Subtab: Clustering ----
        bslib::nav_panel(
          "Clustering",
          bslib::navset_tab(
            id = "clustering_subtab",
            
            ## ---- Subtab: Hierarchical Clustering ----
            bslib::nav_panel(
              "Hierarchical Clustering",
              
              # Static figure from env (gtable)
              h4("Hierarchical clustering heatmap", class = "plot-title full"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotOutput("hc_static_fig_rna", height = "700px"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),

              tags$h4(
                "Interactive Hierarchical Clustering",
                class = "plot-title full",
                style = "text-align:center;margin-top:15px;margin-bottom:10px;"
              ),
              tags$hr(),

              # Interactive heatmap
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("de_heatmap", height = "700px"),
                    type = 4, color = "#2c3e50"
                  )
                ),
                column(width = 1)
              ),

              # ============================================================
              # HEATMAP-ZOOM MULTI-SELECT (revertible) — Hierarchical only
              # Remove this fluidRow block to revert.
              # ============================================================
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  div(
                    style = paste("display:flex; gap:14px; align-items:center;",
                                  "flex-wrap:wrap; padding:8px 0; margin-top:4px;",
                                  "border-top:1px solid #eee;"),
                    uiOutput("hc_zoom_msg_ui"),
                    actionButton(
                      "hc_zoom_open_btn",
                      "Open visible genes in Multi-Gene Explorer",
                      icon  = icon("up-right-from-square"),
                      class = "btn btn-primary btn-sm"
                    )
                  )
                ),
                column(width = 1)
              )
              # ============================================================
              # END HEATMAP-ZOOM MULTI-SELECT
              # ============================================================
            ),

            ## ---- Subtab: Partition Clustering ----
            bslib::nav_panel(
              "Partition Clustering",
              
              # Partition clustering title (dynamic)
              uiOutput("pc_title_ui"),

              # Static figure from env (gtable)
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotOutput("pc_static_fig_rna", height = "700px"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),
              tags$hr(),
              
              # Cluster selector + plots
              selectInput("pc_cluster_pick", "Choose cluster", choices = NULL, width = "300px"),
              uiOutput("pc_cluster_title_ui"),
              uiOutput("pc_cluster_note_ui"),
              
              fluidRow(
                column(
                  width = 9,
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("pc_heatmap_cluster", height = "700px"),
                    type = 4, color = "#2c3e50"
                  )
                ),
                column(
                  width = 3,
                  shinycssloaders::withSpinner(
                    plotOutput("pc_cluster_subgraph", height = "550px"),
                    type = 4, color = "#2c3e50"
                  )
                )
              ),

              # ============================================================
              # HEATMAP-ZOOM MULTI-SELECT (revertible) — Partition clustering
              # Remove this fluidRow block to revert.
              # ============================================================
              fluidRow(
                column(
                  width = 9,
                  div(
                    style = paste("display:flex; gap:14px; align-items:center;",
                                  "flex-wrap:wrap; padding:8px 0; margin-top:4px;",
                                  "border-top:1px solid #eee;"),
                    uiOutput("pc_zoom_msg_ui"),
                    actionButton(
                      "pc_zoom_open_btn",
                      "Open visible genes in Multi-Gene Explorer",
                      icon  = icon("up-right-from-square"),
                      class = "btn btn-primary btn-sm"
                    )
                  )
                )
              ),
              # ============================================================
              # END HEATMAP-ZOOM MULTI-SELECT (Partition)
              # ============================================================

              uiOutput("pc_results_title_ui"),
              div(
                uiOutput("dl_results_de_filtered_ui"),
                style = "margin-bottom: 8px;"
              )
            ),
            
            ## ---- Subtab: Binary Pattern Clustering ----
            bslib::nav_panel(
              "Binary Pattern",
              fluidRow(
                column(4, uiOutput("bp_pattern_ui"))
              ),
              tags$hr(),
              
              uiOutput("bp_static_title_ui"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  plotOutput("bp_static", height = "700px")
                ),
                column(width = 1)
              ),
              tags$hr(),

              uiOutput("bp_interactive_title_ui"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("bp_interactive", height = "700px"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),

              # ============================================================
              # HEATMAP-ZOOM MULTI-SELECT (revertible) — Binary Pattern
              # Remove this fluidRow block to revert.
              # ============================================================
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  div(
                    style = paste("display:flex; gap:14px; align-items:center;",
                                  "flex-wrap:wrap; padding:8px 0; margin-top:4px;",
                                  "border-top:1px solid #eee;"),
                    uiOutput("bp_zoom_msg_ui"),
                    actionButton(
                      "bp_zoom_open_btn",
                      "Open visible genes in Multi-Gene Explorer",
                      icon  = icon("up-right-from-square"),
                      class = "btn btn-primary btn-sm"
                    )
                  )
                ),
                column(width = 1)
              ),
              # ============================================================
              # END HEATMAP-ZOOM MULTI-SELECT (Binary Pattern)
              # ============================================================

              tags$hr(),

              uiOutput("bp_results_title_ui"),
              div(
                uiOutput("dl_results_de_bp_filtered_ui"),
                style = "margin-bottom: 8px;"
              )
            )
          )
        ),
        
        ## ---- Subtab: Enrichment (ORA/GSEA — Stage 4) ----
        ## Replaces the old GO/KEGG "Coming soon" placeholder with the
        ## ORA/GSEA-first structure. Markup + all logic live in functions_shiny.R
        ## (enrichment_panel_ui / wire_enrichment); server side in server.R.
        enrichment_panel_ui("rna_enrich"),
        
        ## ---- Subtab: Keyword Search ----
        bslib::nav_panel(
          "Keyword Search",
          keyword_search_panel_ui(id_prefix = "rna_kws")
        ),

        ## ---- Subtab: Gene Explorer  ----

        bslib::nav_panel(
          "Gene Explorer",
          fluidRow(
            column(
              width = 12,
              
              # Label is rendered ABOVE the flex row (not inside selectizeInput)
              # so the flex row contains only single-line elements: the input
              # box and the buttons. With label-free form-group + zeroed
              # bottom margin + align-items:center, the input and buttons sit
              # on the same horizontal line with no vertical offset.
              tags$label("Select Gene:", `for` = "gene_explorer_gene",
                         style = "font-size: 18px; font-weight: 600; display:block; margin-bottom:5px;"),
              div(
                style = "display:flex; gap:8px; align-items:center; flex-wrap:wrap;",
                htmltools::tagAppendAttributes(
                  selectizeInput(
                    inputId  = "gene_explorer_gene",
                    label    = NULL,
                    choices  = NULL,
                    width    = "500px",
                    multiple = FALSE,
                    options  = list(
                      placeholder = "Type or paste a Gene/Transcript ID...",
                      maxOptions  = 20
                    )
                  ),
                  style = "margin-bottom: 0;"
                ),
                actionButton("gene_explorer_clear_btn", "Clear selection",
                             icon  = icon("trash"),
                             class = "btn btn-default btn-sm")
              ),
              
              br(),
              
              uiOutput("gene_explorer_selected_gene_title"),
              uiOutput("gene_explorer_plot_title"),
              
              tags$hr(style = "margin-top: 5px; margin-bottom: 10px; border-color: #ddd;"),
              
              conditionalPanel(
                condition = "input.gene_explorer_gene != null && input.gene_explorer_gene != ''",
                div(
                  style = "display: flex; justify-content: flex-start; gap: 12px; align-items: center; margin-bottom: 10px;",
                  
                  selectInput(
                    inputId  = "gene_explorer_scale",
                    label    = "Normalization Method:",
                    choices  = NULL,
                    width    = "280px"
                  ),
                  
                  selectInput(
                    inputId  = "gene_explorer_color_by",
                    label    = "Color bars by:",
                    choices  = NULL,
                    width    = "220px"
                  )
                )
              ),
              
              br(),
              
              plotOutput("gene_explorer_counts_plot", height = "400px"),

              br(),

              conditionalPanel(
                condition = "input.gene_explorer_gene != null && input.gene_explorer_gene != ''",
                div(
                  style = "display: flex; justify-content: flex-end; gap: 12px; align-items: center;",
                  downloadButton("download_gene_explorer_plot", "Download Plot")
                )
              ),

              br(),

              # Expression Data — wide-format table for the single selected gene
              conditionalPanel(
                condition = "input.gene_explorer_gene != null && input.gene_explorer_gene != ''",
                div(
                  h4("Expression Data", class = "plot-title full"),
                  DT::DTOutput("gene_explorer_data_table"),
                  br(),
                  div(
                    style = "display: flex; justify-content: flex-end;",
                    downloadButton("download_gene_explorer_data", "Download Data")
                  )
                )
              ),

              br(),

              conditionalPanel(
                condition = "input.gene_explorer_gene != null && input.gene_explorer_gene != ''",
                div(
                  h4("Differential Expression Statistics", class = "plot-title full"),
                  
                  DT::dataTableOutput("gene_explorer_de_table"),
                  
                  div(
                    style = paste(
                      "display:flex;",
                      "justify-content:space-between;",
                      "align-items:flex-start;",
                      "margin-top:10px;",
                      "gap:16px;"
                    ),
                    
                    div(
                      style = "font-size:14px; line-height:1.6;",
                      uiOutput("gene_explorer_de_legend_ui")
                    ),
                    
                    div(
                      style = "display:flex; gap:10px; align-items:flex-end; flex:0 0 auto; min-width:280px;",
                      uiOutput("gene_explorer_cutoffs_ui")
                    )
                  ),
                  
                  tags$hr(),
                  
                  uiOutput("gene_explorer_annot_title_ui"),
                  DT::DTOutput("gene_explorer_annot_dt"),
                  uiOutput("gene_explorer_annot_dl_ui"),

                  br(),

                  # Only when the MultiOmics tab exists (its target)
                  if (SHOW_TAB$multiomics) conditionalPanel(
                    condition = "input.gene_explorer_gene != null && input.gene_explorer_gene != ''",
                    div(
                      style = "margin-top: 10px; margin-bottom: 20px;",
                      actionButton("gene_explorer_multi_btn", "View in MultiOmics Explorer",
                                   class = "btn btn-info btn-sm")
                    )
                  ),

                  br()
                )
              )
            )
          )
        ),

        ## ---- Subtab: Multi-Gene Explorer (RNAseq only) ----

        bslib::nav_panel(
          "Multi-Gene Explorer",
          fluidRow(
            column(
              width = 12,

              h4("Multi-Gene Explorer", class = "plot-title full"),

              # Label rendered ABOVE the flex row so the input box and the
              # action buttons can sit on a single line with center alignment.
              tags$label("Select genes (up to 100):", `for` = "mge_genes_input",
                         style = "font-size:14px; font-weight:600; display:block; margin-top:8px; margin-bottom:5px;"),
              div(
                style = paste("display:flex; gap:12px; align-items:center;",
                              "flex-wrap:wrap;"),
                div(
                  style = "flex:1 1 360px; min-width:300px;",
                  htmltools::tagAppendAttributes(
                    selectizeInput(
                      inputId  = "mge_genes_input",
                      label    = NULL,
                      choices  = NULL,
                      multiple = TRUE,
                      width    = "100%",
                      options  = list(
                        placeholder = "Type or paste Gene/Transcript IDs...",
                        maxOptions  = 50
                      )
                    ),
                    style = "margin-bottom: 0;"
                  )
                ),
                div(
                  style = "display:flex; gap:8px;",
                  actionButton("mge_add_btn",   "Add selected genes",
                               icon = icon("plus"),
                               class = "btn btn-primary btn-sm"),
                  actionButton("mge_clear_btn", "Clear selection",
                               icon = icon("trash"),
                               class = "btn btn-default btn-sm"),
                  # Compact secondary action — opens a modal with the textarea
                  actionButton("mge_paste_open", "Paste feature list",
                               icon  = icon("clipboard-list"),
                               class = "btn btn-default btn-sm")
                )
              ),

              uiOutput("mge_count_ui"),
              uiOutput("mge_summary_ui"),

              # Empty-state gate: hide the Expression Profile / data / DE
              # sections until at least one gene is selected.
              conditionalPanel(
                condition = "output.mge_has_selection == 'yes'",

                br(),

                # 1. Selected-features hierarchical heatmap — shows ALL
                #    selected features, not the paginated subset.
                uiOutput("mge_selected_heatmap_title"),
                uiOutput("mge_selected_heatmap_msg"),
                plotly::plotlyOutput("mge_selected_heatmap", height = "520px"),

                br(),

                # 2. Pagination controls (text + Previous/Next). Drives the
                #    Expression Profile / Data / DE Statistics sections below.
                uiOutput("mge_pagination_ui"),

                # 3. Expression Profile (current page only)
                h4("Expression Profile", class = "plot-title full"),

                div(
                  style = "display:flex; gap:12px; align-items:flex-end; flex-wrap:wrap; margin-bottom:8px;",
                  selectInput("mge_scale",    label = "Normalization Method:",
                              choices = NULL, width = "260px"),
                  selectInput("mge_color_by", label = "Color bars by:",
                              choices = NULL, width = "200px"),
                  selectInput("mge_sort_by",  label = "Sort by:",
                              choices  = c("Selection order"                    = "selection",
                                           "Feature name A–Z"              = "alpha",
                                           "Significance in selected contrast"  = "significance",
                                           "linearFC in selected contrast"      = "linearfc"),
                              selected = "selection",
                              width    = "260px"),
                  uiOutput("mge_sort_contrast_ui")
                ),

                div(
                  style = "overflow:auto; max-height:1200px; border:1px solid #eee; padding:6px; border-radius:6px;",
                  plotOutput("mge_plot", height = "auto")
                ),

                br(),

                div(
                  style = "display:flex; justify-content:flex-end;",
                  downloadButton("download_mge_plot", "Download Plot")
                ),

                br(),

                # 4. Expression Data (current page only)
                h4("Expression Data", class = "plot-title full"),

                DT::DTOutput("mge_data_table"),

                br(),

                div(
                  style = "display:flex; justify-content:flex-end;",
                  downloadButton("download_mge_table", "Download Data")
                ),

                br(),

                # 5. Differential Expression Statistics (current page only)
                h4("Differential Expression Statistics", class = "plot-title full"),

                uiOutput("mge_de_contrast_ui"),

                DT::dataTableOutput("mge_de_table"),

                div(
                  style = paste("display:flex;", "justify-content:space-between;",
                                "align-items:flex-start;", "margin-top:10px;",
                                "gap:16px;"),
                  div(
                    style = "font-size:14px; line-height:1.6;",
                    uiOutput("mge_de_legend_ui")
                  ),
                  div(
                    style = "display:flex; gap:10px; align-items:flex-end; flex:0 0 auto; min-width:280px;",
                    uiOutput("mge_de_cutoffs_ui")
                  )
                ),

                br(),

                # 6. Annotation (full current selection)
                uiOutput("mge_annot_title_ui"),
                DT::DTOutput("mge_annot_dt"),
                uiOutput("mge_annot_dl_ui"),

                br()
              )
            )
          )
        )

      )  # <-- סוגר navset_pill_list של RNAseq
    ),   # <-- סוגר tabPanel("RNAseq", ...)
    
    
    # ---------- tab panel 5 - Proteomics  ----------
    
    # Shown only when data/shiny_payload_proteomics.rds exists (SHOW_TAB in global.R)
    if (SHOW_TAB$proteomics) tabPanel(
      "Proteomics",
      titlePanel("Proteomics Analysis"),
      bslib::navset_pill_list(
        id     = "prot_tabs",
        widths = c(2, 9),
        
        ## ---- Subtab: Exploratory Analyses ----
        bslib::nav_panel(
          "Exploratory Analyses",
          bslib::navset_tab(
            id = "prot_explore_subtab",
            
            ## ---- Subtab: Data distributions ----
            bslib::nav_panel(
              "Data distributions",
              
              h4("Observed vs imputed distributions", class = "plot-title full"),
              shinycssloaders::withSpinner(
                plotOutput("imp_hist_samp_plot", height = "500px"),
                type = 4
              ),
              
              h4("Normalized expression boxplots", class = "plot-title full"),
              shiny::uiOutput("prot_norm_boxplot_controls"),
              shinycssloaders::withSpinner(
                plotOutput("prot_norm_boxplot", height = "550px"),
                type = 4
              ),
              
              # Section 3
              h4("Density plot of normalized intensities", class = "plot-title full"),
              shinycssloaders::withSpinner(
                plotly::plotlyOutput("prot_density_plot", height = "550px"),
                type = 4
              ),
            ),
            
            
            ## ---- Subtab: Data structure ----
            bslib::nav_panel(
              "Data structure",
              
              
              # 3D PCA - iframe output
              h4("3D PCA", class = "plot-title full"),
              shinycssloaders::withSpinner(plotly::plotlyOutput("prot_pca_3d_ui", height = "600px"), type = 4),
              br(), hr(),
              
              
              # 2D PCA - PC1 vs PC2 / PC1 vs PC3
              h4("2D PCA", class = "plot-title full"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  h4("PC1 vs PC2"),
                  shinycssloaders::withSpinner(
                    plotlyOutput("prot_pca_pc1_pc2_plot", width = "100%"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),
              br(),
              
              
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  h4("PC1 vs PC3"),
                  shinycssloaders::withSpinner(
                    plotlyOutput("prot_pca_pc1_pc3_plot", width = "100%"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),
              
              
              br(), hr(),
              
              
              # Sample distance heatmap
              h4("Sample Distance Heatmap", class = "plot-title full"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotlyOutput("prot_sample_corr_heatmap", height = "700px"),
                    type = 4
                  )
                ),
                column(width = 1)
              )
            ),
            
            ## ---- Subtab: DE Protein ----
            bslib::nav_panel(
              "DE Protein",
              fluidRow(
                column(
                  width = 3,
                  tags$h4("Choose contrast:"),
                  selectInput("prot_de_contrast", label = NULL, choices = NULL)
                ),
                column(
                  width = 9,
                  # Volcano plot (DE analysis)
                  h4("Volcano plot", class = "plot-title full"),
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("prot_volcano_plot", height = "520px"),
                    type = 4
                  ),
                  tags$hr(),
                  # MA plot (DE analysis)
                  h4("MA plot", class = "plot-title full"),
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("prot_ma_plot", height = "520px"),
                    type = 4
                  )
                )
              )
            )
          )
        ),
        
        ## ---- Subtab: Stats results ----
        bslib::nav_panel(
          "Stats results",
          # Summary DE statistics table
          h4("Summary DE statistics", class = "plot-title full"),
          tableOutput("prot_DE_genes_stats"),
          tags$hr(),
          
          # Download result tables (Excel)
          h4("Download result tables", class = "plot-title full"),
          uiOutput("prot_results_buttons"),
          tags$hr(),
          
          # Download protein intensity matrices
          h4("Download protein intensity matrices", class = "plot-title full"),
          div(
            style = "margin-bottom: 8px;",
            downloadButton(
              "prot_dl_log2_raw_all_txt",
              "Log2 matrix - All proteins (raw, no imputation)",
              class = "btn btn-success btn-sm",
              style = "margin:4px;"
            )
          ),
          div(
            style = "margin-bottom: 8px;",
            downloadButton(
              "prot_dl_log2_imp1_all_txt",
              "Log2 matrix - All proteins (imputed)",
              class = "btn btn-success btn-sm",
              style = "margin:4px;"
            )
          ),
          div(
            downloadButton(
              "prot_dl_norm_de_txt",
              "Log2 matrix - DE proteins",
              class = "btn btn-success btn-sm",
              style = "margin:4px;"
            )
          )
        ),
        
        ## ---- Subtab: Clustering ----
        bslib::nav_panel(
          "Clustering",
          bslib::navset_tab(
            id = "prot_clustering_subtab",
            
            ## ---- Subtab: Hierarchical Clustering ----
            bslib::nav_panel(
              "Hierarchical Clustering",
              
              # Static figure from env (gtable)
              h4("Hierarchical clustering heatmap", class = "plot-title full"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotOutput("hc_static_fig_prot", height = "700px"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),

              tags$h4(
                "Interactive Hierarchical Clustering",
                class = "plot-title full",
                style = "text-align:center;margin-top:15px;margin-bottom:10px;"
              ),
              tags$hr(),

              # Interactive heatmap
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("prot_de_heatmap", height = "700px"),
                    type = 4, color = "#2c3e50"
                  )
                ),
                column(width = 1)
              ),

              # ============================================================
              # HEATMAP-ZOOM MULTI-SELECT (revertible) — Proteomics HC
              # ============================================================
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  div(
                    style = paste("display:flex; gap:14px; align-items:center;",
                                  "flex-wrap:wrap; padding:8px 0; margin-top:4px;",
                                  "border-top:1px solid #eee;"),
                    uiOutput("prot_hc_zoom_msg_ui"),
                    actionButton(
                      "prot_hc_zoom_open_btn",
                      "Open visible proteins in Multi-Protein Explorer",
                      icon  = icon("up-right-from-square"),
                      class = "btn btn-primary btn-sm"
                    )
                  )
                ),
                column(width = 1)
              )
              # ============================================================
              # END HEATMAP-ZOOM (Proteomics HC)
              # ============================================================
            ),

            ## ---- Subtab: Partition Clustering ----
            bslib::nav_panel(
              "Partition Clustering",
              
              uiOutput("prot_pc_title_ui"),

              # Static figure from env (gtable)
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotOutput("pc_static_fig_prot", height = "700px"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),
              tags$hr(),
              
              selectInput("prot_pc_cluster_pick", "Choose cluster", choices = NULL, width = "300px"),
              uiOutput("prot_pc_cluster_title_ui"),
              uiOutput("prot_pc_cluster_note_ui"),
              
              fluidRow(
                column(
                  width = 9,
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("prot_pc_heatmap_cluster", height = "700px"),
                    type = 4, color = "#2c3e50"
                  )
                ),
                column(
                  width = 3,
                  shinycssloaders::withSpinner(
                    plotOutput("prot_pc_cluster_subgraph", height = "550px"),
                    type = 4, color = "#2c3e50"
                  )
                )
              ),

              # ============================================================
              # HEATMAP-ZOOM MULTI-SELECT (revertible) — Proteomics Partition
              # ============================================================
              fluidRow(
                column(
                  width = 9,
                  div(
                    style = paste("display:flex; gap:14px; align-items:center;",
                                  "flex-wrap:wrap; padding:8px 0; margin-top:4px;",
                                  "border-top:1px solid #eee;"),
                    uiOutput("prot_pc_zoom_msg_ui"),
                    actionButton(
                      "prot_pc_zoom_open_btn",
                      "Open visible proteins in Multi-Protein Explorer",
                      icon  = icon("up-right-from-square"),
                      class = "btn btn-primary btn-sm"
                    )
                  )
                )
              ),
              # ============================================================
              # END HEATMAP-ZOOM (Proteomics Partition)
              # ============================================================

              uiOutput("prot_pc_results_title_ui"),
              div(
                uiOutput("prot_dl_results_de_filtered_ui"),
                style = "margin-bottom: 8px;"
              )
            ),
            
            ## ---- Subtab: Binary Pattern ----
            bslib::nav_panel(
              "Binary Pattern",
              fluidRow(
                column(4, uiOutput("prot_bp_pattern_ui"))
              ),
              tags$hr(),
              
              uiOutput("prot_bp_static_title_ui"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  plotOutput("prot_bp_static", height = "700px")
                ),
                column(width = 1)
              ),
              tags$hr(),

              uiOutput("prot_bp_interactive_title_ui"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("prot_bp_interactive", height = "700px"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),

              # ============================================================
              # HEATMAP-ZOOM MULTI-SELECT (revertible) — Proteomics BP
              # ============================================================
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  div(
                    style = paste("display:flex; gap:14px; align-items:center;",
                                  "flex-wrap:wrap; padding:8px 0; margin-top:4px;",
                                  "border-top:1px solid #eee;"),
                    uiOutput("prot_bp_zoom_msg_ui"),
                    actionButton(
                      "prot_bp_zoom_open_btn",
                      "Open visible proteins in Multi-Protein Explorer",
                      icon  = icon("up-right-from-square"),
                      class = "btn btn-primary btn-sm"
                    )
                  )
                ),
                column(width = 1)
              ),
              # ============================================================
              # END HEATMAP-ZOOM (Proteomics BP)
              # ============================================================

              tags$hr(),

              uiOutput("prot_bp_results_title_ui"),
              div(
                uiOutput("prot_dl_results_de_bp_filtered_ui"),
                style = "margin-bottom: 8px;"
              )
            )
          )
        ),
        
        ## ---- Subtab: Enrichment (generic ORA/GSEA — same code as RNA) ----
        ## Wired to the omics-agnostic enrichment UI (functions_shiny.R:
        ## enrichment_panel_ui / wire_enrichment; server side in server.R).
        ## Proteomics currently has no enrichment payload, so this shows only the
        ## friendly unavailable message until proteomics_env$enrichment (the §14
        ## contract) is added upstream — no UI/server rewrite needed then.
        enrichment_panel_ui("prot_enrich"),
        
        ## ---- Subtab: Keyword Search ----
        bslib::nav_panel(
          "Keyword Search",
          keyword_search_panel_ui(id_prefix = "prot_kws")
        ),

        ## ---- Subtab: Protein Explorer ----
        bslib::nav_panel(
          "Protein Explorer",
          fluidRow(
            column(
              width = 12,
              
              tags$label("Select Protein:", `for` = "prot_explorer_protein",
                         style = "font-size: 18px; font-weight: 600; display:block; margin-bottom:5px;"),
              div(
                style = "display:flex; gap:8px; align-items:center; flex-wrap:wrap;",
                htmltools::tagAppendAttributes(
                  selectizeInput(
                    inputId  = "prot_explorer_protein",
                    label    = NULL,
                    choices  = NULL,
                    width    = "500px",
                    multiple = FALSE,
                    options  = list(
                      placeholder = "Type or paste a Protein ID...",
                      maxOptions  = 20
                    )
                  ),
                  style = "margin-bottom: 0;"
                ),
                actionButton("prot_explorer_clear_btn", "Clear selection",
                             icon  = icon("trash"),
                             class = "btn btn-default btn-sm")
              ),
              
              br(),
              
              uiOutput("prot_explorer_selected_title"),
              uiOutput("prot_explorer_plot_title"),
              
              tags$hr(style = "margin-top: 5px; margin-bottom: 10px; border-color: #ddd;"),
              
              conditionalPanel(
                condition = "input.prot_explorer_protein != null && input.prot_explorer_protein != ''",
                div(
                  style = "display: flex; justify-content: flex-start; gap: 12px; align-items: center; margin-bottom: 10px;",
                  
                  selectInput(
                    inputId  = "prot_explorer_scale",
                    label    = "Missing Value Imputation:",
                    choices  = NULL,
                    width    = "280px"
                  ),
                  
                  selectInput(
                    inputId  = "prot_explorer_color_by",
                    label    = "Color bars by:",
                    choices  = NULL,
                    width    = "220px"
                  )
                )
              ),
              
              br(),
              
              plotOutput("prot_explorer_counts_plot", height = "400px"),

              br(),

              conditionalPanel(
                condition = "input.prot_explorer_protein != null && input.prot_explorer_protein != ''",
                div(
                  style = "display: flex; justify-content: flex-end; gap: 12px; align-items: center;",
                  downloadButton("download_prot_explorer_plot", "Download Plot")
                )
              ),

              br(),

              # Expression Data — wide-format table for the single selected protein
              conditionalPanel(
                condition = "input.prot_explorer_protein != null && input.prot_explorer_protein != ''",
                div(
                  h4("Expression Data", class = "plot-title full"),
                  DT::DTOutput("prot_explorer_data_table"),
                  br(),
                  div(
                    style = "display: flex; justify-content: flex-end;",
                    downloadButton("download_prot_explorer_data", "Download Data")
                  )
                )
              ),

              br(),

              conditionalPanel(
                condition = "input.prot_explorer_protein != null && input.prot_explorer_protein != ''",
                div(
                  h4("Differential Expression Statistics", class = "plot-title full"),
                  
                  DT::dataTableOutput("prot_explorer_de_table"),
                  
                  div(
                    style = paste(
                      "display:flex;",
                      "justify-content:space-between;",
                      "align-items:flex-start;",
                      "margin-top:10px;",
                      "gap:16px;"
                    ),
                    
                    div(
                      style = "font-size:14px; line-height:1.6;",
                      uiOutput("prot_explorer_de_legend_ui")
                    ),
                    
                    div(
                      style = "display:flex; gap:10px; align-items:flex-end; flex:0 0 auto; min-width:280px;",
                      uiOutput("prot_explorer_cutoffs_ui")
                    )
                  ),
                  
                  tags$hr(),

                  uiOutput("prot_explorer_annot_title_ui"),
                  DT::DTOutput("prot_explorer_annot_dt"),
                  uiOutput("prot_explorer_annot_dl_ui"),

                  br(),

                  # Only when the MultiOmics tab exists (its target)
                  if (SHOW_TAB$multiomics) conditionalPanel(
                    condition = "input.prot_explorer_protein != null && input.prot_explorer_protein != ''",
                    div(
                      style = "margin-top: 10px; margin-bottom: 20px;",
                      actionButton("prot_explorer_multi_btn", "View in MultiOmics Explorer",
                                   class = "btn btn-info btn-sm")
                    )
                  ),

                  br()
                )
              )
            )
          )
        ),

        ## ---- Subtab: Multi-Protein Explorer (Proteomics) ----

        bslib::nav_panel(
          "Multi-Protein Explorer",
          fluidRow(
            column(
              width = 12,

              h4("Multi-Protein Explorer", class = "plot-title full"),

              tags$label("Select proteins (up to 100):", `for` = "mpe_genes_input",
                         style = "font-size:14px; font-weight:600; display:block; margin-top:8px; margin-bottom:5px;"),
              div(
                style = paste("display:flex; gap:12px; align-items:center;",
                              "flex-wrap:wrap;"),
                div(
                  style = "flex:1 1 360px; min-width:300px;",
                  htmltools::tagAppendAttributes(
                    selectizeInput(
                      inputId  = "mpe_genes_input",
                      label    = NULL,
                      choices  = NULL,
                      multiple = TRUE,
                      width    = "100%",
                      options  = list(
                        placeholder = "Type or paste Protein IDs...",
                        maxOptions  = 50
                      )
                    ),
                    style = "margin-bottom: 0;"
                  )
                ),
                div(
                  style = "display:flex; gap:8px;",
                  actionButton("mpe_add_btn",   "Add selected proteins",
                               icon = icon("plus"),
                               class = "btn btn-primary btn-sm"),
                  actionButton("mpe_clear_btn", "Clear selection",
                               icon = icon("trash"),
                               class = "btn btn-default btn-sm"),
                  # Compact secondary action — opens a modal with the textarea
                  actionButton("mpe_paste_open", "Paste feature list",
                               icon  = icon("clipboard-list"),
                               class = "btn btn-default btn-sm")
                )
              ),

              uiOutput("mpe_count_ui"),
              uiOutput("mpe_summary_ui"),

              conditionalPanel(
                condition = "output.mpe_has_selection == 'yes'",

                br(),

                # 1. Selected-features hierarchical heatmap — shows ALL
                #    selected features, not the paginated subset.
                uiOutput("mpe_selected_heatmap_title"),
                uiOutput("mpe_selected_heatmap_msg"),
                plotly::plotlyOutput("mpe_selected_heatmap", height = "520px"),

                br(),

                # 2. Pagination controls
                uiOutput("mpe_pagination_ui"),

                # 3. Expression Profile (current page only)
                h4("Expression Profile", class = "plot-title full"),

                div(
                  style = "display:flex; gap:12px; align-items:flex-end; flex-wrap:wrap; margin-bottom:8px;",
                  selectInput("mpe_scale",    label = "Missing Value Imputation:",
                              choices = NULL, width = "260px"),
                  selectInput("mpe_color_by", label = "Color bars by:",
                              choices = NULL, width = "200px"),
                  selectInput("mpe_sort_by",  label = "Sort by:",
                              choices  = c("Selection order"                    = "selection",
                                           "Feature name A–Z"              = "alpha",
                                           "Significance in selected contrast"  = "significance",
                                           "linearFC in selected contrast"      = "linearfc"),
                              selected = "selection",
                              width    = "260px"),
                  uiOutput("mpe_sort_contrast_ui")
                ),

                div(
                  style = "overflow:auto; max-height:1200px; border:1px solid #eee; padding:6px; border-radius:6px;",
                  plotOutput("mpe_plot", height = "auto")
                ),

                br(),

                div(
                  style = "display:flex; justify-content:flex-end;",
                  downloadButton("download_mpe_plot", "Download Plot")
                ),

                br(),

                # 4. Expression Data (current page only)
                h4("Expression Data", class = "plot-title full"),

                DT::DTOutput("mpe_data_table"),

                br(),

                div(
                  style = "display:flex; justify-content:flex-end;",
                  downloadButton("download_mpe_table", "Download Data")
                ),

                br(),

                # 5. Differential Expression Statistics (current page only)
                h4("Differential Expression Statistics", class = "plot-title full"),

                uiOutput("mpe_de_contrast_ui"),

                DT::dataTableOutput("mpe_de_table"),

                div(
                  style = paste("display:flex;", "justify-content:space-between;",
                                "align-items:flex-start;", "margin-top:10px;",
                                "gap:16px;"),
                  div(
                    style = "font-size:14px; line-height:1.6;",
                    uiOutput("mpe_de_legend_ui")
                  ),
                  div(
                    style = "display:flex; gap:10px; align-items:flex-end; flex:0 0 auto; min-width:280px;",
                    uiOutput("mpe_de_cutoffs_ui")
                  )
                ),

                br(),

                # 6. Annotation (full current selection)
                uiOutput("mpe_annot_title_ui"),
                DT::DTOutput("mpe_annot_dt"),
                uiOutput("mpe_annot_dl_ui"),

                br()
              )
            )
          )
        )
      )
    ),


    # ---------- tab panel 6 - Metabolomics ----------

    # Shown only when data/shiny_payload_metabolomics.rds exists (SHOW_TAB in global.R)
    if (SHOW_TAB$metabolomics) tabPanel(
      "Metabolomics",
      titlePanel("Metabolomics Analysis"),
      bslib::navset_pill_list(
        id     = "metab_tabs",
        widths = c(2, 9),

        ## ---- Subtab: Exploratory Analyses ----
        bslib::nav_panel(
          "Exploratory Analyses",
          bslib::navset_tab(
            id = "metab_explore_subtab",

            ## ---- Subtab: Data distributions ----
            bslib::nav_panel(
              "Data distributions",

              h4("Normalized expression boxplots", class = "plot-title full"),
              shiny::uiOutput("metab_norm_boxplot_controls"),
              shinycssloaders::withSpinner(
                plotOutput("metab_norm_boxplot", height = "550px"),
                type = 4
              ),

              h4("Density plot of normalized intensities", class = "plot-title full"),
              shinycssloaders::withSpinner(
                plotly::plotlyOutput("metab_density_plot", height = "550px"),
                type = 4
              )
            ),

            ## ---- Subtab: Data structure ----
            bslib::nav_panel(
              "Data structure",

              h4("3D PCA", class = "plot-title full"),
              shinycssloaders::withSpinner(plotly::plotlyOutput("metab_pca_3d_ui", height = "600px"), type = 4),
              br(), hr(),

              h4("2D PCA", class = "plot-title full"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  h4("PC1 vs PC2"),
                  shinycssloaders::withSpinner(
                    plotlyOutput("metab_pca_pc1_pc2_plot", width = "100%"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),
              br(),

              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  h4("PC1 vs PC3"),
                  shinycssloaders::withSpinner(
                    plotlyOutput("metab_pca_pc1_pc3_plot", width = "100%"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),

              br(), hr(),

              h4("Sample Distance Heatmap", class = "plot-title full"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotlyOutput("metab_sample_corr_heatmap", height = "700px"),
                    type = 4
                  )
                ),
                column(width = 1)
              )
            ),

            ## ---- Subtab: DE Metabolite ----
            bslib::nav_panel(
              "DE Metabolite",
              fluidRow(
                column(
                  width = 3,
                  tags$h4("Choose contrast:"),
                  selectInput("metab_de_contrast", label = NULL, choices = NULL)
                ),
                column(
                  width = 9,
                  h4("Volcano plot", class = "plot-title full"),
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("metab_volcano_plot", height = "520px"),
                    type = 4
                  ),
                  tags$hr(),
                  h4("MA plot", class = "plot-title full"),
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("metab_ma_plot", height = "520px"),
                    type = 4
                  )
                )
              )
            )
          )
        ),

        ## ---- Subtab: Stats results ----
        bslib::nav_panel(
          "Stats results",
          h4("Summary DE statistics", class = "plot-title full"),
          tableOutput("metab_DE_stats"),
          tags$hr(),

          h4("Download result tables", class = "plot-title full"),
          uiOutput("metab_results_buttons"),
          tags$hr(),

          h4("Download metabolite intensity matrices", class = "plot-title full"),
          div(
            style = "margin-bottom: 8px;",
            downloadButton(
              "metab_dl_log2_raw_all_txt",
              "Log2 matrix - All metabolites (raw)",
              class = "btn btn-success btn-sm",
              style = "margin:4px;"
            )
          ),
          div(
            style = "margin-bottom: 8px;",
            downloadButton(
              "metab_dl_log2_norm_all_txt",
              "Log2 matrix - All metabolites (normalized)",
              class = "btn btn-success btn-sm",
              style = "margin:4px;"
            )
          ),
          div(
            downloadButton(
              "metab_dl_norm_de_txt",
              "Log2 matrix - DE metabolites",
              class = "btn btn-success btn-sm",
              style = "margin:4px;"
            )
          )
        ),

        ## ---- Subtab: Clustering ----
        bslib::nav_panel(
          "Clustering",
          bslib::navset_tab(
            id = "metab_clustering_subtab",

            ## ---- Subtab: Hierarchical Clustering ----
            bslib::nav_panel(
              "Hierarchical Clustering",

              h4("Hierarchical clustering heatmap", class = "plot-title full"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotOutput("hc_static_fig_metab", height = "700px"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),

              tags$h4(
                "Interactive Hierarchical Clustering",
                class = "plot-title full",
                style = "text-align:center;margin-top:15px;margin-bottom:10px;"
              ),
              tags$hr(),

              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("metab_de_heatmap", height = "700px"),
                    type = 4, color = "#2c3e50"
                  )
                ),
                column(width = 1)
              ),

              # ============================================================
              # HEATMAP-ZOOM MULTI-SELECT (revertible) — Metabolomics HC
              # ============================================================
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  div(
                    style = paste("display:flex; gap:14px; align-items:center;",
                                  "flex-wrap:wrap; padding:8px 0; margin-top:4px;",
                                  "border-top:1px solid #eee;"),
                    uiOutput("metab_hc_zoom_msg_ui"),
                    actionButton(
                      "metab_hc_zoom_open_btn",
                      "Open visible metabolites in Multi-Metabolite Explorer",
                      icon  = icon("up-right-from-square"),
                      class = "btn btn-primary btn-sm"
                    )
                  )
                ),
                column(width = 1)
              )
              # ============================================================
              # END HEATMAP-ZOOM (Metabolomics HC)
              # ============================================================
            ),

            ## ---- Subtab: Partition Clustering ----
            bslib::nav_panel(
              "Partition Clustering",

              uiOutput("metab_pc_title_ui"),

              # Static figure from env (gtable)
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotOutput("pc_static_fig_metab", height = "700px"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),
              tags$hr(),

              selectInput("metab_pc_cluster_pick", "Choose cluster", choices = NULL, width = "300px"),
              uiOutput("metab_pc_cluster_title_ui"),
              uiOutput("metab_pc_cluster_note_ui"),

              fluidRow(
                column(
                  width = 9,
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("metab_pc_heatmap_cluster", height = "700px"),
                    type = 4, color = "#2c3e50"
                  )
                ),
                column(
                  width = 3,
                  shinycssloaders::withSpinner(
                    plotOutput("metab_pc_cluster_subgraph", height = "550px"),
                    type = 4, color = "#2c3e50"
                  )
                )
              ),

              # ============================================================
              # HEATMAP-ZOOM MULTI-SELECT (revertible) — Metabolomics Partition
              # ============================================================
              fluidRow(
                column(
                  width = 9,
                  div(
                    style = paste("display:flex; gap:14px; align-items:center;",
                                  "flex-wrap:wrap; padding:8px 0; margin-top:4px;",
                                  "border-top:1px solid #eee;"),
                    uiOutput("metab_pc_zoom_msg_ui"),
                    actionButton(
                      "metab_pc_zoom_open_btn",
                      "Open visible metabolites in Multi-Metabolite Explorer",
                      icon  = icon("up-right-from-square"),
                      class = "btn btn-primary btn-sm"
                    )
                  )
                )
              ),
              # ============================================================
              # END HEATMAP-ZOOM (Metabolomics Partition)
              # ============================================================

              uiOutput("metab_pc_results_title_ui"),
              div(
                uiOutput("metab_dl_results_de_filtered_ui"),
                style = "margin-bottom: 8px;"
              )
            ),

            ## ---- Subtab: Binary Pattern ----
            bslib::nav_panel(
              "Binary Pattern",
              fluidRow(
                column(4, uiOutput("metab_bp_pattern_ui"))
              ),
              tags$hr(),

              uiOutput("metab_bp_static_title_ui"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  plotOutput("metab_bp_static", height = "700px")
                ),
                column(width = 1)
              ),
              tags$hr(),

              uiOutput("metab_bp_interactive_title_ui"),
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  shinycssloaders::withSpinner(
                    plotly::plotlyOutput("metab_bp_interactive", height = "700px"),
                    type = 4
                  )
                ),
                column(width = 1)
              ),

              # ============================================================
              # HEATMAP-ZOOM MULTI-SELECT (revertible) — Metabolomics BP
              # ============================================================
              fluidRow(
                column(width = 1),
                column(
                  width = 10,
                  div(
                    style = paste("display:flex; gap:14px; align-items:center;",
                                  "flex-wrap:wrap; padding:8px 0; margin-top:4px;",
                                  "border-top:1px solid #eee;"),
                    uiOutput("metab_bp_zoom_msg_ui"),
                    actionButton(
                      "metab_bp_zoom_open_btn",
                      "Open visible metabolites in Multi-Metabolite Explorer",
                      icon  = icon("up-right-from-square"),
                      class = "btn btn-primary btn-sm"
                    )
                  )
                ),
                column(width = 1)
              ),
              # ============================================================
              # END HEATMAP-ZOOM (Metabolomics BP)
              # ============================================================

              tags$hr(),

              uiOutput("metab_bp_results_title_ui"),
              div(
                uiOutput("metab_dl_results_de_bp_filtered_ui"),
                style = "margin-bottom: 8px;"
              )
            )
          )
        ),

        ## ---- Subtab: Enrichment (generic ORA/GSEA — same code as RNA) ----
        ## Wired to the omics-agnostic enrichment UI (functions_shiny.R:
        ## enrichment_panel_ui / wire_enrichment; server side in server.R).
        ## Metabolomics currently has no enrichment payload, so this shows only
        ## the friendly unavailable message until metabolomics_env$enrichment
        ## (the §14 contract) is added upstream — no UI/server rewrite needed then.
        enrichment_panel_ui("metab_enrich"),

        ## ---- Subtab: Keyword Search ----
        bslib::nav_panel(
          "Keyword Search",
          keyword_search_panel_ui(id_prefix = "metab_kws")
        ),

        ## ---- Subtab: Metabolite Explorer ----
        bslib::nav_panel(
          "Metabolite Explorer",
          fluidRow(
            column(
              width = 12,

              tags$label("Select Metabolite:", `for` = "metab_explorer_protein",
                         style = "font-size: 18px; font-weight: 600; display:block; margin-bottom:5px;"),
              div(
                style = "display:flex; gap:8px; align-items:center; flex-wrap:wrap;",
                htmltools::tagAppendAttributes(
                  selectizeInput(
                    inputId  = "metab_explorer_protein",
                    label    = NULL,
                    choices  = NULL,
                    width    = "500px",
                    multiple = FALSE,
                    options  = list(
                      placeholder = "Type or paste a Metabolite ID...",
                      maxOptions  = 20
                    )
                  ),
                  style = "margin-bottom: 0;"
                ),
                actionButton("metab_explorer_clear_btn", "Clear selection",
                             icon  = icon("trash"),
                             class = "btn btn-default btn-sm")
              ),

              br(),

              uiOutput("metab_explorer_selected_title"),
              uiOutput("metab_explorer_plot_title"),

              tags$hr(style = "margin-top: 5px; margin-bottom: 10px; border-color: #ddd;"),

              conditionalPanel(
                condition = "input.metab_explorer_protein != null && input.metab_explorer_protein != ''",
                div(
                  style = "display: flex; justify-content: flex-start; gap: 12px; align-items: flex-start; margin-bottom: 10px;",
                  # Metabolomics always uses the normalized matrix — the
                  # raw/non-normalized option is no longer exposed. The
                  # wire_protein_explorer call sets fixed_scale = "imp1"
                  # server-side, so no scale selectInput is rendered. We
                  # only display the title + chosen_norm-derived text.
                  # align-items: flex-start keeps "Normalization Method:" on
                  # the same horizontal line as the "Color bars by:" label.
                  div(
                    style = "width: 280px;",
                    tags$label("Normalization Method:",
                               style = "font-weight: bold; display: block; margin-bottom: 5px;"),
                    uiOutput("metab_explorer_norm_text")
                  ),
                  selectInput(
                    inputId  = "metab_explorer_color_by",
                    label    = "Color bars by:",
                    choices  = NULL,
                    width    = "220px"
                  )
                )
              ),

              br(),

              plotOutput("metab_explorer_counts_plot", height = "400px"),

              br(),

              conditionalPanel(
                condition = "input.metab_explorer_protein != null && input.metab_explorer_protein != ''",
                div(
                  style = "display: flex; justify-content: flex-end; gap: 12px; align-items: center;",
                  downloadButton("download_metab_explorer_plot", "Download Plot")
                )
              ),

              br(),

              # Expression Data — wide-format table for the single selected metabolite
              conditionalPanel(
                condition = "input.metab_explorer_protein != null && input.metab_explorer_protein != ''",
                div(
                  h4("Expression Data", class = "plot-title full"),
                  DT::DTOutput("metab_explorer_data_table"),
                  br(),
                  div(
                    style = "display: flex; justify-content: flex-end;",
                    downloadButton("download_metab_explorer_data", "Download Data")
                  )
                )
              ),

              br(),

              conditionalPanel(
                condition = "input.metab_explorer_protein != null && input.metab_explorer_protein != ''",
                div(
                  h4("Differential Expression Statistics", class = "plot-title full"),

                  DT::dataTableOutput("metab_explorer_de_table"),

                  div(
                    style = paste(
                      "display:flex;",
                      "justify-content:space-between;",
                      "align-items:flex-start;",
                      "margin-top:10px;",
                      "gap:16px;"
                    ),
                    div(
                      style = "font-size:14px; line-height:1.6;",
                      uiOutput("metab_explorer_de_legend_ui")
                    ),
                    div(
                      style = "display:flex; gap:10px; align-items:flex-end; flex:0 0 auto; min-width:280px;",
                      uiOutput("metab_explorer_cutoffs_ui")
                    )
                  ),

                  tags$hr(),

                  uiOutput("metab_explorer_annot_title_ui"),
                  DT::DTOutput("metab_explorer_annot_dt"),
                  uiOutput("metab_explorer_annot_dl_ui"),

                  br(), br()
                )
              )
            )
          )
        ),

        ## ---- Subtab: Multi-Metabolite Explorer (Metabolomics) ----

        bslib::nav_panel(
          "Multi-Metabolite Explorer",
          fluidRow(
            column(
              width = 12,

              h4("Multi-Metabolite Explorer", class = "plot-title full"),

              tags$label("Select metabolites (up to 100):", `for` = "mme_genes_input",
                         style = "font-size:14px; font-weight:600; display:block; margin-top:8px; margin-bottom:5px;"),
              div(
                style = paste("display:flex; gap:12px; align-items:center;",
                              "flex-wrap:wrap;"),
                div(
                  style = "flex:1 1 360px; min-width:300px;",
                  htmltools::tagAppendAttributes(
                    selectizeInput(
                      inputId  = "mme_genes_input",
                      label    = NULL,
                      choices  = NULL,
                      multiple = TRUE,
                      width    = "100%",
                      options  = list(
                        placeholder = "Type or paste metabolite IDs...",
                        maxOptions  = 50
                      )
                    ),
                    style = "margin-bottom: 0;"
                  )
                ),
                div(
                  style = "display:flex; gap:8px;",
                  actionButton("mme_add_btn",   "Add selected metabolites",
                               icon = icon("plus"),
                               class = "btn btn-primary btn-sm"),
                  actionButton("mme_clear_btn", "Clear selection",
                               icon = icon("trash"),
                               class = "btn btn-default btn-sm"),
                  # Compact secondary action — opens a modal with the textarea
                  actionButton("mme_paste_open", "Paste feature list",
                               icon  = icon("clipboard-list"),
                               class = "btn btn-default btn-sm")
                )
              ),

              uiOutput("mme_count_ui"),
              uiOutput("mme_summary_ui"),

              conditionalPanel(
                condition = "output.mme_has_selection == 'yes'",

                br(),

                # 1. Selected-features hierarchical heatmap — shows ALL
                #    selected features, not the paginated subset.
                uiOutput("mme_selected_heatmap_title"),
                uiOutput("mme_selected_heatmap_msg"),
                plotly::plotlyOutput("mme_selected_heatmap", height = "520px"),

                br(),

                # 2. Pagination controls
                uiOutput("mme_pagination_ui"),

                # 3. Expression Profile (current page only)
                h4("Expression Profile", class = "plot-title full"),

                div(
                  style = "display:flex; gap:12px; align-items:flex-start; flex-wrap:wrap; margin-bottom:8px;",
                  # Static normalization display (no dropdown). See the
                  # metab single Explorer for the same pattern + rationale.
                  # wire_multi_gene_explorer receives fixed_scale = "imp1"
                  # server-side so no mme_scale selectInput is needed.
                  # align-items: flex-start keeps "Normalization Method:" on
                  # the same horizontal line as the "Color bars by:" label.
                  div(
                    style = "width: 260px;",
                    tags$label("Normalization Method:",
                               style = "font-weight: bold; display: block; margin-bottom: 5px;"),
                    uiOutput("mme_norm_text")
                  ),
                  selectInput("mme_color_by", label = "Color bars by:",
                              choices = NULL, width = "200px"),
                  selectInput("mme_sort_by",  label = "Sort by:",
                              choices  = c("Selection order"                    = "selection",
                                           "Feature name A–Z"              = "alpha",
                                           "Significance in selected contrast"  = "significance",
                                           "linearFC in selected contrast"      = "linearfc"),
                              selected = "selection",
                              width    = "260px"),
                  uiOutput("mme_sort_contrast_ui")
                ),

                div(
                  style = "overflow:auto; max-height:1200px; border:1px solid #eee; padding:6px; border-radius:6px;",
                  plotOutput("mme_plot", height = "auto")
                ),

                br(),

                div(
                  style = "display:flex; justify-content:flex-end;",
                  downloadButton("download_mme_plot", "Download Plot")
                ),

                br(),

                # 4. Expression Data (current page only)
                h4("Expression Data", class = "plot-title full"),

                DT::DTOutput("mme_data_table"),

                br(),

                div(
                  style = "display:flex; justify-content:flex-end;",
                  downloadButton("download_mme_table", "Download Data")
                ),

                br(),

                # 5. Differential Expression Statistics (current page only)
                h4("Differential Expression Statistics", class = "plot-title full"),

                uiOutput("mme_de_contrast_ui"),

                DT::dataTableOutput("mme_de_table"),

                div(
                  style = paste("display:flex;", "justify-content:space-between;",
                                "align-items:flex-start;", "margin-top:10px;",
                                "gap:16px;"),
                  div(
                    style = "font-size:14px; line-height:1.6;",
                    uiOutput("mme_de_legend_ui")
                  ),
                  div(
                    style = "display:flex; gap:10px; align-items:flex-end; flex:0 0 auto; min-width:280px;",
                    uiOutput("mme_de_cutoffs_ui")
                  )
                ),

                br(),

                # 6. Annotation (full current selection)
                uiOutput("mme_annot_title_ui"),
                DT::DTOutput("mme_annot_dt"),
                uiOutput("mme_annot_dl_ui"),

                br()
              )
            )
          )
        )

      )  # end navset_pill_list Metabolomics
    ),  # end tabPanel Metabolomics
    
    # ---------- tab panel 7 - Lipidomics  ----------
    ## tab panel 6 - Lipidomics
    # Shown only when data/shiny_payload_lipidomics.rds exists (SHOW_TAB in global.R)
    if (SHOW_TAB$lipidomics) tabPanel("Lipidomics", "Relevant information has not yet been provided"),
    
    # ---------- tab panel 8 - MultiOmics  ----------
    ## tab panel 7 - MultiOmics
    # Shown only when at least two individual omics payloads exist (SHOW_TAB in global.R)
    if (SHOW_TAB$multiomics) tabPanel(
      HTML("<b><span style='color:red;'>MultiOmics</span></b>"),
      value = "MultiOmics",

      bslib::navset_pill_list(
        id      = "multi_tabs",
        widths  = c(2, 9),

        ## ---- Subtab: MultiOmics Explorer ----
        bslib::nav_panel(
          "MultiOmics Explorer",
          fluidRow(
            column(
              width = 12,

              # Message when multiomics matching data is unavailable
              uiOutput("multi_unavailable_msg"),

              # Search input
              tags$label("Search by Gene or Protein ID:", `for` = "multi_search_query",
                         style = "font-size: 18px; font-weight: 600; display:block; margin-bottom:5px;"),
              div(
                style = "display:flex; gap:8px; align-items:center; flex-wrap:wrap;",
                htmltools::tagAppendAttributes(
                  selectizeInput(
                    inputId  = "multi_search_query",
                    label    = NULL,
                    choices  = NULL,
                    width    = "500px",
                    multiple = FALSE,
                    options  = list(
                      placeholder = "Type a Gene ID or Protein ID...",
                      maxOptions  = 20
                    )
                  ),
                  style = "margin-bottom: 0;"
                ),
                actionButton("multi_search_clear_btn", "Clear selection",
                             icon  = icon("trash"),
                             class = "btn btn-default btn-sm")
              ),

              uiOutput("multi_pair_info"),

              conditionalPanel(
                condition = "input.multi_search_query != null && input.multi_search_query != ''",

                # Hidden per-panel inputs (driven by shared selectors;
                # required by the wire functions that read input$multi_ge_color_by etc.)
                div(style = "display:none;",
                    selectInput("multi_ge_color_by", label = NULL, choices = NULL),
                    selectInput("multi_pe_color_by", label = NULL, choices = NULL)
                ),

                # ============================================================
                # EXPRESSION PROFILE SECTION
                # ============================================================

                # ---- Unified section title ----
                h4("Expression Profile", class = "plot-title full",
                   style = "margin-top:12px; margin-bottom:4px;"),
                tags$hr(style = "margin-top:2px; margin-bottom:10px; border-color:#ddd;"),

                # ---- Shared controls row (before labels and plots) ----
                div(
                  style = "display:flex; flex-wrap:wrap; justify-content:center; align-items:flex-end; gap:20px; margin-bottom:12px;",
                  selectInput(inputId = "multi_shared_color_by",
                              label   = "Color by:",
                              choices = NULL, width = "260px"),
                  selectInput(inputId = "multi_ge_scale",
                              label   = "Gene normalization:",
                              choices = NULL, width = "280px"),
                  selectInput(inputId = "multi_pe_scale",
                              label   = "Protein Missing Value Imputation:",
                              choices = NULL, width = "280px")
                ),

                # ---- Gene (left) / Protein (right) labels above plots ----
                fluidRow(
                  column(6,
                    uiOutput("multi_ge_selected_gene_title"),
                    uiOutput("multi_ge_no_data_msg")
                  ),
                  column(6,
                    uiOutput("multi_pe_selected_title"),
                    uiOutput("multi_pe_no_data_msg")
                  )
                ),

                # ---- Side-by-side plots ----
                fluidRow(
                  column(6, plotOutput("multi_ge_counts_plot", height = "350px")),
                  column(6, plotOutput("multi_pe_counts_plot", height = "350px"))
                ),

                # ---- Download button (right-aligned) ----
                div(
                  style = "display:flex; justify-content:flex-end; margin-top:10px; margin-bottom:6px;",
                  downloadButton("download_multi_combined_plot", "Download Plot")
                ),

                # ============================================================
                # EXPRESSION DATA SECTION (combined RNA + Protein)
                # ------------------------------------------------------------
                # One shared table with an Omics column distinguishing rows.
                # Sample columns are the union of both omics; cells are NA
                # where a sample doesn't exist in the other omics. The DT and
                # download handler are wired in server.R (not via the wire
                # helpers — those still bind multi_ge_data_table /
                # multi_pe_data_table outputs but they're unreferenced and
                # stay dormant since DT::renderDT suspendWhenHidden is on).
                # ============================================================

                h4("Expression Data", class = "plot-title full",
                   style = "margin-top:14px; margin-bottom:4px;"),
                DT::DTOutput("multi_combined_data_table"),
                br(),
                div(
                  style = "display:flex; justify-content:flex-end;",
                  downloadButton("download_multi_combined_data",
                                 "Download Data")
                ),

                # ============================================================
                # DIFFERENTIAL EXPRESSION STATISTICS SECTION
                # ============================================================

                # ---- Unified statistics title ----
                h4("Differential Expression Statistics", class = "plot-title full",
                   style = "margin-top:14px; margin-bottom:4px;"),
                tags$hr(style = "margin-top:2px; margin-bottom:10px; border-color:#ddd;"),

                # ---- Gene (left) / Protein (right) labels above tables ----
                fluidRow(
                  column(6,
                    tags$div(
                      style = "text-align:center; margin-bottom:6px;",
                      tags$span("Gene: ", style = "color:#777; font-size:14px; font-weight:400;"),
                      uiOutput("multi_ge_de_label", inline = TRUE)
                    )
                  ),
                  column(6,
                    tags$div(
                      style = "text-align:center; margin-bottom:6px;",
                      tags$span("Protein: ", style = "color:#777; font-size:14px; font-weight:400;"),
                      uiOutput("multi_pe_de_label", inline = TRUE)
                    )
                  )
                ),

                # ---- DE tables side-by-side ----
                fluidRow(

                  # Gene DE table (left)
                  column(6,
                    DT::dataTableOutput("multi_ge_de_table"),

                    div(
                      style = paste("display:flex;", "justify-content:space-between;",
                                    "align-items:flex-start;", "margin-top:10px;", "gap:16px;"),
                      div(style = "font-size:14px; line-height:1.6;",
                          uiOutput("multi_ge_de_legend_ui")),
                      div(style = "display:flex; gap:10px; align-items:flex-end; flex:0 0 auto; min-width:280px;",
                          uiOutput("multi_ge_cutoffs_ui"))
                    ),

                    br(), br()
                  ),

                  # Protein DE table (right)
                  column(6,
                    DT::dataTableOutput("multi_pe_de_table"),

                    div(
                      style = paste("display:flex;", "justify-content:space-between;",
                                    "align-items:flex-start;", "margin-top:10px;", "gap:16px;"),
                      div(style = "font-size:14px; line-height:1.6;",
                          uiOutput("multi_pe_de_legend_ui")),
                      div(style = "display:flex; gap:10px; align-items:flex-end; flex:0 0 auto; min-width:280px;",
                          uiOutput("multi_pe_cutoffs_ui"))
                    ),

                    br(), br()
                  )

                )  # end DE tables fluidRow
                ,

                # ============================================================
                # ANNOTATION SECTION (gene + protein, side by side)
                # ------------------------------------------------------------
                # One shared "Annotation" title; left = gene annotation
                # (rnaseq feature_annot via multi_ge), right = protein
                # annotation (proteomics feature_annot via multi_pe). Each side
                # renders its own table / fallback message independently, so a
                # missing side does not hide the other. Excel download under
                # each table. Reuses the canonical wire_feature_annotation path.
                # ============================================================

                h4("Annotation", class = "plot-title full",
                   style = "margin-top:14px; margin-bottom:4px;"),
                tags$hr(style = "margin-top:2px; margin-bottom:10px; border-color:#ddd;"),

                fluidRow(
                  # Gene annotation (left)
                  column(6,
                    tags$div(style = "font-weight:600; margin-bottom:6px;",
                             "Gene annotation"),
                    DT::DTOutput("multi_ge_annot_dt"),
                    uiOutput("multi_ge_annot_dl_ui")
                  ),
                  # Protein annotation (right)
                  column(6,
                    tags$div(style = "font-weight:600; margin-bottom:6px;",
                             "Protein annotation"),
                    DT::DTOutput("multi_pe_annot_dt"),
                    uiOutput("multi_pe_annot_dl_ui")
                  )
                ),

                # ============================================================
                # GENE–PROTEIN CORRELATION SECTION
                # ------------------------------------------------------------
                # Interactive Plotly scatter (one point per biological sample,
                # hover shows sample + values) plus a statistics table listing
                # every correlation method that can be computed for the pair.
                # No method dropdown: the scatter is method-agnostic. A friendly
                # message is shown when nothing can be calculated. Methods are
                # detected dynamically via the shared validation path. Isolated
                # to MultiOmics.
                # ============================================================

                h4("Gene–Protein Correlation", class = "plot-title full",
                   style = "margin-top:14px; margin-bottom:4px;"),
                tags$hr(style = "margin-top:2px; margin-bottom:10px; border-color:#ddd;"),

                # Scatter + statistics table when correlation is computable...
                conditionalPanel(
                  condition = "output.multi_corr_available == 'yes'",
                  plotly::plotlyOutput("multi_corr_plot", height = "420px"),
                  tags$div(
                    style = "color:#777; font-size:12px; margin-top:4px; margin-bottom:12px;",
                    paste("The trend line is a simple linear fit shown as a visual",
                          "aid only; it is not specific to Pearson, Spearman, or",
                          "Kendall.")
                  ),
                  h4("Correlation Statistics", class = "plot-title full",
                     style = "margin-top:14px; margin-bottom:4px;"),
                  tags$hr(style = "margin-top:2px; margin-bottom:10px; border-color:#ddd;"),
                  DT::DTOutput("multi_corr_stats_table")
                ),

                # ...otherwise a clear "cannot be calculated" message.
                conditionalPanel(
                  condition = "output.multi_corr_available == 'no'",
                  uiOutput("multi_corr_msg")
                ),

                br(), br()

              )  # end conditionalPanel (multi_search_query)
            )
          )
        )  # end nav_panel MultiOmics Explorer

      )  # end navset_pill_list
    )
  )
)