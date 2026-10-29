# multiomics-shinyapp

A reusable R Shiny application for presenting a researcher's multi-omics analysis
results: RNA-seq, Proteomics and Metabolomics, plus MultiOmics and Lipidomics modules.

> **The MultiOmics and Lipidomics modules are currently under development.**
> Their tabs can appear in the app (see [Which tabs are shown](#which-tabs-are-shown)),
> but they are not complete.

The same code is used for every researcher. Everything researcher-specific lives in
`data/`: the researcher's details, the experiment design table, the analysis payloads
(`.rds`) produced by the upstream pipelines, and optional supporting files.
**No R code changes are needed to create a new researcher instance.**

This README covers running the application on your own computer. Deployment to the
team's Shiny Server is described in a separate guide.

---

## Repository structure

```text
multiomics-shinyapp/
├── global.R              # packages, settings and data loading (shared by ui.R and server.R)
├── ui.R                  # user interface
├── server.R              # server logic (also contains the access check used on the server)
├── functions_shiny.R     # analysis and plotting helpers
├── www/                  # generic images used by the interface
├── data/
│   ├── home_info.example.txt                    # example for data/home_info.txt
│   ├── experiment_design.example.xlsx           # example for data/experiment_design.xlsx
│   ├── Matching_genes_proteins.example.xlsx     # example for data/Matching_genes_proteins.xlsx
│   └── Data_&_Files/
│       └── Data_&_Files_manifest.example.xlsx   # example for Data_&_Files_manifest.xlsx
├── .github/
│   └── CODEOWNERS        # GitHub: whose approval pull requests into main require
├── .gitignore            # keeps researcher files out of Git
└── README.md
```

The app runs in Shiny's two-file mode (`ui.R` + `server.R`, with `global.R` loaded
first). There is intentionally no `app.R`.

### Shared template files vs. local researcher files

| Shared template files (in Git; do not edit per researcher) | Local researcher files (never in Git) |
|---|---|
| `global.R`, `ui.R`, `server.R`, `functions_shiny.R` | `data/home_info.txt` |
| `www/` | `data/experiment_design.xlsx` |
| `data/*.example.*` and `data/Data_&_Files/*.example.xlsx` | `data/Matching_genes_proteins.xlsx` |
| `README.md`, `.gitignore`, `.github/CODEOWNERS` | `data/Data_&_Files/Data_&_Files_manifest.xlsx` and the files it lists |
| | `data/shiny_payload_*.rds` |

The local files are created from the examples (see [Quick start](#quick-start)) and
are ignored by `.gitignore`, so a researcher's information is not picked up by
`git add`. Still, always check `git status` before committing.

---

## Requirements

### R

- Developed and tested with **R 4.5** (R 4.5.2 on Windows).
- The code needs at least **R 4.1** (it uses the native pipe `|>`), but versions older
  than 4.5 have not been tested. Use a current R release; if you must use an older
  one, check that the app works before relying on it.

### R packages

| Group | Packages | Without them |
|---|---|---|
| **Required to start** (loaded by `global.R`) | `shiny`, `shinythemes`, `bslib`, `shinycssloaders`, `DT`, `ggplot2`, `plotly`, `heatmaply`, `pheatmap`, `RColorBrewer`, `dplyr`, `readxl`, `writexl`, `openxlsx`, `bcrypt` | The app does not start |
| **Required for specific plots** (Bioconductor) | `DESeq2` (RNA-seq volcano and MA plots), `limma` (Proteomics and Metabolomics volcano and MA plots) | Those plots show an error |
| **Required for a specific plot** (CRAN) | `ggridges` (enrichment ridge plot) | That plot shows an error |
| **Used directly, installed automatically** | `gridExtra`, `scales`, `tibble`, `jsonlite`, `htmltools`: dependencies of the required packages, listed for completeness | Cannot be missing once the required packages are installed |
| **Optional** | `patchwork` | The enrichment running-score plot is shown without its extra panels |

Install whatever is missing (already installed packages are not reinstalled or
upgraded):

```r
# CRAN packages
pkgs <- c(
  # required to start
  "shiny", "shinythemes", "bslib", "shinycssloaders", "DT",
  "ggplot2", "plotly", "heatmaply", "pheatmap", "RColorBrewer",
  "dplyr", "readxl", "writexl", "openxlsx", "bcrypt",
  # required for a specific plot
  "ggridges",
  # used directly (normally already installed as dependencies)
  "gridExtra", "scales", "tibble", "jsonlite", "htmltools",
  # optional
  "patchwork"
)
install.packages(setdiff(pkgs, rownames(installed.packages())))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(setdiff(c("DESeq2", "limma"), rownames(installed.packages())),
                     update = FALSE)
```

**Platform notes**

- **Windows and macOS:** CRAN and Bioconductor normally provide pre-built packages, so
  no extra tools are needed (the app is developed on Windows; macOS has not been tested).
- **Linux:** packages are usually compiled from source, which needs a compiler and the
  system libraries that the packages declare, among them libcurl, OpenSSL and zlib
  (on Debian/Ubuntu: `sudo apt install r-base-dev libcurl4-openssl-dev libssl-dev zlib1g-dev`).

---

## Quick start

```bash
git clone https://github.com/multiomics-center-Israel/multiomics-shinyapp.git
```

1. **Create the local files from the examples.** The commands never overwrite an
   existing file. Run them from the repository folder:

   ```r
   # R (any OS): returns FALSE instead of overwriting if the file already exists
   # required
   file.copy("data/home_info.example.txt", "data/home_info.txt", overwrite = FALSE)
   file.copy("data/experiment_design.example.xlsx", "data/experiment_design.xlsx", overwrite = FALSE)
   # optional
   file.copy("data/Matching_genes_proteins.example.xlsx",
             "data/Matching_genes_proteins.xlsx", overwrite = FALSE)
   file.copy("data/Data_&_Files/Data_&_Files_manifest.example.xlsx",
             "data/Data_&_Files/Data_&_Files_manifest.xlsx", overwrite = FALSE)
   ```

   ```bash
   # Linux / macOS / Git Bash (-n = do not overwrite)
   cp -n data/home_info.example.txt data/home_info.txt
   cp -n data/experiment_design.example.xlsx data/experiment_design.xlsx
   cp -n data/Matching_genes_proteins.example.xlsx data/Matching_genes_proteins.xlsx                 # optional
   cp -n "data/Data_&_Files/Data_&_Files_manifest.example.xlsx" "data/Data_&_Files/Data_&_Files_manifest.xlsx"  # optional
   ```

2. Edit `data/home_info.txt` with the researcher's details ([format](#home_infotxt-required)).
3. Edit `data/experiment_design.xlsx`: replace the example rows with the researcher's
   samples ([format](#experiment_designxlsx-required)).
4. Copy the researcher's payload files into `data/` using the
   [standard names](#omics-payloads-standard-file-names).
5. Optional:
   - `data/Matching_genes_proteins.xlsx`: replace the example rows with the real
     gene–protein table ([format](#matching_genes_proteinsxlsx-optional)), or delete
     the file if there is none.
   - Data & Files: edit `data/Data_&_Files/Data_&_Files_manifest.xlsx` and put the
     listed files under `data/Data_&_Files/` ([format](#data--files-manifest-optional)),
     or delete the manifest if not needed.
6. Run the app from R / RStudio:

   ```r
   shiny::runApp("path/to/multiomics-shinyapp")
   ```

   The R console should print:

   ```text
   Appstream authentication: OFF (local run)
   Omics tabs shown: rnaseq, proteomics, multiomics      # depends on your payloads
   ```

   Access control (Appstream) applies only on the team's server. It is off when the
   app runs on your own computer.

> **Do not commit researcher information.** The local files, payloads and all other
> researcher files in `data/` are ignored by Git; only the `*.example.*` files are
> tracked. Never force-add ignored files (`git add -f`), and do not put real
> information into the example files.

---

## The `data/` directory

All paths are relative to the application folder and the file names are fixed.

| File | Required | Used for | If missing |
|---|---|---|---|
| `data/home_info.txt` | **Yes** | Home page (researcher, project, description, team) | **The app does not start** |
| `data/experiment_design.xlsx` | **Yes** (for a complete instance) | Experiment design tab (sample table) | The app starts; the tab shows "Experiment design table not available" (or "Experiment design file could not be loaded." if the file is not a valid workbook) |
| `data/shiny_payload_rnaseq.rds` | No | RNAseq tab | Tab hidden |
| `data/shiny_payload_proteomics.rds` | No | Proteomics tab | Tab hidden |
| `data/shiny_payload_metabolomics.rds` | No | Metabolomics tab | Tab hidden |
| `data/shiny_payload_lipidomics.rds` | No | Lipidomics tab (module under development) | Tab hidden |
| `data/shiny_payload_multiomics.rds` | No | Reserved for the MultiOmics module (under development); not used by the current version | No effect |
| `data/Matching_genes_proteins.xlsx` | No | Gene–protein links for the MultiOmics module | The MultiOmics tab shows a "not available" message |
| `data/Data_&_Files/Data_&_Files_manifest.xlsx` | No | Data & Files tab (downloadable files) | "No Data & files are available right now." |
| `data/Data_&_Files/<files listed in the manifest>` | No | Downloads in the Data & Files tab | Entry shown as "Missing file" |

**When changes take effect:** `home_info.txt`, the payloads and
`Matching_genes_proteins.xlsx` are read once when the app starts, so restart the app
after changing them. `experiment_design.xlsx` and the Data & Files manifest are read
each time the app is opened in the browser, so reloading the page is enough.

### `home_info.txt` (required)

Plain text, one field per line, **field name, a TAB character, then the value**:

```text
Researcher<TAB>Researcher Name
Project<TAB>Project title
Description<TAB>Short description of the experiment and its aims.\nSecond paragraph of the description.
Bioinformatics_Team<TAB>vered,michal,oz,shelly,dor
```

- The four field names must stay exactly as shown.
- In `Description`, type the two characters `\n` where a new line should start.
- `Bioinformatics_Team` is a comma-separated list of team-member keys
  (`vered`, `michal`, `oz`, `shelly`, `or`, `dor`), shown in the given order.
  If it is empty, all members are shown.
- Use an editor that keeps the TAB characters (e.g. RStudio, Notepad++).

### `experiment_design.xlsx` (required)

The sample table shown in the **Experiment design** tab. Every researcher instance
should have one; prepare it when setting up a new researcher. Without it the app
still starts, but the tab only shows "Experiment design table not available".

- **Location:** directly inside `data/` (not in a sub-folder), named exactly
  `experiment_design.xlsx`.
- **Worksheet:** the app reads the **first worksheet**; its name does not matter.
- **Layout:** row 1 holds the column names; then one row per sample.
- **Descriptive columns** (for example `SampleID`, `Group`, `Replicate`): any names and
  any number of them. They are displayed as they are, in the same order.
- **Assay columns:** `RNAseq`, `Proteomics`, `Metabolomics`, `Lipidomics`, with
  exactly these names (case-sensitive).
  - Write `YES` (upper case) when the sample was measured with that omics; it is shown
    as a green check mark.
  - Leave the cell empty when it was not measured. Any other value is shown as an
    empty cell.
  - Include only the assay columns you need; any subset works. Each assay column that
    is present is shown with check marks. A misspelled assay column is treated as a
    descriptive column and shows its plain values.
- All sample rows are always shown, without paging.
- Use the same sample names as in the analysis, so the table matches the plots. (The
  app does not cross-check them.)
- Save it as a real Excel workbook (`.xlsx`). If the file cannot be read (for example
  a renamed CSV), the tab shows "Experiment design file could not be loaded.", the
  rest of the app keeps working, and the R console shows the reason. After fixing the
  file, reload the page.

The example file (`experiment_design.example.xlsx`) contains:

| SampleID | Group | Replicate | RNAseq | Proteomics | Metabolomics | Lipidomics |
|---|---|---|---|---|---|---|
| S01 | Control | 1 | YES | YES | YES | |
| S02 | Control | 2 | YES | YES | YES | |
| S03 | Control | 3 | YES | YES | YES | |
| S04 | Control | 4 | YES | YES | | |
| S05 | Treatment | 1 | YES | YES | YES | |
| S06 | Treatment | 2 | YES | YES | YES | |
| S07 | Treatment | 3 | YES | YES | YES | |
| S08 | Treatment | 4 | YES | YES | | |

### Omics payloads: standard file names

Payloads are produced by the upstream multi-omics pipelines and must use exactly these
names (no other names are recognised):

```text
shiny_payload_rnaseq.rds
shiny_payload_proteomics.rds
shiny_payload_metabolomics.rds
shiny_payload_lipidomics.rds
shiny_payload_multiomics.rds
```

### Which tabs are shown

- **Home, Experiment design, Data & Files:** always.
- **RNAseq / Proteomics / Metabolomics / Lipidomics:** when the matching payload file
  exists. Only the file's presence is checked; analyses missing *inside* a payload show
  their own "not available" messages.
- **MultiOmics:** when **at least two** of the four individual payloads exist.
  `shiny_payload_multiomics.rds` is not required and is not counted.

| Payloads present | Omics tabs shown |
|---|---|
| RNA-seq | RNAseq |
| RNA-seq + Proteomics | RNAseq, Proteomics, MultiOmics |
| Proteomics + Lipidomics | Proteomics, Lipidomics, MultiOmics |
| none | (none) |

These rules only decide whether a tab appears. They do not mean that the module
behind it is complete: **the MultiOmics and Lipidomics modules are currently under
development.**

### `Matching_genes_proteins.xlsx` (optional)

A table linking RNA-seq genes to Proteomics proteins, used by the MultiOmics module.
The app works without it.

- **Location:** directly inside `data/`, named exactly `Matching_genes_proteins.xlsx`.
  Create it by copying `data/Matching_genes_proteins.example.xlsx` (Quick start,
  step 1) and replacing the example rows.
- **Worksheet:** the app reads the **first worksheet**; its name does not matter.
- **Columns:** `genes` and `proteins`, with exactly these lower-case names. Other columns
  are not used.
- **Rows:** one gene–protein pair per row.
  - `genes` must use the same IDs as the RNA-seq payload.
  - `proteins` must use the same IDs as the Proteomics payload. A trailing version
    suffix such as `.1` or `.2` is ignored when matching, so `PROT0004.1` in the file
    matches `PROT0004` or `PROT0004.2` in the data.
- **When it is used:** it is read when the app starts, and kept only if the RNA-seq
  and Proteomics payloads are both present and at least one gene and one protein ID in
  the file match them. Otherwise the R console shows a warning, the file is ignored,
  and the MultiOmics tab shows "No relevant multiomics matching information is
  available for the currently loaded datasets."

The example file uses synthetic IDs that do not match any real data:

| genes | proteins |
|---|---|
| EXGENE_0001 | EXPROT0001 |
| EXGENE_0002 | EXPROT0002 |
| EXGENE_0003 | EXPROT0003 |
| EXGENE_0004 | EXPROT0004.1 |
| EXGENE_0005 | EXPROT0005 |
| EXGENE_0006 | EXPROT0006 |

### Data & Files manifest (optional)

`data/Data_&_Files/Data_&_Files_manifest.xlsx` lists the files offered for download in
the **Data & Files** tab. The app reads the **first worksheet**; column names are not
case-sensitive. One row per file:

| Column | Meaning |
|---|---|
| `title` | Link text (rows without a title are skipped) |
| `file_path` | **Where the actual file is**: its path relative to `data/Data_&_Files/`, including the file name and extension (not a descriptive title) |
| `description` | Optional text shown above the link |
| `category` | Optional section heading (rows with the same category are grouped) |
| `updated_at` | Optional date (Excel date or `DD/MM/YYYY`) |

**How `file_path` is resolved:** the app joins `data/Data_&_Files/` and the value of
`file_path`, and offers that file for download under its own file name. Examples:

| File on disk | `file_path` |
|---|---|
| `data/Data_&_Files/methods.pdf` | `methods.pdf` |
| `data/Data_&_Files/reports/analysis_report.html` | `reports/analysis_report.html` |
| `data/Data_&_Files/presentations/project_overview.pptx` | `presentations/project_overview.pptx` |

- Sub-folders are supported. Use forward slashes (`/`); backslashes work only on
  Windows.
- Match the capitalisation of folder and file names exactly. Windows ignores case,
  but other systems do not.
- Keep all files inside `data/Data_&_Files/`. Absolute paths (such as `C:/...`) are not
  supported, and do not use `..` to point outside the folder.
- If `file_path` is empty or does not point to an existing file, the entry is shown in
  grey with "Missing file" and cannot be downloaded.

The example manifest uses the paths above (plus `tables/supplementary_tables.xlsx`).
Its files are not included, so the entries appear as "Missing file" until you add the
files or replace the rows in your local copy.

---

## Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| App fails at startup with `cannot open file 'data/home_info.txt'` | `data/home_info.txt` has not been created. Copy it from `data/home_info.example.txt` (Quick start, step 1) |
| App fails at startup with `there is no package called ...` | Install the missing package (see [Requirements](#requirements)) |
| Home page shows empty researcher/project | Field names in `home_info.txt` changed, or spaces used instead of a TAB |
| Experiment design tab shows "Experiment design table not available" | `data/experiment_design.xlsx` is missing or misnamed |
| Experiment design tab shows "Experiment design file could not be loaded." | `data/experiment_design.xlsx` is not a valid Excel workbook (the R console shows the reason). Save it again as `.xlsx` and reload the page |
| Experiment design tab shows "Experiment design table is empty." | The first worksheet of `data/experiment_design.xlsx` has no content |
| An assay column in the Experiment design table has no check marks | The column name is not exactly `RNAseq`, `Proteomics`, `Metabolomics` or `Lipidomics` (it is then shown with its plain values), or the cells do not contain exactly `YES` |
| An omics tab is missing | Payload file name is not exactly `shiny_payload_<omics>.rds`, or it is not in `data/` |
| MultiOmics tab is missing | Fewer than two of the four individual payloads are present |
| MultiOmics tab shows "No relevant multiomics matching information…" | `Matching_genes_proteins.xlsx` is missing, its columns are not named `genes` and `proteins`, or its IDs do not match the RNA-seq and Proteomics payloads (see the warning in the R console) |
| A plot shows an error mentioning `DESeq2` or `limma` | Install the Bioconductor packages (see [Requirements](#requirements)) |
| Data & Files entry shows "Missing file" | `file_path` does not point to an existing file under `data/Data_&_Files/`. Check folders, spelling, capitalisation and forward slashes |
| Changes in `data/` are not visible | Restart the app (`home_info.txt`, payloads, matching table) or reload the browser page (experiment design, Data & Files manifest) |
