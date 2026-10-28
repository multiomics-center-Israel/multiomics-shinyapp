# multiomics-shinyapp

A reusable R Shiny application for presenting multi-omics analysis results
(RNA-seq, Proteomics, Metabolomics, Lipidomics, MultiOmics) to a researcher.

The same code is used for every researcher. Everything researcher-specific lives in
`data/`: a short text file with the researcher's details, the analysis payloads
(`.rds`) produced by the upstream pipelines, and optional supporting files.
**No R code changes are needed to create a new researcher instance.**

---

## Repository structure

```text
multiomics-shinyapp/
├── global.R              # packages, settings, data loading (shared by ui.R and server.R)
├── ui.R                  # user interface
├── server.R              # server logic + Appstream authentication
├── functions_shiny.R     # analysis and plotting helpers
├── www/                  # generic images used by the interface
├── data/
│   ├── home_info.example.txt                     # example (tracked)
│   └── Data_&_Files/
│       └── Data_&_Files_manifest.example.xlsx    # example (tracked)
├── .gitignore
└── README.md
```

The app runs in Shiny's two-file mode (`ui.R` + `server.R`, with `global.R` loaded
first). There is intentionally no `app.R`.

### Shared template files vs. local researcher files

| Shared template files (in Git — do not edit per researcher) | Local researcher files (never in Git) |
|---|---|
| `global.R`, `ui.R`, `server.R`, `functions_shiny.R` | `data/home_info.txt` (created from the example) |
| `www/` | `data/Data_&_Files/Data_&_Files_manifest.xlsx` (created from the example) |
| `data/home_info.example.txt` | `data/shiny_payload_*.rds` |
| `data/Data_&_Files/Data_&_Files_manifest.example.xlsx` | `data/experiment_design.xlsx`, `data/Matching_genes_proteins.xlsx` |
| `README.md`, `.gitignore` | Files in `data/Data_&_Files/` listed in the manifest |

Everything in the right-hand column is ignored by `.gitignore`, so a researcher's
information is not picked up by `git add`. Still, always check `git status` before
committing.

---

## Requirements

- **R ≥ 4.1** (developed and tested with R 4.5).
- **CRAN packages:**

  ```r
  install.packages(c(
    "shiny", "shinythemes", "shinycssloaders", "bslib", "DT", "htmltools",
    "ggplot2", "plotly", "heatmaply", "pheatmap", "RColorBrewer", "viridisLite",
    "scales", "ggridges", "gridExtra", "patchwork",
    "dplyr", "stringr", "tibble", "jsonlite",
    "readxl", "writexl", "openxlsx",
    "bcrypt"
  ))
  ```

  (`patchwork` is optional: without it, one enrichment plot is shown without its
  extra panels.)

- **Bioconductor packages** — needed only if the payloads contain these model objects
  (DESeq2 results for RNA-seq, limma fits):

  ```r
  install.packages("BiocManager")
  BiocManager::install(c("DESeq2", "limma"))
  ```

---

## Quick start

```bash
git clone https://github.com/multiomics-center-Israel/multiomics-shinyapp.git
```

1. Create the two local runtime files from the examples. The commands below never
   overwrite an existing file. Run them from the repository folder:

   ```r
   # R (any OS) — returns FALSE instead of overwriting if the file already exists
   file.copy("data/home_info.example.txt", "data/home_info.txt", overwrite = FALSE)
   file.copy("data/Data_&_Files/Data_&_Files_manifest.example.xlsx",
             "data/Data_&_Files/Data_&_Files_manifest.xlsx", overwrite = FALSE)
   ```

   ```bash
   # Linux / macOS / Git Bash (-n = do not overwrite)
   cp -n data/home_info.example.txt data/home_info.txt
   cp -n "data/Data_&_Files/Data_&_Files_manifest.example.xlsx" "data/Data_&_Files/Data_&_Files_manifest.xlsx"
   ```

2. Edit `data/home_info.txt` with the researcher's details (see below).
3. Copy the researcher's payload files into `data/` using the standard names (below).
4. Optionally add the supporting files, edit `data/Data_&_Files/Data_&_Files_manifest.xlsx`
   and put the listed files next to it (or delete the manifest if not needed).
5. Run locally from R / RStudio:

   ```r
   shiny::runApp("path/to/multiomics-shinyapp")
   ```

   The R console should print:

   ```text
   Appstream authentication: OFF (local run)
   Omics tabs shown: rnaseq, proteomics, multiomics      # depends on your payloads
   ```

> **Do not commit researcher information.** The runtime files, payloads and all other
> researcher files in `data/` are ignored by Git; only the `*.example.*` files are
> tracked. Never force-add ignored files (`git add -f`), and do not put real
> information into the example files.

---

## The `data/` directory

All paths are relative to the application folder and the file names are fixed.

| File | Required | Used for | If missing |
|---|---|---|---|
| `data/home_info.txt` (copy of `home_info.example.txt`) | **Yes** | Home page (researcher, project, description, team) | **The app does not start** |
| `data/shiny_payload_rnaseq.rds` | No | RNA-seq tab | Tab hidden |
| `data/shiny_payload_proteomics.rds` | No | Proteomics tab | Tab hidden |
| `data/shiny_payload_metabolomics.rds` | No | Metabolomics tab | Tab hidden |
| `data/shiny_payload_lipidomics.rds` | No | Lipidomics tab (currently a placeholder page) | Tab hidden |
| `data/shiny_payload_multiomics.rds` | No | Reserved for a future MultiOmics payload (not read yet) | No effect |
| `data/experiment_design.xlsx` | No | Experiment design tab (sample table) | "Experiment design table not available" |
| `data/Matching_genes_proteins.xlsx` | No | MultiOmics gene–protein view | MultiOmics shows a "not available" message |
| `data/Data_&_Files/Data_&_Files_manifest.xlsx` (copy of the `.example.xlsx`) | No | Data & Files tab (downloadable documents) | "No Data & files are available right now." |
| `data/Data_&_Files/<files listed in the manifest>` | No | Downloads in the Data & Files tab | Entry shown as "Missing file" |

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
  (`vered`, `michal`, `oz`, `shelly`, `dor`), shown in the given order.
  If it is empty, all members are shown.
- Use an editor that keeps the TAB characters (e.g. RStudio, Notepad++).

### Omics payloads — standard file names

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

- **Home, Experiment design, Data & Files** — always.
- **RNAseq / Proteomics / Metabolomics / Lipidomics** — when the matching payload file
  exists. Only the file's presence is checked; analyses missing *inside* a payload show
  their own "not available" messages.
- **MultiOmics** — when **at least two** of the four individual payloads exist.
  `shiny_payload_multiomics.rds` is not required and is not counted.

| Payloads present | Omics tabs shown |
|---|---|
| RNA-seq | RNAseq |
| RNA-seq + Proteomics | RNAseq, Proteomics, MultiOmics |
| Proteomics + Lipidomics | Proteomics, Lipidomics, MultiOmics |
| none | (none) |

### Data & Files manifest (optional)

`data/Data_&_Files/Data_&_Files_manifest.xlsx`, sheet `documentation_manifest`,
one row per downloadable file:

| Column | Meaning |
|---|---|
| `title` | Link text (rows without a title are skipped) |
| `file_path` | File name inside `data/Data_&_Files/` |
| `description` | Optional text shown above the link |
| `category` | Optional section heading (rows with the same category are grouped) |
| `updated_at` | Optional date (Excel date or `DD/MM/YYYY`) |

The example rows are placeholders; their files are not included, so they appear as
"Missing file" until you add the files or replace the rows in your local copy.

### Supporting Excel files (optional)

- `experiment_design.xlsx` — sample table shown as-is. If it has the columns
  `RNAseq`, `Proteomics`, `Metabolomics`, `Lipidomics`, cells containing `YES` are shown
  as green check marks.
- `Matching_genes_proteins.xlsx` — columns `genes` and `proteins`; used by MultiOmics
  only when its IDs match the loaded RNA-seq and Proteomics data.

---

## Deployment on Shiny Server and Appstream authentication

1. Copy the whole application folder (code + `www/` + the researcher's `data/`) to
   `/srv/shiny-server/<Researcher>/<App>/`.
   **Use only ASCII letters, digits, `_` and `-` in `<Researcher>` and `<App>`**
   (no spaces, Hebrew or other characters) — otherwise authentication cannot match the
   application path.
2. Make sure the Shiny Server user can read the folder, then restart the app
   (`touch restart.txt` inside the folder). Payloads are loaded once per app process,
   so restart again after changing files in `data/`.
3. Give access through the existing Appstream portal: add the app URL
   `http://<server>:<port>/<Researcher>/<App>/` to each authorised user's Appstream
   permissions, in the existing format.

How authentication works:

- On Shiny Server the app requires an Appstream login. It detects its own path
  (`<Researcher>/<App>`) from its folder, so no settings need to be edited.
- A visitor without valid Appstream credentials for this app sees "Access denied" and
  never receives any researcher content (the Home details are sent only after login).
  Other users are not affected.
- When run locally (RStudio / `shiny::runApp`) authentication is off.

The server log shows `Appstream authentication: ON for <Researcher>/<App> (...)` at
startup and one line per denied access with the reason.

---

## Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| App fails at startup with `cannot open file 'data/home_info.txt'` | `data/home_info.txt` has not been created — copy it from `data/home_info.example.txt` (Quick start, step 1) |
| Home page shows empty researcher/project | Field names in `home_info.txt` changed, or spaces used instead of a TAB |
| An omics tab is missing | Payload file name is not exactly `shiny_payload_<omics>.rds`, or it is not in `data/` |
| MultiOmics tab is missing | Fewer than two of the four individual payloads are present |
| `there is no package called ...` | Install the package (see Requirements) |
| Error mentioning `DESeq2` or `limma` | Install the Bioconductor packages |
| Data & Files entry shows "Missing file" | The file named in `file_path` is not in `data/Data_&_Files/` |
| On the server everyone gets "Access denied" | Folder name contains spaces/non-ASCII characters, the app URL is not in the user's Appstream permissions, or the user did not open the app through Appstream |
| Changes in `data/` are not visible on the server | Restart the app: `touch restart.txt` |
