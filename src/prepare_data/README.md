


run_spliser.py --> aggregate_sse.py --> make_splice_table.py --> ssu_to_splice_table.py --> 



# `run_spliser.py`

## Overview
The `run_spliser.py` script automates the SpliSER pipeline, a tool for analyzing splice site expression using RNA-Seq data. It processes BAM files to extract junction information, runs SpliSER on individual samples.

---

## Features
1. **Junction Extraction**:
   - Uses `regtools` to extract junction data from BAM files.
2. **SpliSER Analysis**:
   - Runs the SpliSER process for each sample to calculate splice site expression.
3. **Combination**:
   - Merges SpliSER outputs into a consolidated file.
4. **Differential Analysis**:
   - Runs DiffSpliSER for detecting differential splice site usage.
5. **Parallel Processing**:
   - Processes multiple BAM files in parallel for efficiency.
6. **Logging**:
   - Logs pipeline progress and errors to a file for easy troubleshooting.

---

## Input Requirements
1. **Configuration File**:
   - A JSON file (`spliser_config.json`) defining the pipeline parameters.
2. **BAM Files**:
   - A TSV file (`bam_list.tsv`) listing the paths to BAM files.
3. **Gene Annotation File**:
   - A GTF file (e.g., `gencode.v46.primary_assembly.annotation.gtf`) for transcript annotation.
4. **SpliSER Python Script**:
   - The main SpliSER script (`spliser.py`) to process, combine, and output results.

---

## Output Files
1. **Junction Files**:
   - Created by `regtools`, stored in `junction_files_dir`.
2. **SpliSER Process Outputs**:
   - Per-sample SpliSER results stored in `process_output_dir`.
3. **Combined Results**:
   - Merged SpliSER results saved in `output_combined`.
4. **DiffSpliSER Results**:
   - Differential splice site expression outputs stored in `output_diffspliser`.
5. **Samples File**:
   - Tracks processed samples, outputs, and BAM file paths.

---

## Usage

### Step 1: Prepare Input Files
1. Create a JSON configuration file (`spliser_config.json`) with all necessary paths and parameters.
2. Ensure `bam_list.tsv` lists the paths to the BAM files.
3. Verify that the GTF file and `spliser.py` script are accessible.

### Step 2: Execute the Script
Run the script with the configuration file:
```bash
python run_spliser.py --config spliser_config.json
```



```markdown
## Configuration File: `spliser_config.json`

The configuration file defines paths, parameters, and behavior of the SpliSER pipeline. Below is an explanation of its fields.

### Fields and Descriptions

1. **`bam_file_list`**:
   - Path to the TSV file listing BAM files to process.
   - Example: `"./bam_list.tsv"`

2. **`spliser_script`**:
   - Path to the SpliSER Python script.
   - Example: `"../SpliSER/SpliSER_v0_1_8.py"`

3. **Directories**:
   - **`spliser_dir`**: Main output directory for pipeline results.
   - **`junction_files_dir`**: Directory to store junction files from `regtools`.
   - **`process_output_dir`**: Directory for SpliSER outputs.
   - **`samples_file`**: Path to a TSV file tracking processed samples.
   - **`output_combined`**: Directory for combined results.
   - **`output_diffspliser`**: Directory for DiffSpliSER outputs.

4. **`gtf_file`**:
   - Path to the GTF annotation file.
   - Example: `"../GRCh38/gencode.v46.primary_assembly.annotation.gtf"`

5. **`n_threads`**:
   - Number of threads for parallel processing.
   - Example: `60`

6. **Processing Flags**:
   - **`skip_processing`**: If `true`, skips BAM processing and only runs the combine/output steps.
   - **`only_process`**: If `true`, runs only the processing step and skips combine/output steps.

7. **Command Arguments**:
   - **`regtools_args`**:
     - Arguments for `regtools` junction extraction.
     - Example: `["-a", "8", "-m", "50", "-M", "500000", "-s", "RF"]`
   - **`spliser_process_args`**:
     - Arguments for the SpliSER `process` step.
     - Example: `["--isStranded", "-s", "rf", "-m", "500000", "--beta2Cryptic"]`
   - **`spliser_combine_args`**:
     - Arguments for the SpliSER `combine` step.
     - Example: `["--isStranded", "-s", "rf", "--n_jobs", "30"]`
   - **`spliser_output_args`**:
     - Arguments for the SpliSER `output` step.
     - Example: `[]`

---

### Example Configuration
```json
{
    "bam_file_list": "./bam_list.tsv",
    "spliser_script": "../SpliSER/SpliSER_v0_1_8.py",
    "spliser_dir": "./output",
    "junction_files_dir": "./output/junction_files",
    "samples_file": "./output/samples_file.tsv",
    "process_output_dir": "./output/process",
    "output_combined": "./output/spliser_combined",
    "output_diffspliser": "./output/spliser_output",
    "gtf_file": "../GRCh38/gencode.v46.primary_assembly.annotation.gtf",
    "n_threads": 60,
    "skip_processing": false,
    "only_process": true,
    "regtools_args": ["-a", "8", "-m", "50", "-M", "500000", "-s", "RF"],
    "spliser_process_args": ["--isStranded", "-s", "rf", "-m", "500000", "--beta2Cryptic"],
    "spliser_combine_args": ["--isStranded", "-s", "rf", "--n_jobs", "30"],
    "spliser_output_args": []
}
```



# `aggregate_sse.py`

This script processes SpliSER output `.tsv` files to aggregate Splice Site Efficiency (SSE) values for splice sites, annotate sites as exon starts or ends, and filter out sites that do not meet specific criteria. It is designed to handle parallel processing of multiple files for efficient data aggregation.

## Overview

In each SpliSER `.tsv` file, the script extracts the following details for each splice site:
- **Region** (chromosome),
- **Site** (position),
- **Strand**,
- **Partners** (splice partners),
- **SSE** (Splice Site Efficiency).

### Key Steps:
1. **Annotation**: Each site is annotated as:
   - **Exon end**: If all partners are downstream.
   - **Exon start**: If all partners are upstream.
   - **Skipped**: If the site has both upstream and downstream partners (e.g., indicative of circular RNA).
2. **Parallel Processing**: Each `.tsv` file is processed independently using multiple threads.
3. **Aggregation**: The script combines the results from all files, calculates the average SSE for each site across samples, and applies filters based on user-defined thresholds.
4. **Output**:
   - Aggregated SSE values.
   - Skipped sites with both upstream and downstream partners.
5. **Statistics**: Provides statistics on sites filtered by SSE and occurrence thresholds.

## Input Arguments

### Required Arguments:
- `--input_dir`: Path to the directory containing SpliSER `.tsv` files.
- `--output_file`: Path to the output file for aggregated SSE values.
- `--skipped_output_file`: Path to the output file for skipped sites (both upstream and downstream partners).
- `--log_file`: Path to the log file for entries with SSE = 0.

### Optional Arguments:
- `--n_threads`: Number of threads for parallel processing (default: 4).
- `--min_alpha_count`: Minimum alpha count required to include a site (default: 5).
- `--min_file_count`: Minimum number of files in which a site must appear to be included in the aggregated results (default: 2).

## Outputs

1. **Aggregated SSE File**:
   - Contains sites that passed the filters with average SSE values.
   - Format: Tab-separated file with columns:
     - `Region`, `Site`, `Strand`, `Average_SSE`, `Annotation`.
2. **Skipped Sites File**:
   - Lists sites excluded due to having both upstream and downstream partners.
   - Format: Tab-separated file with columns:
     - `Region`, `Site`, `Strand`, `File`, `Reason`.
3. **Log File**:
   - Records sites with SSE = 0 for further investigation.

## Example Usage

### Example 1: Processing K562 Data
```
python aggregate_sse.py \
    --input_dir ../../k562 \
    --output_file k562_filt_aggregated_ssu.tsv \
    --skipped_output_file k562_up_and_down_partners.tsv \
    --log_file zero_sse_log.txt \
    --n_threads 56 \
    --min_alpha_count 5 \
    --min_file_count 2
```

# `make_splice_table.py`

This script processes a GTF file to extract detailed transcript and exon information, identify primary transcripts, annotate splice sites, and determine paralog status for genes. It generates a splice site table and a gene-to-transcript mapping file, providing comprehensive data for downstream splicing analysis.

## Overview

The script performs the following steps:
1. **Read and Parse GTF File**:
   - Extract attributes such as `gene_id`, `transcript_id`, `gene_name`, `gene_type`, and `tag`.
   - Filter for protein-coding genes and transcripts.

2. **Select Primary Transcripts**:
   - Prioritizes transcripts using tags like `MANE_Select`, `Ensembl_canonical`, and `basic`.
   - Retains one primary transcript per gene.

3. **Process Exons**:
   - Identifies and processes exon start and end sites for primary transcripts.
   - Excludes single-exon transcripts where start and end sites match transcript boundaries.

4. **Determine Paralog Status**:
   - Checks each gene against a local paralog file (`paralogs_GRCh38.txt`) to determine if it has a paralog.

5. **Filter by Chromosome**:
   - Includes only standard chromosomes (`chr1`-`chr22`, `chrX`, and `chrY`).

6. **Generate Outputs**:
   - Saves a splice site table with detailed information about transcripts and exons.
   - Saves a gene-to-transcript mapping file.

## Input Arguments

### Required Arguments:
- `-g`, `--gtf_file`: Path to the gzipped GTF file (e.g., `gencode.v29.primary_assembly.annotation.gtf.gz`).
- `-o`, `--output_file`: Path to save the splice site table.
- `-m`, `--gene_transcript_output_file`: Path to save the gene-to-transcript mapping file.

## Outputs

1. **Splice Site Table**:
   - Format: Tab-separated file with columns:
     - `gene_id`: Gene identifier.
     - `paralog_status`: `1` if the gene has a paralog, `0` otherwise.
     - `seqname`: Chromosome name.
     - `strand`: Strand (`+` or `-`).
     - `transcript_start`: Start position of the transcript.
     - `transcript_end`: End position of the transcript.
     - `exon_end_sites`: Comma-separated exon end positions.
     - `exon_start_sites`: Comma-separated exon start positions.

2. **Gene-to-Transcript Mapping File**:
   - Format: Tab-separated file with columns:
     - `gene_id`: Gene identifier.
     - `transcript_id`: Transcript identifier.
     - `gene_name`: Gene name.
     - `selected_tag`: Tag used to select the primary transcript (e.g., `MANE_Select`, `basic`).

## Example Usage

```
python make_splice_table.py \
    -g ../../GRCh38/gencode.v29.primary_assembly.annotation.gtf.gz \
    -o ../../splice_tables/GRCh38_v29_splice_table.txt \
    -m gene_2_tx.tsv
```



# `ssu_to_splice_table.py`

This script integrates aggregated Splice Site Efficiency (SSE) data into a splice table, matching sites based on chromosome, strand, and transcript start and end coordinates. It updates the splice table by appending SSE values for exon start and end positions, allowing for enhanced analysis of RNA sequencing data.

## Overview

The script processes each row in the splice table by:
1. Extracting the chromosome, strand, transcription start, and transcription end sites.
2. Filtering the aggregated SSE file for rows that match the same chromosome and strand, and are within the transcript's start and end positions.
3. Identifying and appending:
   - Exon end positions and their corresponding SSE values.
   - Exon start positions and their corresponding SSE values.
4. Optionally filtering to retain only exon positions that match annotated positions in the splice table.

## Input Arguments

### Required Arguments:
- `--splice-table`: Path to the splice table file, which must contain 8 columns:
  1. **Transcript ID**: Identifier for the transcript.
  2. **Paralog Status**: `0` for no paralog, `1` for paralog.
  3. **Chromosome**: Chromosome name (e.g., `chr1`).
  4. **Strand**: Strand (`+` or `-`).
  5. **Transcription Start Site**: Start site of transcription (1-based).
  6. **Transcription End Site**: End site of transcription (1-based).
  7. **Exon End Sites**: Comma-separated exon end positions (1-based, trailing comma).
  8. **Exon Start Sites**: Comma-separated exon start positions (1-based, trailing comma).

- `--aggregated-sse`: Path to the `.tsv` file created by `aggregate_sse.py`, containing:
  - Chromosome, strand, site positions, annotations, and average SSE values.

- `--output`: Path to save the updated splice table.

### Optional Arguments:
- `--annotated-only`: If set, filters exon positions and SSE values to include only those that match annotated positions in the splice table.

## Output

The script produces an updated splice table that includes four new columns appended to the input splice table:
1. **Exon End Positions**: Comma-separated positions of exon ends that match aggregated SSE data.
2. **Exon End SSE Values**: Corresponding SSE values for the exon end positions.
3. **Exon Start Positions**: Comma-separated positions of exon starts that match aggregated SSE data.
4. **Exon Start SSE Values**: Corresponding SSE values for the exon start positions.

## Example Usage

### Without Filtering to Annotated Sites
```bash
python ssu_to_splice_table.py \
    -s ../../splice_tables/GRCh38_v29_splice_table.txt \
    -a ./hepg2_filt_aggregated_ssu.tsv \
    -o ../../splice_tables/hepg2_filt_rna_seq_splice_table.txt




    