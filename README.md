# Bacterial Genomics Pipeline

Snakemake pipeline for whole-genome sequencing of clinical bacterial isolates — trimming, assembly, species identification, resistance, typing, and SNP phylogenomics.

---

## Requirements

- Linux (tested on Ubuntu/Fedora)
- [Pixi](https://pixi.sh) installed (`curl -fsSL https://pixi.sh/install.sh | bash`)

---

## Installation

```bash
git clone https://github.com/OddAlexander/Genomics-Workflow.git ~/genomics
cd ~/genomics
pixi install
```

### Databases

**GAMBIT**
```bash
pixi run --environment identification gambit-db-genomes download -d /databases/gambit_db/
```

**Kraken2** — recommended: PlusPF 8 GB (~8 GB, requires ≥16 GB RAM):
```bash
mkdir -p /databases/kraken2_db_light
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_8gb_20240904.tar.gz -P /databases/kraken2_db_light/
tar -xzf /databases/kraken2_db_light/k2_pluspf_8gb_20240904.tar.gz -C /databases/kraken2_db_light/
```

**skani** (~3 GB sketch, ~25 GB temporary during build):
```bash
sudo mkdir -p /databases/fastani_db /databases/skani_db && sudo chown $USER /databases/fastani_db /databases/skani_db
pixi run --environment identification datasets download genome taxon bacteria \
    --reference --assembly-source refseq --assembly-level complete \
    --include genome --filename /databases/fastani_db/bacteria_refs.zip
unzip /databases/fastani_db/bacteria_refs.zip -d /databases/fastani_db/
find /databases/fastani_db/ncbi_dataset/data -name "*.fna" > /databases/fastani_db/reference_list.txt
pixi run --environment identification skani sketch \
    -l /databases/fastani_db/reference_list.txt -o /databases/skani_db/bacteria -t 16
rm -rf /databases/fastani_db/ncbi_dataset /databases/fastani_db/bacteria_refs.zip
```

**MOB-suite** (~1.5 GB): `pixi run --environment mobsuite mob_init`

**AMRFinder**: `pixi run --environment amrfinder4 amrfinder --update`

---

## Usage

### Input format

```
data/
└── 26-03-2026/
    ├── 001k/
    │   ├── 001k_R1.fastq.gz
    │   └── 001k_R2.fastq.gz
    └── 002k/
        ├── 002k_R1.fastq.gz
        └── 002k_R2.fastq.gz
```

Filenames can be anything, but must end with `_R1.fastq.gz` / `_R2.fastq.gz`.

### Main pipeline (`Snakefile`)

```bash
# Dry run
pixi run snakemake -n --cores 16

# One sample / one date / multiple samples / all
pixi run snakemake --cores 16 --resources mem_mb=60000 --config samples=001k
pixi run snakemake --cores 16 --resources mem_mb=60000 --config samples=19-03-2026
pixi run snakemake --cores 16 --resources mem_mb=60000 --config "samples=[001k,002k]"
pixi run snakemake --cores 16 --resources mem_mb=60000

# One specific file
pixi run snakemake results/26-03-2026/001k/MLST/mlst.tsv
```

> **`mem_mb`** should be set to total RAM minus ~4 GB (e.g. `12000` on a 16 GB machine). Kraken2, Shovill, and skani reserve 12, 12, and 4 GB respectively.

### Phylogenomics pipeline (`Snakefile_phylo`)

Runs independently of the main pipeline — takes raw reads directly as input.

```bash
# With an external reference (.gbk from Prokka — recommended)
pixi run snakemake -s Snakefile_phylo --cores 16 \
    --config ref=/databases/prokka_refs/19-03-2026_005a/19-03-2026_005a.gbk

# Without reference — first sample is auto-annotated with Prokka
pixi run snakemake -s Snakefile_phylo --cores 16

# Filter by date or sample ID
pixi run snakemake -s Snakefile_phylo --cores 16 \
    --config ref=... samples=19-03-2026
pixi run snakemake -s Snakefile_phylo --cores 16 \
    --config ref=... "samples=[19-03-2026/005a,19-03-2026/007b]"

# Custom run name (default: date+time)
pixi run snakemake -s Snakefile_phylo --cores 16 \
    --config ref=... run_name=salmonella_outbreak

# SNP threshold for outbreak cluster coloring in the report
pixi run snakemake -s Snakefile_phylo --cores 16 \
    --config ref=... snp_threshold=20
```

Output is written to `results_phylo/<run_name>/` (default: `DD-MM-YYYY_HHMM`):

```
results_phylo/07-05-2026_1621/
├── <date>/<sample-id>/
│   ├── Trimmed/      fastp
│   ├── QC/           FastQC + fastp JSON + HTML
│   ├── Assembly/     Shovill
│   ├── QUAST/
│   ├── ID_Skani/     species identification (species.txt)
│   ├── Prokka/       only for auto-reference (first sample)
│   └── MultiQC/
├── Snippy/
│   └── <sample-id>/  snippy SNP calling against reference
├── Core/             snippy-core — core.full.aln
├── Gubbins/          recombination-filtered alignment + tree (hybrid: FastTree/RAxML)
│                     gubbins.recombination_predictions.gff
├── SNP_Dists/        snp_dists.tsv (raw core genome) + snp_dists_gubbins.tsv
├── IQtree/           iqtree.treefile (GTR+G, UFBoot 1000; skipped if < 3 taxa)
├── MultiQC_run/      multiqc_report.html (aggregated QC for the full run)
└── phylo_report.html sample overview, reference, SNP matrices, recombination summary,
                      tree download, and tool versions
```

The phylo report SNP matrices support threshold-based cluster coloring when `snp_threshold` is set: cells ≤ threshold are shaded green (same cluster), cells > threshold are shaded red (different cluster).

#### Creating custom references

```bash
scripts/pixi/pixi_annotate_prokka.sh 19-03-2026
scripts/pixi/pixi_annotate_prokka.sh 19-03-2026/005a --genus Staphylococcus --species aureus
```

---

## Pipeline overview

### Main pipeline

```
FASTQ
  │
  ▼
fastp → FastQC
  │
  ├── Kraken2 + Bracken
  ▼
Shovill → QUAST
  │
  ├── GAMBIT → species.txt  ← checkpoint (controls DAG routing)
  ├── skani  (only on uncertain GAMBIT match)
  ├── MLST
  ├── AMRFinder
  ├── MOB-suite
  ├── MEfinder
  ├── [S. aureus]          → StaphScope
  ├── [Klebsiella/E. coli] → Kleborate
  ├── [S. pyogenes]        → emmtyper
  ├── [P. aeruginosa]      → Pasty
  ├── [Salmonella]         → SeqSero2
  └── [H. influenzae]      → hicap
  │
  ▼
MultiQC → HTML report
```

### Phylogenomics pipeline

```
FASTQ (raw reads)
  │
  ▼
fastp → FastQC → Shovill → QUAST → skani
  │                │
  │                └── Prokka  (auto-reference only)
  ▼
Snippy (per sample, against reference)
  │
  ▼
snippy-core → core.full.aln
  │                │
  │                └── snp-dists  → SNP distance matrix (raw)
  ▼
Gubbins (recombination removal, hybrid tree builder)
  │                │
  │                └── snp-dists  → SNP distance matrix (filtered)
  └── IQ-TREE2  → ML tree  (skipped if < 3 taxa)
  │
  ▼
MultiQC + HTML report (phylo_report.html)
```

---

## Species-specific typing

| Species | Tool | Output |
|---------|------|--------|
| *S. aureus* | StaphScope | spa type, SCCmec, MRSA/MSSA, virulence, lineage |
| *Klebsiella* spp. + *E. coli* | Kleborate v3 | K/O type, ESBL, carbapenemase, virulence |
| *S. pyogenes* | emmtyper | emm type |
| *P. aeruginosa* | Pasty | O-antigen serotype |
| *Salmonella* spp. | SeqSero2 | O and H antigen serotyping |
| *H. influenzae* | hicap | Capsule type (a–f / non-typeable) |
| All | MEfinder | Insertion sequences and transposons |
| All | MOB-suite | Plasmid typing (relaxase, conjugation) |

ECTyper (*E. coli* O:H serotyping) requires a MASH sketch from Zenodo:

```bash
mkdir -p databases/ectyper
wget -O databases/ectyper/EnteroRef_GTDBSketch_20231003_V2.msh \
    "https://zenodo.org/records/13969103/files/EnteroRef_GTDBSketch_20231003_V2.msh?download=1"
```

---

## Configuration (`config.yaml`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `data_dir` | `data/` | Directory containing FASTQ files |
| `results_dir` | `results/` | Output directory |
| `threads` | `8` | Threads per rule |
| `kraken2_db` | `/databases/kraken2_db_light/` | Kraken2 database |
| `kraken2_mem_mb` | `12000` | MB RAM per Kraken2 job |
| `gambit_db` | `/databases/gambit_db/` | GAMBIT database |
| `skani_db` | `/databases/skani_db/bacteria` | skani sketch database |
| `skani_threshold` | `0.3` | GAMBIT `closest.distance` above this triggers skani |
| `skani_mem_mb` | `4000` | MB RAM per skani job |
| `skani_threads` | `8` | Threads per skani job |
| `shovill_mem_mb` | `12000` | MB RAM per Shovill job |

**Phylogenomics-only options** (passed via `--config`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ref` | *(first sample)* | Path to reference `.gbk` (Prokka-annotated) |
| `run_name` | `DD-MM-YYYY_HHMM` | Output subdirectory name under `results_phylo/` |
| `snp_threshold` | `0` (disabled) | SNP cutoff for cluster coloring in the report |
| `mincov` | `10` | Snippy minimum read coverage |
| `minfrac` | `0.9` | Snippy minimum allele frequency |

---

## Troubleshooting

```bash
# View log
cat results/26-03-2026/001k/logs/amrfinder.log

# Force rerun of one rule
pixi run snakemake --cores 16 --forcerun amrfinder

# Resume an interrupted run
pixi run snakemake --cores 16 --resources mem_mb=60000 --rerun-incomplete
```

---

## Tools

| Tool | Version | Purpose |
|------|---------|---------|
| fastp | ≥0.23 | Read trimming |
| FastQC + MultiQC | ≥0.12 | Sequence quality |
| Shovill | ≥1.1 | Assembly (SPAdes) |
| QUAST | ≥5.2 | Assembly statistics |
| Kraken2 + Bracken | ≥2.1 | Species identification (reads) |
| GAMBIT | ≥0.5 | Species identification (assembly) — controls DAG routing |
| skani | ≥0.2 | ANI-based species ID when GAMBIT match is uncertain |
| MLST | ≥2.23 | Sequence typing (PubMLST) |
| AMRFinder | v4.x | Resistance and virulence genes |
| MOB-suite | ≥3.1 | Plasmid typing |
| MobileElementFinder | ≥1.0 | MGE detection |
| StaphScope | ≥1.0 | *S. aureus*: spa, SCCmec, MRSA, virulence, lineage |
| emmtyper | ≥0.2 | *S. pyogenes* emm typing |
| Kleborate | ≥3.0 | *Klebsiella* / *E. coli* typing |
| Pasty | ≥2.2 | *P. aeruginosa* O-antigen |
| SeqSero2 | ≥1.3 | *Salmonella* serotyping |
| hicap | ≥1.0 | *H. influenzae* capsule type |
| Snippy + snippy-core | ≥4.6 | SNP calling and core alignment |
| Gubbins | ≥3.3 | Recombination removal |
| IQ-TREE2 | ≥2.2 | ML phylogenetic tree |
| snp-dists | ≥0.8 | Pairwise SNP distances |
| Prokka | ≥1.14 | Annotation for custom reference genomes |
| chewBBACA | ≥3.3 | cgMLST allele calling |
