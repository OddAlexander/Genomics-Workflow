# Bakteriell genomikk-pipeline

Snakemake-pipeline for helgenomsekvensering av kliniske bakterieisolater — trimming, assembly, artsidentifikasjon, resistens, typing og SNP-fylogeni.

---

## Krav

- Linux (testet på Ubuntu/Fedora)
- [Pixi](https://pixi.sh) installert (`curl -fsSL https://pixi.sh/install.sh | bash`)

---

## Installasjon

```bash
git clone https://github.com/OddAlexander/Genomics-Workflow.git ~/genomics
cd ~/genomics
pixi install
```

### Databaser

**GAMBIT**
```bash
pixi run --environment identification gambit-db-genomes download -d /databases/gambit_db/
```

**Kraken2** — anbefalt: PlusPF 8 GB (~8 GB, krever ≥16 GB RAM):
```bash
mkdir -p /databases/kraken2_db_light
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_8gb_20240904.tar.gz -P /databases/kraken2_db_light/
tar -xzf /databases/kraken2_db_light/k2_pluspf_8gb_20240904.tar.gz -C /databases/kraken2_db_light/
```

**skani** (~3 GB ferdig skisse, ~25 GB midlertidig under bygging):
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

## Bruk

### Inputformat

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

Filnavn kan være hva som helst, men må slutte på `_R1.fastq.gz` / `_R2.fastq.gz`.

### Hovedpipeline (`Snakefile`)

```bash
# Tørrkjøring
pixi run snakemake -n --cores 16

# Én prøve / én dato / flere prøver / alle
pixi run snakemake --cores 16 --resources mem_mb=60000 --config samples=001k
pixi run snakemake --cores 16 --resources mem_mb=60000 --config samples=19-03-2026
pixi run snakemake --cores 16 --resources mem_mb=60000 --config "samples=[001k,002k]"
pixi run snakemake --cores 16 --resources mem_mb=60000

# Én spesifikk fil
pixi run snakemake results/26-03-2026/001k/MLST/mlst.tsv
```

> **`mem_mb`** bør settes til total RAM minus ~4 GB (f.eks. `12000` på 16 GB maskin). Kraken2, Shovill og skani reserverer henholdsvis 12, 12 og 4 GB.

### Slektskapsanalyse (`Snakefile_phylo`)

Kjøres uavhengig av hovedpipelinen — tar råreads direkte som input.

```bash
# Med ekstern referanse (.gbk fra Prokka — anbefalt)
pixi run snakemake -s Snakefile_phylo --cores 16 \
    --config ref=/databases/prokka_refs/19-03-2026_005a/19-03-2026_005a.gbk

# Uten referanse — første sample prokka-annoteres automatisk
pixi run snakemake -s Snakefile_phylo --cores 16

# Filtrer på dato eller prøve-ID
pixi run snakemake -s Snakefile_phylo --cores 16 \
    --config ref=... samples=19-03-2026
pixi run snakemake -s Snakefile_phylo --cores 16 \
    --config ref=... "samples=[19-03-2026/005a,19-03-2026/007b]"

# Egendefinert mappenavn (standard: dato+klokkeslett)
pixi run snakemake -s Snakefile_phylo --cores 16 \
    --config ref=... run_name=salmonella_utbrudd
```

Output lagres i `results_phylo/<run_name>/` (standard: `DD-MM-YYYY_HHMM`):

```
results_phylo/07-05-2026_1621/
├── <dato>/<prøve-id>/
│   ├── Trimmed/      fastp
│   ├── QC/           FastQC + fastp JSON + HTML
│   ├── Assembly/     Shovill
│   ├── QUAST/
│   ├── ID_Skani/     artsidentifikasjon (species.txt)
│   ├── Prokka/       kun ved auto-referanse (første sample)
│   └── MultiQC/
├── Snippy/
│   └── <prøve-id>/   snippy SNP-kalling mot referanse
├── Core/             snippy-core — core.aln, core.full.aln
├── Gubbins/          rekombinasjonsfiltert alignment + tre
├── SNP_Dists/        snp_dists.tsv
├── IQtree/           iqtree.treefile  (GTR+G, UFBoot 1000)
├── RAxML/            raxml.raxml.bestTree  (GTR+G, 100 bootstrap)
├── MultiQC_run/      multiqc_report.html  (samlet QC for hele kjøringen)
└── phylo_report.html prøveoversikt, SNP-matrise og fylogenetisk tre
```

#### Lage egne referanser

```bash
scripts/pixi/pixi_annotate_prokka.sh 19-03-2026
scripts/pixi/pixi_annotate_prokka.sh 19-03-2026/005a --genus Staphylococcus --species aureus
```

---

## Pipeline-flyt

### Hovedpipeline

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
  ├── GAMBIT → species.txt  ← sjekkpunkt (styrer DAG-routing)
  ├── skani  (kun ved usikker GAMBIT-match)
  ├── MLST
  ├── AMRFinder
  ├── MOB-suite
  ├── MEfinder
  ├── [S. aureus]         → StaphScope
  ├── [Klebsiella/E.coli] → Kleborate
  ├── [S. pyogenes]       → emmtyper
  ├── [P. aeruginosa]     → Pasty
  ├── [Salmonella]        → SeqSero2
  └── [H. influenzae]     → hicap
  │
  ▼
MultiQC → HTML-rapport
```

### Slektskapsanalyse

```
FASTQ (råreads)
  │
  ▼
fastp → FastQC → Shovill → QUAST → skani
  │                │
  │                └── Prokka  (kun ved auto-ref)
  ▼
Snippy (per prøve, mot referanse)
  │
  ▼
snippy-core → core.full.aln
  │
  ▼
Gubbins (rekombinasjonsfjerning)
  │
  ├── snp-dists  → SNP-avstandsmatrise
  ├── IQ-TREE2   → ML-tre
  └── RAxML-NG   → ML-tre (alternativt)
  │
  ▼
MultiQC + HTML-rapport (phylo_report.html)
```

---

## Artsspesifikk typing

| Art | Verktøy | Output |
|-----|---------|--------|
| *S. aureus* | StaphScope | spa-type, SCCmec, MRSA/MSSA, virulens, lineage |
| *Klebsiella* spp. + *E. coli* | Kleborate v3 | K/O-type, ESBL, karbapenemase, virulens |
| *S. pyogenes* | emmtyper | emm-type |
| *P. aeruginosa* | Pasty | O-antigen serotype |
| *Salmonella* spp. | SeqSero2 | O- og H-antigen serotyping |
| *H. influenzae* | hicap | Kapseltype (a–f / ikke-typbar) |
| Alle | MEfinder | Insersjonssekvenser og transposoner |
| Alle | MOB-suite | Plasmidtyping (relaxase, konjugasjon) |

ECTyper (E. coli O:H-serotyping) krever en MASH-skisse fra Zenodo:

```bash
mkdir -p databases/ectyper
wget -O databases/ectyper/EnteroRef_GTDBSketch_20231003_V2.msh \
    "https://zenodo.org/records/13969103/files/EnteroRef_GTDBSketch_20231003_V2.msh?download=1"
```

---

## Konfigurasjon (`config.yaml`)

| Parameter | Standard | Beskrivelse |
|-----------|----------|-------------|
| `data_dir` | `data/` | Mappe med FASTQ-filer |
| `results_dir` | `results/` | Utdatamappe |
| `threads` | `8` | Tråder per regel |
| `kraken2_db` | `/databases/kraken2_db_light/` | Kraken2-database |
| `kraken2_mem_mb` | `12000` | MB RAM per Kraken2-jobb |
| `gambit_db` | `/databases/gambit_db/` | GAMBIT-database |
| `skani_db` | `/databases/skani_db/bacteria` | skani-skissebase |
| `skani_threshold` | `0.3` | GAMBIT `closest.distance` over denne utløser skani |
| `skani_mem_mb` | `4000` | MB RAM per skani-jobb |
| `skani_threads` | `8` | Tråder per skani-jobb |
| `shovill_mem_mb` | `12000` | MB RAM per Shovill-jobb |

---

## Feilsøking

```bash
# Se logg
cat results/26-03-2026/001k/logs/amrfinder.log

# Tving omkjøring av én regel
pixi run snakemake --cores 16 --forcerun amrfinder

# Gjenoppta avbrutt kjøring
pixi run snakemake --cores 16 --resources mem_mb=60000 --rerun-incomplete
```

---

## Verktøy

| Verktøy | Versjon | Formål |
|---------|---------|--------|
| fastp | ≥0.23 | Trimming |
| FastQC + MultiQC | ≥0.12 | Sekvenskvalitet |
| Shovill | ≥1.1 | Assembly (SPAdes) |
| QUAST | ≥5.2 | Assemblystatistikk |
| Kraken2 + Bracken | ≥2.1 | Artsidentifikasjon (reads) |
| GAMBIT | ≥0.5 | Artsidentifikasjon (assembly) — styrer DAG |
| skani | ≥0.2 | ANI-basert artsidentifikasjon ved usikker GAMBIT-match |
| MLST | ≥2.23 | Sekvenstyping (PubMLST) |
| AMRFinder | v4.x | Resistens- og virulensgener |
| MOB-suite | ≥3.1 | Plasmidtyping |
| MobileElementFinder | ≥1.0 | MGE-deteksjon |
| StaphScope | ≥1.0 | *S. aureus*: spa, SCCmec, MRSA, virulens, lineage |
| emmtyper | ≥0.2 | *S. pyogenes* emm-typing |
| Kleborate | ≥3.0 | *Klebsiella* / *E. coli* typing |
| Pasty | ≥2.2 | *P. aeruginosa* O-antigen |
| SeqSero2 | ≥1.3 | *Salmonella* serotyping |
| hicap | ≥1.0 | *H. influenzae* kapseltype |
| Snippy + snippy-core | ≥4.6 | SNP-kalling og core-alignment |
| Gubbins | ≥3.3 | Rekombinasjonsfjerning |
| IQ-TREE2 | ≥2.2 | ML-fylogenetisk tre |
| RAxML-NG | ≥1.2 | ML-fylogenetisk tre (alternativt) |
| snp-dists | ≥0.8 | Parvise SNP-avstander |
| Prokka | ≥1.14 | Annotasjon for egne referansegenomer |
| chewBBACA | ≥3.3 | cgMLST allelkalling |
