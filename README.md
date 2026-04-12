# Bakteriell genomikk-pipeline

Snakemake-basert pipeline for helgenomsekvensering av kliniske bakterieisolater.

---

## Innhold

```
genomics/
├── Snakefile                  # Hovedpipeline (per prøve)
├── slektskap_Snakefile        # Slektskapsanalyse/utbruddspipeline
├── config.yaml                # Konfigurasjon
├── pixi.toml                  # Pixi-miljødefinisjon
├── report_template.html       # HTML-rapportmal
├── scripts/
│   └── make_report.py         # Genererer HTML-rapport per prøve
└── results/                   # Analyseresultater (ikke i git)
    └── <dato>/<prøve-id>/
        ├── Trimmed/           # Trimma reads (fastp)
        ├── Assembly/          # Genommontering (Shovill/SPAdes)
        ├── QC/                # FastQC + fastp JSON
        ├── MultiQC/           # MultiQC-rapport
        ├── QUAST/             # Assemblystatistikk
        ├── ID_Kraken2/        # Kraken2 + Bracken artsidentifikasjon
        ├── ID_GAMBIT/         # GAMBIT artsidentifikasjon (på assembly)
        ├── ID_Skani/          # skani (kjøres kun ved usikker GAMBIT-match)
        ├── species.txt        # GAMBIT artsidentifikasjon (styrer DAG-routing)
        ├── MLST/              # Sekvenstyping (PubMLST)
        ├── AMRFinder/         # Resistens- og virulensgen
        ├── MOBSuite/          # Plasmidtyping
        ├── MEfinder/          # Mobile genetiske element
        ├── SpaTyper/          # S. aureus spa-typing
        ├── SCCmec/            # S. aureus SCCmec-kassetttyping
        ├── AgrVATE/           # S. aureus agr-typing
        ├── EmmTyper/          # S. pyogenes emm-typing
        ├── Kleborate/         # Klebsiella/E. coli K/O-type, ESBL, virulens
        ├── Pasty/             # P. aeruginosa O-antigen serotyping
        ├── SeqSero2/          # Salmonella serotyping
        ├── Hicap/             # H. influenzae kapseltype
        ├── pipeline_summary.html  # HTML-rapport (sluttprodukt)
        └── logs/
    └── pipeline_run_<YYYYMMDD_HHMM>.tsv  # Kjøringslogg (én fil per kjøring)
```

---

## Forutsetninger

- Linux (testet på Ubuntu/Fedora)
- [Pixi](https://pixi.sh) installert
- Masse diskplass

---

## Installasjon

### 1. Klon repositoriet

```bash
git clone https://github.com/OddAlexander/Genomics-Workflow.git ~/genomics
cd ~/genomics
```

### 2. Installer Pixi

```bash
curl -fsSL https://pixi.sh/install.sh | bash
source ~/.bashrc
```

### 3. Installer alle pixi-miljøer

```bash
pixi install
```

### 4. Last ned databaser

#### GAMBIT

```bash
pixi run --environment identification gambit-db-genomes download -d /databases/gambit_db/
```

#### Kraken2

Anbefalt: **PlusPF 8 GB** (~8 GB, passer i RAM på maskiner med ≥16 GB):

```bash
mkdir -p /databases/kraken2_db
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_8gb_20240904.tar.gz -P /databases/kraken2_db/
tar -xzf /databases/kraken2_db/k2_pluspf_8gb_20240904.tar.gz -C /databases/kraken2_db/
```

Alternativ: **PlusPF full** (~80 GB, krever ≥80 GB RAM for full ytelse):

```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240904.tar.gz -P /databases/kraken2_db/
tar -xzf /databases/kraken2_db/k2_pluspf_20240904.tar.gz -C /databases/kraken2_db/
```

#### skani (~3 GB skisse, ~6 500 komplette RefSeq-referansegenomer)

```bash
# Last ned genomer fra NCBI (~25 GB midlertidig diskplass)
sudo mkdir -p /databases/fastani_db /databases/skani_db && sudo chown $USER /databases/fastani_db /databases/skani_db
pixi run --environment identification datasets download genome taxon bacteria \
    --reference --assembly-source refseq --assembly-level complete \
    --include genome --filename /databases/fastani_db/bacteria_refs.zip
unzip /databases/fastani_db/bacteria_refs.zip -d /databases/fastani_db/
find /databases/fastani_db/ncbi_dataset/data -name "*.fna" > /databases/fastani_db/reference_list.txt

# Bygg skani-skissebase (~3 GB, beholder bare denne)
pixi run --environment identification skani sketch \
    -l /databases/fastani_db/reference_list.txt \
    -o /databases/skani_db/bacteria -t 16

# Slett rågenomene etter at skissen er ferdig (~25 GB frigjøres)
rm -rf /databases/fastani_db/ncbi_dataset /databases/fastani_db/bacteria_refs.zip
```

#### MOB-suite (~1.5 GB)

```bash
pixi run --environment mobsuite mob_init
```

#### AMRFinder

```bash
pixi run --environment amrfinder4 amrfinder --update
```

---

## Bruk

### Hovedpipeline

```bash
cd ~/genomics

# Tørrkjøring
pixi run snakemake -n --cores 16

# Én prøve
pixi run snakemake --cores 16 --resources mem_mb=60000 --config samples=001k

# Alle prøver fra én dato
pixi run snakemake --cores 16 --resources mem_mb=60000 --config samples=19-03-2026

# Flere spesifikke prøver
pixi run snakemake --cores 16 --resources mem_mb=60000 --config "samples=[001k,002k]"

# Alle prøver (auto-deteksjon)
pixi run snakemake --cores 16 --resources mem_mb=60000

# Kjør én spesifikk regel for én prøve
pixi run snakemake results/26-03-2026/001k/MLST/mlst.tsv
```

> **`mem_mb`** skaleres etter tilgjengelig RAM. Standard jobbreservasjoner: Kraken2=12 GB, Shovill=12 GB, skani=4 GB. Sett `mem_mb` til ca. total RAM minus 4 GB (f.eks. `12000` på en 16 GB maskin).

**Forventet mappestruktur for input:**

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

Filnavnene kan hete hva som helst så lenge de slutter på `_R1.fastq.gz` / `_R2.fastq.gz`.

### Slektskapsanalyse (utbrudd) — IKKE FERDIG

Kjøres etter hovedpipelinen når du mistenker utbrudd:

```bash
pixi run snakemake --snakefile slektskap_Snakefile --cores 16 \
  --config \
    samples="[001k,002k,003k,004k]" \
    species="Enterococcus_faecium" \
    outbreak_name="VRE_Post4_Mars2026"
```

---

## Pipeline-flyt

```
FASTQ-filer
    │
    ▼
fastp (trimming) → FastQC
    │
    ├── Kraken2 + Bracken (artsidentifikasjon, parallelt)
    │
    ▼
Shovill (Assembly) → QUAST
    │
    ├── GAMBIT (artsidentifikasjon på assembly) → species.txt
    │       └── SJEKKPUNKT: styrer DAG-routing
    │
    ├── skani (kun hvis GAMBIT closest.distance > skani_threshold)
    ├── MLST
    ├── AMRFinder
    ├── MOB-suite (plasmid)
    ├── MEfinder (MGE)
    ├── [S. aureus]        → spaTyper, SCCmec, AgrVATE
    ├── [Klebsiella/E.coli]→ Kleborate
    ├── [S. pyogenes]      → emmtyper
    ├── [P. aeruginosa]    → Pasty
    ├── [Salmonella]       → SeqSero2
    └── [H. influenzae]    → hicap
    │
    ▼
MultiQC
    │
    ▼
HTML-rapport (pipeline_summary.html)
```

---

## Artsspesifikk typing

| Art | Verktøy | Resultat |
|-----|---------|---------|
| *S. aureus* | spaTyper, SCCmec, AgrVATE | spa-type, SCCmec-type (MRSA/MSSA), agr-gruppe |
| *Klebsiella* spp. + *E. coli* | Kleborate v3 | K-type, O-type, ESBL, karbapenemase, virulens |
| *S. pyogenes* | emmtyper | emm-type |
| *P. aeruginosa* | Pasty | O-antigen serotype |
| *Salmonella* spp. | SeqSero2 | Serotyping (O- og H-antigen) |
| *H. influenzae* | hicap | Kapseltype (a–f, ikke-typbar) |
| Alle | MEfinder | Insertionssekvenser og transposoner |
| Alle | MOB-suite | Plasmidtyping (relaxase, konjugasjon) |

ECTyper (E. coli O:H-serotyping) krever MASH-skissen fra Zenodo. Last ned med:

```bash
mkdir -p databases/ectyper
wget -O databases/ectyper/EnteroRef_GTDBSketch_20231003_V2.msh \
    "https://zenodo.org/records/13969103/files/EnteroRef_GTDBSketch_20231003_V2.msh?download=1"
```

---

## Konfigurasjon (config.yaml)

| Parameter | Standard | Beskrivelse |
|-----------|---------|-------------|
| `data_dir` | `data/` | Mappe med FASTQ-filer |
| `results_dir` | `results/` | Utdatamappe |
| `threads` | `16` | Antall tråder per regel |
| `kraken2_db` | `/databases/kraken2_db/` | Sti til Kraken2-database |
| `kraken2_mem_mb` | `25000` | MB RAM per Kraken2-jobb |
| `gambit_db` | `/databases/gambit_db/` | Sti til GAMBIT-database |
| `skani_db` | `/databases/skani_db/bacteria` | Sti til skani-skissebase (katalog) |
| `skani_threshold` | `0.3` | GAMBIT `closest.distance` over denne utløser skani |
| `skani_mem_mb` | `4000` | MB RAM per skani-jobb |
| `skani_threads` | `8` | Antall tråder per skani-jobb |
| `shovill_mem_mb` | `12000` | MB RAM per Shovill/SPAdes-jobb (SPAdes begrenses til dette / 1024 GB) |

---

## Feilsøking

```bash
# Sjekk logg for en spesifikk regel
cat results/26-03-2026/001k/logs/amrfinder.log

# Tving omkjøring av én regel
pixi run snakemake --cores 16 --forcerun amrfinder

# Tving full omkjøring av én prøve
pixi run snakemake --cores 16 --resources mem_mb=60000 \
  results/26-03-2026/001k/pipeline_summary.html --forceall

# Gjenoppta avbrutt kjøring
pixi run snakemake --cores 16 --resources mem_mb=60000 --rerun-incomplete
```

---

## Verktøy

| Verktøy | Versjon | Formål |
|---------|---------|--------|
| fastp | ≥0.23 | Trimming og QC |
| FastQC + MultiQC | ≥0.12 | Sekvenskvalitet |
| Shovill | ≥1.1 | Genommontering (SPAdes) |
| QUAST | ≥5.2 | Assemblystatistikk |
| Kraken2 + Bracken | ≥2.1 | Artsidentifikasjon og renhetskontroll |
| GAMBIT | ≥0.5 | Artsidentifikasjon på assembly (SJEKKPUNKT — styrer DAG) |
| skani | ≥0.2 | Rask ANI-beregning ved usikker GAMBIT-match |
| MLST | ≥2.23 | Sekvenstyping (PubMLST) |
| AMRFinder | v4.x | Resistens- og virulensgen-deteksjon |
| MOB-suite | ≥3.1 | Plasmidtyping |
| MobileElementFinder | ≥1.0 | MGE-deteksjon (CGE) |
| spaTyper | ≥0.3 | *S. aureus* spa-typing |
| SCCmec | ≥1.0 | MRSA-kassetttyping |
| AgrVATE | ≥1.0 | *S. aureus* agr-typing |
| emmtyper | ≥0.2 | *S. pyogenes* emm-typing |
| Kleborate | ≥3.0 | *Klebsiella* K/O-type, ESBL, virulens |
| Pasty | ≥2.2 | *P. aeruginosa* O-antigen serotyping |
| SeqSero2 | ≥1.3 | *Salmonella* serotyping |
| hicap | ≥1.0 | *H. influenzae* kapseltype |
| Snippy | ≥4.6 | SNP-kalling mot referansegenom |
| IQ-TREE | ≥2.2 | Fylogenetisk tre (ML) |
| snp-dists | ≥0.8 | SNP-avstandsmatrise |
| chewBBACA | ≥3.3 | cgMLST allelkalling |
