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
├── scripts/                   # Hjelpeskript
└── results/                   # Analyseresultater (ikke i git)
    └── <dato>/<prøve-id>/
        ├── Trimmed/           # Trimma reads (fastp)
        ├── Assembly/          # Genommontering (Shovill/SPAdes)
        ├── QC/                # FastQC + fastp JSON
        ├── MultiQC/           # MultiQC-rapport
        ├── QUAST/             # Assemblystatistikk
        ├── ID_Kraken2/        # Kraken2 + Bracken renhetskontroll
        ├── ID_GAMBIT/         # GAMBIT artsidentifikasjon
        ├── species.txt        # Identifisert art
        ├── MLST/              # Sekvenstyping (PubMLST)
        ├── AMRFinder/         # Resistens- og virulensgenar
        ├── MOBSuite/          # Plasmidtyping
        ├── MEfinder/          # Mobile genetiske element
        ├── SpaTyper/          # S. aureus spa-typing
        ├── SCCmec/            # S. aureus SCCmec-kassetttyping
        ├── AgrVATE/           # S. aureus agr-typing
        ├── EmmTyper/          # S. pyogenes emm-typing
        ├── ECTyper/           # E. coli O:H-serotyping
        ├── Kleborate/         # Klebsiella K/O-type, ESBL, virulens
        └── logs/
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

Dette installerer alle verktøy, inkludert MobileElementFinder via PyPI.

### 4. Last ned GAMBIT-database

```bash
pixi run --environment identification gambit-db-genomes download -d /databases/gambit_db/

```
Ev. manuell nedlasting fra : https://github.com/jlumpe/gambit

### 5. Initialiser MOB-suite databaser (~1.5 GB)

```bash
pixi run --environment mobsuite mob_init
```

### 6. Oppdater AMRFinder-database

```bash
pixi run --environment amrfinder4 amrfinder --update
```

---

## Databaser

### Kraken2 (PlusPF, ~80 GB)

```bash
mkdir -p /databases/kraken2_db
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240904.tar.gz -P /databases/kraken2_db/
tar -xzf /databases/kraken2_db/k2_pluspf_20240904.tar.gz -C /databases/kraken2_db/
```

---

## Bruk

### Hovedpipeline

```bash
cd ~/genomics

# Tørrkjøring
pixi run snakemake -n --cores 8

# Én prøve
pixi run snakemake --cores 8 --config samples=001k

# Flere prøver
pixi run snakemake --cores 8 --config "samples=[001k,002k]"

# Alle prøver (auto-deteksjon)
pixi run snakemake --cores 8
```

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

### Slektskapsanalyse (utbrudd) -- IKKE FERDIG

Kjøres etter hovedpipelinen når du mistenker utbrudd:

```bash
pixi run snakemake --snakefile slektskap_Snakefile --cores 8 \
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
    ├── Kraken2 + Bracken (renhetskontroll, parallelt)
    │
    ▼
Shovill (Assembly) → QUAST
    │
    ▼
GAMBIT
    │
    SJEKKPUNKT: art identifisert
    │
    ├── MLST
    ├── AMRFinder
    ├── MOB-suite (plasmid)
    ├── MEfinder (MGE)
    ├── [S. aureus]   → spaTyper, SCCmec, AgrVATE
    ├── [Klebsiella]  → Kleborate
    ├── [E. coli]     → ECTyper
    └── [S. pyogenes] → emmtyper
    │
    ▼
MultiQC
```

---

## Artsspesifikk typing

| Art | Verktøy | Resultat |
|-----|---------|---------|
| *S. aureus* | spaTyper, SCCmec, AgrVATE | spa-type, SCCmec-type, agr-gruppe |
| *Klebsiella* spp. | Kleborate v3 | K-type, O-type, ESBL, karbapenemase, virulens |
| *E. coli / Shigella* | ECTyper | O:H-serotype |
| *S. pyogenes* | emmtyper | emm-type |
| Alle | MEfinder | Insertionssekvenser og transposoner |
| Alle | MOB-suite | Plasmidtyping (relaxase, konjugasjon) |

---

## Konfigurasjon (config.yaml)

| Parameter | Standard | Beskrivelse |
|-----------|---------|-------------|
| `data_dir` | `data/` | Mappe med FASTQ-filer |
| `results_dir` | `results/` | Utdatamappe |
| `threads` | `8` | Antall tråder per regel |
| `kraken2_db` | `/databases/kraken2_db/` | Sti til Kraken2-database |
| `gambit_db`  | `/databases/gambit_db/`  | Sti til GAMBIT-database |

---

## Feilsøking

```bash
# Sjekk logg for ein spesifikk regel
cat results/26-03-2026/001k/logs/amrfinder.log

# Tving omkjøring av én regel
pixi run snakemake --cores 8 --forcerun amrfinder

# Tving full omkjøring av én prøve
pixi run snakemake --cores 8 results/26-03-2026/001k/.done --forceall
```

---

## Verktøy

| Verktøy | Versjon | Formål |
|---------|---------|--------|
| fastp | ≥0.23 | Trimming og QC |
| FastQC + MultiQC | ≥0.12 | Sekvenskvalitet |
| Shovill | ≥1.1 | Genommontering (SPAdes) |
| QUAST | ≥5.2 | Assemblystatistikk |
| Kraken2 + Bracken | ≥2.1 | Renhetskontroll av reads |
| GAMBIT | ≥0.5 | Artsidentifikasjon (på assembly) |
| MLST | ≥2.23 | Sekvenstyping (PubMLST) |
| AMRFinder | v4.x | Resistens- og virulensgen-deteksjon |
| MOB-suite | ≥3.1 | Plasmidtyping |
| MobileElementFinder | ≥1.0 | MGE-deteksjon (CGE) |
| spaTyper | ≥0.3 | *S. aureus* spa-typing |
| SCCmec | ≥1.0 | MRSA-kassetttyping |
| AgrVATE | ≥1.0 | *S. aureus* agr-typing |
| ECTyper | ≥2.0 | *E. coli* O:H-serotyping |
| emmtyper | ≥0.2 | *S. pyogenes* emm-typing |
| Kleborate | ≥3.0 | *Klebsiella* K/O-type, ESBL, virulens |
| Snippy | ≥4.6 | SNP-kalling mot referansegenom |
| IQ-TREE | ≥2.2 | Fylogenetisk tre (ML) |
| snp-dists | ≥0.8 | SNP-avstandsmatrise |
| chewBBACA | ≥3.3 | cgMLST allelkalling |
