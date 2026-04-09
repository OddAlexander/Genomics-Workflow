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
├── report_template.html   what too    # HTML-rapportmal
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
        ├── early_species.txt  # Mash Screen artsidentifikasjon (tidlig, styrerer DAG)
        ├── species.txt        # GAMBIT artsidentifikasjon
        ├── MLST/              # Sekvenstyping (PubMLST)
        ├── AMRFinder/         # Resistens- og virulensgen
        ├── MOBSuite/          # Plasmidtyping
        ├── MEfinder/          # Mobile genetiske element
        ├── SpaTyper/          # S. aureus spa-typing
        ├── SCCmec/            # S. aureus SCCmec-kassetttyping
        ├── AgrVATE/           # S. aureus agr-typing
        ├── EmmTyper/          # S. pyogenes emm-typing
        ├── ECTyper/           # E. coli O:H-serotyping
        ├── Kleborate/         # Klebsiella K/O-type, ESBL, virulens
        ├── Pasty/             # P. aeruginosa O-antigen serotyping
        ├── SeqSero2/          # Salmonella serotyping
        ├── Hicap/             # H. influenzae kapseltype
        ├── Report/            # HTML-rapport (sluttprodukt)
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

### 4. Last ned GAMBIT-database

```bash
pixi run --environment identification gambit-db-genomes download -d /databases/gambit_db/
```

### 5. Last ned Mash Screen-database (~3 GB)

```bash
mkdir -p /databases/mash_db
wget -P /databases/mash_db \
  https://gembox.cbcb.umd.edu/mash/refseq.genomes%2Bplasmid.k21s1000.msh
```

### 6. Initialiser MOB-suite databaser (~1.5 GB)

```bash
pixi run --environment mobsuite mob_init
```

### 7. Oppdater AMRFinder-database

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

# Alle prøver fra én dato
pixi run snakemake --cores 8 --config samples=19-03-2026

# Flere spesifikke prøver
pixi run snakemake --cores 8 --config "samples=[001k,002k]"

# Alle prøver (auto-deteksjon)
pixi run snakemake --cores 8

# Kjør én spesifikk regel for én prøve
pixi run snakemake results/26-03-2026/001k/MLST/mlst.tsv
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
    ├── Mash Screen (artsidentifikasjon på reads)
    │       └── SJEKKPUNKT: early_species.txt (styrer DAG-routing)
    │
    ├── Kraken2 + Bracken (artsidentifikasjon, parallelt)
    │
    ▼
Shovill (Assembly) → QUAST
    │
    ├── GAMBIT (artsidentifikasjon på assembly) → species.txt
    │
    ├── MLST
    ├── AMRFinder
    ├── MOB-suite (plasmid)
    ├── MEfinder (MGE)
    ├── [S. aureus]        → spaTyper, SCCmec, AgrVATE
    ├── [Klebsiella/E.coli]→ Kleborate
    ├── [E. coli/Shigella] → ECTyper
    ├── [S. pyogenes]      → emmtyper
    ├── [P. aeruginosa]    → Pasty
    ├── [Salmonella]       → SeqSero2
    └── [H. influenzae]    → hicap
    │
    ▼
MultiQC
    │
    ▼
HTML-rapport (Report/report.html)
    └── Artsidentifikasjon (Mash/Bracken/GAMBIT konkordans + renhet)
        AMR, plasmider, MGE, ST-type, artsspesifikk typing
```

---

## Artsspesifikk typing

| Art | Verktøy | Resultat |
|-----|---------|---------|
| *S. aureus* | spaTyper, SCCmec, AgrVATE | spa-type, SCCmec-type (MRSA/MSSA), agr-gruppe |
| *Klebsiella* spp. + *E. coli* | Kleborate v3 | K-type, O-type, ESBL, karbapenemase, virulens |
| *E. coli / Shigella* | ECTyper | O:H-serotype |
| *S. pyogenes* | emmtyper | emm-type |
| *P. aeruginosa* | Pasty | O-antigen serotype |
| *Salmonella* spp. | SeqSero2 | Serotyping (O- og H-antigen) |
| *H. influenzae* | hicap | Kapseltype (a–f, ikke-typbar) |
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
| `gambit_db` | `/databases/gambit_db/` | Sti til GAMBIT-database |
| `mash_db` | `/databases/mash_db/refseq.genomes+plasmid.k21s1000.msh` | Sti til Mash-database |

---

## Feilsøking

```bash
# Sjekk logg for en spesifikk regel
cat results/26-03-2026/001k/logs/amrfinder.log

# Tving omkjøring av én regel
pixi run snakemake --cores 8 --forcerun amrfinder

# Tving full omkjøring av én prøve
pixi run snakemake --cores 8 results/26-03-2026/001k/Report/report.html --forceall
```

---

## Verktøy

| Verktøy | Versjon | Formål |
|---------|---------|--------|
| fastp | ≥0.23 | Trimming og QC |
| FastQC + MultiQC | ≥0.12 | Sekvenskvalitet |
| Shovill | ≥1.1 | Genommontering (SPAdes) |
| QUAST | ≥5.2 | Assemblystatistikk |
| Mash Screen | ≥2.3 | Tidleg artsidentifikasjon på reads (styrer DAG) |
| Kraken2 + Bracken | ≥2.1 | Artsidentifikasjon og renhetskontroll |
| GAMBIT | ≥0.5 | Artsidentifikasjon på assembly |
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
| Pasty | ≥2.2 | *P. aeruginosa* O-antigen serotyping |
| SeqSero2 | ≥1.3 | *Salmonella* serotyping |
| hicap | ≥1.0 | *H. influenzae* kapseltype |
| Snippy | ≥4.6 | SNP-kalling mot referansegenom |
| IQ-TREE | ≥2.2 | Fylogenetisk tre (ML) |
| snp-dists | ≥0.8 | SNP-avstandsmatrise |
| chewBBACA | ≥3.3 | cgMLST allelkalling |
