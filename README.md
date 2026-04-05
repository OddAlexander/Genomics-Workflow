# Bakteriell genomikk-pipeline

Snakemake-basert pipeline for helgenomsekvensering av kliniske bakterieisolater.  

---

## Innhold

```
genomics/
├── Snakefile                  # Hovudpipeline (per prøve)
├── slektskap_Snakefile        # Slektskapsanalyse/utbruddspipeline
├── cgmlst_Snakefile           # cgMLST-analyse med opplastingslenker
├── config.yaml                # Konfigurasjon
├── report_template.html       # HTML-rapportmal (per prøve)
├── pixi.toml                  # Pixi-miljødefinisjon
├── results/                   # Analyseresultater (ikke i git)
│   └── <prøve-id>/
│       ├── pipeline_sammendrag.html
│       ├── pipeline_sammendrag.json
│       ├── Assembly/
│       ├── MLST/
│       ├── AMR/
│       ├── Plasmids/
│       ├── MGE/
│       ├── SCCmec/
│       ├── SpaType/
│       ├── Kleborate/
│       ├── ECTyper/
│       ├── EmmType/
│       ├── ID_Kraken2/
│       ├── ID_FastANI/
│       ├── QC/
│       └── logs/
├── slektskapsanalyser/        # Output fra utbruddsanalyser (ikke i git)
│   └── <utbruddsnavn>/
│       ├── rapport/utbrudd_rapport.html
│       ├── rapport/snp_matrix.tsv
│       └── rapport/outbreak_tree.treefile
└── cgmlst_analyser/           # Output fra cgMLST-analyser (ikke i git)
    └── <analyse_navn>/
        ├── allelkalling/results_alleles.tsv
        └── rapport/cgmlst_rapport.html
```

---

## Forutsetninger

- Linux (testet på Ubuntu/Fedora)
- [Pixi](https://pixi.sh) installert
- Masse diskplass :P

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
cd ~/genomics
pixi install
```

### 4. Installer MobileElementFinder

```bash
pixi run --environment mefinder -- install-mefinder
```

### 5. Initialiser MOB-suite databaser (~1.5 GB)

```bash
pixi run --environment mobsuite -- mob_init
```

### 6. Oppdater AMRFinder-database

```bash
pixi run --environment amrfinder4 -- amrfinder --update
```

---

## Databaser

### Kraken2 (PlusPF, ~80 GB)

```bash
mkdir -p /databases/kraken2_db
# Last ned fra https://benlangmead.github.io/aws-indexes/k2
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240904.tar.gz -P /databases/kraken2_db/
tar -xzf /databases/kraken2_db/k2_pluspf_20240904.tar.gz -C /databases/kraken2_db/
echo 'export KRAKEN2_DB=/databases/kraken2_db' >> ~/.bashrc
```

### FastANI (kliniske referansegenomer)

```bash
mkdir -p /databases/fastANI_refs
cd /databases/fastANI_refs
pixi run -- datasets download genome taxon \
    "Enterococcus faecium" "Enterococcus faecalis" \
    "Klebsiella pneumoniae" "Klebsiella oxytoca" \
    "Staphylococcus aureus" "Escherichia coli" \
    "Pseudomonas aeruginosa" "Acinetobacter baumannii" \
    "Enterobacter cloacae" "Streptococcus pyogenes" \
    "Streptococcus pneumoniae" "Clostridioides difficile" \
    "Salmonella enterica" "Campylobacter jejuni" \
    "Haemophilus influenzae" "Neisseria meningitidis" \
    "Proteus mirabilis" \
    --reference --include genome \
    --filename /databases/fastANI_refs/clinical_refs.zip
unzip clinical_refs.zip
sudo chmod -R 644 /databases/fastANI_refs/ncbi_dataset/data/
sudo chmod -R 755 /databases/fastANI_refs/ncbi_dataset/
find /databases/fastANI_refs/ncbi_dataset/data -name "*.fna" | sort > /databases/fastANI_refs/reference_list.txt
echo 'export FASTANI_REFS=/databases/fastANI_refs/reference_list.txt' >> ~/.bashrc
```

### Snippy-referansegenomer

```bash
mkdir -p /databases/snippy_refs
cd /databases/snippy_refs
# Eksempel: E. faecium
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/250/945/GCF_000250945.1_Aus0004_v1/GCF_000250945.1_Aus0004_v1_genomic.gbff.gz"
gunzip GCF_000250945.1_Aus0004_v1_genomic.gbff.gz
mv GCF_000250945.1_Aus0004_v1_genomic.gbff efaecium_Aus0004.gbff
echo 'export SNIPPY_REF=/databases/snippy_refs/efaecium_Aus0004.gbff' >> ~/.bashrc
source ~/.bashrc
```

---

## Bruk

### Hovudpipeline

```bash
cd ~/genomics

# Tørkjøring
pixi run -- snakemake -n --cores 8

# Én prøve
pixi run -- snakemake --cores 8 --config samples=001k

# Alle prøver i én dato-run
pixi run -- snakemake --cores 8 --config run_date=26-03-2026

# Flere datoer
pixi run -- snakemake --cores 8 --config "run_date=[26-03-2026,05-03-2026]"

# Alle prøver (auto-deteksjon av dato-mapper)
pixi run -- snakemake --cores 8
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

> Filnavn må matche prøve-ID. Bruk rename-kommande ved behov:
> ```bash
> for d in data/*/*/; do
>     s=$(basename "$d")
>     r1=$(ls "$d"*_R1*.fastq.gz 2>/dev/null | head -1)
>     r2=$(ls "$d"*_R2*.fastq.gz 2>/dev/null | head -1)
>     [[ -n "$r1" ]] && mv "$r1" "${d}${s}_R1.fastq.gz"
>     [[ -n "$r2" ]] && mv "$r2" "${d}${s}_R2.fastq.gz"
> done
> ```

### Slektskapsanalyse (utbrudd) - IKKE FERDIG

Kjøres etter hovedpipelinen når du mistenker utbrudd:

```bash
pixi run -- snakemake --snakefile slektskap_Snakefile --cores 8 \
  --config \
    samples="[001k,002k,003k,004k]" \
    species="Enterococcus_faecium" \
    outbreak_name="VRE_Post4_Mars2026"
```

Med CT/HC-numre fra Enterobase/cgmlst.org:
```bash
pixi run -- snakemake --snakefile slektskap_Snakefile --cores 8 \
  --config \
    samples="[001k,002k,003k]" \
    species="Enterococcus_faecium" \
    outbreak_name="VRE_Post4_Mars2026" \
    'cluster_types={"001k":"CT5","002k":"CT5","003k":"CT12"}'
```

### cgMLST-analyse (IKKE FERDIG)

Allelkalling og lokal klyngeanalyse. Produserer opplastingsfil for offisiell CT/HC-typing:

```bash
# Kun allelkalling (for nettoppslag)
pixi run -- snakemake --snakefile cgmlst_Snakefile --cores 8 \
  --config \
    samples="[001k,002k,003k]" \
    species="Enterococcus_faecium" \
    analyse_navn="VRE_cgMLST_Mars2026" \
  --until allelkalling

# Full analyse med rapport
pixi run -- snakemake --snakefile cgmlst_Snakefile --cores 8 \
  --config \
    samples="[001k,002k,003k]" \
    species="Enterococcus_faecium" \
    analyse_navn="VRE_cgMLST_Mars2026"
```

Last opp `cgmlst_analyser/<navn>/allelkalling/results_alleles.tsv` til:
- **E. faecium / E. coli**: https://enterobase.warwick.ac.uk
- **K. pneumoniae / S. aureus**: https://www.cgmlst.org

---

## Pipeline-flyt

```
FASTQ-filer
    │
    ▼
fastp → FastQC → Shovill → QUAST/CheckM → MultiQC
                               │
                    Kraken2 + Bracken
                    FastANI
                               │
                    SJEKKPUNKT: art identifisert
                               │
    ┌──────────┬───────────────┼──────────────┬────────────┬──────────┐
    ▼          ▼               ▼              ▼            ▼          ▼
  MLST    AMRFinder        MOB-suite    MobileElement   [S. aureus] [Klebsiella]
          (v4, +org)     (plasmid-typ)   Finder (CGE)   spaTyper    Kleborate
                               │                        SCCmec      (K/O-type,
                               ▼                        AgrVATE      ESBL/Carb)
                     AMR-plasmid-kobling          [E. coli]  [S. pyogenes]
                               │                 ECTyper     emmtyper
                               ▼
                    pipeline_sammendrag.html
```

---

## Artsspesifikk typing

| Art | Verktøy | Resultat i rapport |
|-----|---------|-------------------|
| *S. aureus* | spaTyper, SCCmec, AgrVATE | MRSA/MSSA, SCCmec-type, spa-type (CC) |
| *Klebsiella* spp. | Kleborate v3 | K-type, O-type, ESBL, karbapenemase, virulensscoring |
| *E. coli / Shigella* | ECTyper, Kleborate | O:H-serotype |
| *S. pyogenes* | emmtyper | emm-type |
| Alle | MobileElementFinder (CGE) | Insertionssekvenser og transposoner |

---

## Konfigurasjon (config.yaml)

| Parameter | Standard | Beskrivelse |
|-----------|---------|-------------|
| `data_dir` | `data/` | Mappe med FASTQ-filer |
| `results_dir` | `results/` | Utdatamappe |
| `threads` | `8` | Antall tråder per regel |
| `kraken2_db` | `$KRAKEN2_DB` | Sti til Kraken2-database |
| `fastani_refs` | `$FASTANI_REFS` | Sti til FastANI-referanseliste |

---

## Feilsøking

```bash
# Sjekk logg for en spesifikk regel
cat results/001k/logs/amrfinder.log

# Tving omkjøring av én regel
pixi run -- snakemake --cores 8 --config samples=001k --forcerun amrfinder

# Tving full omkjøring av én prøve
pixi run -- snakemake --cores 8 --config samples=001k --forceall
```

---

## Verktøy og databaser

| Verktøy | Versjon | Formål |
|---------|---------|--------|
| fastp | ≥0.23 | Trimming og QC |
| FastQC + MultiQC | ≥0.12 | Sekvenskvalitet |
| Shovill | ≥1.1 | Genommontering (SPAdes) |
| QUAST | ≥5.2 | Assemblystatistikk |
| CheckM | ≥1.2 | Komplethet og kontaminering |
| Kraken2 + Bracken | ≥2.1 | Taksonomisk klassifisering |
| FastANI | ≥1.34 | ANI-artsbekreftelse |
| MLST | ≥2.23 | Sekvenstyping (PubMLST) |
| AMRFinder | v4.x | Resistens- og virulensgen-deteksjon |
| MOB-suite | ≥3.1 | Plasmidtyping |
| MobileElementFinder | 1.1.2 | MGE-deteksjon (CGE) |
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
| cgmlst-dists | ≥0.4 | cgMLST-avstandsmatrise |


