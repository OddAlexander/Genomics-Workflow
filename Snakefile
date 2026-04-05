# Snakefile -- Bakteriell genomikk-pipeline (Optimalisert versjon)

import os, csv, json, re

from datetime import datetime

from pathlib import Path

 

# --- Farger og Logging ---

CYAN='\033[0;36m'; RESET='\033[0m'; BOLD='\033[1m'

def logm(m):    print(f"{CYAN}[{datetime.now().strftime('%H:%M:%S')}]{RESET} {m}")

def wlog(lf, st, msg, tool="", sample=""):

    os.makedirs(os.path.dirname(lf), exist_ok=True)

    with open(lf,'a') as f: f.write(f"[{datetime.now().isoformat()}] {st} | tool={tool} sample={sample}\n  {msg}\n\n")

# --- Konfigurasjon ---

configfile: "config.yaml"

DATA_DIR    = config.get("data_dir", "data/")

RESULTS_DIR = config.get("results_dir", "results/")

THREADS     = config.get("threads", 8)

KRAKEN2_DB  = config.get("kraken2_db", "/run/media/odd/Database/kraken2_db/")

# --- Artskart  ---

MLST_SCHEME = {
    "Enterococcus_faecium":"efaecium","Enterococcus_faecalis":"efaecalis","Enterococcus_hirae":"ehirae","Enterococcus_mundtii":"emundtii",
    "Klebsiella_pneumoniae":"klebsiella","Klebsiella_variicola":"klebsiella","Klebsiella_quasipneumoniae":"klebsiella","Klebsiella_quasivariicola":"klebsiella","Klebsiella_africana":"klebsiella",
    "Klebsiella_oxytoca":"koxytoca","Klebsiella_michiganensis":"koxytoca","Klebsiella_grimontii":"koxytoca",
    "Klebsiella_aerogenes":"kaerogenes","Enterobacter_cloacae":"ecloacae","Enterobacter_hormaechei":"ecloacae","Enterobacter_asburiae":"ecloacae",
    "Staphylococcus_aureus":"saureus","Staphylococcus_epidermidis":"sepidermidis","Staphylococcus_haemolyticus":"shaemolyticus","Staphylococcus_capitis":"scapitis",
    "Streptococcus_pneumoniae":"spneumoniae","Streptococcus_pyogenes":"spyogenes","Streptococcus_agalactiae":"sagalactiae","Streptococcus_dysgalactiae":"sdysgalactiae",
    "Escherichia_coli":"ecoli","Shigella_sonnei":"ecoli","Shigella_flexneri":"ecoli",
    "Pseudomonas_aeruginosa":"paeruginosa","Acinetobacter_baumannii":"abaumannii","Acinetobacter_pittii":"abaumannii","Acinetobacter_nosocomialis":"abaumannii",
    "Salmonella_enterica":"salmonella","Campylobacter_jejuni":"campylobacter","Campylobacter_coli":"campylobacter",
    "Clostridioides_difficile":"cdifficile","Haemophilus_influenzae":"hinfluenzae","Stenotrophomonas_maltophilia":"smaltophilia",
    "Neisseria_gonorrhoeae":"neisseria","Neisseria_meningitidis":"neisseria",
}
KLEB_PRESET = {
    "Klebsiella_pneumoniae":"kpsc","Klebsiella_variicola":"kpsc","Klebsiella_quasipneumoniae":"kpsc","Klebsiella_quasivariicola":"kpsc","Klebsiella_africana":"kpsc",
    "Klebsiella_oxytoca":"kosc","Klebsiella_michiganensis":"kosc","Klebsiella_grimontii":"kosc",
    "Escherichia_coli":"escherichia","Shigella_sonnei":"escherichia","Shigella_flexneri":"escherichia",
}
AMR_ORG = {
    "Enterococcus_faecium":"Enterococcus_faecium","Enterococcus_faecalis":"Enterococcus_faecalis",
    "Klebsiella_pneumoniae":"Klebsiella_pneumoniae","Klebsiella_variicola":"Klebsiella_pneumoniae","Klebsiella_quasipneumoniae":"Klebsiella_pneumoniae",
    "Klebsiella_oxytoca":"Klebsiella_oxytoca","Klebsiella_michiganensis":"Klebsiella_oxytoca","Klebsiella_aerogenes":"Klebsiella_pneumoniae",
    "Staphylococcus_aureus":"Staphylococcus_aureus","Streptococcus_pneumoniae":"Streptococcus_pneumoniae",
    "Streptococcus_pyogenes":"Streptococcus_pyogenes","Streptococcus_agalactiae":"Streptococcus_agalactiae",
    "Escherichia_coli":"Escherichia","Shigella_sonnei":"Escherichia","Shigella_flexneri":"Escherichia",
    "Pseudomonas_aeruginosa":"Pseudomonas_aeruginosa","Acinetobacter_baumannii":"Acinetobacter_baumannii",
    "Salmonella_enterica":"Salmonella","Campylobacter_jejuni":"Campylobacter","Campylobacter_coli":"Campylobacter",
    "Clostridioides_difficile":"Clostridioides_difficile","Haemophilus_influenzae":"Haemophilus_influenzae",
    "Neisseria_gonorrhoeae":"Neisseria_gonorrhoeae","Neisseria_meningitidis":"Neisseria_meningitidis",
}
SA  = {"Staphylococcus_aureus"}
EC  = {"Escherichia_coli","Shigella_sonnei","Shigella_flexneri","Shigella_dysenteriae","Shigella_boydii"}
GAS = {"Streptococcus_pyogenes"}

 
# 2. Python-hjelpefunksjoner

def get_species_param(wildcards, tool_dict):
    species_file = os.path.join(RESULTS_DIR, wildcards.sample, "species.txt")
    
    # 1. Håndter dry-run eller hvis species.txt ikke er laget ennå
    if not os.path.exists(species_file):
        return tool_dict.get("default", "all")

    # 2. Les arten fra fila
    with open(species_file, "r") as f:
        species = f.read().strip().replace(" ", "_")

    # 3. Slå opp i ordboken som sendes inn (f.eks. KLEB_PRESET)
    return tool_dict.get(species, tool_dict.get("default", "all"))

# --- Sample Deteksjon ---
# --- 1. Finn alle filer og lagre stier med en gang ---
# Dette finner sample-navnet og den unike delen av filnavnet ({any})
SAMPLES, ANY = glob_wildcards(os.path.join(DATA_DIR, "{sample}/{any}_R1.fastq.gz"))

# Lag en enkel ordbok for lynraskt oppslag: { 'sample1': 'full/sti/til/R1.fastq.gz' }
READS_R1 = {s: os.path.join(DATA_DIR, s, f"{a}_R1.fastq.gz") for s, a in zip(SAMPLES, ANY)}

print(f"Fant {len(SAMPLES)} samples.")

# --- Rules ---
rule all:
    input:
        expand(os.path.join(RESULTS_DIR, "{sample}/pipeline_sammendrag.html"), sample=SAMPLES)

rule fastp:
    input:
        r1 = lambda wc: READS_R1[wc.sample],
        r2 = lambda wc: READS_R1[wc.sample].replace("_R1", "_R2")
    output:
        r1   = os.path.join(RESULTS_DIR, "{sample}/Trimmed/{sample}_R1.fastq.gz"),
        r2   = os.path.join(RESULTS_DIR, "{sample}/Trimmed/{sample}_R2.fastq.gz"),
        json = os.path.join(RESULTS_DIR, "{sample}/QC/fastp.json")
    threads: 4
    shell:
        "pixi run fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -j {output.json} --thread {threads}"

rule assemble:
    input:
        r1 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/{sample}_R1.fastq.gz"),
        r2 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/{sample}_R2.fastq.gz")
    output:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    params:
        outdir = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "Assembly")
    threads: THREADS
    shell:
        "pixi run shovill --R1 {input.r1} --R2 {input.r2} --outdir {params.outdir} --cpus {threads} --force"

rule identify_species:
    input:
        r1 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/{sample}_R1.fastq.gz"),
        r2 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/{sample}_R2.fastq.gz")
    output:
        sp  = os.path.join(RESULTS_DIR, "{sample}/species.txt"),
        rep = os.path.join(RESULTS_DIR, "{sample}/ID_Kraken2/kraken2_report.txt")
    shell:
        """
        pixi run --environment identification kraken2 --db {KRAKEN2_DB} --paired --report {output.rep} {input.r1} {input.r2} > /dev/null
        grep -P '\\tS\\t' {output.rep} | sort -nrk3 | head -n1 | cut -f6 | sed 's/^[[:space:]]*//' > {output.sp}
        """

rule fastqc:
    input:
        r1 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/{sample}_R1.fastq.gz"),
        r2 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/{sample}_R2.fastq.gz")
    output:
        os.path.join(RESULTS_DIR, "{sample}/QC/fastqc/{sample}_R1_fastqc.html")
    params:
        outdir = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "QC/fastqc")
    threads: 2
    shell:
        "pixi run fastqc {input.r1} {input.r2} -o {params.outdir} --threads {threads}"

rule quast:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/QUAST/report.html")
    params:
        outdir = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "QUAST")
    shell:
        "pixi run quast {input.fa} -o {params.outdir}"

rule multiqc:
    input:
        fastp  = os.path.join(RESULTS_DIR, "{sample}/QC/fastp.json"),
        fastqc = os.path.join(RESULTS_DIR, "{sample}/QC/fastqc/{sample}_R1_fastqc.html"),
        quast  = os.path.join(RESULTS_DIR, "{sample}/QUAST/report.html")
    output:
        os.path.join(RESULTS_DIR, "{sample}/MultiQC/multiqc_report.html")
    params:
        qcdir   = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "QC"),
        quastdir = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "QUAST"),
        outdir  = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "MultiQC")
    shell:
        "pixi run multiqc {params.qcdir} {params.quastdir} -o {params.outdir}"

rule kleborate:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa"),
        sp = os.path.join(RESULTS_DIR, "{sample}/species.txt")
    output:
        os.path.join(RESULTS_DIR, "{sample}/Kleborate/{sample}_kleborate.txt")
    params:
        preset = lambda wc: get_species_param(wc, KLEB_PRESET)
    shell:
        "pixi run kleborate --data ./kleborate_data -a {input.fa} -o {output} --species {params.preset}"

rule report:
    input:
        kleb    = os.path.join(RESULTS_DIR, "{sample}/Kleborate/{sample}_kleborate.txt"),
        quast   = os.path.join(RESULTS_DIR, "{sample}/QUAST/report.html"),
        multiqc = os.path.join(RESULTS_DIR, "{sample}/MultiQC/multiqc_report.html")
    output:
        os.path.join(RESULTS_DIR, "{sample}/pipeline_sammendrag.html")
    shell:
        "echo 'Lager rapport for {wildcards.sample}' > {output}"