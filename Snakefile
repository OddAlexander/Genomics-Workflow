# Snakefile -- Bakteriell genomikk-pipeline
import os
from pathlib import Path

# --- Konfigurasjon ---
configfile: "config.yaml"

DATA_DIR    = config.get("data_dir", "data/")
RESULTS_DIR = config.get("results_dir", "results/")
THREADS     = config.get("threads", 8)
KRAKEN2_DB  = config.get("kraken2_db", "/databases/kraken2_db_mini/")

# --- Artsgrupper og skjema ---
SA  = {"Staphylococcus_aureus"}
EC  = {"Escherichia_coli", "Shigella_sonnei", "Shigella_flexneri", "Shigella_dysenteriae", "Shigella_boydii"}
GAS = {"Streptococcus_pyogenes"}

KLEB_PRESET = {
    "Klebsiella_pneumoniae":"kpsc", "Klebsiella_variicola":"kpsc", "Klebsiella_quasipneumoniae":"kpsc",
    "Klebsiella_quasivariicola":"kpsc", "Klebsiella_africana":"kpsc",
    "Klebsiella_oxytoca":"kosc", "Klebsiella_michiganensis":"kosc", "Klebsiella_grimontii":"kosc",
    "Escherichia_coli":"escherichia", "Shigella_sonnei":"escherichia", "Shigella_flexneri":"escherichia",
}

MLST_SCHEME = {
    "Enterococcus_faecium":"efaecium", "Enterococcus_faecalis":"efaecalis", "Enterococcus_hirae":"ehirae", "Enterococcus_mundtii":"emundtii",
    "Klebsiella_pneumoniae":"klebsiella", "Klebsiella_variicola":"klebsiella", "Klebsiella_quasipneumoniae":"klebsiella", "Klebsiella_quasivariicola":"klebsiella", "Klebsiella_africana":"klebsiella",
    "Klebsiella_oxytoca":"koxytoca", "Klebsiella_michiganensis":"koxytoca", "Klebsiella_grimontii":"koxytoca",
    "Klebsiella_aerogenes":"kaerogenes", "Enterobacter_cloacae":"ecloacae", "Enterobacter_hormaechei":"ecloacae", "Enterobacter_asburiae":"ecloacae",
    "Staphylococcus_aureus":"saureus", "Staphylococcus_epidermidis":"sepidermidis", "Staphylococcus_haemolyticus":"shaemolyticus", "Staphylococcus_capitis":"scapitis",
    "Streptococcus_pneumoniae":"spneumoniae", "Streptococcus_pyogenes":"spyogenes", "Streptococcus_agalactiae":"sagalactiae", "Streptococcus_dysgalactiae":"sdysgalactiae",
    "Escherichia_coli":"ecoli", "Shigella_sonnei":"ecoli", "Shigella_flexneri":"ecoli",
    "Pseudomonas_aeruginosa":"paeruginosa", "Acinetobacter_baumannii":"abaumannii", "Acinetobacter_pittii":"abaumannii", "Acinetobacter_nosocomialis":"abaumannii",
    "Salmonella_enterica":"salmonella", "Campylobacter_jejuni":"campylobacter", "Campylobacter_coli":"campylobacter",
    "Clostridioides_difficile":"cdifficile", "Haemophilus_influenzae":"hinfluenzae", "Stenotrophomonas_maltophilia":"smaltophilia",
    "Neisseria_gonorrhoeae":"neisseria", "Neisseria_meningitidis":"neisseria",
}

AMR_ORG = {
    "Enterococcus_faecium":"Enterococcus_faecium", "Enterococcus_faecalis":"Enterococcus_faecalis",
    "Klebsiella_pneumoniae":"Klebsiella_pneumoniae", "Klebsiella_variicola":"Klebsiella_pneumoniae", "Klebsiella_quasipneumoniae":"Klebsiella_pneumoniae",
    "Klebsiella_oxytoca":"Klebsiella_oxytoca", "Klebsiella_michiganensis":"Klebsiella_oxytoca", "Klebsiella_aerogenes":"Klebsiella_pneumoniae",
    "Staphylococcus_aureus":"Staphylococcus_aureus", "Streptococcus_pneumoniae":"Streptococcus_pneumoniae",
    "Streptococcus_pyogenes":"Streptococcus_pyogenes", "Streptococcus_agalactiae":"Streptococcus_agalactiae",
    "Escherichia_coli":"Escherichia", "Shigella_sonnei":"Escherichia", "Shigella_flexneri":"Escherichia",
    "Pseudomonas_aeruginosa":"Pseudomonas_aeruginosa", "Acinetobacter_baumannii":"Acinetobacter_baumannii",
    "Salmonella_enterica":"Salmonella", "Campylobacter_jejuni":"Campylobacter", "Campylobacter_coli":"Campylobacter",
    "Clostridioides_difficile":"Clostridioides_difficile", "Haemophilus_influenzae":"Haemophilus_influenzae",
    "Neisseria_gonorrhoeae":"Neisseria_gonorrhoeae", "Neisseria_meningitidis":"Neisseria_meningitidis",
}

# --- Verktøy per artsgruppe ---
ALWAYS_TOOLS = [
    "QUAST/report.html",
    "MultiQC/multiqc_report.html",
    "MLST/{sample}_mlst.tsv",
    "AMRFinder/{sample}_amrfinder.tsv",
]

SPECIES_TOOLS = [
    (KLEB_PRESET.keys(), ["Kleborate/kleborate_output.tsv"]),
    (SA,                  ["SpaTyper/spatyper.txt", "SCCmec/sccmec.txt", "AgrVATE/agrvate.txt"]),
    (GAS,                 ["EmmTyper/emmtyper.txt"]),
]

# --- Hjelpefunksjoner ---
def read_species(sample):
    sp = os.path.join(RESULTS_DIR, sample, "species.txt")
    return open(sp).read().strip().replace(" ", "_") if os.path.exists(sp) else ""

def get_all_outputs(wildcards):
    checkpoints.identify_species.get(sample=wildcards.sample)
    sp    = read_species(wildcards.sample)
    s     = wildcards.sample
    tools = list(ALWAYS_TOOLS)
    for species_set, outputs in SPECIES_TOOLS:
        if sp in species_set:
            tools += outputs
    return [os.path.join(RESULTS_DIR, s, t.format(sample=s)) for t in tools]

# --- Sample Deteksjon ---
SAMPLES, ANY = glob_wildcards(os.path.join(DATA_DIR, "{sample}/{any}_R1.fastq.gz"))
READS_DICT = {s: os.path.join(DATA_DIR, s, f"{a}") for s, a in zip(SAMPLES, ANY)}
print(f"Fant {len(SAMPLES)} samples.")

# --- Rules ---
rule all:
    input:
        expand(os.path.join(RESULTS_DIR, "{sample}/.done"), sample=SAMPLES)

rule done:
    input:
        get_all_outputs
    output:
        touch(os.path.join(RESULTS_DIR, "{sample}/.done"))

rule fastp:
    input:
        r1 = lambda wc: f"{READS_DICT[wc.sample]}_R1.fastq.gz",
        r2 = lambda wc: f"{READS_DICT[wc.sample]}_R2.fastq.gz"
    output:
        r1   = os.path.join(RESULTS_DIR, "{sample}/Trimmed/{sample}_R1.fastq.gz"),
        r2   = os.path.join(RESULTS_DIR, "{sample}/Trimmed/{sample}_R2.fastq.gz"),
        json = os.path.join(RESULTS_DIR, "{sample}/QC/fastp.json")
    threads: THREADS
    shell:
        "pixi run fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -j {output.json} --thread {threads}"

checkpoint identify_species:
    input:
        r1 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/{sample}_R1.fastq.gz"),
        r2 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/{sample}_R2.fastq.gz")
    output:
        sp      = os.path.join(RESULTS_DIR, "{sample}/species.txt"),
        rep     = os.path.join(RESULTS_DIR, "{sample}/ID_Kraken2/kraken2_report.txt"),
        bracken = os.path.join(RESULTS_DIR, "{sample}/ID_Kraken2/bracken_species.txt")
    shell:
        """
        pixi run --environment identification kraken2 --db {KRAKEN2_DB} --memory-mapping --paired --report {output.rep} {input.r1} {input.r2} > /dev/null
        pixi run --environment identification bracken -d {KRAKEN2_DB} -i {output.rep} -o {output.bracken} -l S
        RESULT=$(tail -n +2 {output.bracken} | sort -t$'\\t' -k6 -rn | head -n1 | cut -f1 | sed 's/ /_/g')
        echo "${{RESULT:-Unknown}}" > {output.sp}
        """

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

rule mlst:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa"),
        sp = os.path.join(RESULTS_DIR, "{sample}/species.txt")
    output:
        os.path.join(RESULTS_DIR, "{sample}/MLST/{sample}_mlst.tsv")
    params:
        scheme = lambda wc: f"--scheme {MLST_SCHEME[read_species(wc.sample)]}" if read_species(wc.sample) in MLST_SCHEME else ""
    shell:
        "pixi run mlst {params.scheme} {input.fa} > {output}"

rule kleborate:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa"),
        sp = os.path.join(RESULTS_DIR, "{sample}/species.txt")
    output:
        os.path.join(RESULTS_DIR, "{sample}/Kleborate/kleborate_output.tsv")
    params:
        outdir = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "Kleborate"),
        preset = lambda wc: KLEB_PRESET.get(read_species(wc.sample), "kpsc")
    shell:
        """
        pixi run --environment identification kleborate -a {input.fa} -o {params.outdir} -p {params.preset}
        mv $(ls {params.outdir}/*_output.txt | grep -v hAMRonization) {output}
        """

rule spatyper:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/SpaTyper/spatyper.txt")
    shell:
        "pixi run spaTyper -f {input.fa} -o {output}"

rule sccmec:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/SCCmec/sccmec.txt")
    shell:
        "pixi run sccmec --input {input.fa} --output {output}"

rule agrvate:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/AgrVATE/agrvate.txt")
    shell:
        "pixi run agrvate -i {input.fa} -o {output}"

rule emmtyper:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/EmmTyper/emmtyper.txt")
    shell:
        "pixi run emmtyper {input.fa} > {output}"

rule amrfinder:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa"),
        sp = os.path.join(RESULTS_DIR, "{sample}/species.txt")
    output:
        os.path.join(RESULTS_DIR, "{sample}/AMRFinder/{sample}_amrfinder.tsv")
    params:
        org = lambda wc: f"--organism {AMR_ORG[read_species(wc.sample)]}" if read_species(wc.sample) in AMR_ORG else ""
    shell:
        "pixi run --environment amrfinder4 amrfinder --nucleotide {input.fa} {params.org} --output {output} --threads {THREADS} --plus"

rule multiqc:
    input:
        fastp   = os.path.join(RESULTS_DIR, "{sample}/QC/fastp.json"),
        fastqc  = os.path.join(RESULTS_DIR, "{sample}/QC/fastqc/{sample}_R1_fastqc.html"),
        quast   = os.path.join(RESULTS_DIR, "{sample}/QUAST/report.html")
    output:
        os.path.join(RESULTS_DIR, "{sample}/MultiQC/multiqc_report.html")
    params:
        qcdir    = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "QC"),
        quastdir = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "QUAST"),
        outdir   = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "MultiQC")
    shell:
        "pixi run multiqc {params.qcdir} {params.quastdir} -o {params.outdir}"

