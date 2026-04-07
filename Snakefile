# Snakefile -- Bakteriell genomikk-pipeline
import os

# --- Konfigurasjon ---
configfile: "config.yaml"

DATA_DIR    = config.get("data_dir", "data/")
RESULTS_DIR = config.get("results_dir", "results/")
THREADS     = config.get("threads", 8)
KRAKEN2_DB  = config.get("kraken2_db", "/databases/kraken2_db/")

# --- Artsgrupper og skjemaer ---
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
    "MLST/mlst.tsv",
    "AMRFinder/amrfinder.tsv",
    "MOBSuite/mobtyper.tsv",
    "MEfinder/mefinder.tsv",
]

SPECIES_TOOLS = [
    (KLEB_PRESET.keys(), ["Kleborate/kleborate_output.tsv"]),
    (SA,                  ["SpaTyper/spatyper.txt", "SCCmec/sccmec.txt", "AgrVATE/agrvate.txt"]),
    (GAS,                 ["EmmTyper/emmtyper.txt"]),
    (EC,                  ["ECTyper/ectyper.tsv"]),
]

# --- Hjelpefunksjoner ---
def read_species(sample):
    sp_path = checkpoints.identify_species.get(sample=sample).output.sp
    with open(sp_path, "r") as f:
        species = f.read().strip().replace(" ", "_")
    return species if species else "unknown_species"

def get_all_outputs(wildcards):
    sp    = read_species(wildcards.sample)
    s     = wildcards.sample
    tools = list(ALWAYS_TOOLS)
    for species_set, outputs in SPECIES_TOOLS:
        if sp in species_set:
            tools += outputs
    return [os.path.join(RESULTS_DIR, s, t) for t in tools]

def get_multiqc_inputs(wildcards):
    return [f for f in get_all_outputs(wildcards) if "multiqc_report.html" not in f]

# --- Sample Deteksjon ---
# {sample} fanger hele den relative stien under DATA_DIR, f.eks. "19-03-2026/005a"
# slik at mappestrukturen i results/ speiler data/
ALL_SAMPLES, ANY = glob_wildcards(os.path.join(DATA_DIR, "{sample}/{any}_R1.fastq.gz"))
READS_DICT = {s: os.path.join(DATA_DIR, s, a) for s, a in zip(ALL_SAMPLES, ANY)}
assert len(READS_DICT) == len(ALL_SAMPLES), \
    f"Duplikate sample-ID-er oppdaget: {len(ALL_SAMPLES)} filer, men bare {len(READS_DICT)} unike navn"

_filter = config.get("samples", None)
if _filter:
    _filter = _filter if isinstance(_filter, list) else [_filter]
    SAMPLES = [s for s in ALL_SAMPLES if any(s.endswith(f) for f in _filter)]
    if not SAMPLES:
        raise ValueError(f"Ingen samples funnet for filter: {_filter}")
else:
    SAMPLES = ALL_SAMPLES
print(f"Kjører {len(SAMPLES)}/{len(ALL_SAMPLES)} samples.")

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
        r1   = os.path.join(RESULTS_DIR, "{sample}/Trimmed/R1.fastq.gz"),
        r2   = os.path.join(RESULTS_DIR, "{sample}/Trimmed/R2.fastq.gz"),
        json = os.path.join(RESULTS_DIR, "{sample}/QC/fastp.json")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/fastp.log")
    threads: THREADS
    shell:
        "pixi run fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -j {output.json} --thread {threads} > {log} 2>&1"

checkpoint identify_species:
    input:
        r1 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/R1.fastq.gz"),
        r2 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/R2.fastq.gz")
    output:
        sp      = os.path.join(RESULTS_DIR, "{sample}/species.txt"),
        rep     = os.path.join(RESULTS_DIR, "{sample}/ID_Kraken2/kraken2_report.txt"),
        bracken = os.path.join(RESULTS_DIR, "{sample}/ID_Kraken2/bracken_species.txt")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/identify_species.log")
    params:
        db = KRAKEN2_DB
    shell:
        """
        pixi run --environment identification kraken2 --db {params.db} --memory-mapping --paired --report {output.rep} {input.r1} {input.r2} > /dev/null 2>> {log}
        pixi run --environment identification bracken -d {params.db} -i {output.rep} -o {output.bracken} -l S >> {log} 2>&1
        RESULT=$(tail -n +2 {output.bracken} | sort -t$'\\t' -k6 -rn | head -n1 | cut -f1 | sed 's/ /_/g')
        echo "${{RESULT:-Unknown}}" > {output.sp}
        echo "Identifisert art: ${{RESULT:-Unknown}}" >> {log}
        """

rule assemble:
    input:
        r1 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/R1.fastq.gz"),
        r2 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/R2.fastq.gz")
    output:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/assemble.log")
    params:
        outdir = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "Assembly")
    threads: THREADS
    shell:
        "pixi run shovill --R1 {input.r1} --R2 {input.r2} --outdir {params.outdir} --cpus {threads} --force > {log} 2>&1"

rule fastqc:
    input:
        r1 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/R1.fastq.gz"),
        r2 = os.path.join(RESULTS_DIR, "{sample}/Trimmed/R2.fastq.gz")
    output:
        os.path.join(RESULTS_DIR, "{sample}/QC/fastqc/R1_fastqc.html")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/fastqc.log")
    params:
        outdir = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "QC/fastqc")
    threads: 2
    shell:
        "pixi run fastqc {input.r1} {input.r2} -o {params.outdir} --threads {threads} > {log} 2>&1"

rule quast:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/QUAST/report.html")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/quast.log")
    params:
        outdir = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "QUAST")
    shell:
        "pixi run quast {input.fa} -o {params.outdir} > {log} 2>&1"

rule mlst:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa"),
        sp = os.path.join(RESULTS_DIR, "{sample}/species.txt")
    output:
        os.path.join(RESULTS_DIR, "{sample}/MLST/mlst.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/mlst.log")
    params:
        scheme = lambda wc: f"--scheme {MLST_SCHEME[sp]}" if (sp := read_species(wc.sample)) in MLST_SCHEME else ""
    shell:
        "pixi run mlst {params.scheme} {input.fa} > {output} 2> {log}"

rule kleborate:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa"),
        sp = os.path.join(RESULTS_DIR, "{sample}/species.txt")
    output:
        os.path.join(RESULTS_DIR, "{sample}/Kleborate/kleborate_output.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/kleborate.log")
    params:
        outdir = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "Kleborate"),
        preset = lambda wc: KLEB_PRESET.get(read_species(wc.sample), "kpsc")
    shell:
        """
        pixi run --environment identification kleborate -a {input.fa} -o {params.outdir} -p {params.preset} > {log} 2>&1
        find {params.outdir} -maxdepth 1 -name '*_output.txt' ! -name '*hAMRonization*' -exec mv {{}} {output} \\;
        """

rule spatyper:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/SpaTyper/spatyper.txt")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/spatyper.log")
    shell:
        "pixi run spaTyper -f {input.fa} -o {output} > {log} 2>&1"

rule sccmec:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/SCCmec/sccmec.txt")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/sccmec.log")
    shell:
        "pixi run sccmec --input {input.fa} --output {output} > {log} 2>&1"

rule agrvate:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/AgrVATE/agrvate.txt")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/agrvate.log")
    shell:
        "pixi run agrvate -i {input.fa} -o {output} > {log} 2>&1"

rule emmtyper:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/EmmTyper/emmtyper.txt")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/emmtyper.log")
    shell:
        "pixi run emmtyper {input.fa} > {output} 2> {log}"

rule amrfinder:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa"),
        sp = os.path.join(RESULTS_DIR, "{sample}/species.txt")
    output:
        os.path.join(RESULTS_DIR, "{sample}/AMRFinder/amrfinder.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/amrfinder.log")
    params:
        org = lambda wc: f"--organism {AMR_ORG[sp]}" if (sp := read_species(wc.sample)) in AMR_ORG else ""
    threads: THREADS
    shell:
        "pixi run --environment amrfinder4 amrfinder --nucleotide {input.fa} {params.org} --output {output} --threads {threads} --plus > {log} 2>&1"

rule ectyper:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/ECTyper/ectyper.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/ectyper.log")
    params:
        outdir = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "ECTyper")
    shell:
        "pixi run ectyper -i {input.fa} -o {params.outdir} > {log} 2>&1 && mv {params.outdir}/output.tsv {output}"

rule mob_typer:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/MOBSuite/mobtyper.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/mob_typer.log")
    shell:
        "pixi run --environment mobsuite mob_typer --infile {input.fa} --out_file {output} > {log} 2>&1"

rule mefinder:
    input:
        fa = os.path.join(RESULTS_DIR, "{sample}/Assembly/contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "{sample}/MEfinder/mefinder.tsv")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/mefinder.log")
    params:
        outdir = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "MEfinder")
    shell:
        "pixi run --environment mefinder mefinder find --contig {input.fa} --output {params.outdir}/mefinder > {log} 2>&1 && mv {params.outdir}/mefinder.results.tsv {output}"

rule multiqc:
    input:
        fastp   = os.path.join(RESULTS_DIR, "{sample}/QC/fastp.json"),
        fastqc  = os.path.join(RESULTS_DIR, "{sample}/QC/fastqc/R1_fastqc.html"),
        dynamic = get_multiqc_inputs
    output:
        os.path.join(RESULTS_DIR, "{sample}/MultiQC/multiqc_report.html")
    log:
        os.path.join(RESULTS_DIR, "{sample}/logs/multiqc.log")
    params:
        sample_dir = lambda wc: os.path.join(RESULTS_DIR, wc.sample),
        outdir     = lambda wc: os.path.join(RESULTS_DIR, wc.sample, "MultiQC")
    shell:
        "pixi run multiqc {params.sample_dir} -o {params.outdir} > {log} 2>&1"
