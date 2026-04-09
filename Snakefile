# Snakefile -- Bakteriell genomikk-pipeline

# --- Konfigurasjon ---
configfile: "config.yaml"

DATA_DIR    = config.get("data_dir", "data/").rstrip("/")
RESULTS_DIR = config.get("results_dir", "results/").rstrip("/")
THREADS     = config.get("threads", 8)
KRAKEN2_DB  = config.get("kraken2_db", "/databases/kraken2_db_mini/")
KRAKEN2_MEM = config.get("kraken2_mem_mb", 10000)  # MB RAM reservert per Kraken2-jobb -- begrenser antall samtidige jobber via --resources mem_mb=<tilgjengelig RAM>
GAMBIT_DB   = config.get("gambit_db", "/databases/gambit_db/")
MASH_DB     = config.get("mash_db", "/databases/mash_db/refseq.msh")

# --- Artsgrupper og skjemaer ---
SA  = {"Staphylococcus_aureus"}
EC  = {"Escherichia_coli", "Shigella_sonnei", "Shigella_flexneri", "Shigella_dysenteriae", "Shigella_boydii"}
GAS = {"Streptococcus_pyogenes"}
PA  = {"Pseudomonas_aeruginosa"}
SAL = {"Salmonella_enterica"}
HI  = {"Haemophilus_influenzae"}

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
    "ID_Kraken2/kraken2_report.txt",
    "ID_Kraken2/bracken_species.txt",
    "ID_GAMBIT/gambit.csv",
]

SPECIES_TOOLS = [
    (set(KLEB_PRESET.keys()), ["Kleborate/kleborate_output.tsv"]),
    (SA,                  ["SpaTyper/spatyper.txt", "SCCmec/sccmec.txt", "AgrVATE/agrvate.txt"]),
    (GAS,                 ["EmmTyper/emmtyper.txt"]),
    (EC,                  ["ECTyper/ectyper.tsv"]),
    (PA,                  ["Pasty/pasty.tsv"]),
    (SAL,                 ["SeqSero2/seqsero2.tsv"]),
    (HI,                  ["Hicap/hicap.tsv"]),
]

# --- Hjelpefunksjoner ---
def read_species(sample):
    sp_path = checkpoints.identify_species_early.get(sample=sample).output.sp
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
    return [f"{RESULTS_DIR}/{s}/{t}" for t in tools]

# --- Sample Deteksjon ---
# {sample} fanger hele den relative stien under DATA_DIR, f.eks. "19-03-2026/005a" slik at mappestrukturen i results/ speiler data/
ALL_SAMPLES, ANY = glob_wildcards(f"{DATA_DIR}/{{sample}}/{{any}}_R1.fastq.gz")
READS_DICT = {s: f"{DATA_DIR}/{s}/{a}" for s, a in zip(ALL_SAMPLES, ANY)}
assert len(set(ALL_SAMPLES)) == len(ALL_SAMPLES), \
    f"Duplikate sample-ID-er oppdaget: {len(ALL_SAMPLES) - len(set(ALL_SAMPLES))} duplikater funnet"

_filter = config.get("samples", None)
if _filter:
    _filter = _filter if isinstance(_filter, list) else [_filter]
    SAMPLES = [s for s in ALL_SAMPLES if any(f in s for f in _filter)]
    if not SAMPLES:
        raise ValueError(f"Ingen samples funnet for filter: {_filter}")
else:
    SAMPLES = ALL_SAMPLES
print(f"Kjører {len(SAMPLES)}/{len(ALL_SAMPLES)} samples.")

# --- Rules ---
rule all:
    input:
        expand(f"{RESULTS_DIR}/{{sample}}/pipeline_summary.html", sample=SAMPLES)

rule report:
    input:
        get_all_outputs
    output:
        f"{RESULTS_DIR}/{{sample}}/pipeline_summary.html"
    params:
        template   = "report_template.html",
        results    = RESULTS_DIR,
        kraken2_db = KRAKEN2_DB
    shell:
        "python scripts/make_report.py --sample {wildcards.sample} --results-dir {params.results} --template {params.template} --output {output} --kraken2-db {params.kraken2_db}"

rule fastp:
    input:
        r1 = lambda wc: f"{READS_DICT[wc.sample]}_R1.fastq.gz",
        r2 = lambda wc: f"{READS_DICT[wc.sample]}_R2.fastq.gz"
    output:
        r1   = f"{RESULTS_DIR}/{{sample}}/Trimmed/R1.fastq.gz",
        r2   = f"{RESULTS_DIR}/{{sample}}/Trimmed/R2.fastq.gz",
        json = f"{RESULTS_DIR}/{{sample}}/QC/fastp.json"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/fastp.log"
    threads: THREADS
    shell:
        "pixi run fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -j {output.json} --thread {threads} 2>&1 | tee {log}"

checkpoint identify_species_early:
    input:
        r1 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R1.fastq.gz",
        r2 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R2.fastq.gz"
    output:
        sp = f"{RESULTS_DIR}/{{sample}}/early_species.txt"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/identify_species_early.log"
    params:
        db = MASH_DB
    threads: THREADS
    shell:
        # Extract genus+species from top mash screen hit. 
        """
        pixi run --environment identification mash screen -w -p {threads} {params.db} {input.r1} {input.r2} \
            > {output.sp}.mash_raw 2>{log}
        sort -grk1 {output.sp}.mash_raw \
            | awk -F'\\t' 'NR==1{{print $6}}' \
            | grep -oP '[A-Z][a-z]+ [a-z]+' \
            | awk 'NR==1{{gsub(/ /,"_"); print}}' > {output.sp}
        rm {output.sp}.mash_raw
        """

rule kraken2_qc:
    input:
        r1 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R1.fastq.gz",
        r2 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R2.fastq.gz"
    output:
        rep     = f"{RESULTS_DIR}/{{sample}}/ID_Kraken2/kraken2_report.txt",
        bracken = f"{RESULTS_DIR}/{{sample}}/ID_Kraken2/bracken_species.txt"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/kraken2_qc.log"
    params:
        db = KRAKEN2_DB
    threads: THREADS
    resources:
        mem_mb = KRAKEN2_MEM
    shell:
        """
        pixi run --environment identification kraken2 --db {params.db} \
            --threads {threads} --paired --report {output.rep} {input.r1} {input.r2} 2>&1 | tee {log}
        pixi run --environment identification bracken -d {params.db} \
            -i {output.rep} -o {output.bracken} -l S 2>&1 | tee -a {log}
        """

rule identify_species:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        sp  = f"{RESULTS_DIR}/{{sample}}/species.txt",
        csv = f"{RESULTS_DIR}/{{sample}}/ID_GAMBIT/gambit.csv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/identify_species.log"
    params:
        db = GAMBIT_DB
    shell:
        """
        pixi run --environment identification gambit -d {params.db} query \
            --output {output.csv} --outfmt csv {input.fa} 2>&1 | tee {log}
        awk -F',' 'NR==2{{print $2}}' {output.csv} | sed 's/ /_/g' > {output.sp}
        """

rule assemble:
    input:
        r1 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R1.fastq.gz",
        r2 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R2.fastq.gz"
    output:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/assemble.log"
    params:
        outdir = lambda wc: f"{RESULTS_DIR}/{wc.sample}/Assembly"
    threads: THREADS
    shell:
        "pixi run shovill --R1 {input.r1} --R2 {input.r2} --outdir {params.outdir} --cpus {threads} --force 2>&1 | tee {log}"

rule fastqc:
    input:
        r1 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R1.fastq.gz",
        r2 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R2.fastq.gz"
    output:
        r1 = f"{RESULTS_DIR}/{{sample}}/QC/fastqc/R1_fastqc.html",
        r2 = f"{RESULTS_DIR}/{{sample}}/QC/fastqc/R2_fastqc.html"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/fastqc.log"
    params:
        outdir = lambda wc: f"{RESULTS_DIR}/{wc.sample}/QC/fastqc"
    threads: 2
    shell:
        "pixi run fastqc {input.r1} {input.r2} -o {params.outdir} --threads {threads} 2>&1 | tee {log}"

rule quast:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/QUAST/report.html"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/quast.log"
    params:
        outdir = lambda wc: f"{RESULTS_DIR}/{wc.sample}/QUAST"
    shell:
        "pixi run quast {input.fa} -o {params.outdir} 2>&1 | tee {log}"

rule mlst:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa",
        sp = f"{RESULTS_DIR}/{{sample}}/early_species.txt"
    output:
        f"{RESULTS_DIR}/{{sample}}/MLST/mlst.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/mlst.log"
    run:
        with open(input.sp) as f:
            sp = f.read().strip()
        scheme = f"--scheme {MLST_SCHEME[sp]}" if sp in MLST_SCHEME else ""
        shell("pixi run mlst " + scheme + " {input.fa} > {output} 2>{log}")

rule kleborate:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa",
        sp = f"{RESULTS_DIR}/{{sample}}/species.txt"
    output:
        f"{RESULTS_DIR}/{{sample}}/Kleborate/kleborate_output.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/kleborate.log"
    params:
        outdir = lambda wc: f"{RESULTS_DIR}/{wc.sample}/Kleborate"
    run:
        with open(input.sp) as f:
            sp = f.read().strip()
        preset = KLEB_PRESET.get(sp, "kpsc")
        shell("""
        pixi run --environment identification kleborate -a {input.fa} -o {params.outdir} -p """ + preset + """ 2>&1 | tee {log}
        find {params.outdir} -maxdepth 1 -name '*_output.txt' ! -name '*hAMRonization*' -exec mv {{}} {output} \\;
        """)

rule spatyper:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/SpaTyper/spatyper.txt"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/spatyper.log"
    shell:
        "pixi run spaTyper -f {input.fa} --output {output} 2>&1 | tee {log}"

rule sccmec:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/SCCmec/sccmec.txt"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/sccmec.log"
    params:
        outdir = lambda wc: f"{RESULTS_DIR}/{wc.sample}/SCCmec"
    shell:
        """
        pixi run sccmec --input {input.fa} --outdir {params.outdir} 2>&1 | tee {log}
        mv {params.outdir}/sccmec.tsv {output}
        """

rule agrvate:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/AgrVATE/agrvate.txt"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/agrvate.log"
    resources:
        agrvate_jobs = 1 #bruker samme midlertidig fil, som kan skape konflikter ved parallell kjøring
    shell:
        "pixi run agrvate -i {input.fa} -m -f 2>&1 | tee {log} && mv contigs-results/contigs-summary.tab {output}"

rule emmtyper:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/EmmTyper/emmtyper.txt"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/emmtyper.log"
    shell:
        "pixi run emmtyper {input.fa} > {output} 2>{log}"

rule amrfinder:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa",
        sp = f"{RESULTS_DIR}/{{sample}}/early_species.txt"
    output:
        f"{RESULTS_DIR}/{{sample}}/AMRFinder/amrfinder.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/amrfinder.log"
    threads: THREADS
    run:
        with open(input.sp) as f:
            sp = f.read().strip()
        org = f"--organism {AMR_ORG[sp]}" if sp in AMR_ORG else ""
        shell("pixi run --environment amrfinder4 amrfinder --nucleotide {input.fa} " + org + " --output {output} --threads {threads} --plus 2>&1 | tee {log}")

rule ectyper:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/ECTyper/ectyper.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/ectyper.log"
    params:
        outdir = lambda wc: f"{RESULTS_DIR}/{wc.sample}/ECTyper"
    shell:
        "pixi run ectyper -i {input.fa} -o {params.outdir} 2>&1 | tee {log} && mv {params.outdir}/output.tsv {output}"

rule mob_typer:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/MOBSuite/mobtyper.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/mob_typer.log"
    shell:
        "pixi run --environment mobsuite mob_typer --infile {input.fa} --out_file {output} 2>&1 | tee {log}"

rule mefinder:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/MEfinder/mefinder.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/mefinder.log"
    params:
        outdir = lambda wc: f"{RESULTS_DIR}/{wc.sample}/MEfinder"
    shell:
        # Shovill legger til metadata i FASTA-headers som mefinder ikke takler - strippes med sed.
        # --temp-dir per sample unngår konflikter mellom parallelle jobber på /tmp/mge_finder.
        """
        sed '/^>/s/ .*//' {input.fa} > {params.outdir}/contigs.fa
        pixi run --environment mefinder mefinder find {params.outdir}/mefinder \
            -c {params.outdir}/contigs.fa --temp-dir {params.outdir}/tmp 2>&1 | tee {log}
        mv {params.outdir}/mefinder.csv {output}
        """

rule seqsero2:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/SeqSero2/seqsero2.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/seqsero2.log"
    params:
        outdir  = lambda wc: f"{RESULTS_DIR}/{wc.sample}/SeqSero2",
        fa_dir  = lambda wc: f"{RESULTS_DIR}/{wc.sample}/Assembly"
    shell:
        # SeqSero2 strips path from input — must run from the assembly directory.
        """
        cd {params.fa_dir}
        pixi run --environment seqsero2 SeqSero2_package.py -t 4 -m k \
            -i contigs.fa -d {params.outdir} 2>&1 | tee {log}
        mv {params.outdir}/SeqSero_result.tsv {output}
        """

rule hicap:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/Hicap/hicap.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/hicap.log"
    params:
        outdir = lambda wc: f"{RESULTS_DIR}/{wc.sample}/Hicap"
    shell:
        """
        pixi run --environment hicap hicap -q {input.fa} -o {params.outdir} 2>&1 | tee {log}
        mv {params.outdir}/*.tsv {output}
        """

rule pasty:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/Pasty/pasty.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/pasty.log"
    params:
        outdir = lambda wc: f"{RESULTS_DIR}/{wc.sample}/Pasty"
    shell:
        "pixi run pasty -i {input.fa} -o {params.outdir} -p pasty --force 2>&1 | tee {log}"

rule multiqc:
    input:
        fastp   = f"{RESULTS_DIR}/{{sample}}/QC/fastp.json",
        fastqc  = f"{RESULTS_DIR}/{{sample}}/QC/fastqc/R1_fastqc.html",
        dynamic = lambda wc: [f for f in get_all_outputs(wc) if "multiqc_report.html" not in f]
    output:
        f"{RESULTS_DIR}/{{sample}}/MultiQC/multiqc_report.html"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/multiqc.log"
    params:
        sample_dir = lambda wc: f"{RESULTS_DIR}/{wc.sample}",
        outdir     = lambda wc: f"{RESULTS_DIR}/{wc.sample}/MultiQC"
    shell:
        "pixi run multiqc {params.sample_dir} -o {params.outdir} 2>&1 | tee {log}"
