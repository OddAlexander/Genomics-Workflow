# Snakefile -- bacterial genomics pipeline

# --- Configuration ---
configfile: "config.yaml"
include: "rules/common.smk"
from pathlib import Path

DATA_DIR    = config.get("data_dir", "data/").rstrip("/")
RESULTS_DIR = config.get("results_dir", "results/").rstrip("/")
THREADS     = config.get("threads", 8)
KRAKEN2_DB       = config.get("kraken2_db",     "/databases/kraken2_db_pluspf_16GB/")
KRAKEN2_MEM      = config.get("kraken2_mem_mb", 58000)  # MB RAM reserved per job -- caps parallel jobs via --resources mem_mb=<available RAM>
SHOVILL_MEM      = config.get("shovill_mem_mb",  12000)  # MB RAM reserved per Shovill/SPAdes job (SPAdes --ram is set to this / 1024)
SKANI_MEM        = config.get("skani_mem_mb", 10000)    # MB RAM reserved per skani job; sized for the full-DB fallback peak (~8 GB for GTDB r226) since any sample may trigger it. Small-DB searches use far less.
SKANI_THREADS    = config.get("skani_threads", 8)
CHECKM_MEM       = config.get("checkm_mem_mb", 22000)   # MB RAM per CheckM job (lineage_wf --reduced_tree, pplacer peaks ~14-18 GB; 22 GB allows 1 job per 30 GB host)
SKANI_DB         = config.get("skani_db",         "/databases/skani_db/bacteria")       # ANI confirmation in report
SKANI_DB_FULL    = config.get("skani_db_full",    "/databases/skani_db_full")           # fallback when normal skani ANI is low (dir with markers.bin/sketches.db)
SKANI_ANI_MIN    = config.get("skani_ani_min", 95)                                      # top-hit ANI below this triggers full-DB fallback
GAMBIT_DB        = config.get("gambit_db",        "/databases/gambit_db")              # Primary species ID (drives DAG routing)
ECTYPER_MASH     = config.get("ectyper_mash",     "/databases/ectyper/EnteroRef_GTDBSketch_20231003_V2.msh")
PLASMIDFINDER_DB = config.get("plasmidfinder_db", "/databases/plasmidfinder_db/")
LRE_DIR          = config.get("lrefinder_dir",    f"{workflow.basedir}/.lrefinder")   # cloned once via scripts/fetch_lrefinder.sh
GBS_SBG_DIR      = config.get("gbssbg_dir",       f"{workflow.basedir}/.gbssbg")      # cloned once via scripts/fetch_gbssbg.sh
RUN_LOG          = f"{RESULTS_DIR}/run_log.tsv"

# --- Species groups and schemes ---
SA  = {"Staphylococcus_aureus"}
EC  = {"Escherichia_coli", "Shigella_sonnei", "Shigella_flexneri", "Shigella_dysenteriae", "Shigella_boydii"}
GAS = {"Streptococcus_pyogenes"}
PA  = {"Pseudomonas_aeruginosa"}
SAL = {"Salmonella_enterica"}
HI  = {"Haemophilus_influenzae"}
EFM = {"Enterococcus_faecium", "Enterococcus_faecalis"}      # LRE-Finder targets
GBS = {"Streptococcus_agalactiae"}                           # GBS-SBG serotyping

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

# --- Tools per species group ---
ALWAYS_TOOLS = [
    "QC/fastp.json",
    "QC/self_coverage.tsv",
    "QUAST/report.html",
    "CheckM/quality.tsv",
    "MultiQC/multiqc_report.html",
    "MLST/mlst.tsv",
    "AMRFinder/amrfinder.tsv",
    "MOBSuite/mobtyper.tsv",
    "PlasmidFinder/results_tab.tsv",
    "MEfinder/mefinder.tsv",
    "ID_Kraken2/kraken2_report.txt",
    "ID_Kraken2/bracken_species.txt",
    "ID_Skani/skani.tsv",
    "ID_Gambit/gambit.csv",
]

SPECIES_TOOLS = [
    (set(KLEB_PRESET.keys()), ["Kleborate/kleborate_output.tsv"]),
    (SA,                  ["Staphscope/Staphscope_final_report/staphscope_comprehensive_report.tsv"]),
    (GAS,                 ["EmmTyper/emmtyper.txt"]),
    (EC,                  ["ECTyper/ectyper.tsv"]),
    (PA,                  ["Pasty/pasty.tsv"]),
    (SAL,                 ["SeqSero2/seqsero2.tsv"]),
    (HI,                  ["Hicap/hicap.tsv"]),
    (EFM,                 ["LRE-Finder/lre.res"]),
    (GBS,                 ["GBSSeroTyper/serotype.tsv"]),
]

# --- Helper functions ---
def read_species(sample):
    sp_path = checkpoints.identify_species.get(sample=sample).output.sp
    with open(sp_path, "r") as f:
        species = f.read().strip().replace(" ", "_")
    return species if species else "unknown_species"

def species_flag(mapping, fmt="{}"):
    """Make a params function: reads input.sp, looks it up in mapping, formats with fmt.
    Returns '' on miss so the flag is silently dropped from the shell command."""
    def fn(wildcards, input):
        sp = Path(input.sp).read_text().strip()
        return fmt.format(mapping[sp]) if sp in mapping else ""
    return fn

def get_all_outputs(wildcards):
    sp    = read_species(wildcards.sample)
    s     = wildcards.sample
    tools = list(ALWAYS_TOOLS)
    for species_set, outputs in SPECIES_TOOLS:
        if sp in species_set:
            tools += outputs
    return [f"{RESULTS_DIR}/{s}/{t}" for t in tools]

# --- Sample detection ---
# {sample} captures the full relative path under DATA_DIR, e.g. "19-03-2026/005a",
# so the results/ tree mirrors the data/ tree.
ALL_SAMPLES, ANY = glob_wildcards(f"{DATA_DIR}/{{sample}}/{{any}}_R1.fastq.gz")
READS_DICT = {s: f"{DATA_DIR}/{s}/{a}" for s, a in zip(ALL_SAMPLES, ANY)}
assert_unique_samples(ALL_SAMPLES)
SAMPLES = filter_samples(ALL_SAMPLES, config.get("samples"))
print(f"Running {len(SAMPLES)}/{len(ALL_SAMPLES)} samples.")

# --- Rules ---
rule all:
    input:
        [f"{RESULTS_DIR}/{s}/{Path(s).name}_pipeline_summary.html" for s in SAMPLES],
        [f"{RESULTS_DIR}/{s}/{Path(s).name}_pipeline_summary.pdf"  for s in SAMPLES],

rule report:
    input:
        get_all_outputs
    output:
        html = f"{RESULTS_DIR}/{{sample}}/{{sname}}_pipeline_summary.html",
        pdf  = f"{RESULTS_DIR}/{{sample}}/{{sname}}_pipeline_summary.pdf",
    wildcard_constraints:
        sname = r"[^/]+"
    params:
        template   = "scripts/templates/report_template.html",
        results    = RESULTS_DIR,
        kraken2_db = KRAKEN2_DB,
        skani_db   = SKANI_DB,
        gambit_db  = GAMBIT_DB,
        log_path   = RUN_LOG
    shell:
        "python scripts/make_report.py --sample {wildcards.sample} --results-dir {params.results} --template {params.template} --output {output.html} --pdf {output.pdf} --kraken2-db {params.kraken2_db} --skani-db {params.skani_db} --gambit-db {params.gambit_db} --log-path {params.log_path}"

rule fastp:
    input:
        r1 = lambda wc: f"{READS_DICT[wc.sample]}_R1.fastq.gz",
        r2 = lambda wc: f"{READS_DICT[wc.sample]}_R2.fastq.gz"
    output:
        r1   = f"{RESULTS_DIR}/{{sample}}/Trimmed/R1.fastq.gz",
        r2   = f"{RESULTS_DIR}/{{sample}}/Trimmed/R2.fastq.gz",
        json = f"{RESULTS_DIR}/{{sample}}/QC/fastp.json",
        html = f"{RESULTS_DIR}/{{sample}}/QC/fastp.html"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/fastp.log"
    threads: THREADS
    shell:
        "pixi run fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -j {output.json} -h {output.html} --thread {threads} 2>&1 | tee {log}"

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
        mem_mb     = KRAKEN2_MEM
    shell:
        """
        pixi run --environment identification kraken2 --db {params.db} \
            --threads {threads} --paired --memory-mapping --report {output.rep} {input.r1} {input.r2} 2>&1 | tee {log}
        pixi run --environment identification bracken -d {params.db} \
            -i {output.rep} -o {output.bracken} -l S 2>&1 | tee -a {log}
        """

checkpoint identify_species:
    # GAMBIT (primary) + SKANI (report) species identification.
    # GAMBIT predicted.name -> species.txt, which drives DAG routing 
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        sp         = f"{RESULTS_DIR}/{{sample}}/species.txt",
        gambit_csv = f"{RESULTS_DIR}/{{sample}}/ID_Gambit/gambit.csv",
        tsv        = f"{RESULTS_DIR}/{{sample}}/ID_Skani/skani.tsv",
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/identify_species.log"
    params:
        gambit_db     = GAMBIT_DB,
        skani_db      = SKANI_DB,
        skani_db_full = SKANI_DB_FULL,
        skani_ani_min = SKANI_ANI_MIN,
    threads: SKANI_THREADS
    resources:
        mem_mb = SKANI_MEM
    retries: 2
    shell:
        """
        pixi run --environment identification gambit -d {params.gambit_db} query \
            -o {output.gambit_csv} -f csv {input.fa} 2>&1 | tee {log}
        pixi run --environment identification skani search \
            -q {input.fa} -d {params.skani_db} -o {output.tsv} -t {threads} 2>&1 | tee -a {log}
        ani=$(awk -F'\t' 'NR==2{{print $3}}' {output.tsv})
        if [ -z "$ani" ] || awk "BEGIN{{exit !($ani < {params.skani_ani_min})}}"; then
            printf '\033[1;93m================================================================\033[0m\n' >&2
            printf '\033[1;91m  Normal species ID failed (skani top-hit ANI=%s < %s).\033[0m\n' "${{ani:-none}}" "{params.skani_ani_min}" >&2
            printf '\033[1;93m  Running full skani search against %s ...\033[0m\n' "{params.skani_db_full}" >&2
            printf '\033[1;93m================================================================\033[0m\n' >&2
            pixi run --environment identification skani search \
                -q {input.fa} -d {params.skani_db_full} -o {output.tsv} -t {threads} 2>&1 | tee -a {log}
        fi
        sp=$(awk -F',' 'NR==2{{
            gsub(/"/, "", $2);
            n = split($2, a, " ");
            if (n >= 2) {{ print a[1] "_" a[2]; exit }}
        }}' {output.gambit_csv})
        echo "${{sp:-unknown_species}}" > {output.sp}
        """

rule assemble:
    input:
        r1 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R1.fastq.gz",
        r2 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R2.fastq.gz"
    output:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/assemble.log"
    threads: THREADS
    resources:
        mem_mb = SHOVILL_MEM
    shell:
        "pixi run shovill --R1 {input.r1} --R2 {input.r2} --outdir $(dirname {output.fa}) --cpus {threads} --ram $(( {resources.mem_mb} / 1024 )) --force 2>&1 | tee {log}"

rule fastqc:
    input:
        r1 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R1.fastq.gz",
        r2 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R2.fastq.gz"
    output:
        r1 = f"{RESULTS_DIR}/{{sample}}/QC/fastqc/R1_fastqc.html",
        r2 = f"{RESULTS_DIR}/{{sample}}/QC/fastqc/R2_fastqc.html"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/fastqc.log"
    threads: 2
    shell:
        "pixi run fastqc {input.r1} {input.r2} -o $(dirname {output.r1}) --threads {threads} 2>&1 | tee {log}"

rule quast:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/QUAST/report.html"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/quast.log"
    threads: THREADS
    shell:
        "pixi run quast {input.fa} -o $(dirname {output}) --threads {threads} 2>&1 | tee {log}"

rule self_coverage:
    # Map trimmed reads back to the sample's own Shovill
    # assembly with bwa-mem2 and summarise with `samtools coverage`. Because the
    # assembly was built from these exact reads, mapping rate is ~98%+ and the
    # resulting depth/breadth reflect the *actual* sequencing yield -- distinct
    # from the assembly-length estimate (over-counts non-aligning reads)
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa",
        r1 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R1.fastq.gz",
        r2 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R2.fastq.gz"
    output:
        cov = f"{RESULTS_DIR}/{{sample}}/QC/self_coverage.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/self_coverage.log"
    threads: THREADS
    shell:
        """
        work=$(mktemp -d)
        cp {input.fa} "$work/asm.fa"
        pixi run bwa-mem2 index "$work/asm.fa" 2>&1 | tee {log}
        pixi run bwa-mem2 mem -t {threads} "$work/asm.fa" {input.r1} {input.r2} 2>>{log} \
          | pixi run samtools sort -@ {threads} -O bam -o "$work/aln.bam" - 2>>{log}
        pixi run samtools coverage "$work/aln.bam" > {output.cov} 2>>{log}
        rm -rf "$work"
        """

rule checkm:
    # Requires CheckM data root set once: `pixi run --environment checkm checkm data setRoot /databases/checkm_data`.
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        tsv = f"{RESULTS_DIR}/{{sample}}/CheckM/quality.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/checkm.log"
    threads: THREADS
    resources:
        mem_mb = CHECKM_MEM
    shell:
        """
        outdir=$(dirname {output.tsv})
        rm -rf "$outdir/input" "$outdir/bins" "$outdir/storage"
        mkdir -p "$outdir/input"
        ln -sf "$(readlink -f {input.fa})" "$outdir/input/contigs.fasta"
        pixi run --environment checkm checkm lineage_wf \
            "$outdir/input" "$outdir" -t {threads} -x fasta \
            --reduced_tree --tab_table -f {output.tsv} 2>&1 | tee {log}
        rm -rf "$outdir/input" "$outdir/bins" "$outdir/storage" "$outdir/lineage.ms"
        """

rule mlst:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa",
        sp = f"{RESULTS_DIR}/{{sample}}/species.txt"
    output:
        f"{RESULTS_DIR}/{{sample}}/MLST/mlst.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/mlst.log"
    params:
        scheme = species_flag(MLST_SCHEME, "--scheme {}")
    shell:
        "pixi run mlst {params.scheme} {input.fa} > {output} 2>{log}"

rule kleborate:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa",
        sp = f"{RESULTS_DIR}/{{sample}}/species.txt"
    output:
        f"{RESULTS_DIR}/{{sample}}/Kleborate/kleborate_output.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/kleborate.log"
    params:
        preset = lambda wc, input: KLEB_PRESET.get(Path(input.sp).read_text().strip(), "kpsc")
    shell:
        """
        outdir=$(dirname {output})
        pixi run --environment identification kleborate -a {input.fa} -o "$outdir" -p {params.preset} 2>&1 | tee {log}
        find "$outdir" -maxdepth 1 -name '*_output.txt' ! -name '*hAMRonization*' -exec mv {{}} {output} \\;
        """

rule staphscope:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/Staphscope/Staphscope_final_report/staphscope_comprehensive_report.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/staphscope.log"
    params:
        outdir = lambda wc: f"{RESULTS_DIR}/{wc.sample}/Staphscope"
    threads: THREADS
    shell:
        "pixi run --environment staphscope staphscope -i {input.fa} -o {params.outdir} -t {threads} 2>&1 | tee {log}"

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
        sp = f"{RESULTS_DIR}/{{sample}}/species.txt"
    output:
        f"{RESULTS_DIR}/{{sample}}/AMRFinder/amrfinder.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/amrfinder.log"
    params:
        organism = species_flag(AMR_ORG, "--organism {}")
    threads: THREADS
    shell:
        "pixi run --environment amrfinder4 amrfinder --nucleotide {input.fa} {params.organism} --output {output} --threads {threads} --plus 2>&1 | tee {log}"

rule ectyper:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa",
    output:
        f"{RESULTS_DIR}/{{sample}}/ECTyper/ectyper.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/ectyper.log"
    params:
        mash = ECTYPER_MASH,
    shell:
        """
        outdir=$(dirname {output})
        ref_flag=""
        [ -f {params.mash} ] && [ "$(stat -c%s {params.mash})" -gt 1000000 ] && ref_flag="-r {params.mash}"
        pixi run ectyper -i {input.fa} -o "$outdir" $ref_flag 2>&1 | tee {log}
        mv "$outdir/output.tsv" {output}
        """

rule mob_typer:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/MOBSuite/mobtyper.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/mob_typer.log"
    shell:
        "pixi run --environment mobsuite mob_typer --infile {input.fa} --out_file {output} 2>&1 | tee {log}"

rule plasmidfinder:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/PlasmidFinder/results_tab.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/plasmidfinder.log"
    params:
        db = PLASMIDFINDER_DB
    retries: 2
    shell:
        "pixi run --environment plasmidfinder plasmidfinder.py -i {input.fa} -o $(dirname {output}) -p {params.db} -x 2>&1 | tee {log}"

rule mefinder_db:
    # mefinder's default --db-path is /tmp/mge_finder/database, shared across all
    # parallel jobs -- two jobs racing on makeblastdb corrupt each other's DB.
    # Build it once into a stable per-results location and point every job there.
    output:
        marker = f"{RESULTS_DIR}/.mefinder_db/mge_records.nhr"
    log:
        f"{RESULTS_DIR}/.mefinder_db/build.log"
    params:
        db_dir = f"{RESULTS_DIR}/.mefinder_db"
    shell:
        """
        mkdir -p {params.db_dir}
        FNA=$(pixi run --environment mefinder python -c 'import mgedb, os; print(os.path.join(os.path.dirname(mgedb.__file__), "data/sequences.d/mge_records.fna"))')
        pixi run --environment mefinder makeblastdb -in "$FNA" -dbtype nucl -title mge_records -out {params.db_dir}/mge_records 2>&1 | tee {log}
        """


rule mefinder:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa",
        db = f"{RESULTS_DIR}/.mefinder_db/mge_records.nhr"
    output:
        f"{RESULTS_DIR}/{{sample}}/MEfinder/mefinder.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/mefinder.log"
    params:
        outdir = lambda wc: f"{RESULTS_DIR}/{wc.sample}/MEfinder",
        db_dir = f"{RESULTS_DIR}/.mefinder_db"
    shell:
        """
        sed '/^>/s/ .*//' {input.fa} > {params.outdir}/contigs.fa
        pixi run --environment mefinder mefinder find {params.outdir}/mefinder \
            -c {params.outdir}/contigs.fa --db-path {params.db_dir} \
            --temp-dir {params.outdir}/tmp 2>&1 | tee {log}
        rm -f {params.outdir}/contigs.fa
        rm -rf {params.outdir}/tmp
        [ -f {output} ] || mv {params.outdir}/mefinder.csv {output}
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
    shell:
        """
        outdir=$(dirname {output})
        pixi run --environment hicap hicap -q {input.fa} -o "$outdir" 2>&1 | tee {log}
        mv "$outdir/contigs.tsv" {output}
        """

rule pasty:
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        f"{RESULTS_DIR}/{{sample}}/Pasty/pasty.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/pasty.log"
    shell:
        "pixi run pasty -i {input.fa} -o $(dirname {output}) -p pasty --force 2>&1 | tee {log}"

rule lrefinder:
    input:
        r1 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R1.fastq.gz",
        r2 = f"{RESULTS_DIR}/{{sample}}/Trimmed/R2.fastq.gz",
    output:
        res = f"{RESULTS_DIR}/{{sample}}/LRE-Finder/lre.res",
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/lrefinder.log",
    threads: THREADS
    shell:
        # `-1t1` enables one-to-one read-template matching (recommended for mixed
        # 23S copies); `-cge` selects the CGE-style result formatting.
        """
        if [ ! -f {LRE_DIR}/LRE-Finder.py ]; then
            echo "ERROR: LRE-Finder not installed at {LRE_DIR}/" >&2
            echo "       Run: scripts/fetch_lrefinder.sh" >&2
            exit 1
        fi
        outdir=$(dirname {output.res})
        mkdir -p "$outdir"
        pixi run --environment lrefinder python {LRE_DIR}/LRE-Finder.py \
            -ipe {input.r1} {input.r2} \
            -o "$outdir/lre" \
            -t_db {LRE_DIR}/elmDB/elm \
            -1t1 -cge 2>&1 | tee {log}
        """

rule gbsserotyper:
    # Perl + BLAST+-based; not on conda. Install once: scripts/fetch_gbssbg.sh
    # Output TSV: Name / Serotype / Uncertainty.
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        tsv = f"{RESULTS_DIR}/{{sample}}/GBSSeroTyper/serotype.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/gbsserotyper.log"
    params:
        dir  = GBS_SBG_DIR,
        name = lambda wc: wc.sample.replace("/", "_"),
    shell:
        """
        if [ ! -f {params.dir}/GBS-SBG.pl ]; then
            echo "ERROR: GBS-SBG not installed at {params.dir}/" >&2
            echo "       Run: scripts/fetch_gbssbg.sh" >&2
            exit 1
        fi
        mkdir -p "$(dirname {output.tsv})"
        pixi run --environment identification perl {params.dir}/GBS-SBG.pl \
            {input.fa} -name {params.name} -best \
            -ref {params.dir}/GBS-SBG.fasta 2>&1 | tee {log} > {output.tsv}
        """

rule multiqc:
    input:
        fastp   = f"{RESULTS_DIR}/{{sample}}/QC/fastp.json",
        fastqc  = f"{RESULTS_DIR}/{{sample}}/QC/fastqc/R1_fastqc.html",       
        dynamic = lambda wc: [f for f in get_all_outputs(wc) if "multiqc_report.html" not in f]
    output:
        f"{RESULTS_DIR}/{{sample}}/MultiQC/multiqc_report.html"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/multiqc.log"
    shell:
        "pixi run multiqc $(dirname $(dirname {output})) -o $(dirname {output}) --force 2>&1 | tee {log}"
