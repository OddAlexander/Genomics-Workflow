# Snakefile -- bacterial genomics pipeline

# --- Configuration ---
configfile: "config.yaml"
include: "rules/common.smk"
from pathlib import Path

DATA_DIR    = config.get("data_dir", "data/").rstrip("/")
RESULTS_DIR = config.get("results_dir", "results/").rstrip("/")
THREADS     = config.get("threads", 8)
KRAKEN2_DB       = config.get("kraken2_db", "/databases/kraken2_db_light/")
KRAKEN2_MEM      = config.get("kraken2_mem_mb", 12000)  # MB RAM reserved per job -- caps parallel jobs via --resources mem_mb=<available RAM>
SHOVILL_MEM      = config.get("shovill_mem_mb",  12000)  # MB RAM reserved per Shovill/SPAdes job (SPAdes --ram is set to this / 1024)
SKANI_MEM        = config.get("skani_mem_mb", 4000)     # MB RAM reserved per skani job (much lower than fastANI)
SKANI_THREADS    = config.get("skani_threads", 8)
CHECKM_MEM       = config.get("checkm_mem_mb", 22000)   # MB RAM per CheckM job (lineage_wf --reduced_tree, pplacer peaks ~14-18 GB; 22 GB allows 1 job per 30 GB host)
SKANI_DB         = config.get("skani_db", "/databases/skani_db/bacteria")     # Path to skani sketch -- used for species ID (replaces GAMBIT)
PLASMIDFINDER_DB = config.get("plasmidfinder_db", "/databases/plasmidfinder_db/")
LRE_DIR          = config.get("lrefinder_dir", f"{workflow.basedir}/.lrefinder")  # cloned at runtime, same pattern as ReporTree
RUN_LOG          = f"{RESULTS_DIR}/run_log.tsv"

# --- Species groups and schemes ---
SA  = {"Staphylococcus_aureus"}
EC  = {"Escherichia_coli", "Shigella_sonnei", "Shigella_flexneri", "Shigella_dysenteriae", "Shigella_boydii"}
GAS = {"Streptococcus_pyogenes"}
PA  = {"Pseudomonas_aeruginosa"}
SAL = {"Salmonella_enterica"}
HI  = {"Haemophilus_influenzae"}
EFM = {"Enterococcus_faecium", "Enterococcus_faecalis"}      # LRE-Finder targets

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
        expand(f"{RESULTS_DIR}/{{sample}}/pipeline_summary.html", sample=SAMPLES)

rule report:
    input:
        get_all_outputs
    output:
        f"{RESULTS_DIR}/{{sample}}/pipeline_summary.html"
    params:
        template   = "scripts/templates/report_template.html",
        results    = RESULTS_DIR,
        kraken2_db = KRAKEN2_DB,
        log_path   = RUN_LOG
    shell:
        "python scripts/make_report.py --sample {wildcards.sample} --results-dir {params.results} --template {params.template} --output {output} --kraken2-db {params.kraken2_db} --log-path {params.log_path}"

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
            --threads {threads} --paired --report {output.rep} {input.r1} {input.r2} 2>&1 | tee {log}
        pixi run --environment identification bracken -d {params.db} \
            -i {output.rep} -o {output.bracken} -l S 2>&1 | tee -a {log}
        """

checkpoint identify_species:
    # SKANI species ID -- writes both the canonical skani.tsv and species.txt
    # (Genus_species form, first two whitespace-separated tokens of Ref_name).
    # This is a Snakemake checkpoint because species.txt drives DAG routing for
    # the species-specific rules (Kleborate, StaphScope, AMRFinder, ...).
    input:
        fa = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa"
    output:
        sp  = f"{RESULTS_DIR}/{{sample}}/species.txt",
        tsv = f"{RESULTS_DIR}/{{sample}}/ID_Skani/skani.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/identify_species.log"
    params:
        db = SKANI_DB
    threads: SKANI_THREADS
    resources:
        mem_mb = SKANI_MEM
    retries: 2
    shell:
        # Ref_name (column 6) looks like "NC_016845.1 Klebsiella pneumoniae subsp. ..."
        # so the first token is the accession, not the genus. Walk tokens and pick
        # the first "Genus species" pair (Capital+lowercase, lowercase) -> Genus_species.
        # Falls through to "unknown_species" for accession-only / 'sp.' / phage hits.
        """
        pixi run --environment identification skani search \
            -q {input.fa} -d {params.db} -o {output.tsv} -t {threads} 2>&1 | tee {log}
        awk -F'\\t' 'NR==2{{
            n = split($6, a, " ");
            for (i = 1; i < n; i++) {{
                if (a[i] ~ /^[A-Z][a-z]+$/ && a[i+1] ~ /^[a-z]+$/) {{
                    print a[i] "_" a[i+1]; exit;
                }}
            }}
            print "unknown_species";
        }}' {output.tsv} > {output.sp}
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
    # SeqSphere-style depth: map trimmed reads back to the sample's own Shovill
    # assembly with bwa-mem2 and summarise with `samtools coverage`. Because the
    # assembly was built from these exact reads, mapping rate is ~98%+ and the
    # resulting depth/breadth reflect the *actual* sequencing yield -- distinct
    # from the assembly-length estimate (over-counts non-aligning reads) and
    # from the varcall reference-based depth (under-counts when reference differs).
    # All intermediates land in a tmp dir and are wiped at the end -- only the
    # samtools coverage TSV is kept.
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
    # CheckM expects a directory of FASTAs, one per genome/bin. Symlink the
    # assembly under a private input/ so the run produces a single-row report.
    # --reduced_tree drops memory from ~40 GB to ~14 GB for marginal accuracy loss.
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

ECTYPER_MASH = config.get("ectyper_mash",
              "/databases/ectyper/EnteroRef_GTDBSketch_20231003_V2.msh")

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
        pixi run --environment mefinder bash -c '
            FNA=$(python -c "import mgedb, os; print(os.path.join(os.path.dirname(mgedb.__file__), \"data/sequences.d/mge_records.fna\"))")
            makeblastdb -in "$FNA" -dbtype nucl -title mge_records -out {params.db_dir}/mge_records
        ' 2>&1 | tee {log}
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
        # Shovill appends metadata to FASTA headers that mefinder can't parse -- strip with sed.
        # --db-path points at a pre-built shared BLAST DB to avoid the race on /tmp/mge_finder.
        # --temp-dir per sample isolates the per-job working directory.
        # mefinder versions differ in output format: older writes mefinder.csv, newer mefinder.tsv.
        # The prefix passed to `mefinder find` is {params.outdir}/mefinder, so the .tsv path equals {output}.
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
    # Detects linezolid resistance markers in E. faecium / E. faecalis directly
    # from reads (KMA k-mer alignment, no assembly): cfr/optrA/poxtA acquired
    # genes + 23S rRNA G2576T/G2505A + L3/L4 mutations, with mosaicism %.
    # LRE-Finder is not on conda/PyPI; install it once with `scripts/fetch_lrefinder.sh`
    # (clones the repo to .lrefinder/ at the project root). The rule fails loudly
    # with that hint if the script hasn't been run.
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

rule multiqc:
    input:
        fastp   = f"{RESULTS_DIR}/{{sample}}/QC/fastp.json",
        fastqc  = f"{RESULTS_DIR}/{{sample}}/QC/fastqc/R1_fastqc.html",
        # lambda triggers checkpoint resolution via get_all_outputs → read_species → checkpoints.identify_species.get(...)
        dynamic = lambda wc: [f for f in get_all_outputs(wc) if "multiqc_report.html" not in f]
    output:
        f"{RESULTS_DIR}/{{sample}}/MultiQC/multiqc_report.html"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/multiqc.log"
    shell:
        "pixi run multiqc $(dirname $(dirname {output})) -o $(dirname {output}) --force 2>&1 | tee {log}"
