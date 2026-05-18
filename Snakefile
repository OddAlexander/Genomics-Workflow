# Snakefile -- bacterial genomics pipeline

# --- Configuration ---
configfile: "config.yaml"
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
GAMBIT_DB        = config.get("gambit_db", "/databases/gambit_db/")
SKANI_DB         = config.get("skani_db", "/databases/skani_db/bacteria")     # Path to skani sketch (directory) -- empty = skani is skipped
SKANI_THRESHOLD  = config.get("skani_threshold", 0.3)      # skani runs when GAMBIT closest.distance exceeds this
PLASMIDFINDER_DB = config.get("plasmidfinder_db", "/databases/plasmidfinder_db/")
RUN_LOG          = f"{RESULTS_DIR}/run_log.tsv"

# --- Species groups and schemes ---
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

# --- Tools per species group ---
ALWAYS_TOOLS = [
    "QC/fastp.json",
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
    "ID_GAMBIT/gambit.csv",
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
assert len(set(ALL_SAMPLES)) == len(ALL_SAMPLES), \
    f"Duplicate sample IDs: {[s for s in set(ALL_SAMPLES) if ALL_SAMPLES.count(s) > 1]}"

_filter = config.get("samples", None)
if _filter:
    _filter = _filter if isinstance(_filter, list) else [_filter]
    SAMPLES = [s for s in ALL_SAMPLES if any(f in s for f in _filter)]
    if not SAMPLES:
        raise ValueError(f"No samples found for filter: {_filter}")
else:
    SAMPLES = ALL_SAMPLES
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
        template   = "report_template.html",
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

rule skani:
    input:
        fa     = f"{RESULTS_DIR}/{{sample}}/Assembly/contigs.fa",
        gambit = f"{RESULTS_DIR}/{{sample}}/ID_GAMBIT/gambit.csv"
    output:
        f"{RESULTS_DIR}/{{sample}}/ID_Skani/skani.tsv"
    log:
        f"{RESULTS_DIR}/{{sample}}/logs/skani.log"
    params:
        db        = SKANI_DB,
        threshold = SKANI_THRESHOLD
    threads: SKANI_THREADS
    resources:
        mem_mb = SKANI_MEM
    retries: 2
    run:
        import csv
        row  = next(csv.DictReader(open(input.gambit)))
        dist = float(row.get("closest.distance") or "inf")
        if params.db and dist > params.threshold:
            shell("pixi run --environment identification skani search -q {input.fa} -d {params.db}"
                  " -o {output} -t {threads} 2>&1 | tee {log}")
        else:
            reason = "no database" if not params.db else f"GAMBIT confident (dist={dist:.3f})"
            shell(f"echo -e 'status\\treason\\nskipped\\t{reason}' > {{output}}")

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
        "pixi run multiqc $(dirname $(dirname {output})) -o $(dirname {output}) 2>&1 | tee {log}"
