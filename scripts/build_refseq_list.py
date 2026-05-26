#!/usr/bin/env python3
"""Scan a refseq directory and emit refseq_list.json for Snakefile_varcall.

Each top-level subdir of --refseq-root is treated as one species reference,
laid out as written by `datasets download genome accession ...`:

    <refseq-root>/<shortname>/ncbi_dataset/data/data_summary.tsv
    <refseq-root>/<shortname>/ncbi_dataset/data/<accession>/*_genomic.fna
    <refseq-root>/<shortname>/ncbi_dataset/data/<accession>/genomic.gff   (optional)

Species name comes from data_summary.tsv column 1 ("Organism Scientific Name");
the first two words form the underscored "Genus_species" key. If multiple
assembly accessions exist (e.g. GCA_ + GCF_), GCF_ wins. If no GFF is present,
that entry's "gff" field is set to "-" so bcftools csq is skipped at runtime.
"""
import argparse, csv, json, re, sys
from pathlib import Path


SPECIES_RE = re.compile(r"^([A-Z][a-z]+)\s+([a-z]+)")


def parse_species(ref_dir):
    tsv = ref_dir / "ncbi_dataset" / "data" / "data_summary.tsv"
    if not tsv.exists():
        return None
    rows = list(csv.DictReader(tsv.read_text().splitlines(), delimiter="\t"))
    if not rows:
        return None
    m = SPECIES_RE.match(rows[0].get("Organism Scientific Name", ""))
    return f"{m.group(1)}_{m.group(2)}" if m else None


def find_assembly_dir(ref_dir):
    """Pick the preferred assembly. GCF_ (RefSeq) wins over GCA_ (GenBank);
    within a class, highest version (last by name sort) wins."""
    data = ref_dir / "ncbi_dataset" / "data"
    if not data.is_dir():
        return None
    accs = [d for d in data.iterdir() if d.is_dir()]
    gcf = sorted((d for d in accs if d.name.startswith("GCF_")), key=lambda d: d.name)
    gca = sorted((d for d in accs if d.name.startswith("GCA_")), key=lambda d: d.name)
    if gcf:
        return gcf[-1]
    if gca:
        return gca[-1]
    return None


def find_fasta(asm_dir):
    matches = sorted(asm_dir.glob("*_genomic.fna"))
    return matches[0] if matches else None


def find_gff(asm_dir):
    for pattern in ("*.gff3", "*.gff", "genomic.gff", "genomic.gff3"):
        m = sorted(asm_dir.glob(pattern))
        if m:
            return m[0]
    return None


def build_index(root):
    entries, notes = {}, []
    for ref_dir in sorted(root.iterdir()):
        if not ref_dir.is_dir():
            continue
        species = parse_species(ref_dir)
        if not species:
            notes.append(f"SKIP {ref_dir.name}: no parseable species in data_summary.tsv")
            continue
        asm = find_assembly_dir(ref_dir)
        if not asm:
            notes.append(f"SKIP {ref_dir.name}: no GCF_/GCA_ assembly directory")
            continue
        fasta = find_fasta(asm)
        if not fasta:
            notes.append(f"SKIP {ref_dir.name}: no *_genomic.fna under {asm}")
            continue
        gff = find_gff(asm)
        if species in entries:
            notes.append(f"SKIP {ref_dir.name}: duplicate species {species} (already from another dir)")
            continue
        entries[species] = {
            "fasta": str(fasta),
            "gff":   str(gff) if gff else "-",
        }
    return entries, notes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    # Default output sits next to the Snakefiles (project root, one level up from scripts/).
    default_output = Path(__file__).resolve().parent.parent / "refseq_list.json"
    ap.add_argument("--refseq-root", default="/databases/refseq",
                    help="Root containing per-species refseq subdirs (default: %(default)s)")
    ap.add_argument("--output", default=str(default_output),
                    help="Output JSON path (use '-' for stdout, default: %(default)s)")
    args = ap.parse_args()

    root = Path(args.refseq_root)
    if not root.is_dir():
        sys.exit(f"refseq root not found: {root}")

    entries, notes = build_index(root)

    payload = json.dumps(entries, indent=2, ensure_ascii=False) + "\n"
    if args.output == "-":
        sys.stdout.write(payload)
    else:
        Path(args.output).write_text(payload)
        print(f"Wrote {len(entries)} entries to {args.output}", file=sys.stderr)

    for n in notes:
        print(n, file=sys.stderr)

    no_gff = sorted(sp for sp, e in entries.items() if e["gff"] == "-")
    if no_gff:
        print(f"\n{len(no_gff)} species have no GFF (bcftools csq will be skipped):",
              file=sys.stderr)
        for sp in no_gff:
            print(f"  - {sp}", file=sys.stderr)
        print(
            "\nTo enable consequence annotation, drop a GFF3 alongside each *_genomic.fna\n"
            "  e.g. `datasets download genome accession <ACC> --include gff3`\n"
            "and re-run this script.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
