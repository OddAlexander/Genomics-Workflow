# Shared helpers for all four pipelines (main / phylo / cgmlst / varcall).
# Each Snakefile defines its own glob_wildcards() (the input pattern differs --
# main/phylo/varcall glob from data/, cgmlst globs assembled contigs from
# results/) and then calls these helpers.


def assert_unique_samples(all_samples):
    """Raise if glob_wildcards picked up duplicate sample IDs -- shouldn't
    happen with a well-formed data/ tree, and silently failing here would let
    Snakemake later rebuild whichever rule's output came last."""
    seq = list(all_samples)
    if len(set(seq)) != len(seq):
        dups = sorted({s for s in seq if seq.count(s) > 1})
        raise ValueError(f"Duplicate sample IDs: {dups}")


def filter_samples(all_samples, samples_filter):
    """Filter ALL_SAMPLES by the optional `--config samples=...` value.

    Preserves the user-given order so SAMPLES[0] is the first match of token #1
    (downstream rules rely on this -- the phylo auto-reference picker uses
    SAMPLES[0] as its Bakta seed). Deduplicates: a sample matched by more than
    one token only appears once.

    Matching is exact against any path component (date or leaf), or the full
    path -- never a substring match, so '3w' will not match '3wous'.

    samples_filter may be a string or list:
      samples=001k             -> sample whose leaf name is exactly '001k'
      samples=19-03-2026       -> all samples under that date directory
      samples=[3q,3o,3n]       -> these three, in that order
      None / falsy             -> all samples (original glob order)
    """
    if not samples_filter:
        return list(all_samples)
    tokens = samples_filter if isinstance(samples_filter, list) else [samples_filter]
    seen, filtered = set(), []
    for tok in tokens:
        tok = str(tok)
        for s in all_samples:
            parts = s.split("/")
            if (s == tok or tok in parts) and s not in seen:
                filtered.append(s)
                seen.add(s)
    if not filtered:
        raise ValueError(f"No samples found for filter: {samples_filter}")
    return filtered
