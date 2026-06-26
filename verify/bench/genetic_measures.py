#!/usr/bin/env python3
"""
Genetic-measures scorer for MARS benchmarking.

Reads a multiple sequence alignment (PHYLIP sequential or FASTA) and computes
the standard genetic-measure suite used in the MARS paper (Ayad & Pissis 2017):

  L           alignment length (columns)
  PM_sites    polymorphic (segregating) columns: >=2 distinct non-gap chars
  indel_cols  columns containing at least one gap
  AVPD        average pairwise distance = mean over all unordered pairs of the
              number of columns at which the two sequences differ. This is the
              headline metric; calibrated to reproduce the paper's Table 5/6.
  ti/tv/subs  pairwise transition / transversion / substitution counts summed
              over all unordered pairs (both-non-gap, differing columns).

Gap handling for AVPD is configurable; the default ('pairdiff') counts a column
as a difference whenever the two characters differ (gap vs base = differ;
gap vs gap = same). This matches the paper's convention (see CALIBRATION notes).

Usage:
  python3 genetic_measures.py MSA.phylip
  python3 genetic_measures.py --format fasta msa.fasta
  python3 genetic_measures.py --all-variants msa.phy   # print candidate AVPDs
"""
import argparse
import sys
import itertools
import numpy as np

GAP = "-"

# transitions: purine<->purine (A<->G) or pyrimidine<->pyrimidine (C<->T)
_PUR = frozenset("AG")
_PYR = frozenset("CT")


def read_phylip(path):
    """Interleaved PHYLIP: first line '<nseq> <ncol>'; then blocks of nseq
    lines. The first block carries a name token at the start of each line;
    subsequent blocks are sequence-only. Sequence characters are grouped in
    space-separated 10-mer blocks."""
    with open(path) as fh:
        header = fh.readline().split()
        nseq, ncol = int(header[0]), int(header[1])
        lines = [ln.strip() for ln in fh if ln.strip()]
    order = []
    chunks = [[] for _ in range(nseq)]
    for k, ln in enumerate(lines):
        if k < nseq:
            parts = ln.split(None, 1)
            order.append(parts[0])
            body = parts[1] if len(parts) > 1 else ""
        else:
            body = ln
        chunks[k % nseq].append(body.replace(" ", "").replace("\t", ""))
    seqs = {}
    for name, c in zip(order, chunks):
        seqs[name] = "".join(c).upper()
    return order, seqs


def read_fasta(path):
    order, seqs, name, buf = [], {}, None, []
    for ln in open(path):
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf).upper()
            name = ln[1:].split()[0]
            order.append(name)
            buf = []
        else:
            buf.append(ln)
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return order, seqs


def to_matrix(order, seqs):
    n = len(order)
    L = len(seqs[order[0]])
    M = np.empty((n, L), dtype="U1")
    for r, s in enumerate(order):
        sq = seqs[s]
        if len(sq) != L:
            sys.exit("unequal length: %s has %d != %d" % (s, len(sq), L))
        for c in range(L):
            M[r, c] = sq[c] if c < len(sq) else GAP
    return M


def avpd_pairdiff(M):
    """gap counts as difference; gap-gap = same."""
    n = M.shape[0]
    tot = 0
    for i, j in itertools.combinations(range(n), 2):
        tot += int(np.sum(M[i] != M[j]))
    return tot / (n * (n - 1) / 2)


def avpd_variants(M):
    """Return dict of candidate AVPD conventions for calibration."""
    n, L = M.shape
    npairs = n * (n - 1) // 2
    incgap = nogap_mm = nogap_cols = 0  # accumulators
    for i, j in itertools.combinations(range(n), 2):
        a, b = M[i], M[j]
        both_nongap = (a != GAP) & (b != GAP)
        incgap += int(np.sum(a != b))
        nogap_mm += int(np.sum(both_nongap & (a != b)))
        nogap_cols += int(np.sum(both_nongap))
    return {
        "raw_incgap": incgap / npairs,            # mean raw mismatches, gap=diff
        "raw_nogap": nogap_mm / npairs,           # mean raw mismatches, gaps excluded
        "pdist_nogap": (nogap_mm / nogap_cols) if nogap_cols else float("nan"),
        "pdist_incgap": (incgap / (npairs * L)) if L else float("nan"),
        "pdist_nogap_x1000": 1000.0 * (nogap_mm / nogap_cols) if nogap_cols else float("nan"),
    }


def titv_subs(M):
    n = M.shape[0]
    ti = tv = 0
    for i, j in itertools.combinations(range(n), 2):
        a, b = M[i], M[j]
        both = (a != GAP) & (b != GAP) & (a != b)
        for c in np.where(both)[0]:
            x, y = a[c], b[c]
            if ({x, y} <= _PUR) or ({x, y} <= _PYR):
                ti += 1
            else:
                tv += 1
    return ti, tv, ti + tv


def pm_and_indel(M):
    n, L = M.shape
    pm = indel = 0
    for c in range(L):
        col = M[:, c]
        has_gap = np.any(col == GAP)
        if has_gap:
            indel += 1
        nongap = col[col != GAP]
        if nongap.size >= 2 and len(set(nongap.tolist())) >= 2:
            pm += 1
        elif nongap.size < 2 and not has_gap:
            pass
    return pm, indel


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("msa", help="alignment file")
    ap.add_argument("--format", choices=["phylip", "fasta", "auto"], default="auto")
    ap.add_argument("--all-variants", action="store_true",
                    help="print all candidate AVPD conventions (for calibration)")
    args = ap.parse_args()

    fmt = args.format
    if fmt == "auto":
        with open(args.msa) as fh:
            first = fh.readline().split()
        fmt = "phylip" if (len(first) == 2 and first[0].isdigit()) else "fasta"
    order, seqs = (read_phylip if fmt == "phylip" else read_fasta)(args.msa)
    M = to_matrix(order, seqs)
    n, L = M.shape

    if args.all_variants:
        v = avpd_variants(M)
        for k in ("raw_incgap", "raw_nogap", "pdist_nogap", "pdist_incgap",
                  "pdist_nogap_x1000"):
            print("%-18s %.4f" % (k, v[k]))
        return

    avpd = avpd_pairdiff(M)
    pm, indel = pm_and_indel(M)
    ti, tv, subs = titv_subs(M)
    print("n=%d L=%d PM=%d indel_cols=%d AVPD=%.4f ti=%d tv=%d subs=%d"
          % (n, L, pm, indel, avpd, ti, tv, subs))


if __name__ == "__main__":
    main()
